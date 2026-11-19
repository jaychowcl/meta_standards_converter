# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Bounded public-response capture at the existing HTTP adapter seam.

No production imports, network shortcuts, pacing changes or credentials in the
trace. Bodies are the original application bytes (after HTTP content decoding).
Repeated phases reuse an earlier successful body, explicitly marked in traces.
"""
from collections import defaultdict
from contextlib import contextmanager
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import signal
import time
from unittest.mock import patch

import requests
from tests.support.recorded_replay import request_identity


class CaptureLimit(BaseException):
    """Stop the entire input, including provider retry loops, at its budget."""


def utc():
    return datetime.now(timezone.utc).isoformat()


def write_json(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False, default=str) + '\n')


class Capture:
    def __init__(self, root, accession, revision, seconds=300, requests_limit=200):
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)
        self.accession, self.revision = accession, revision
        self.seconds, self.requests_limit = seconds, requests_limit
        self.started = time.monotonic()
        self.captured_at = utc()
        self.physical_requests = 0
        self.phase = None
        self.events, self.files, self.cache = [], {}, {}
        self.original_send = requests.adapters.HTTPAdapter.send

    def remaining(self):
        return max(0, self.seconds - (time.monotonic() - self.started))

    def check(self):
        if self.remaining() <= 0 or self.physical_requests >= self.requests_limit:
            raise CaptureLimit('aggregate per-input 300-second/200-request ceiling reached')

    def response(self, request, attempt):
        response = requests.Response()
        response.status_code = attempt['status']
        response.headers.update(attempt['headers'])
        response.encoding = attempt['encoding']
        response._content = (self.root / self.files[attempt['response']]['path']).read_bytes()
        response._content_consumed = True
        response.url, response.request = request.url, request
        return response

    def send(self, adapter, request, **kwargs):
        self.check()
        identity = request_identity(request)
        key = json.dumps(identity, sort_keys=True)
        event = {'phase': self.phase, 'identity': identity, 'started_at': utc(), 'physical': True}
        self.events.append(event)
        cached = self.cache.get(key)
        if cached and cached[0] != self.phase:
            event.update(physical=False, attempt=dict(cached[1]), reused_from_phase=cached[0], finished_at=utc())
            return self.response(request, event['attempt'])
        self.physical_requests += 1
        try:
            response = self.original_send(adapter, request, **kwargs)
        except (requests.Timeout, requests.ConnectionError) as exc:
            event.update(attempt={'exception': 'timeout' if isinstance(exc, requests.Timeout) else 'connection-error', 'repeat': 1}, finished_at=utc())
            self.flush()
            raise
        event['attempt'] = {'repeat': 1, 'status': response.status_code, 'encoding': response.encoding,
                            'headers': {k: v for k, v in response.headers.items() if k.lower() in
                                        {'content-type', 'content-length', 'etag', 'last-modified', 'retry-after', 'location'}}}
        original_iter = response.iter_content
        original_close = response.close
        parts = []
        finished = False

        def finalize(complete):
            nonlocal finished
            if finished:
                return
            finished = True
            content = b''.join(parts)
            digest = hashlib.sha256(content).hexdigest()
            name = f'bodies/{digest}.body'
            path = self.root / name
            path.parent.mkdir(exist_ok=True)
            if not path.exists():
                path.write_bytes(content)
            self.files[digest] = {'path': name, 'sha256': digest, 'bytes': len(content)}
            event['attempt']['response'] = digest
            event.update(finished_at=utc(), complete_body=complete)
            if complete and response.status_code == 200:
                self.cache[key] = (self.phase, dict(event['attempt']))
            self.flush()

        def iterate(chunk_size=1, decode_unicode=False):
            # MSC reads byte chunks. Do not silently change a caller requesting
            # decoded text; requests' normal .text path calls .content first.
            if decode_unicode:
                raise ValueError('audit capture requires byte iteration')
            try:
                for chunk in original_iter(chunk_size=chunk_size, decode_unicode=False):
                    parts.append(chunk)
                    yield chunk
                finalize(True)
            except BaseException:
                finalize(False)
                raise

        def close():
            finalize(response._content_consumed)
            original_close()

        response.iter_content, response.close = iterate, close
        return response

    def flush(self):
        write_json(self.root / 'trace.json', {'captured_at': self.captured_at,
                   'physical_requests': self.physical_requests, 'events': self.events})

    def manifest(self, phase):
        requests_by_key = defaultdict(list)
        identities = {}
        incomplete = []
        for event in self.events:
            if event['phase'] != phase:
                continue
            key = json.dumps(event['identity'], sort_keys=True)
            identities[key] = event['identity']
            attempt = event.get('attempt')
            if not attempt or ('exception' not in attempt and
                               ('response' not in attempt or event.get('complete_body') is False)):
                incomplete.append(event['identity'])
                continue
            requests_by_key[key].append(attempt)
        doc = {'schema_version': '1.0', 'cohort': self.accession + '/' + phase,
               'source_revision': self.revision, 'captured_at': self.captured_at,
               'evidence_kind': 'recorded-http', 'transport': 'http-adapter',
               'files': self.files, 'requests': [{'identity': identities[key], 'attempts': attempts}
                   for key, attempts in requests_by_key.items()], 'incomplete_requests': incomplete}
        write_json(self.root / (phase + '.manifest.json'), doc)

    @contextmanager
    def active(self):
        old_handler = signal.getsignal(signal.SIGALRM)
        def alarm(signum, frame):
            raise CaptureLimit('aggregate per-input wall-clock ceiling reached')
        signal.signal(signal.SIGALRM, alarm)
        signal.setitimer(signal.ITIMER_REAL, self.remaining())
        try:
            with patch.object(requests.adapters.HTTPAdapter, 'send', self.send_adapter):
                yield self
        finally:
            signal.setitimer(signal.ITIMER_REAL, 0)
            signal.signal(signal.SIGALRM, old_handler)
            self.flush()

    @property
    def send_adapter(self):
        def send(adapter, request, **kwargs):
            return self.send(adapter, request, **kwargs)
        return send
