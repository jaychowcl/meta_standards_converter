# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Wall-clock bounds and physical-attempt accounting for metadata streams."""
from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass, field
import socket
from threading import Lock, Thread, Event
import time

import requests


class MetadataResponseTimeout(requests.Timeout):
    """The complete metadata response exceeded its wall-clock allowance."""


class BatchPayloadFailure(Exception):
    """A batch payload should be retried with fewer explicitly requested IDs."""


class BatchRequestLimit(Exception):
    """An original XML batch exhausted its shared recovery allowance."""


class _UnchangedRequestLimit(BatchRequestLimit):
    """The unchanged wire request has already used its three retries."""


@dataclass
class AttemptAllowance:
    maximum: int = 8
    used: int = 0
    requests: dict = field(default_factory=dict, repr=False)

    def charge(self, identity=None):
        if self.used >= self.maximum:
            raise BatchRequestLimit('XML batch physical request allowance exhausted')
        if identity is not None:
            if self.requests.get(identity, 0) >= 4:
                raise _UnchangedRequestLimit('Unchanged XML request allowance exhausted')
            self.requests[identity] = self.requests.get(identity, 0) + 1
        self.used += 1


attempt_allowance = ContextVar('msc_xml_attempt_allowance', default=None)
split_payload = ContextVar('msc_split_xml_payload', default=False)


def charge_attempt(url=None, params=None):
    allowance = attempt_allowance.get()
    if allowance is not None:
        # Request identity is private and never emitted; params may contain a key.
        allowance.charge((url, repr(params)) if url is not None else None)


def _socket(response):
    """Only real sockets may be interrupted; injected clients need no socket."""
    raw = getattr(response, 'raw', None)
    connection = getattr(raw, '_connection', None)
    candidate = getattr(connection, 'sock', None)
    if isinstance(candidate, socket.socket):
        return candidate
    # HTTP/1.0 and Connection: close detach the connection before body reading.
    fp = getattr(getattr(raw, '_fp', None), 'fp', None)
    candidate = getattr(getattr(fp, 'raw', None), '_sock', None)
    return candidate if isinstance(candidate, socket.socket) else None


@contextmanager
def response_deadline(response, seconds):
    """Interrupt a trickling socket, including waits inside iter_content()."""
    from ..metadata.retrieval_budget import _active, _check_cancelled, _cancelled
    budget = _active.get()
    started = getattr(response, '_msc_request_started_at', None)
    started = started if isinstance(started, (int, float)) else time.monotonic()
    deadline = started + seconds
    if budget is not None and budget.deadline is not None:
        deadline = min(deadline, time.monotonic() + budget.remaining())
    lock, state = Lock(), {'active': True, 'expired': False}
    sock = _socket(response)
    cancelled, finished = _cancelled.get(), Event()

    def expire():
        with lock:
            if not state['active']:
                return
            state['expired'] = True
            if sock is not None:
                try:
                    sock.shutdown(socket.SHUT_RDWR)
                except OSError:
                    pass

    def check():
        _check_cancelled()
        if budget is not None:
            budget.check()
        if state['expired'] or time.monotonic() >= deadline:
            raise MetadataResponseTimeout('metadata response deadline exceeded')

    def monitor():
        while not finished.wait(min(.1, max(0, deadline-time.monotonic()))):
            if time.monotonic() >= deadline or (cancelled is not None and cancelled.is_set()):
                expire()
                return
    watcher = Thread(target=monitor, name='msc-metadata-deadline', daemon=True)
    watcher.start()
    try:
        check()
        try:
            yield check
        except Exception:
            check()
            raise
        check()
    finally:
        with lock:
            state['active'] = False
        finished.set()
        watcher.join()
