# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Bounded new-live integrity audit, with isolated repeated route measurements.

Production stays at 2368e51. Historical guards and evidence are never rewritten.
Capture timings include HTTP evidence recording; payload writing is separate.
Offline measurements suppress transport pacing only, and never fetch a fallback.
"""
import argparse
from collections import Counter
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import resource
import statistics
import subprocess
import sys
import time
from unittest.mock import patch

REPO = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPO), str(REPO / 'src')]
from scripts.audit_authorship import policy, runtime, save_result
from scripts.audit_authorship_report import read
from scripts.audit_authorship_candidate import verify_seal, sha256
from scripts.audit_additional_conversions import clean_revision, initialize, reference_inputs
from scripts.audit_recorded_provider import primitive
from tests.support.authorship_capture import Capture, CaptureLimit, utc, write_json
from tests.support.recorded_replay import ReplayTransport

PRODUCTION = '2368e51f0367a8eccf490e1c633e2638b63ed6ff'
MODES = ('off', 'curators', 'standard')
COHORT = {'schema_version': '1.0', 'production_revision': PRODUCTION,
    'mirror_accession': 'SRP012607', 'inputs': [
        {'provider': 'geo', 'accession': 'GSE11817', 'reference': 'E-GEOD-11817', 'evidence': 'live'},
        {'provider': 'geo', 'accession': 'GSE26884', 'reference': 'E-GEOD-26884', 'evidence': 'live'},
        {'provider': 'geo', 'accession': 'GSE37703', 'reference': 'E-GEOD-37703', 'evidence': 'live'},
        {'provider': 'sra', 'accession': 'SRP012607', 'reference': 'E-GEOD-37703', 'evidence': 'live'},
        {'provider': 'ena', 'accession': 'SRP012607', 'reference': 'E-GEOD-37703', 'evidence': 'live'}]}


class Events:
    def __init__(self):
        self.events = []

    def emit(self, event_type, **values):
        self.events.append({'event_type': event_type, **values})


def identity(spec):
    return spec['provider'] + '-' + spec['accession']


def payload_hash(payload):
    return hashlib.sha256(json.dumps(primitive(payload), sort_keys=True, ensure_ascii=False).encode()).hexdigest()


def request_counts(events):
    physical = [e for e in events if e.get('physical')]
    identities = Counter(json.dumps(e['identity'], sort_keys=True) for e in physical)
    return {'physical_requests': len(physical), 'reused_responses': len(events) - len(physical),
        'repeated_physical_identities': sum(max(0, n - 1) for n in identities.values()),
        'response_statuses': dict(Counter(str(e.get('attempt', {}).get('status',
            e.get('attempt', {}).get('exception', 'unfinished'))) for e in physical))}


def summary(records):
    values = [r['wall_seconds'] for r in records]
    return {'repetitions': len(values), 'median_seconds': statistics.median(values),
        'min_seconds': min(values), 'max_seconds': max(values),
        'median_cpu_seconds': statistics.median(r['cpu_seconds'] for r in records),
        'peak_rss_bytes': max(r['peak_rss_bytes'] for r in records),
        'payloads_identical': len({r['payload_sha256'] for r in records}) == 1,
        'outcomes_identical': len({json.dumps(r['items'], sort_keys=True) for r in records}) == 1}


def capture(root, spec, revision):
    from meta_standards_converter.converters import Converter
    from meta_standards_converter.conversion_logging import operation
    from meta_standards_converter.helpers.request_helper import HostRequestGate
    case = root / 'cases' / identity(spec)
    if case.exists():
        raise FileExistsError(case)
    cap = Capture(case / 'http', spec['accession'], revision)
    outcomes = dict(spec, tested_revision=revision, started_at=utc(),
        limits={'seconds': 300, 'physical_requests': 200}, phases={})
    original_wait = HostRequestGate.wait
    waits = []
    def wait(gate, *args, **kwargs):
        started = time.perf_counter()
        try:
            return original_wait(gate, *args, **kwargs)
        finally:
            waits.append(time.perf_counter() - started)
    with cap.active(), patch.object(HostRequestGate, 'wait', wait):
        for mode in MODES:
            cap.phase = mode
            emitter = Events(); started = time.perf_counter(); cpu = time.process_time()
            start_events, start_waits = len(cap.events), len(waits)
            detail = {'policy': asdict(policy(mode)), 'timing_kind': 'live acquisition with explicit prior-phase response reuse',
                'capture_io_included_in_call': True, 'artifact_write_seconds': 0.0}
            outcomes['phases'][mode] = detail
            try:
                cap.check()
                with operation('live_integrity_capture', heartbeat=False, emitter=emitter):
                    result = Converter().convert(spec['accession'], in_type=spec['provider'], out_type='json',
                        preparation=policy(mode), runtime=runtime(cap.remaining(), 200-cap.physical_requests))
                detail['conversion_call_seconds'] = time.perf_counter() - started
                detail['cpu_seconds'] = time.process_time() - cpu
                write_start = time.perf_counter()
                detail['result'] = save_result(result, case / mode / 'captured-json')
                detail['artifact_write_seconds'] = time.perf_counter() - write_start
            except BaseException as exc:
                if isinstance(exc, (KeyboardInterrupt, SystemExit)):
                    raise
                detail.setdefault('conversion_call_seconds', time.perf_counter() - started)
                detail.setdefault('cpu_seconds', time.process_time() - cpu)
                detail['limitation'] = type(exc).__name__ + ': ' + str(exc)
            finally:
                detail['elapsed_seconds'] = time.perf_counter() - started
                detail['gate_wait_seconds'] = sum(waits[start_waits:])
                detail['request_counts'] = request_counts(cap.events[start_events:])
                detail['peak_rss_bytes'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024
                detail['events'] = emitter.events
                cap.manifest(mode)
                write_json(case / 'acquisition.json', outcomes)
            print(identity(spec) + '/' + mode + ': ' + json.dumps({k: detail.get(k) for k in
                ('conversion_call_seconds', 'artifact_write_seconds', 'request_counts', 'limitation')}, default=str), flush=True)
            if cap.remaining() <= 0 or cap.physical_requests >= 200:
                break
    outcomes.update(finished_at=utc(), elapsed_seconds=time.monotonic()-cap.started,
        physical_requests=cap.physical_requests, request_counts=request_counts(cap.events))
    write_json(case / 'acquisition.json', outcomes)


def measure(root, spec, mode, route, revision):
    """One process per cell: one warm-up and three measured conversions."""
    case = root / 'cases' / identity(spec); stage = case / mode
    target = stage / ('timing-' + route + '.json')
    if target.exists():
        raise FileExistsError(target)
    acquisition = read(case / 'acquisition.json')
    phase = acquisition['phases'].get(mode, {})
    limit = phase.get('limitation')
    manifest = case / 'http' / (mode + '.manifest.json')
    if not phase or limit or not manifest.exists():
        write_json(target, {'unavailable': limit or 'Capture phase outside aggregate input budget'})
        return
    value, input_type, output_type, selected = spec['accession'], spec['provider'], route, policy(mode)
    if route == 'saved-0':
        source = stage / 'source-json/payload-0.json'
        if not source.exists() or not read(source):
            write_json(target, {'unavailable': 'Fresh JSON route has no usable payload'})
            return
        value = stage / 'saved-0.json'; write_json(value, read(source)[0])
        input_type, output_type, selected = 'json', 'magetab', policy('off')
        manifest = root / 'empty.manifest.json'
    with ReplayTransport(root / 'empty.manifest.json').offline():
        from meta_standards_converter.converters import Converter
        from meta_standards_converter.conversion_logging import operation
        from meta_standards_converter.helpers.request_helper import HostRequestGate
        records, first = [], None
        for number in range(4):
            replay = ReplayTransport(manifest); emitter = Events()
            wall_start, cpu_start = time.perf_counter(), time.process_time()
            error, result = None, None
            try:
                with replay.offline(), patch.object(HostRequestGate, 'wait', return_value=0.0):
                    with operation('live_integrity_replay', heartbeat=False, emitter=emitter):
                        result = Converter().convert(value, in_type=input_type, out_type=output_type,
                            preparation=selected, runtime=runtime())
                wall, cpu = time.perf_counter()-wall_start, time.process_time()-cpu_start
                replay.assert_complete()
            except BaseException as exc:
                if isinstance(exc, (KeyboardInterrupt, SystemExit)):
                    raise
                wall, cpu = time.perf_counter()-wall_start, time.process_time()-cpu_start
                error = type(exc).__name__ + ': ' + str(exc)
            record = {'iteration': number, 'warmup': number == 0, 'wall_seconds': wall,
                'cpu_seconds': cpu, 'peak_rss_bytes': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
                'rss_scope': 'process lifetime including imports and warm-up', 'error': error,
                'status': result.status if result else None,
                'items': [{k: getattr(i, k) for k in ('execution', 'completeness', 'validation')} for i in result.items] if result else [],
                'payload_sha256': payload_hash([i.payload for i in result.items]) if result else None,
                'events': emitter.events, 'request_trace': replay.events, 'replay_errors': replay.errors}
            records.append(record)
            if number == 1 and result is not None and error is None:
                directory = stage / ({'json': 'source-json', 'magetab': 'replayed-magetab', 'saved-0': 'saved-0-magetab'}[route])
                write_start = time.perf_counter(); first = save_result(result, directory)
                record['artifact_write_seconds'] = time.perf_counter()-write_start
                if route in ('json', 'magetab'):
                    write_json(stage / ('trace-' + route + '.json'), {'events': replay.events, 'errors': replay.errors})
        detail = {'case': identity(spec), 'mode': mode, 'route': route, 'revision': revision,
            'live_requests': 0, 'network': 'strict replay; no live fallback',
            'pacing': 'transport gate suppressed for offline compute measurement only',
            'timed_scope': 'Converter call including event collection, excluding artifact writes and manifest setup',
            'manifest_sha256': sha256(manifest), 'records': records,
            'summary': summary(records[1:]), 'result': first}
        write_json(target, detail)
        print(identity(spec) + '/' + mode + '/' + route + ': ' + json.dumps(detail['summary']), flush=True)


def initialize_bundle(root, revision):
    initialize(root, COHORT, revision)
    parent = root.parent
    seals = {name: verify_seal(parent / name) for name in
        ('run-20261119', 'fix-validation', 'additional-conversions', 'metadata-integrity-validation')}
    environment = read(root / 'environment.json')
    environment.update(historical_seals_before=seals, source_tree=subprocess.check_output(
        ['git', 'rev-parse', 'HEAD:src'], cwd=REPO, text=True).strip(),
        protected_refs_before=subprocess.check_output(['git', 'show-ref'], cwd=REPO, text=True).splitlines())
    write_json(root / 'environment.json', environment)


def collect_measurements(root):
    rows = []
    for spec in COHORT['inputs']:
        case = root / 'cases' / identity(spec); outcomes = {}
        for mode in MODES:
            outcomes[mode] = {}
            for route in ('json', 'magetab', 'saved-0'):
                path = case / mode / ('timing-' + route + '.json')
                detail = read(path) if path.exists() else {'unavailable': 'Measurement not executed'}
                rows.append(dict(case=identity(spec), mode=mode, route=route, **{k: detail[k]
                    for k in ('summary', 'unavailable') if k in detail}))
                if detail.get('result'):
                    outcomes[mode][route] = detail['result']
                else:
                    outcomes[mode][route + '_limitation'] = detail.get('unavailable') or [r['error'] for r in detail.get('records', [])]
        write_json(case / 'replay.json', outcomes)
    write_json(root / 'runtime-summary.json', rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('init', 'capture', 'references', 'measure', 'collect'))
    parser.add_argument('--bundle', type=Path, required=True)
    parser.add_argument('--revision', required=True)
    parser.add_argument('--case')
    parser.add_argument('--mode', choices=MODES)
    parser.add_argument('--route', choices=('json', 'magetab', 'saved-0'))
    args = parser.parse_args()
    revision = clean_revision(args.revision, PRODUCTION)
    root = args.bundle.resolve()
    if (root / 'SHA256SUMS').exists():
        raise ValueError('Sealed evidence is immutable')
    selected = [s for s in COHORT['inputs'] if identity(s) == args.case]
    if args.action in ('capture', 'measure') and len(selected) != 1:
        raise ValueError('Exactly one explicit case is required')
    if args.action == 'init':
        initialize_bundle(root, revision)
    elif args.action == 'capture':
        capture(root, selected[0], revision)
    elif args.action == 'references':
        reference_inputs(root, COHORT, revision)
    elif args.action == 'measure':
        if not args.mode or not args.route:
            raise ValueError('Explicit mode and route required')
        measure(root, selected[0], args.mode, args.route, revision)
    else:
        collect_measurements(root)


if __name__ == '__main__':
    main()
