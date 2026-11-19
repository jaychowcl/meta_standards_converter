# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Explicit additional cohort; bounded capture and independent offline assessment.

Never changes production, historical evidence, or baseline audit guards. Native
snapshots are parser evidence, not HTTP acquisition. Historical golden values
are retained for provenance but are never used as a scientific oracle here.
"""
import argparse
from collections import Counter
import copy
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from unittest.mock import patch

REPO = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPO), str(REPO / 'src')]
from scripts.audit_authorship import capture_case, policy, runtime, save_result
from scripts.audit_authorship_candidate import sha256, seal_result
from scripts.audit_authorship_report import build, read
from scripts.audit_recorded_provider import build_services, primitive
from tests.support.authorship_capture import Capture, utc, write_json
from tests.support.recorded_replay import ReplayTransport
from tests.support.authorship_fidelity import table, idf_comparison, scoped_sdrf_cells
from tests.support.deposited_fidelity import without_added_protocol_placeholders, without_derived_files, assert_deposited_sdrf

COHORT = REPO / 'tests/fixtures/recorded_provider/additional/cohort.json'
MODES = ('off', 'curators', 'standard')


def clean_revision(expected, production):
    def git(*args): return subprocess.check_output(['git', *args], cwd=REPO, text=True).strip()
    actual = git('rev-parse', 'HEAD')
    if actual != git('rev-parse', '--verify', expected + '^{commit}'):
        raise ValueError('Expected revision differs from tested HEAD')
    if git('status', '--porcelain', '--untracked-files=all'):
        raise ValueError('Commit audit code and tests before running')
    if git('diff', production, '--', 'src', 'pyproject.toml'):
        raise ValueError('Production code differs from approved revision')
    return actual


def validate_cohort(cohort):
    if cohort.get('schema_version') != '1.0': raise ValueError('Unknown cohort schema')
    seen = set()
    for spec in cohort['inputs']:
        key = (spec['provider'], spec['accession'])
        if key in seen: raise ValueError('Duplicate input')
        seen.add(key)
        if spec['provider'] not in {'geo', 'sra', 'ena'} or spec['evidence'] not in {'live', 'parser-merge-replay'}:
            raise ValueError('Unsupported route')
    return cohort


def initialize(root, cohort, revision):
    if root.exists(): raise FileExistsError('Use a fresh bundle')
    root.mkdir(parents=True)
    write_json(root / 'cohort.json', cohort)
    import meta_standards_converter
    write_json(root / 'environment.json', {'tested_revision': revision,
        'production_revision': cohort['production_revision'], 'production_tree': subprocess.check_output(
            ['git', 'rev-parse', 'HEAD:src'], cwd=REPO, text=True).strip(),
        'python': sys.executable, 'version': sys.version, 'msc_origin': str(meta_standards_converter.__file__),
        'created_at': utc(), 'limits': {'seconds_per_input': 300, 'requests_per_input': 200},
        'settings': {'presets': MODES, 'network_workers': 1, 'expand_studies': False,
                     'saved_json': 'fresh JSON, enrichment off, empty network manifest'}})
    (root / 'cases').mkdir()
    write_json(root / 'empty.manifest.json', {'schema_version': '1.0', 'cohort': 'additional-offline',
        'source_revision': revision, 'captured_at': utc(), 'evidence_kind': 'recorded-http', 'files': {}, 'requests': []})


def reference_inputs(root, cohort, revision):
    from meta_standards_converter.sources.magetab import AEWebFetcher
    refs = root / 'references'; refs.mkdir(exist_ok=True)
    for accession in dict.fromkeys(spec['reference'] for spec in cohort['inputs']):
        dest = root / 'reference-acquisition' / accession
        if (dest / 'result.json').exists(): raise FileExistsError(dest)
        if accession.startswith('E-MTAB-'):
            for kind in ('idf', 'sdrf'):
                shutil.copyfile(REPO / 'tests/fixtures/fidelity' / f'{accession}.{kind}.txt', refs / f'{accession}.{kind}.txt')
            result = {'kind': 'frozen-deposited-reference', 'live_requests': 0}
        else:
            cap = Capture(dest / 'http', accession, revision); cap.phase = 'reference'
            try:
                with cap.active():
                    fetched = AEWebFetcher().resolve(accession)
                    for item in (fetched.idf, *fetched.sdrfs):
                        # Decoded UTF-8 metadata is saved; the HTTP capture retains original response bytes.
                        (refs / item.name).write_text(item.text)
                result = {'kind': 'live-deposited-reference', 'live_requests': cap.physical_requests}
            except BaseException as exc:
                if isinstance(exc, (KeyboardInterrupt, SystemExit)): raise
                result = {'limitation': type(exc).__name__ + ': ' + str(exc), 'live_requests': cap.physical_requests}
            finally: cap.manifest('reference')
        result['files'] = {p.name: sha256(p) for p in refs.glob(accession + '.*.txt')}
        write_json(dest / 'result.json', result)
        print('Reference ' + accession + ': ' + json.dumps(result), flush=True)


def frozen_input(root, spec, archive, revision):
    case = root / 'cases' / (spec['provider'] + '-' + spec['accession'])
    if case.exists(): raise FileExistsError(case)
    panel = REPO / 'tests/fixtures/recorded_provider/full-panel'
    snapshots = []
    for mode in MODES:
        origin = panel / (spec['reference'] + ('.native' if mode == 'off' else '.standard') + '.manifest.json')
        doc = read(origin)
        write_json(case / 'original-manifests' / (mode + '.json'), doc)
        for file_id, item in doc['files'].items():
            rel = item['path']
            source = archive / rel.removeprefix('.dev/comms/') if rel.startswith('.dev/comms/') else REPO / rel
            if sha256(source) != item['sha256']: raise ValueError('Historical evidence hash mismatch: ' + str(source))
            target = case / 'http' / 'bodies' / (item['sha256'] + source.suffix)
            target.parent.mkdir(parents=True, exist_ok=True)
            if not target.exists(): shutil.copyfile(source, target)
            item['path'] = str(target.relative_to(case / 'http'))
            if file_id == 'records' and mode == 'off':
                snapshots.append({'url': 'frozen-parser-record:' + spec['accession'], 'body': str(target.relative_to(case)),
                                  'sha256': item['sha256'], 'original_path': str(source)})
        doc.pop('external_evidence', None)
        write_json(case / 'http' / (mode + '.manifest.json'), doc)
    write_json(case / 'source-snapshots.json', snapshots)
    write_json(case / 'acquisition.json', dict(spec, physical_requests=0, elapsed_seconds=0,
        tested_revision=revision, phases={mode: {'evidence': 'parser-merge-replay'} for mode in MODES}))


def consumed_manifest(replay):
    """Derive strict request closure from a successful offline selection probe."""
    replay.assert_no_errors()
    doc = copy.deepcopy(replay.manifest)
    selected = []
    for request in doc['requests']:
        # Match the public request identity without relying on private canonical serialization.
        count = sum(event.get('identity') == request['identity'] for event in replay.events)
        attempts = []
        for attempt in request['attempts']:
            take = min(count, attempt['repeat'])
            if take: attempts.append(dict(attempt, repeat=take))
            count -= take
        if count: raise ValueError('Selected requests exceed frozen attempts')
        if attempts: selected.append(dict(request, attempts=attempts))
    doc['requests'] = selected
    doc['selection_policy'] = 'Observed offline requests only; unused original responses retained; no live fallback'
    return doc


def replay_input(root, spec):
    from meta_standards_converter.converters import Converter
    from meta_standards_converter.helpers.request_helper import HostRequestGate
    case = root / 'cases' / (spec['provider'] + '-' + spec['accession'])
    acquisition = read(case / 'acquisition.json'); results = {}
    frozen = spec['evidence'] == 'parser-merge-replay'
    for mode in MODES:
        stage = case / mode; phase = {}; results[mode] = phase
        manifest = case / 'http' / (mode + '.manifest.json')
        unavailable = acquisition['phases'].get(mode, {}).get('limitation')
        if not manifest.exists() or unavailable:
            for route in ('json', 'magetab', 'saved-0'):
                phase[route + '_limitation'] = 'Acquisition unavailable: ' + (unavailable or 'phase outside input budget')
            write_json(case / 'replay.json', results)
            continue
        for output in ('json', 'magetab'):
            try:
                replay = ReplayTransport(manifest)
                with replay.offline(), patch.object(HostRequestGate, 'wait', return_value=0.0):
                    services = build_services(replay, replay.manifest['scenarios'][0]['services']) if frozen else {}
                    outcome = Converter(services=services).convert(spec['accession'], in_type=spec['provider'],
                        out_type=output, preparation=policy(mode), runtime=runtime())
                replay.assert_no_errors()
                if frozen:
                    selected = case / 'http' / f'{mode}-{output}-selected.manifest.json'
                    write_json(selected, consumed_manifest(replay))
                    write_json(stage / f'selection-{output}.json', {'events': replay.events, 'errors': replay.errors})
                    replay = ReplayTransport(selected)
                    with replay.offline(), patch.object(HostRequestGate, 'wait', return_value=0.0):
                        services = build_services(replay, replay.manifest['scenarios'][0]['services'])
                        checked = Converter(services=services).convert(spec['accession'], in_type=spec['provider'],
                            out_type=output, preparation=policy(mode), runtime=runtime())
                    if primitive(outcome.items[0].payload) != primitive(checked.items[0].payload):
                        raise AssertionError('Strict replay differs from frozen request-selection probe')
                    outcome = checked
                replay.assert_complete()
                phase[output] = save_result(outcome, stage / ('replayed-' + output))
                if output == 'json':
                    previous = stage / 'source-json'
                    if previous.exists(): previous.rename(stage / 'captured-json')
                    save_result(outcome, previous)
            except Exception as exc:
                phase[output + '_limitation'] = type(exc).__name__ + ': ' + str(exc)
            finally:
                if 'replay' in locals(): write_json(stage / ('trace-' + output + '.json'), {'events': replay.events, 'errors': replay.errors})
        payload = stage / 'source-json/payload-0.json'
        if payload.exists():
            for index, package in enumerate(read(payload) or []):
                saved = stage / f'saved-{index}.json'; write_json(saved, package)
                try:
                    replay = ReplayTransport(root / 'empty.manifest.json')
                    with replay.offline():
                        outcome = Converter().convert(saved, out_type='magetab', preparation=policy('off'), runtime=runtime())
                    replay.assert_complete()
                    phase[f'saved-{index}'] = save_result(outcome, stage / f'saved-{index}-magetab')
                except Exception as exc: phase[f'saved-{index}_limitation'] = type(exc).__name__ + ': ' + str(exc)
        write_json(case / 'replay.json', results)
        print(case.name + '/' + mode + ': ' + json.dumps({k: v.get('status') if isinstance(v, dict) else v for k, v in phase.items()}), flush=True)


def deposited_inputs(root, cohort):
    from meta_standards_converter.converters import Converter
    reports = []
    for accession in dict.fromkeys(spec['reference'] for spec in cohort['inputs']):
        source = root / 'references' / (accession + '.idf.txt')
        reference = {kind: table(root / 'references' / f'{accession}.{kind}.txt') for kind in ('idf', 'sdrf')}
        dest = root / 'deposited' / accession; report = {'accession': accession, 'routes': {}}
        isolation = ReplayTransport(root / 'empty.manifest.json')
        with isolation.offline():
            loaded = Converter().convert(source, out_type='json', preparation=policy('off'))
            save_result(loaded, dest / 'ae2json')
            saved = dest / 'saved.json'; write_json(saved, primitive(loaded.items[0].payload[0]))
            for route, value in [('direct', source), ('via-json', saved)]:
                result = Converter().convert(value, out_type='magetab', preparation=policy('off'))
                save_result(result, dest / route)
                document = primitive(result.items[0].payload[0])
                detail = {k: getattr(result.items[0], k) for k in ('execution', 'validation', 'completeness')}
                detail['status'] = result.status
                detail['idf'] = idf_comparison(reference['idf'], document, documented=True, source_sdrf=reference['sdrf'])
                actual = next(r[1] for r in document if r[0] == 'SDRF File')
                projected = without_derived_files(reference['sdrf'])
                try:
                    actual = without_added_protocol_placeholders(reference['idf'], document)
                    detail['sdrf'] = {'verdict': 'preserved', **assert_deposited_sdrf(projected, actual, accession=accession)}
                except AssertionError as exc: detail['sdrf'] = {'verdict': 'discrepancy', 'reason': str(exc)}
                for key, rows in [('source', reference['sdrf']), ('projected', projected), ('generated', actual)]:
                    detail['sdrf'][key + '_populated_cells'] = len(scoped_sdrf_cells(rows, accession)[1])
                report['routes'][route] = detail
        isolation.assert_complete(); reports.append(report)
        write_json(dest / 'comparison.json', report)
        write_json(root / 'deposited-results.json', reports)
        print('Deposited ' + accession, flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['init', 'capture', 'references', 'frozen', 'replay', 'deposited', 'report'])
    parser.add_argument('--bundle', type=Path, required=True)
    parser.add_argument('--revision', required=True)
    parser.add_argument('--cohort', type=Path, default=COHORT)
    parser.add_argument('--archive-root', type=Path)
    parser.add_argument('--case')
    args = parser.parse_args()
    cohort = validate_cohort(read(args.cohort)); revision = clean_revision(args.revision, cohort['production_revision'])
    root = args.bundle.resolve()
    if (root / 'SHA256SUMS').exists(): raise ValueError('Sealed output is immutable')
    if args.action == 'init': initialize(root, cohort, revision); return
    if read(root / 'cohort.json') != cohort: raise ValueError('Cohort differs from initialized bundle')
    selected = [s for s in cohort['inputs'] if not args.case or s['provider'] + '-' + s['accession'] == args.case]
    if args.case and not selected: raise ValueError('Unknown case')
    write_json(root / 'actions' / (utc().replace(':', '-') + '.json'), {'action': args.action, 'case': args.case, 'tested_revision': revision})
    if args.action == 'references': reference_inputs(root, cohort, revision)
    elif args.action == 'deposited': deposited_inputs(root, cohort)
    elif args.action == 'report': build(root, cohort=cohort)
    else:
        for spec in selected:
            if args.action == 'capture' and spec['evidence'] == 'live':
                capture_case(root, spec['provider'], spec['accession'], spec['reference'], revision=revision)
            elif args.action == 'frozen' and spec['evidence'] == 'parser-merge-replay':
                if not args.archive_root: raise ValueError('Explicit archive root required')
                frozen_input(root, spec, args.archive_root, revision)
            elif args.action == 'replay': replay_input(root, spec)


if __name__ == '__main__': main()
