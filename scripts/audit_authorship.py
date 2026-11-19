# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Reproducible MSC author/contact and complete MAGE-TAB audit at a fixed SHA.

Capture is opt-in, bounded to 300 seconds/200 physical requests per input across
all presets. The default replay never falls back to live network access.
"""
import argparse
from dataclasses import asdict
import json
import hashlib
from pathlib import Path
import shutil
import subprocess
import sys
import time

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / 'src'))
from scripts.audit_recorded_provider import primitive
from tests.support.authorship_capture import Capture, CaptureLimit, utc, write_json
from tests.support.recorded_replay import ReplayTransport

BASELINE = '942ef678e9625ed61b93438c02fd65b58079f466'
CASES = [('geo', 'GSE100', 'E-GEOD-100'), ('geo', 'GSE60450', 'E-GEOD-60450'),
         ('geo', 'GSE52564', 'E-GEOD-52564'), ('ae', 'E-MTAB-6486', 'E-MTAB-6486'),
         ('sra', 'SRP045534', 'E-GEOD-60450'), ('ena', 'ERP000263', 'E-MTAB-308'),
         ('sra', 'SRP002056', None), ('ena', 'SRP002056', None), ('ena', 'DRP000001', None)]

SETTINGS = {'resource_profile': 'standard', 'resource_overrides': {'network_workers': 1},
            'expand_studies': False, 'input_options': 'MSC defaults (GEO remove_empty=True)',
            'output_options': 'MSC defaults; automatic platform selection',
            'metadata_timeout_seconds': 'remaining aggregate 300-second capture budget',
            'metadata_max_requests': '200 minus physical requests already used by this input',
            'saved_json_export': 'enrichment off, expansion/linked metadata/peer provider disabled, empty replay manifest'}


def policy(mode):
    from meta_standards_converter.converters import PreparationPolicy
    return PreparationPolicy(enrichment=mode, expand_studies=False,
                             linked_metadata=mode == 'standard', peer_provider=mode == 'standard')


def runtime(seconds=300, requests_limit=200):
    from meta_standards_converter.converters import RuntimeOptions
    return RuntimeOptions(metadata_timeout_seconds=max(.1, seconds), metadata_max_requests=max(1, requests_limit),
                          resource_overrides={'network_workers': 1})


def save_result(result, directory):
    directory.mkdir(parents=True, exist_ok=True)
    write_json(directory / 'result.json', result.to_dict())
    for index, item in enumerate(result.items):
        write_json(directory / f'payload-{index}.json', primitive(item.payload))
        if item.payload and isinstance(item.payload, list) and isinstance(item.payload[0], list):
            from meta_standards_converter.magetab.writer import MAGETabWriter
            for number, document in enumerate(item.payload):
                if any(isinstance(row, list) and row and row[0] == 'SDRF File' for row in document):
                    MAGETabWriter().write(document, str(directory / f'magetab-{index}-{number}'))
    return result.to_dict()


def capture_case(root, provider, accession, reference, *, revision=BASELINE):
    from meta_standards_converter.converters import Converter
    dest = root / 'cases' / (provider + '-' + accession)
    if (dest / 'acquisition.json').exists():
        raise FileExistsError(f'Frozen capture already exists: {dest}')
    cap = Capture(dest / 'http', accession, revision)
    outcomes = {'provider': provider, 'accession': accession, 'reference': reference, 'tested_revision': revision,
                'started_at': utc(), 'limits': {'seconds': 300, 'physical_requests': 200},
                'effective_settings': SETTINGS, 'phases': {}}
    with cap.active():
        for mode in ('off', 'curators', 'standard'):
            cap.phase = mode
            started = time.monotonic()
            try:
                cap.check()
                result = Converter().convert(accession, out_type='json', in_type='ae_accession' if provider == 'ae' else provider,
                    preparation=policy(mode), runtime=runtime(cap.remaining(), 200 - cap.physical_requests))
                outcomes['phases'][mode] = {'result': save_result(result, dest / mode / 'source-json'),
                    'elapsed_seconds': time.monotonic() - started, 'policy': asdict(policy(mode))}
            except CaptureLimit as exc:
                outcomes['phases'][mode] = {'limitation': str(exc), 'elapsed_seconds': time.monotonic() - started}
                break
            except Exception as exc:
                outcomes['phases'][mode] = {'limitation': type(exc).__name__ + ': ' + str(exc),
                    'elapsed_seconds': time.monotonic() - started}
            finally:
                cap.manifest(mode)
                write_json(dest / 'acquisition.json', outcomes)
    outcomes.update(finished_at=utc(), physical_requests=cap.physical_requests,
                    elapsed_seconds=time.monotonic() - cap.started)
    write_json(dest / 'acquisition.json', outcomes)
    print(json.dumps({'case': dest.name, 'requests': cap.physical_requests,
                      'seconds': round(outcomes['elapsed_seconds'], 2),
                      'phases': {key: value.get('result', {}).get('status', value.get('limitation'))
                                 for key, value in outcomes['phases'].items()}}), flush=True)


def replay_case(dest, *, fresh_json=False):
    from meta_standards_converter.converters import Converter
    acquisition = json.loads((dest / 'acquisition.json').read_text())
    outcomes = {}
    for mode, phase in acquisition['phases'].items():
        manifest = dest / 'http' / (mode + '.manifest.json')
        stage = dest / mode
        result_info = {}
        for output in ('json', 'magetab'):
            replay = ReplayTransport(manifest)
            try:
                with replay.offline():
                    result = Converter().convert(acquisition['accession'], out_type=output,
                        in_type='ae_accession' if acquisition['provider'] == 'ae' else acquisition['provider'], preparation=policy(mode), runtime=runtime())
                result_info[output] = save_result(result, stage / ('replayed-' + output))
                replay.assert_complete()
                if fresh_json and output == 'json':
                    save_result(result, stage / 'source-json')
            except BaseException as exc:
                if isinstance(exc, (KeyboardInterrupt, SystemExit)): raise
                result_info[output + '_limitation'] = type(exc).__name__ + ': ' + str(exc)
            write_json(stage / ('trace-' + output + '.json'), {'events': replay.events, 'errors': replay.errors})
        payload = stage / 'source-json/payload-0.json'
        if payload.exists():
            for index, package in enumerate(json.loads(payload.read_text()) or []):
                saved = stage / f'saved-{index}.json'
                write_json(saved, package)
                # Empty manifest: no DNS or HTTP is allowed during saved JSON export.
                empty = dest / 'http/empty.manifest.json'
                if not empty.exists():
                    doc = json.loads(manifest.read_text())
                    doc.update(files={}, requests=[])
                    write_json(empty, doc)
                isolation = ReplayTransport(empty)
                try:
                    with isolation.offline():
                        converted = Converter().convert(saved, out_type='magetab', preparation=policy('off'), runtime=runtime())
                    result_info[f'saved-{index}'] = save_result(converted, stage / f'saved-{index}-magetab')
                    isolation.assert_complete()
                except Exception as exc:
                    result_info[f'saved-{index}_limitation'] = type(exc).__name__ + ': ' + str(exc)
        outcomes[mode] = result_info
    write_json(dest / 'replay.json', outcomes)
    print('Replayed ' + dest.name, flush=True)


def deposited(root, selected=None, *, sealed_references=None, revision=BASELINE):
    from meta_standards_converter.converters import Converter
    from tests.support.deposited_fidelity import DEPOSITED_PANEL, load_deposited_reference, without_derived_files, without_added_protocol_placeholders, assert_deposited_sdrf
    from tests.support.authorship_fidelity import idf_comparison, table, scoped_sdrf_cells
    references = root / 'references'
    references.mkdir(exist_ok=True)
    original = REPO / 'tests/fixtures/fidelity'
    if sealed_references is None:
        shutil.copyfile(original / 'provenance.json', references / 'original-provenance.json')
    manifest = root / 'empty.manifest.json'
    write_json(manifest, {'schema_version': '1.0', 'cohort': 'deposited', 'source_revision': revision,
                         'captured_at': utc(), 'evidence_kind': 'recorded-http', 'files': {}, 'requests': []})
    reports = json.loads((root / 'deposited-results.json').read_text()) if selected and (root / 'deposited-results.json').exists() else []
    for accession in [*DEPOSITED_PANEL, 'E-MTAB-6486']:
        if selected and accession not in selected: continue
        if sealed_references is not None:
            directory = sealed_references
            reference = {kind: table(directory / f'{accession}.{kind}.txt') for kind in ('idf', 'sdrf')}
        elif accession in DEPOSITED_PANEL:
            reference = load_deposited_reference(original, accession)
            directory = original
        else:
            directory = REPO / 'tests/fixtures/studies/E-MTAB-6486/inputs'
            reference = {kind: table(directory / f'{accession}.{kind}.txt') for kind in ('idf', 'sdrf')}
            shutil.copyfile(directory.parent / 'provenance.json', references / 'E-MTAB-6486.provenance.json')
        source = references / (accession + '.idf.txt')
        for kind in ('idf', 'sdrf'):
            shutil.copyfile(directory / f'{accession}.{kind}.txt', references / f'{accession}.{kind}.txt')
        dest = root / 'deposited' / accession
        report = {'accession': accession, 'reference_sha256': {
            kind: hashlib.sha256((references / f'{accession}.{kind}.txt').read_bytes()).hexdigest()
            for kind in ('idf', 'sdrf')}, 'routes': {}}
        isolation = ReplayTransport(manifest)
        with isolation.offline():
            loaded = Converter().convert(source, out_type='json', preparation=policy('off'))
            save_result(loaded, dest / 'ae2json')
            package = dest / 'saved.json'
            if loaded.items and loaded.items[0].payload:
                write_json(package, primitive(loaded.items[0].payload[0]))
            for name, input_value in [('direct', source), ('via-json', package)]:
                outcome = Converter().convert(input_value, out_type='magetab', preparation=policy('off'))
                save_result(outcome, dest / name)
                detail = {'status': outcome.status, 'execution': outcome.items[0].execution,
                          'validation': outcome.items[0].validation, 'completeness': outcome.items[0].completeness}
                if outcome.items[0].payload:
                    document = primitive(outcome.items[0].payload[0])
                    sdrf = without_added_protocol_placeholders(reference['idf'], document)
                    detail['idf'] = idf_comparison(reference['idf'], document, documented=True, source_sdrf=reference['sdrf'])
                    projected = without_derived_files(reference['sdrf'])
                    try:
                        detail['sdrf'] = {'verdict': 'preserved', **assert_deposited_sdrf(projected, sdrf, accession=accession)}
                    except AssertionError as exc:
                        detail['sdrf'] = {'verdict': 'discrepancy', 'reason': str(exc)}
                    detail['sdrf']['source_populated_cells'] = len(scoped_sdrf_cells(reference['sdrf'], accession)[1])
                    detail['sdrf']['projected_populated_cells'] = len(scoped_sdrf_cells(projected, accession)[1])
                    detail['sdrf']['generated_populated_cells'] = len(scoped_sdrf_cells(sdrf, accession)[1])
                report['routes'][name] = detail
        isolation.assert_complete()
        write_json(dest / 'comparison.json', report)
        reports = [previous for previous in reports if previous['accession'] != accession] + [report]
        print('Deposited ' + accession, flush=True)
    write_json(root / 'deposited-results.json', reports)


def frozen_unicode(root, *, sealed=False):
    from pytest import MonkeyPatch
    from meta_standards_converter.converters import Converter
    from tests.support.contracts import GEO, install_replay
    from tests.support.authorship_sources import citation_facts
    import xml.etree.ElementTree as ET
    dest = root / 'frozen-GSE328265'
    for directory in ('inputs', 'provider_responses'):
        shutil.copytree(GEO / directory, dest / directory, dirs_exist_ok=True)
    shutil.copyfile(GEO / ('original-provenance.json' if sealed else 'provenance.json'), dest / 'original-provenance.json')
    raw = GEO / 'provider_responses/pubmed.xml'
    doc = {'kind': 'xml', 'xml': ET.fromstring(raw.read_bytes()), 'url': 'frozen esummary 42129775',
           'sha256': hashlib.sha256(raw.read_bytes()).hexdigest(), 'body': 'provider_responses/pubmed.xml',
           'member': None, 'reference_assisted': False}
    expected = citation_facts([doc])['42129775'][0]
    reports = {}
    for mode in ('off', 'curators', 'standard'):
        reports[mode] = {}
        for output in ('json', 'magetab'):
            isolation = ReplayTransport(root / 'empty.manifest.json')
            with isolation.offline(), MonkeyPatch.context() as patches:
                calls = install_replay(patches)
                result = Converter().convert('GSE328265', in_type='geo', out_type=output, preparation=policy(mode),
                    output_options={'platform_handler': 'single_cell_sequencing'} if output == 'magetab' else None)
                save_result(result, dest / mode / output)
            isolation.assert_complete()
            reports[mode][output] = {'status': result.status, 'requests': calls}
            if output == 'json' and result.items[0].payload:
                saved = dest / mode / 'saved.json'
                write_json(saved, primitive(result.items[0].payload[0]))
                with ReplayTransport(root / 'empty.manifest.json').offline():
                    converted = Converter().convert(saved, out_type='magetab', preparation=policy('off'),
                        output_options={'platform_handler': 'single_cell_sequencing'})
                save_result(converted, dest / mode / 'via-json')
                reports[mode]['via-json'] = {'status': converted.status}
            if output == 'magetab' and result.items[0].payload:
                rows = {row[0]: row[1:] for row in result.items[0].payload[0]}
                authors = rows.get('Publication Author List', [])
                reports[mode]['authors'] = authors
                reports[mode]['author_identity_order_unicode'] = (authors == [expected['author_list']]) if mode != 'off' else 'not hydrated by policy'
    write_json(dest / 'assessment.json', {'scope': 'Frozen one-sample contract; not a full 31-sample live study',
        'expected_from_raw_pubmed': expected, 'modes': reports})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['capture', 'replay', 'deposited', 'frozen', 'report', 'candidate', 'candidate-review'])
    parser.add_argument('--bundle', type=Path, required=True)
    parser.add_argument('--source-bundle', type=Path, help='Sealed baseline input; required for candidate')
    parser.add_argument('--candidate-revision', help='Expected conversion revision; must equal HEAD for a new candidate run')
    parser.add_argument('--candidate-bundle', type=Path, help='Sealed conversion outputs for analysis-only candidate-review')
    parser.add_argument('--case', action='append', help='provider-accession; omit for the declared cohort')
    args = parser.parse_args()
    if args.action in {'candidate', 'candidate-review'}:
        if not args.source_bundle or not args.candidate_revision or args.case:
            parser.error('candidate requires --source-bundle and --candidate-revision and always replays the full cohort')
        from scripts.audit_authorship_candidate import run, review
        if args.action == 'candidate-review':
            if not args.candidate_bundle:
                parser.error('candidate-review requires --candidate-bundle')
            review(args.source_bundle, args.candidate_bundle, args.bundle, args.candidate_revision)
        else:
            if args.candidate_bundle:
                parser.error('--candidate-bundle is only valid for candidate-review')
            run(args.source_bundle, args.bundle, args.candidate_revision)
        return
    if args.source_bundle or args.candidate_revision or args.candidate_bundle:
        parser.error('candidate options are not allowed for baseline actions')
    # A test-only audit commit may advance HEAD; production source must still
    # match the declared baseline before any evidence can be labelled with it.
    changed = subprocess.run(['git', 'diff', '--quiet', BASELINE, '--', 'src', 'pyproject.toml'], cwd=REPO).returncode
    untracked = subprocess.check_output(['git', 'ls-files', '--others', '--exclude-standard', '--', 'src'], cwd=REPO, text=True)
    if changed or untracked:
        raise RuntimeError('Production source differs from audit baseline 942ef67; refusing misleading evidence labels.')
    args.bundle.mkdir(parents=True, exist_ok=True)
    if args.action == 'capture':
        import meta_standards_converter
        environment = {'baseline': BASELINE, 'head': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=REPO, text=True).strip(),
                       'python': sys.executable, 'version': sys.version, 'msc_origin': str(meta_standards_converter.__file__),
                       'captured_at': utc(), 'cohort': CASES}
        write_json(args.bundle / 'environment.json', environment)
        for provider, accession, reference in CASES:
            if args.case and provider + '-' + accession not in args.case: continue
            capture_case(args.bundle, provider, accession, reference)
    elif args.action == 'deposited':
        deposited(args.bundle, args.case)
    elif args.action == 'frozen':
        frozen_unicode(args.bundle)
    elif args.action == 'report':
        from scripts.audit_authorship_report import build
        build(args.bundle)
    else:
        for dest in sorted((args.bundle / 'cases').iterdir()):
            if not (dest / 'acquisition.json').exists() or (args.case and dest.name not in args.case): continue
            replay_case(dest)


if __name__ == '__main__':
    main()
