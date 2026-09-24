# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================

"""Full-panel coverage cannot be established by self-reported index counts."""
import hashlib
import json
from pathlib import Path

import pytest

from scripts import audit_recorded_provider as runner


# Independent approved roster, deliberately not imported from the runner.
NATIVE = ['E-GEOD-48034', 'E-GEOD-48035', 'E-GEOD-52564', 'E-GEOD-60450',
          'E-MTAB-13131', 'E-MTAB-13662', 'E-MTAB-14416', 'E-MTAB-14560',
          'E-MTAB-14566', 'E-MTAB-14960', 'E-MTAB-16192', 'E-MTAB-16253',
          'E-MTAB-16856', 'E-MTAB-308', 'E-MTAB-5481', 'E-MTAB-5485']
ENA = ['E-GEOD-48034', 'E-GEOD-60450', 'E-MTAB-14560', 'E-MTAB-14566',
       'E-MTAB-16253', 'E-MTAB-308']
ROSTER = [(s, r) for s in NATIVE for r in ['reference', 'native', 'linked', 'standard']]
ROSTER += [(s, 'ena') for s in ENA]
ROSTER += [(s, r) for s in ['E-MTAB-16847', 'GSE100'] for r in ['reference', 'direct']]
FORMATS = ['json', 'magetab', 'tsv']
ACCESSIONS = dict(zip(NATIVE, ['SRP026125', 'SRP026126', 'SRP033200', 'SRP045534',
    'ERP148790', 'ERP156296', 'ERP163837', 'ERP165490', 'ERP165595', 'ERP170611',
    'ERP185409', 'ERP185509', 'ERP191728', 'ERP000263', 'ERP022096', 'ERP022298']))


def valid_manifest(study, route, cohort):
    ae = 'E-GEOD-100' if study == 'GSE100' else study
    options = {'enrichment': 'off', 'options': {'expand_studies': False}}
    source = {'kind': 'file', 'file': 'idf'}
    services = {}
    intent = 'deposited-off'
    if route in ['native', 'ena', 'standard', 'linked']:
        source = {'kind': 'accession', 'value': ACCESSIONS[study]}
        options['in_type'] = 'ena' if route == 'ena' else 'sra'
        intent = 'parser-native-off'
    if route == 'standard':
        options['enrichment'] = 'standard'
        intent = 'production-standard'
    if route == 'linked':
        options.pop('enrichment')
        options['input_options'] = {'include_peer': False, 'enrich_from_geo_ae': True}
        services = {'enricher': {'kind': 'no-publication'}, 'ae2json': {'kind': 'local-ae'}}
        intent = 'linked-only'
    if study == 'GSE100' and route == 'direct': source = {'kind': 'file', 'file': 'geo'}
    return {'cohort': cohort, 'files': {'idf': {'path': ae + '.idf.txt'},
            'sdrf': {'path': ae + '.sdrf.txt'}, 'geo': {'path': 'tests/fixtures/studies/GSE100/inputs/geo.xml'}},
            'scenarios': [{'id': f'{study}/{route}', 'formats': FORMATS,
             'classification': 'complete-captured-coverage', 'route_intent': intent,
             'source': source, 'services': services, 'options': options,
             'expected': {fmt: {} for fmt in FORMATS}}]}


def save(path, value):
    path.write_text(json.dumps(value))
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.fixture
def panel(tmp_path, monkeypatch):
    root = tmp_path / 'panel'
    root.mkdir()
    rows = []
    for study, route in ROSTER:
        name = f'{study}.{route}.json'
        cohort = f'{study}-{route}-test'
        manifest = valid_manifest(study, route, cohort)
        rows.append({'study': study, 'route': route, 'manifest_id': cohort,
                     'manifest': name, 'manifest_sha256': save(root / name, manifest)})
    index = {'schema_version': '1.0', 'panel': 'workflow-provenance-18-studies',
             'expected_studies': 18, 'expected_routes': 74, 'expected_facade_calls': 222,
             'routes': rows}
    path = root / 'index.json'
    save(path, index)
    calls = []
    def fake_audit(manifest, *, evidence_root=None, artifacts_dir=None):
        calls.append((Path(manifest), evidence_root, artifacts_dir))
        if artifacts_dir is not None:
            Path(artifacts_dir).mkdir(parents=True)
            (Path(artifacts_dir) / 'payload.json').write_text('{}')
        return {'classification': 'complete-captured-coverage', 'errors': [],
                'scenarios': [{'calls': [{'format': fmt, 'status': 'passed'} for fmt in FORMATS]}]}
    monkeypatch.setattr(runner, 'audit', fake_audit)
    return path, index, calls


def execute(panel, **kwargs):
    path, index, calls = panel
    save(path, index)
    return runner.audit_panel(path, evidence_root=path.parent, **kwargs)


def test_complete_explicit_roster_executes_74_routes_222_calls(panel, tmp_path):
    report = execute(panel, artifacts_dir=tmp_path / 'artifacts')
    assert report['classification'] == 'complete-captured-coverage'
    assert len(panel[2]) == 74
    assert report['passed_facade_calls'] == 222
    assert len({str(c[2]) for c in panel[2]}) == 74
    assert all(c[0].is_relative_to(panel[0].parent) for c in panel[2])


def test_omitted_route_and_tampered_declared_counts_cannot_pass(panel):
    panel[1]['routes'].pop()
    panel[1]['expected_routes'] = 73
    panel[1]['expected_facade_calls'] = 219
    report = execute(panel)
    assert report['classification'] == 'harness-error'
    assert len(report['routes']) == 74
    assert len(panel[2]) == 73  # Available routes still execute.
    assert any('roster' in e or 'counts' in e for e in report['errors'])


def test_same_counts_wrong_roster_fails(panel):
    panel[1]['routes'][0]['study'] = 'E-MTAB-UNAPPROVED'
    report = execute(panel)
    assert report['classification'] == 'harness-error'
    assert len(panel[2]) == 73


@pytest.mark.parametrize('key', ['manifest_id', 'manifest'])
def test_duplicate_manifest_identity_or_path_is_rejected(panel, key):
    panel[1]['routes'][1][key] = panel[1]['routes'][0][key]
    report = execute(panel)
    assert report['classification'] == 'harness-error'
    assert len(panel[2]) < 74


def test_duplicate_study_route_is_not_counted_twice(panel):
    panel[1]['routes'][1].update(study=panel[1]['routes'][0]['study'], route=panel[1]['routes'][0]['route'])
    report = execute(panel)
    assert report['classification'] == 'harness-error'
    assert len(panel[2]) < 74


@pytest.mark.parametrize('change', ['hash', 'missing', 'formats', 'scenario_id', 'package_count', 'id'])
def test_one_bad_manifest_does_not_skip_other_available_routes(panel, change):
    row = panel[1]['routes'][0]
    p = panel[0].parent / row['manifest']
    doc = json.loads(p.read_text())
    if change == 'hash': row['manifest_sha256'] = '0' * 64
    elif change == 'missing': p.unlink()
    else:
        if change == 'formats': doc['scenarios'][0]['formats'] = ['json']
        elif change == 'scenario_id': doc['scenarios'][0]['id'] = 'wrong/reference'
        elif change == 'package_count': doc['scenarios'].append(doc['scenarios'][0])
        elif change == 'id': doc['cohort'] = 'wrong-cohort'
        row['manifest_sha256'] = save(p, doc)
    report = execute(panel)
    assert report['classification'] == 'harness-error'
    assert len(panel[2]) == 73
    assert sum(r['classification'] == 'harness-error' for r in report['routes']) == 1


def test_manifest_path_cannot_escape_index_root(panel, tmp_path):
    panel[1]['routes'][0]['manifest'] = '../outside.json'
    report = execute(panel)
    assert report['classification'] == 'harness-error'
    assert len(panel[2]) == 73


def test_existing_artifact_directory_is_not_reused(panel, tmp_path):
    out = tmp_path / 'artifacts'
    out.mkdir()
    (out / 'sentinel').write_text('preserve')
    report = execute(panel, artifacts_dir=out)
    assert report['classification'] == 'harness-error'
    assert not panel[2]
    assert (out / 'sentinel').read_text() == 'preserve'


def test_json_only_supplemental_cohort_cannot_fill_missing_route(panel):
    removed = panel[1]['routes'].pop()
    panel[1]['fresh_cohorts'] = [removed | {'formats': ['json'], 'classification': 'complete-captured-coverage'}]
    report = execute(panel)
    assert report['classification'] == 'harness-error'
    assert report['passed_facade_calls'] == 219


def test_actual_returned_call_coverage_is_checked(panel, monkeypatch):
    def omitted(*args, **kwargs):
        return {'classification': 'complete-captured-coverage', 'scenarios': [{'calls': [{'format': 'json', 'status': 'passed'}]}]}
    monkeypatch.setattr(runner, 'audit', omitted)
    report = execute(panel)
    assert report['classification'] == 'harness-error'
    assert report['passed_facade_calls'] == 74


def test_explicit_external_root_required(panel):
    report = runner.audit_panel(panel[0])
    assert report['classification'] == 'harness-error'
    assert not panel[2]


def test_unmet_audit_route_retains_failure_and_runs_remaining_routes(panel, monkeypatch):
    original = runner.audit
    failed_manifest = panel[1]['routes'][0]['manifest']
    def missing(manifest, **kwargs):
        result = original(manifest, **kwargs)
        if Path(manifest).name == failed_manifest:
            return {'classification': 'harness-error', 'errors': ['unmet captured request'], 'scenarios': []}
        return result
    monkeypatch.setattr(runner, 'audit', missing)
    report = execute(panel)
    assert report['classification'] == 'harness-error'
    assert len(panel[2]) == 74
    assert report['passed_facade_calls'] == 219
    failed = [r for r in report['routes'] if r['classification'] == 'harness-error']
    assert len(failed) == 1
    assert failed[0]['audit']['errors'] == ['unmet captured request']


def test_changed_manifest_during_execution_is_rejected(panel, monkeypatch):
    original = runner.audit
    def changing(manifest, **kwargs):
        result = original(manifest, **kwargs)
        p = Path(manifest)
        p.write_text(p.read_text() + ' ')
        return result
    monkeypatch.setattr(runner, 'audit', changing)
    report = execute(panel)
    assert report['classification'] == 'harness-error'
    assert all('changed during route audit' in r['errors'][0] for r in report['routes'])


def test_panel_cli_dispatch_and_failure_exit(panel, monkeypatch, capsys):
    observed = []
    def failed(path, **kwargs):
        observed.append((path, kwargs))
        return {'classification': 'harness-error', 'errors': ['missing route']}
    monkeypatch.setattr(runner, 'audit_panel', failed)
    assert runner.main(['--panel', str(panel[0]), '--evidence-root', str(panel[0].parent)]) == 2
    assert observed[0][0] == panel[0]
    assert 'missing route' in capsys.readouterr().out


@pytest.mark.fake_process
def test_optimized_python_cannot_accept_panel(panel):
    import subprocess
    import sys
    result = subprocess.run([sys.executable, '-O', str(Path(runner.__file__)), '--panel',
                             str(panel[0]), '--evidence-root', str(panel[0].parent)],
                            capture_output=True, text=True, check=False)
    assert result.returncode == 2
    report = json.loads(result.stdout)
    assert report['classification'] == 'harness-error'
    assert 'optimized Python' in report['errors'][0]


def test_real_reduced_route_executes_while_73_unavailable_routes_are_reported(tmp_path):
    evidence = Path(__file__).parent / 'fixtures' / 'recorded_provider'
    manifest = json.loads((evidence / 'success.manifest.json').read_text())
    manifest['cohort'] = 'panel-real-route'
    manifest['scenarios'][0]['id'] = 'E-MTAB-16847/reference'
    name = 'real-route.json'
    row = {'study': 'E-MTAB-16847', 'route': 'reference', 'manifest_id': manifest['cohort'],
           'manifest': name, 'manifest_sha256': save(tmp_path / name, manifest)}
    index = {'schema_version': '1.0', 'panel': 'workflow-provenance-18-studies',
             'expected_studies': 18, 'expected_routes': 74, 'expected_facade_calls': 222,
             'routes': [row]}
    path = tmp_path / 'index.json'
    save(path, index)
    report = runner.audit_panel(path, evidence_root=evidence)
    assert report['classification'] == 'harness-error'
    assert len(report['routes']) == 74
    assert report['passed_facade_calls'] == 3
    passed = [r for r in report['routes'] if r['classification'] == 'complete-captured-coverage']
    assert len(passed) == 1
    assert passed[0]['audit']['transport']['unexpected_misses'] == 0


@pytest.mark.parametrize('change', ['loaded_export', 'accession', 'provider', 'disable_peer',
                                   'disable_linked', 'expand', 'archive_idf', 'array_reference', 'geo_direct'])
def test_equal_roster_cannot_substitute_wrong_scientific_route(panel, change):
    study, route = ('E-MTAB-14416', 'standard')
    if change == 'archive_idf': route = 'reference'
    if change == 'array_reference': study, route = 'GSE100', 'reference'
    if change == 'geo_direct': study, route = 'GSE100', 'direct'
    row = next(r for r in panel[1]['routes'] if (r['study'], r['route']) == (study, route))
    path = panel[0].parent / row['manifest']
    doc = json.loads(path.read_text())
    scenario = doc['scenarios'][0]
    if change == 'loaded_export':
        scenario['route_intent'] = 'loaded-export-off'
        scenario['options']['enrichment'] = 'off'
        scenario['source'] = {'kind': 'json', 'file': 'golden'}
    elif change == 'accession': scenario['source']['value'] = 'ERP999999'
    elif change == 'provider': scenario['options']['in_type'] = 'ena'
    elif change in ['disable_peer', 'disable_linked']:
        scenario['options']['input_options'] = {'include_peer' if change == 'disable_peer' else 'enrich_from_geo_ae': False}
    elif change == 'expand': scenario['options']['options']['expand_studies'] = True
    elif change in ['archive_idf', 'array_reference']: doc['files']['idf']['path'] = 'E-MTAB-OTHER.idf.txt'
    elif change == 'geo_direct': scenario['source'] = {'kind': 'file', 'file': 'idf'}
    row['manifest_sha256'] = save(path, doc)
    report = execute(panel)
    assert report['classification'] == 'harness-error'
    assert len(panel[2]) == 73
    assert sum(r['classification'] == 'harness-error' for r in report['routes']) == 1
