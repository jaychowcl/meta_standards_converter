# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Audit boundaries: independent evidence, explicit revision and strict replay."""
import copy
import hashlib
from types import SimpleNamespace

import pytest

from scripts.audit_additional_conversions import COHORT, clean_revision, consumed_manifest, initialize, read, validate_cohort
from tests.support.authorship_capture import write_json
from tests.support.authorship_sources import documents, native_coverage


def test_explicit_cohort_has_seven_inputs_and_five_studies():
    cohort = validate_cohort(read(COHORT))
    assert len(cohort['inputs']) == 7
    assert len({s['reference'] for s in cohort['inputs']}) == 5
    assert sum(s['evidence'] == 'live' for s in cohort['inputs']) == 5
    assert sum(s['evidence'] == 'parser-merge-replay' for s in cohort['inputs']) == 2


@pytest.mark.parametrize('mutation', ['duplicate', 'provider', 'evidence', 'schema'])
def test_invalid_cohort_rejected(mutation):
    cohort = copy.deepcopy(read(COHORT))
    if mutation == 'duplicate': cohort['inputs'].append(cohort['inputs'][0])
    elif mutation == 'schema': cohort['schema_version'] = 'unknown'
    else: cohort['inputs'][0][mutation] = 'unknown'
    with pytest.raises(ValueError): validate_cohort(cohort)


@pytest.mark.parametrize('problem', ['revision', 'dirty', 'production', None])
def test_tested_revision_and_production_guard(monkeypatch, problem):
    def git(args, **kwargs):
        if args[1] == 'status': return ' M file' if problem == 'dirty' else ''
        if args[1] == 'diff': return 'changed source' if problem == 'production' else ''
        return 'wrong' if '--verify' in args and problem == 'revision' else 'actual'
    monkeypatch.setattr('scripts.audit_additional_conversions.subprocess.check_output', git)
    if problem:
        with pytest.raises(ValueError): clean_revision('actual', 'production')
    else: assert clean_revision('actual', 'production') == 'actual'


def test_initialization_cannot_overwrite_existing_evidence(tmp_path):
    (tmp_path / 'evidence').write_text('original')
    with pytest.raises(FileExistsError): initialize(tmp_path, {}, 'unused')
    assert (tmp_path / 'evidence').read_text() == 'original'


def test_request_selection_preserves_attempt_order_and_removes_unused():
    identity = {'method': 'GET', 'url': 'https://example.org/a', 'body_sha256': '0' * 64}
    manifest = {'requests': [{'identity': identity, 'attempts': [
        {'exception': 'timeout', 'repeat': 1}, {'status': 200, 'response': 'body', 'repeat': 3}]},
        {'identity': dict(identity, url='https://example.org/unused'), 'attempts': [{'status': 200, 'response': 'unused', 'repeat': 1}]}]}
    replay = SimpleNamespace(manifest=manifest, events=[{'identity': identity}] * 2, assert_no_errors=lambda: None)
    selected = consumed_manifest(replay)
    assert selected['requests'] == [{'identity': identity, 'attempts': [
        {'exception': 'timeout', 'repeat': 1}, {'status': 200, 'response': 'body', 'repeat': 1}]}]
    assert manifest['requests'][0]['attempts'][1]['repeat'] == 3


def test_selected_evidence_excludes_unused_reference_and_preserves_native_snapshot(tmp_path):
    # Write JSON through the helper to avoid relying on string escape spelling.
    write_json(tmp_path / 'snapshot.json', {'xml': ['<ROOT><SAMPLE accession="ERS1"/><RUN accession="ERR1"/></ROOT>']})
    digest = hashlib.sha256((tmp_path / 'snapshot.json').read_bytes()).hexdigest()
    write_json(tmp_path / 'source-snapshots.json', [{'body': 'snapshot.json', 'sha256': digest, 'url': 'frozen-parser-record:ERP1'}])
    write_json(tmp_path / 'http/curators.manifest.json', {'requests': [{'invalid': 'unconsumed reference'}]})
    write_json(tmp_path / 'http/curators-json-selected.manifest.json', {'requests': [], 'files': {}})
    docs = documents(tmp_path, 'curators')
    assert len(docs) == 1 and docs[0]['reference_assisted'] is False
    coverage = native_coverage(docs, {'sample': [{'sra_run': [{'run': 'ERR1'}], 'accession': [{'value': 'ERS1'}]}]})
    assert coverage['source_runs'] == ['ERR1']
    assert coverage['source_samples'] == ['ERS1']
    (tmp_path / 'snapshot.json').write_text('{}')
    with pytest.raises(AssertionError, match='corrupt native snapshot'): documents(tmp_path, 'curators')


@pytest.mark.parametrize('swap', [False, True])
def test_geo_channel_scalar_labels_keep_channel_ownership(swap):
    from tests.test_authorship_sources import document
    from tests.support.authorship_sources import geo_checks
    from tests.support.authorship_fidelity import EQUIVALENT, LOSS
    docs = [document('<MINiML><Sample iid="GSM1"><Channel><Label>Cy3</Label></Channel>'
                     '<Channel><Label>Cy5</Label></Channel></Sample></MINiML>')]
    labels = ['Cy5', 'Cy3'] if swap else ['Cy3', 'Cy5']
    package = {'series': {'iid': 'GSE1'}, 'sample': [{'iid': 'GSM1', 'channel': [{'label': v} for v in labels]}]}
    checks = geo_checks(docs, package)
    assert len(checks) == 2
    assert {r['classification'] for r in checks} == {LOSS if swap else EQUIVALENT}


@pytest.mark.parametrize('original,canonical', [
    ('Characteristics [site]', 'Characteristics[site]'),
    ('Characteristics [cell line]', 'Characteristics[cell line]'),
    ('Comment [Sample_characteristics]', 'Comment[Sample_characteristics]'),
    ('Comment [ArrayExpress FTP file]', 'Comment[ArrayExpress FTP file]'),
    ('Comment [ercc_mix]', 'Comment[ercc_mix]'),
    ('Factor Value [SITE]', 'Factor Value[SITE]'),
    ('FactorValue [TISSUE ORIGIN]', 'Factor Value[TISSUE ORIGIN]'),
    ('FactorValue [treatment]', 'Factor Value[treatment]'),
])
def test_observed_deposited_header_spellings_preserve_values_and_qualifiers(original, canonical):
    from tests.support.deposited_fidelity import assert_deposited_sdrf
    source = [['Source Name', original, 'Term Source REF'], ['sample', 'original', 'EFO']]
    target = [['Source Name', canonical, 'Term Source REF'], ['sample', 'original', 'EFO']]
    assert_deposited_sdrf(source, target, accession='audit')
    target[1][1] = 'different'
    with pytest.raises(AssertionError): assert_deposited_sdrf(source, target, accession='audit')
    target[1][1] = 'original'; target[1][2] = 'wrong vocabulary'
    with pytest.raises(AssertionError): assert_deposited_sdrf(source, target, accession='audit')


def test_geo_contributor_scope_uses_source_relationships_not_generated_presence():
    from tests.test_authorship_sources import document
    from tests.support.authorship_sources import geo_checks
    from tests.support.authorship_fidelity import LOSS
    docs = [document('<MINiML><Series iid="GSE1"><Sample-Ref ref="GSM1"/></Series>'
        '<Sample iid="GSM1"><Contact-Ref ref="owned"/></Sample>'
        '<Contributor iid="owned"><Email>owned@example.org</Email></Contributor>'
        '<Contributor iid="orphan"><Email>platform@example.org</Email></Contributor></MINiML>')]
    checks = geo_checks(docs, {'series': {'iid': 'GSE1'}, 'sample': [{'iid': 'GSM1'}], 'contributor': []})
    contacts = [c for c in checks if c['scope'].startswith('contact/')]
    assert len(contacts) == 1
    assert contacts[0]['scope'] == 'contact/owned'
    assert contacts[0]['classification'] == LOSS


@pytest.mark.parametrize('value,removed', [('  ', True), ('EFO', False), ('0', False), (' EFO ', False)])
def test_deposited_empty_qualifier_projection_never_removes_populated_values(value, removed):
    from scripts.audit_additional_review import blank_projection
    source = [['Source Name', 'Protocol REF', 'Term Source REF'], ['s', 'p', value]]
    projected, changes = blank_projection(source, 'E-GEOD-31356')
    assert (len(projected[0]) == 2) is removed
    assert bool(changes) is removed
    if not removed: assert projected[1][-1] == value
    assert source[1][-1] == value


def test_blank_characteristic_projection_keeps_populated_owned_comments():
    from scripts.audit_additional_review import blank_projection
    source = [['Source Name', 'Characteristics [site]', 'Comment [Sample_characteristics]'], ['s', '  ', 'supplied text']]
    projected, changes = blank_projection(source, 'E-GEOD-13398')
    assert projected[1] == ['s', '', 'supplied text']
    assert len(changes) == 1
    assert blank_projection(source, 'unrelated')[0] == source
