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
