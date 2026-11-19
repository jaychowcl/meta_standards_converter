# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Candidate replay retains historical revision guards and sealed inputs."""
import json
from pathlib import Path
import pytest

from scripts.audit_additional_conversions import candidate_preflight, COHORT, clean_revision
from scripts.audit_authorship_candidate import seal_result


@pytest.mark.parametrize('problem', ['', 'dirty', 'revision', 'baseline', 'origin', 'overlap', 'exists'])
def test_candidate_requires_clean_revision_origin_and_separate_sealed_inputs(tmp_path, monkeypatch, problem):
    source = tmp_path / 'source'
    source.mkdir()
    cohort = json.loads(COHORT.read_text())
    (source / 'cohort.json').write_text(json.dumps(cohort))
    (source / 'environment.json').write_text(json.dumps({
        'tested_revision': 'historical', 'production_revision': 'wrong' if problem == 'baseline' else cohort['production_revision']}))
    seal_result(source)
    def git(args, **kwargs):
        if 'status' in args:
            return ' M changed' if problem == 'dirty' else ''
        return 'different' if '--verify' in args and problem == 'revision' else 'candidate'
    monkeypatch.setattr('scripts.audit_additional_conversions.subprocess.check_output', git)
    if problem == 'origin':
        monkeypatch.setattr('meta_standards_converter.__file__', '/tmp/unrelated-msc/__init__.py')
    target = source / 'nested' if problem == 'overlap' else tmp_path / 'candidate'
    if problem == 'exists':
        target.mkdir()
    if problem:
        with pytest.raises((ValueError, FileExistsError)):
            candidate_preflight(source, target, 'candidate')
    else:
        _, _, actual, evidence = candidate_preflight(source, target, 'candidate')
        assert actual == cohort and evidence['tested_revision'] == 'candidate'
        assert evidence['source_seal_before']['files'] == 2 and not target.exists()


def test_historical_guard_still_rejects_modified_production(monkeypatch):
    def git(args, **kwargs):
        if 'diff' in args: return 'production changed'
        if 'status' in args: return ''
        return 'baseline'
    monkeypatch.setattr('scripts.audit_additional_conversions.subprocess.check_output', git)
    with pytest.raises(ValueError, match='Production code differs'):
        clean_revision('baseline', 'historical-production')
