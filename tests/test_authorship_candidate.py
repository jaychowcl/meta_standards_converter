# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Candidate replay cannot overwrite history, mislabel revisions or reuse stale JSON."""
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from scripts.audit_authorship import BASELINE, replay_case
from scripts.audit_authorship_candidate import candidate_preflight, verify_seal, seal_result, legacy_idf_comparison
from tests.support.authorship_fidelity import LOSS
from scripts.audit_authorship_report import reference_assistance
from tests.support.authorship_capture import write_json


def sealed(tmp_path):
    root = tmp_path / 'sealed'
    root.mkdir()
    write_json(root / 'environment.json', {'baseline': BASELINE})
    seal_result(root)
    with pytest.raises(FileExistsError):
        seal_result(root)
    return root


@pytest.mark.parametrize('change', ['corruption', 'extra', 'escape', 'absolute'])
def test_seal_rejects_corruption_or_unsealed_bytes(tmp_path, change):
    root = sealed(tmp_path)
    assert verify_seal(root)['files'] == 1
    if change == 'corruption': (root / 'environment.json').write_text('{}')
    elif change == 'extra': (root / 'extra.json').write_text('{}')
    else:
        manifest = root / 'SHA256SUMS'
        manifest.write_text(manifest.read_text().replace('environment.json', '../environment.json' if change == 'escape' else '/environment.json'))
    with pytest.raises(ValueError): verify_seal(root)


@pytest.mark.parametrize('output', ['same', 'inside', 'parent', 'existing'])
def test_output_guard_precedes_any_write(tmp_path, output):
    source = sealed(tmp_path)
    destination = {'same': source, 'inside': source / 'new', 'parent': tmp_path, 'existing': tmp_path / 'exists'}[output]
    if output == 'existing': destination.mkdir()
    before = (source / 'SHA256SUMS').read_bytes()
    with pytest.raises((ValueError, FileExistsError)):
        candidate_preflight(source, destination, 'unused')
    assert (source / 'SHA256SUMS').read_bytes() == before


@pytest.mark.parametrize('problem', ['revision', 'dirty', None])
def test_candidate_requires_the_actual_clean_revision(tmp_path, monkeypatch, problem):
    source = sealed(tmp_path)
    def git(args, **kwargs):
        if args[1] == 'status': return ' M src/a.py\n' if problem == 'dirty' else ''
        if '--verify' in args: return ('other' if problem == 'revision' else 'candidate') + '\n'
        return 'candidate\n'
    monkeypatch.setattr('scripts.audit_authorship_candidate.subprocess.check_output', git)
    if problem:
        with pytest.raises(ValueError): candidate_preflight(source, tmp_path / 'new', 'candidate')
    else:
        _, output, evidence = candidate_preflight(source, tmp_path / 'new', 'candidate')
        assert evidence['tested_revision'] == 'candidate'
        assert evidence['acquisition_revision'] == BASELINE
        assert not output.exists()


def test_fresh_saved_route_uses_new_source_json_and_retains_older_capture(tmp_path, monkeypatch):
    dest = tmp_path / 'geo-GSE1'
    write_json(dest / 'acquisition.json', {'accession': 'GSE1', 'provider': 'geo', 'phases': {'off': {}}})
    write_json(dest / 'http/off.manifest.json', {'schema_version': '1.0', 'cohort': 'test', 'source_revision': BASELINE,
        'captured_at': 'frozen', 'evidence_kind': 'recorded-http', 'files': {}, 'requests': []})
    old = dest / 'off/baseline-source-json/payload-0.json'
    write_json(old, [{'supplied': 'older'}])
    exported = []
    class Result:
        def __init__(self, package): self.items = [SimpleNamespace(payload=[package])]
        def to_dict(self): return {'status': 'success', 'items': []}
    def convert(self, value, *, out_type, **kwargs):
        if isinstance(value, Path):
            package = json.loads(value.read_text())
            exported.append(package)
        else: package = {'supplied': 'fresh'}
        return Result(package)
    monkeypatch.setattr('meta_standards_converter.converters.Converter.convert', convert)
    replay_case(dest, fresh_json=True)
    assert exported == [{'supplied': 'fresh'}]
    assert json.loads(old.read_text()) == [{'supplied': 'older'}]
    assert json.loads((dest / 'off/saved-0.json').read_text()) == {'supplied': 'fresh'}


def test_reference_assistance_uses_actual_replay_not_original_capture(tmp_path):
    write_json(tmp_path / 'http/trace.json', {'events': [
        {'phase': 'off', 'identity': {'url': 'https://www.ebi.ac.uk/biostudies/E-GEOD-1'}, 'attempt': {'response': 'hash'}}]})
    write_json(tmp_path / 'off/trace-json.json', {'events': [], 'errors': []})
    evidence = reference_assistance(tmp_path, 'off', None)
    assert evidence['reference_assisted'] is False
    assert evidence['trace'] == 'off/trace-json.json'


@pytest.mark.parametrize('mutation', [None, 'role', 'phone', 'name', 'address', 'other-contact'])
def test_legacy_review_only_accepts_the_exact_source_supported_contact_completion(mutation):
    full = {'Person First Name': 'Gordon', 'Person Mid Initials': 'K', 'Person Last Name': 'Smyth',
        'Person Email': 'smyth@wehi.edu.au', 'Person Phone': 'supplied phone', 'Person Roles': None,
        'Person Affiliation': 'Smyth, Bioinformatics, Walter and Eliza Hall Institute of Medical Research',
        'Person Address': 'Smyth, Bioinformatics, Walter and Eliza Hall Institute of Medical Research, 1G Royal Pde, Parkville, Vic, 3052, Australia'}
    short = {**full, 'Person Roles': 'submitter',
        'Person Affiliation': 'Bioinformatics, Walter and Eliza Hall Institute of Medical Research',
        'Person Address': 'Bioinformatics, Walter and Eliza Hall Institute of Medical Research, 1G Royal Pde, Parkville, Vic, Australia'}
    other = {'Person First Name': 'Another', 'Person Last Name': 'Contact'}
    new = [{**full, 'Person Roles': 'submitter'}, other.copy()]
    if mutation == 'role': new[0]['Person Roles'] = None
    elif mutation == 'phone': new[0]['Person Phone'] = 'wrong phone'
    elif mutation == 'name': new[0]['Person First Name'] = 'Another'
    elif mutation == 'address': new[0]['Person Address'] = full['Person Address'].replace('3052', '9999')
    elif mutation == 'other-contact': new.pop()
    def rows(people): return [[field, *(p.get(field) for p in people)] for field in full]
    comparisons = legacy_idf_comparison(rows([full, short, other]), rows(new), 'sra-SRP045534', 'standard')
    losses = [r for r in comparisons if r['classification'] == LOSS]
    assert bool(losses) == bool(mutation)
    if not mutation:
        assert any(r.get('generated_identity') for r in comparisons)
        assert any(r['classification'] == LOSS for r in legacy_idf_comparison(rows([full, short]), rows(new), 'sra-SRP002056', 'standard'))
