# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Read a sealed authorship audit and evaluate a clean, committed candidate offline."""
from collections import Counter
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

from scripts.audit_authorship import BASELINE, CASES, REPO, deposited, frozen_unicode, policy, replay_case, save_result
from scripts.audit_authorship_report import build, read
from tests.support.authorship_capture import utc, write_json
from tests.support.authorship_fidelity import EQUIVALENT, LOSS, idf_comparison
from tests.support.deposited_fidelity import assert_deposited_sdrf, without_added_protocol_placeholders
from tests.support.recorded_replay import ReplayTransport


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def seal_result(root):
    """Seal only a newly completed result directory; never refresh an existing seal."""
    manifest = Path(root) / 'SHA256SUMS'
    files = [p for p in sorted(Path(root).rglob('*')) if p.is_file()]
    with manifest.open('x') as stream:
        for path in files:
            stream.write(sha256(path) + '  ' + str(path.relative_to(root)) + '\n')


def verify_seal(root):
    """Check every declared byte and reject paths escaping the sealed root."""
    root = Path(root).resolve()
    manifest = root / 'SHA256SUMS'
    seen = set()
    for line in manifest.read_text().splitlines():
        digest, name = line.split('  ', 1)
        relative = Path(name)
        path = root / relative
        if relative.is_absolute() or '..' in relative.parts or not path.resolve().is_relative_to(root):
            raise ValueError('Unsafe sealed path: ' + name)
        if name in seen or len(digest) != 64 or sha256(path) != digest:
            raise ValueError('Sealed evidence hash mismatch: ' + name)
        seen.add(name)
    actual = {str(p.relative_to(root)) for p in root.rglob('*') if p.is_file()}
    if actual != seen | {'SHA256SUMS'}:
        raise ValueError('Sealed file inventory differs from SHA256SUMS')
    return {'files': len(seen), 'sha256_manifest': sha256(manifest)}


def _pytest_version():
    from importlib.metadata import PackageNotFoundError, version
    try:
        return version('pytest')
    except PackageNotFoundError:
        return 'unavailable'


def candidate_preflight(source, output, expected_revision):
    source, output = Path(source).resolve(), Path(output).resolve()
    if source == output or source.is_relative_to(output) or output.is_relative_to(source):
        raise ValueError('Candidate output must be separate from the sealed input')
    if output.exists():
        raise FileExistsError('Candidate output already exists; use a fresh result directory')
    revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=REPO, text=True).strip()
    expected = subprocess.check_output(['git', 'rev-parse', '--verify', expected_revision + '^{commit}'], cwd=REPO, text=True).strip()
    if revision != expected:
        raise ValueError('Candidate revision does not match the actual tested HEAD')
    dirty = subprocess.check_output(['git', 'status', '--porcelain', '--untracked-files=all'], cwd=REPO, text=True)
    if dirty:
        raise ValueError('Commit the candidate source, tests and runner before replay')
    import meta_standards_converter
    origin = Path(meta_standards_converter.__file__).resolve()
    if not origin.is_relative_to((REPO / 'src').resolve()):
        raise ValueError('Imported MSC code is outside the candidate checkout')
    seal = verify_seal(source)
    if read(source / 'environment.json')['baseline'] != BASELINE:
        raise ValueError('Expected the original 942ef67 acquisition baseline')
    return source, output, {'tested_revision': revision, 'acquisition_revision': BASELINE,
        'source_bundle': str(source), 'source_seal_before': seal, 'started_at': utc(),
        'python': sys.executable, 'python_version': sys.version, 'pytest_version': _pytest_version(), 'msc_origin': str(origin), 'live_requests': 0,
        'settings': 'Frozen request manifests; candidate JSON generated before JSON-to-MAGE-TAB; legacy JSON exported separately.'}


def legacy_idf_comparison(previous, current, case, mode):
    """Keep every old field, with the one independently verified contact join explicit."""
    rows = idf_comparison(previous, current, documented=True)
    if (case, mode) != ('sra-SRP045534', 'standard'):
        return rows
    def matched(document):
        fields = [r for r in document if r[0].startswith('Person ')]
        people = [dict(zip([r[0] for r in fields], values)) for values in zip(*(r[1:] for r in fields))]
        return [p for p in people if p.get('Person Email') == 'smyth@wehi.edu.au']
    old, new = matched(previous), matched(current)
    full = 'Smyth, Bioinformatics, Walter and Eliza Hall Institute of Medical Research, 1G Royal Pde, Parkville, Vic, 3052, Australia'
    short = 'Bioinformatics, Walter and Eliza Hall Institute of Medical Research, 1G Royal Pde, Parkville, Vic, Australia'
    if len(old) != 2 or len(new) != 1 or [p.get('Person Address') for p in old] != [full, short]:
        return rows
    if old[0].get('Person Affiliation') != 'Smyth, Bioinformatics, Walter and Eliza Hall Institute of Medical Research' \
            or old[1].get('Person Affiliation') != 'Walter and Eliza Hall Institute of Medical Research':
        return rows
    for field in old[0].keys() | old[1].keys():
        if field in {'Person Address', 'Person Affiliation'}:
            continue
        left, right = old[0].get(field), old[1].get(field)
        if left and right and ' '.join(str(left).split()).casefold() != ' '.join(str(right).split()).casefold():
            return rows
    expected = {k: v or old[1].get(k) for k, v in old[0].items()}
    if new[0] != expected:
        return rows  # Never excuse missing roles, identifiers or other contact facts.
    identity = ['Gordon', 'K', 'Smyth', 'smyth@wehi.edu.au']
    for row in rows:
        scope, bound, occurrence = row['identity']
        if scope != 'contact' or json.loads(bound) != identity:
            continue
        field = row['field']
        if occurrence == 2 and row['generated'] is None and row['source'] == old[1].get(field):
            row.update(classification=EQUIVALENT, generated=expected.get(field), generated_identity=[scope, bound, 1],
                reason='Verified GEO/AE shared-email contact consolidation: exact fuller source address/affiliation and every supplied remaining field retained on the surviving contact.')
        elif occurrence == 1 and row['source'] is None and row['generated'] == old[1].get(field):
            row.update(classification=EQUIVALENT,
                reason='Supplied field completed from the second email-matched contact; original contributors retained in JSON.')
    return rows


def legacy_exports(source, root):
    """Compare older JSON exports to their frozen outputs, retaining every difference.

    The old JSON is the input, not a proxy for information already lost before it
    was saved. Raw-source checks are performed separately on the fresh route.
    """
    from meta_standards_converter.converters import Converter
    results = []
    for provider, accession, reference in CASES:
        name = provider + '-' + accession
        for mode in ('off', 'curators', 'standard'):
            original = source / 'cases' / name / mode
            dest = root / 'cases' / name / mode
            saved = dest / 'legacy-saved-0.json'
            shutil.copyfile(original / 'saved-0.json', saved)
            before = sha256(saved)
            isolation = ReplayTransport(root / 'empty.manifest.json')
            with isolation.offline():
                converted = Converter().convert(saved, out_type='magetab', preparation=policy('off'))
            isolation.assert_complete()
            outcome = save_result(converted, dest / 'legacy-saved-0-magetab')
            assert sha256(saved) == before == sha256(original / 'saved-0.json')
            previous = read(original / 'saved-0-magetab/payload-0.json')[0]
            current = read(dest / 'legacy-saved-0-magetab/payload-0.json')[0]
            differences = legacy_idf_comparison(previous, current, name, mode)
            old_sdrf = next(r[1] for r in previous if r[0] == 'SDRF File')
            try:
                new_sdrf = without_added_protocol_placeholders(previous, current, allow_mixed_rows=True)
                from tests.support.integrity_legacy import reviewed_legacy_tables
                definitions = dict(zip(next(r[1:] for r in previous if r[0] == 'Protocol Name'),
                                       next(r[1:] for r in previous if r[0] == 'Protocol Description')))
                old_sdrf, new_sdrf, reviewed = reviewed_legacy_tables(read(saved), old_sdrf, new_sdrf, definitions)
                sdrf = {'verdict': 'preserved', **assert_deposited_sdrf(old_sdrf, new_sdrf, accession=reference or accession)}
                sdrf['reviewed_transformations'] = reviewed
            except (AssertionError, KeyError, StopIteration) as exc:
                sdrf = {'verdict': 'discrepancy', 'reason': str(exc)}
            record = {'case': name, 'mode': mode, 'input_sha256': before, 'input_unchanged': True,
                'result': outcome, 'idf': differences, 'sdrf': sdrf,
                'scope': 'All fields already represented in the frozen JSON export; pre-existing omissions cannot be recovered without source evidence.'}
            write_json(dest / 'legacy-comparison.json', record)
            results.append(record)
    write_json(root / 'legacy-results.json', results)
    return results


def reference_corrections(source, root):
    """Create labelled derivatives from independent frozen PubMed evidence only."""
    from tests.support.authorship_sources import citation_facts, documents
    dest = root / 'reference-corrections'
    dest.mkdir()
    records = []
    for accession, case, pmid, expected, action in (
        ('E-GEOD-100', 'geo-GSE100', '12524544', '10.1038/ng1080', 'corrected-derivative'),
        ('E-MTAB-6486', 'ae-E-MTAB-6486', '29925636', '10.1126/scitranslmed.aao2565', 'curator-review-required'),
    ):
        evidence = citation_facts(documents(source / 'cases' / case, 'standard'))[pmid]
        supported = [p for p in evidence if p.get('doi') == expected]
        if not supported:
            raise AssertionError('Frozen PubMed does not support ' + expected)
        original = source / 'references' / (accession + '.idf.txt')
        content = original.read_bytes()
        line = next(line for line in content.splitlines(keepends=True) if line.startswith(b'Publication DOI\t'))
        values = line.rstrip(b'\r\n').split(b'\t')[1:]
        populated = [v for v in values if v]
        assert len(populated) == 1
        supplied = populated[0]
        record = {'accession': accession, 'status': action, 'pmid': pmid, 'original_doi': supplied.decode(),
            'pubmed_supported_candidate': expected, 'original_file': str(original),
            'original_sha256': sha256(original), 'evidence': supported, 'created_at': utc()}
        if action == 'corrected-derivative':
            assert supplied.decode() == pmid
            corrected_line = line.replace(b'\t' + supplied, b'\t' + expected.encode(), 1)
            corrected = dest / (accession + '.corrected.idf.txt')
            corrected.write_bytes(content.replace(line, corrected_line, 1))
            record.update(corrected_file=corrected.name, corrected_sha256=sha256(corrected),
                change='Only the Publication DOI cell is corrected; the sealed source remains unchanged.')
            shutil.copyfile(source / 'references' / (accession + '.sdrf.txt'), dest / (accession + '.sdrf.txt'))
        else:
            record['change'] = 'No replacement: the converter continues to preserve the supplied placeholder pending curator review.'
        write_json(dest / (accession + '.json'), record)
        records.append(record)
    write_json(dest / 'records.json', records)


def contact_acceptance(source, root):
    """Check the full frozen study, including contributor and sample-reference retention."""
    relative = Path('cases/sra-SRP045534/standard')
    old = read(source / relative / 'source-json/payload-0.json')[0]
    new = read(root / relative / 'source-json/payload-0.json')[0]
    def refs(package):
        return [(s.get('iid'), s.get('contact_ref'), s.get('contributor_ref')) for s in package['sample']]
    def contacts(path):
        document = read(path)[0]
        rows = [r for r in document if r[0].startswith('Person ')]
        return [dict(zip([r[0] for r in rows], values)) for values in zip(*(r[1:] for r in rows))
                if dict(zip([r[0] for r in rows], values)).get('Person Last Name') == 'Smyth']
    previous = contacts(source / relative / 'saved-0-magetab/payload-0.json')
    outputs = {route: contacts(root / relative / route / 'payload-0.json')
               for route in ('replayed-magetab', 'saved-0-magetab', 'legacy-saved-0-magetab')}
    checks = {'baseline_four': len(previous) == 4,
        'json_contributors_preserved': old['contributor'] == new['contributor'],
        'sample_references_preserved': refs(old) == refs(new)}
    for route, people in outputs.items():
        matched = [p for p in people if p.get('Person Email') == 'smyth@wehi.edu.au']
        checks[route] = len(people) == 3 and len(matched) == 1 and '3052' in matched[0]['Person Address'] \
            and matched[0]['Person Roles'] == 'submitter' and sum(not p.get('Person Email') for p in people) == 2
    record = {'baseline': previous, 'candidate': outputs, 'checks': checks}
    write_json(root / 'contact-acceptance.json', record)
    return all(checks.values())


def validate(summary, legacy, frozen):
    """Gate execution, independent source preservation and route agreement."""
    routes = summary['live_routes']
    checks = {
        '81_fresh_routes': routes['observed'] == 81 and routes['execution'] == {'succeeded': 81} and routes['validation'] == {'valid': 81},
        '19_deposited_sdrfs': summary['curated_input_sdrf_preserved'] == 19,
        '19_deposited_idfs': summary['curated_input_idf_source_verdicts'].get(LOSS, 0) == 0,
        'independent_source_no_loss': summary['source_checks_off_and_curators'].get(LOSS, 0) == 0,
        '27_direct_fresh_json_pairs': all(m.get('direct_vs_saved_json', {}).get('sdrf_verdict') == 'equivalent'
            and not m['direct_vs_saved_json']['idf_discrepancies'] and not m['replay_limitations']
            for c in summary['cases'] for m in c['modes'].values()),
        '27_older_json_exports': len(legacy) == 27 and all(r['input_unchanged'] and r['result']['items'][0]['execution'] == 'succeeded'
            and r['result']['items'][0]['validation'] == 'valid' for r in legacy),
        '27_older_json_metadata_preserved': len(legacy) == 27 and all(r['sdrf']['verdict'] == 'preserved'
            and not any(f['classification'] == LOSS for f in r['idf']) for r in legacy),
        'frozen_unicode': all(v['author_identity_order_unicode'] is True for k, v in frozen['modes'].items() if k != 'off'),
    }
    return checks


def review(source, candidate, output, expected_revision):
    """Reassess sealed outputs without reconversion; record the analysis revision separately."""
    source, candidate, output = (Path(p).resolve() for p in (source, candidate, output))
    if any(output == p or output.is_relative_to(p) or p.is_relative_to(output) for p in (source, candidate)):
        raise ValueError('Review output must be separate from both sealed evidence bundles')
    if output.exists():
        raise FileExistsError('Use a fresh review directory')
    environment = read(candidate / 'candidate-environment.json')
    expected = subprocess.check_output(['git', 'rev-parse', '--verify', expected_revision + '^{commit}'], cwd=REPO, text=True).strip()
    if expected != environment['tested_revision']:
        raise ValueError('Review revision does not match the recorded conversion revision')
    if subprocess.check_output(['git', 'status', '--porcelain', '--untracked-files=all'], cwd=REPO, text=True):
        raise ValueError('Commit the analysis code before reviewing sealed outputs')
    record = {'conversion_revision': expected,
        'analysis_revision': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=REPO, text=True).strip(),
        'source_bundle': str(source), 'candidate_bundle': str(candidate),
        'source_seal': verify_seal(source), 'candidate_seal': verify_seal(candidate),
        'reviewed_at': utc(), 'conversion_rerun': False, 'results': []}
    assert record['source_seal'] == environment['source_seal_before'], 'Review must use the exact original acquisition bundle'
    for saved in read(candidate / 'legacy-results.json'):
        name, mode = saved['case'], saved['mode']
        relative = Path('cases') / name / mode
        previous_path = source / relative / 'saved-0-magetab/payload-0.json'
        current_path = candidate / relative / 'legacy-saved-0-magetab/payload-0.json'
        previous, current = read(previous_path)[0], read(current_path)[0]
        original_sdrf = next(r[1] for r in previous if r[0] == 'SDRF File')
        reference = next(ref or accession for provider, accession, ref in CASES if provider + '-' + accession == name)
        sdrf = assert_deposited_sdrf(original_sdrf, without_added_protocol_placeholders(previous, current, allow_mixed_rows=True), accession=reference)
        fields = legacy_idf_comparison(previous, current, name, mode)
        assert not any(r['classification'] == LOSS for r in fields), (name, mode)
        record['results'].append({'case': name, 'mode': mode, 'input_unchanged': saved['input_unchanged'],
            'baseline_output_sha256': sha256(previous_path), 'candidate_output_sha256': sha256(current_path),
            'sdrf': {'verdict': 'preserved', **sdrf}, 'idf': fields})
    assert len(record['results']) == 27
    assert verify_seal(source) == record['source_seal'] and verify_seal(candidate) == record['candidate_seal']
    output.mkdir(parents=True)
    write_json(output / 'legacy-metadata-review.json', record)
    seal_result(output)
    return record


def write_report(root, summary, environment, checks, legacy):
    modes = [m for c in summary['cases'] for m in c['modes'].values()]
    lines = ['# MSC authorship and metadata fix validation', '',
        f"Tested candidate: `{environment['tested_revision']}`. Frozen acquisition: `{BASELINE}`. No live acquisition.", '',
        'The original audit remains sealed. This result uses freshly generated candidate JSON for all JSON-to-MAGE-TAB routes; older JSON is assessed separately.', '',
        '| Validation | Result |', '| --- | --- |']
    lines.extend(f'| {key} | {"PASS" if value else "FAIL"} |' for key, value in checks.items())
    lines.extend(['', f"Independent off/curators source checks ({summary['source_check_denominator']}): `{summary['source_checks_off_and_curators']}`.",
        f"Curated-input IDF: {summary['curated_input_idf_source_verdicts'].get(EQUIVALENT, 0)}/{summary['curated_input_idf_source_populated_cells']} populated source cells preserved/equivalent.",
        f"Curated-input SDRF: {summary['curated_input_sdrf_populated_cells']:,} eligible source cells; {summary['curated_input_sdrf_excluded_cells']:,} derived-file cells excluded by the existing export policy.", '',
        '## Issues and fixes', '',
        '- Native PubMed initials/suffix loss → use CollectiveName or ForeName/Initials + LastName + Suffix, preserving order and Unicode.',
        '- Untagged GEO characteristics lost → normalize scalars as unspecified before harmonization, retaining occurrences and channel scope.',
        '- Repeated RelatedExperiment rows lost → import every populated cell; JSON retains occurrences and IDF exports unique targets in first-seen order.',
        '- Duplicate compatible contacts → consolidate only supported compatible identities, compare structured components before rendering, and retain fuller details; conflicting identifiers/addresses and ambiguous chains remain separate.',
        '- E-GEOD-100 PMID in DOI field → labelled corrected derivative only, with frozen PubMed evidence and both hashes.',
        '- E-MTAB-6486 placeholder DOI → curator-review record only; source precedence remains unchanged.', '',
        '## Interpretation and limits', '',
        'Source fidelity, curated-reference agreement and curated-input preservation remain separate in summary.json and field-discrepancies.jsonl. '
        'Reference-only annotations do not inflate converter-failure counts. Standard results consuming the comparison reference are explicitly reference-assisted; '
        'the actual candidate replay request traces, not merely the acquisition trace, establish this.',
        'ERP000263 standard remains partial on a scientific source conflict. DRP000001 has no verified curated counterpart. '
        'The frozen GSE328265 input contains one sample. Independent native checks cover the declared scopes; unscoped literal matches and ambiguous graph joins remain unresolved.',
        'The stable main base already adds validated blank required-protocol definitions/columns. These are projected only after strict checks of type, ontology, blank description, position and absence of qualifiers; existing source protocols are never removed.',
        'legacy-results.json compares all older-JSON exports with their frozen output and records intentional changes separately. Previously discarded source information cannot be recovered by exporting older JSON.', '',
        f"Candidate JSON differs from baseline capture in {sum(m.get('candidate_json_equals_baseline_capture') is False for m in modes)}/27 cases/presets; this is not a replay mismatch.",
        f"Older-JSON SDRF comparisons: {dict(Counter(r['sdrf']['verdict'] for r in legacy))}.",
        f"Validation interpreter `{environment['python']}`, Python {environment['python_version']}, installed pytest {environment['pytest_version']}. "
        'No integration pins, original fixtures, main integration or remote publication are part of this candidate.', '',
        '## Reproduce', '', '```sh',
        f'PYTHONDONTWRITEBYTECODE=1 {sys.executable} scripts/audit_authorship.py candidate \\',
        f'  --source-bundle {environment["source_bundle"]} \\',
        f'  --candidate-revision {environment["tested_revision"]} --bundle /path/to/fresh-results', '```', ''])
    (root / 'REPORT.md').write_text('\n'.join(lines))


def run(source, output, expected_revision):
    source, root, environment = candidate_preflight(source, output, expected_revision)
    root.mkdir(parents=True)
    write_json(root / 'candidate-environment.json', environment)
    shutil.copytree(source / 'references', root / 'references')
    for provider, accession, _ in CASES:
        name = provider + '-' + accession
        original, dest = source / 'cases' / name, root / 'cases' / name
        shutil.copytree(original / 'http', dest / 'http')
        shutil.copyfile(original / 'acquisition.json', dest / 'acquisition.json')
        for mode in ('off', 'curators', 'standard'):
            shutil.copytree(original / mode / 'source-json', dest / mode / 'baseline-source-json')
        replay_case(dest, fresh_json=True)
    deposited(root, sealed_references=source / 'references', revision=environment['tested_revision'])
    # The contract helper reads its resources through this explicitly selected root.
    from pytest import MonkeyPatch
    from tests.support import contracts
    with MonkeyPatch.context() as patches:
        patches.setattr(contracts, 'GEO', source / 'frozen-GSE328265')
        frozen_unicode(root, sealed=True)
    legacy = legacy_exports(source, root)
    reference_corrections(source, root)
    summary = build(root, candidate=True)
    checks = validate(summary, legacy, read(root / 'frozen-GSE328265/assessment.json'))
    checks['SRP045534_contacts_four_to_three'] = contact_acceptance(source, root)
    assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=REPO, text=True).strip() == environment['tested_revision'], 'Candidate HEAD changed during replay'
    assert not subprocess.check_output(['git', 'status', '--porcelain', '--untracked-files=all'], cwd=REPO, text=True), 'Candidate checkout changed during replay'
    environment.update(source_seal_after=verify_seal(source), finished_at=utc())
    assert environment['source_seal_after'] == environment['source_seal_before']
    write_json(root / 'candidate-environment.json', environment)
    write_json(root / 'candidate-checks.json', checks)
    write_report(root, summary, environment, checks, legacy)
    seal_result(root)
    if not all(checks.values()):
        raise AssertionError('Candidate validation failed: ' + ', '.join(k for k, v in checks.items() if not v))
