# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Independently review intentional changes in the older additional JSON exports."""
from copy import deepcopy
import json
from pathlib import Path
import traceback
from tests.support.integrity_legacy import reviewed_legacy_tables
from tests.support.deposited_fidelity import without_added_protocol_placeholders, assert_deposited_sdrf, NODE_HEADERS
from tests.support.authorship_fidelity import idf_comparison, LOSS
from scripts.audit_authorship_candidate import verify_seal, sha256


def placeholder_names(package, old, new):
    original = deepcopy(old)
    def definitions(document):
        fields = {r[0]: r[1:] for r in document if r[0].startswith('Protocol ')}
        return {name: {key: values[i] for key, values in fields.items() if key != 'Protocol Name'}
                for i, name in enumerate(fields['Protocol Name'])}
    left, right = definitions(old), definitions(new)
    authored = {p['name'] for p in package['series'].get('protocols', [])}
    mapping = {}
    for name in left.keys() - right.keys():
        facts = left[name]
        assert name not in authored
        expected = {'sample collection protocol': 'EFO_0005518', 'nucleic acid sequencing protocol': 'EFO_0004170'}
        assert facts['Protocol Type'] in expected and facts['Protocol Term Source REF'] == 'EFO'
        assert facts['Protocol Term Accession Number'] == expected[facts['Protocol Type']]
        assert all(not value for key, value in facts.items() if key not in {'Protocol Type', 'Protocol Term Source REF', 'Protocol Term Accession Number'})
        matches = [n for n in right.keys() - left.keys() if right[n] == facts]
        assert len(matches) == 1
        mapping[name] = matches[0]
    for row in original:
        if row[0] == 'Protocol Name': row[1:] = [mapping.get(n, n) for n in row[1:]]
        if row[0] == 'SDRF File':
            for i, header in enumerate(row[1][0]):
                if header == 'Protocol REF':
                    for cells in row[1][1:]: cells[i] = mapping.get(cells[i], cells[i])
    return original, mapping


def added_cels(package, old, new):
    new = deepcopy(new)
    if 'Array Data File' in old[0] or 'Array Data File' not in new[0]: return new, []
    indices = [i for i, h in enumerate(new[0]) if h == 'Array Data File']
    assert len(indices) == 1
    i = indices[0]
    end = next((j for j in range(i + 1, len(new[0])) if new[0][j] in NODE_HEADERS), len(new[0]))
    assert set(new[0][i + 1:end]) <= {'Comment[ArrayExpress FTP file]', 'Comment[File URI]'}
    samples = {s['iid']: s for s in package['sample']}
    changes = []
    for row in new[1:]:
        sample = samples[row[new[0].index('Source Name')]]
        expected = [f['value'].strip(' \t\r\n') for f in sample.get('supplementary_data', []) if f.get('type', '').casefold() == 'cel']
        assert expected == [row[i]]
        assert all(value == row[i] for value in row[i + 1:end])
        changes.append({'sample': sample['iid'], 'restored_CEL': row[i]})
    return [[value for j, value in enumerate(row) if not i <= j < end] for row in new], changes


def completed_native_files(package, old, new):
    """Verify each added field against an exact source sample/run/file URI."""
    new = deepcopy(new)
    added = set(new[0]) - set(old[0])
    checked = {'Comment[FASTQ_MD5]': 'md5', 'Comment[FASTQ_BYTES]': 'bytes'}
    indices = [i for i, header in enumerate(new[0]) if header in added & checked.keys()]
    samples = {s['iid']: s for s in package['sample']}
    changes = []
    for row in new[1:]:
        for index in indices:
            if not row[index]: continue
            sample = samples[row[new[0].index('Source Name')]]
            runs = {row[j] for j, h in enumerate(new[0]) if h in {'Comment[ENA_RUN]', 'Comment[SRA_RUN]'}}
            uri_index = max(j for j in range(index) if new[0][j] == 'Comment[FASTQ_URI]')
            uri = row[uri_index]
            files = [f for run in sample.get('sra_run', []) if run['run'] in runs
                     for f in [*run.get('files', []), *run.get('fastq_files', [])] if f.get('uri') == uri]
            formats = {f['format'].lower() for f in files if f.get('format')}
            assert files and formats and formats <= {'fastq', 'fq', 'fastq.gz', 'fq.gz'}
            field = checked[new[0][index]]
            assert {f[field] for f in files if f.get(field)} == {row[index]}
            changes.append({'sample': sample['iid'], 'runs': sorted(runs), 'uri': uri,
                            'field': field, 'source_value': row[index]})
    return [[v for i, v in enumerate(row) if i not in indices] for row in new], changes


def review(source, root, output):
    source, root, output = (Path(p).resolve() for p in (source, root, output))
    if output.exists() or any(output.is_relative_to(p) for p in (source, root)):
        raise ValueError('Use a fresh review output outside sealed inputs')
    seals = {str(p): verify_seal(p) for p in (source, root)}
    results = []
    for saved in root.glob('cases/*/*/legacy-saved-0.json'):
        relative = saved.parent.relative_to(root)
        case, mode = relative.parts[1:]
        original = json.loads((source / relative / 'saved-0-magetab/payload-0.json').read_text())[0]
        new = json.loads((saved.parent / 'legacy-saved-0-magetab/payload-0.json').read_text())[0]
        package = json.loads(saved.read_text())
        report = {'case': case, 'mode': mode, 'input_sha256': sha256(saved),
                  'historical_output_sha256': sha256(source / relative / 'saved-0-magetab/payload-0.json'),
                  'candidate_output_sha256': sha256(saved.parent / 'legacy-saved-0-magetab/payload-0.json')}
        try:
            old, names = placeholder_names(package, original, new)
            report['placeholder_names'] = names
            report['idf'] = idf_comparison(old, new, documented=True)
            assert not any(f['classification'] == LOSS for f in report['idf']), 'IDF source loss'
            a = next(row[1] for row in old if row[0] == 'SDRF File')
            b = without_added_protocol_placeholders(old, new, allow_mixed_rows=True)
            b, cels = added_cels(package, a, b)
            report['restored_CELs'] = cels
            definitions = dict(zip(next(row[1:] for row in old if row[0] == 'Protocol Name'), next(row[1:] for row in old if row[0] == 'Protocol Description')))
            if case.startswith('geo-'):
                a, b, changes = reviewed_legacy_tables(package, a, b, definitions)
            else:
                b, changes = completed_native_files(package, a, b)
            report['reviewed_transformations'] = changes
            report['stats'] = assert_deposited_sdrf(a, b, accession=case.split('-', 1)[1].replace('GSE', 'E-GEOD-'))
            report['passed'] = True
        except Exception:
            report.update(passed=False, error=traceback.format_exc())
            print(case, mode, report['error'][-1000:])
        results.append(report)
    assert all(verify_seal(p) == seals[str(p)] for p in (source, root))
    with output.open('x') as stream: json.dump({'results': results, 'input_seals': seals}, stream, indent=2)
    print('passed', sum(r['passed'] for r in results), 'of', len(results))
    return results


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('source', 'root', 'output'): parser.add_argument('--' + name, type=Path, required=True)
    if not all(r['passed'] for r in review(**vars(parser.parse_args()))): raise SystemExit(1)
