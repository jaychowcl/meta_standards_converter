# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Independent, occurrence-aware checks for the concrete cohort repair targets."""
import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path

from scripts.audit_authorship_candidate import verify_seal, sha256
from tests.support.authorship_sources import documents, provenance


def read(path): return json.loads(Path(path).read_text())
def normalized(value): return ' '.join(str(value or '').split())
def nodes(data, kind):
    return [s for p in data['series']['assay_paths'] for s in p['steps'] if s['kind'] == kind]
def comments(node): return {c['name'].upper(): c['value'] for c in node.get('comments', [])}


def check(source, candidate, output):
    source, candidate, output = map(Path, (source, candidate, output))
    assert not output.exists()
    seals = {str(p): verify_seal(p) for p in (source, candidate)}
    records = []
    def record(issue, mode, entity, field, expected, actual, evidence, reference=None, reason=None):
        records.append({'issue': issue, 'mode': mode, 'entity': entity, 'field': field,
                        'source': expected, 'generated': actual, 'reference': reference,
                        'evidence': evidence, 'reason': reason, 'passed': expected == actual})
    modes = ('off', 'curators', 'standard')
    summary = read(candidate / 'summary.json')
    baseline = read(source / 'summary.json')
    record('blank-characteristic-owners', 'curated-input', 'five deposited references', 'eligible preserved cells',
           3194, summary['curated_input_sdrf_cell_review'].get('equivalent-or-documented-transformation'),
           [{'path': str(source / 'references'), 'baseline_preserved': baseline['curated_input_sdrf_cell_review'].get('equivalent-or-documented-transformation')}],
           reason='All 3194 eligible cells are preserved, including the 48 previously lost comments; 212 derived-file cells remain excluded by policy.')
    for case in ('geo-GSE13398', 'geo-GSE31356', 'geo-GSE52778'):
        docs = documents(source / 'cases' / case, 'off')
        raw = {s.get('iid'): (s, d) for d in docs if d['kind'] == 'xml' for s in d['xml'].iter('Sample')}
        for mode in modes:
            data = read(candidate / 'cases' / case / mode / 'saved-0.json')
            for sid, (sample, doc) in raw.items():
                evidence = [provenance(doc, f'Sample[@iid={sid}]/Channel')]
                paths = [p for p in data['series']['assay_paths'] if any(s.get('sample_ref') == sid for s in p['steps'])]
                if case == 'geo-GSE13398':
                    record('distinct-channel-materials', mode, sid, 'source/extract/assay counts', [2, 2, 1],
                           [len({s['name'] for p in paths for s in p['steps'] if s['kind'] == k}) for k in ('source', 'extract', 'assay')], evidence,
                           reason='Distinct sample/channel occurrences retain a shared assay; no donor identity is inferred.')
                elif case == 'geo-GSE31356':
                    expected = [n.text.strip() for n in sample.findall('Supplementary-Data') if n.get('type') == 'CEL']
                    actual = [s.get('link', {}).get('value', s['name']) for p in paths for s in p['steps'] if s['kind'] == 'array_data_file']
                    record('cel-links', mode, sid, 'raw CEL URI', expected, actual, evidence,
                           reason='Surrounding source URI whitespace is normalized only in the interpreted/export view.')
                else:
                    growth = normalized(sample.findtext('Channel/Growth-Protocol'))
                    treatment = normalized(sample.findtext('Channel/Treatment-Protocol'))
                    definitions = {p['name']: normalized(p.get('description')) for p in data['series']['protocols']}
                    sequences = [[definitions.get(s['protocol_ref']) for s in p['steps'] if s['kind'] == 'protocol_application'] for p in paths]
                    record('growth-before-treatment', mode, sid, 'every generated branch preserves chronology', True,
                           bool(sequences) and all(growth in s and treatment in s and s.index(growth) < s.index(treatment) for s in sequences), evidence)
    # Native original-file metadata is read directly from frozen SRA XML.
    for case in ('sra-SRP033351', 'sra-ERP148790', 'sra-ERP185509'):
        docs = documents(source / 'cases' / case, 'off')
        originals = {}
        for doc in docs:
            if doc['kind'] != 'xml': continue
            for run in doc['xml'].iter('RUN'):
                for file in run.iter('SRAFile'):
                    if file.get('supertype') == 'Original':
                        originals[run.get('accession'), file.get('filename')] = (file.attrib, provenance(doc, f'RUN[@accession={run.get("accession")}]/SRAFile[@filename={file.get("filename")}]'))
        assert originals
        for mode in modes:
            saved = candidate / 'cases' / case / mode / 'saved-0.json'
            if not saved.exists(): continue
            data = read(saved)
            for (run, name), (raw, evidence) in originals.items():
                supplied = [f for sample in data['sample'] for r in sample.get('sra_run', []) if r['run'] == run for f in r.get('files', [])
                            if f.get('filename') == name and f.get('role') == 'Original']
                record('submitted-provenance', mode, run + '/' + name, 'original filename and checksum', True,
                       any(f.get('md5') == raw.get('md5') for f in supplied), [evidence])
                if case != 'sra-ERP148790': continue
                matched = [n for p in data['series']['assay_paths'] if any(s['kind'] == 'scan' and s['name'] == run for s in p['steps'])
                           for n in p['steps'] if n['kind'] == 'array_data_file' and n['name'] == name]
                record('fastq-checksum-size', mode, run + '/' + name, 'every compatible occurrence has source checksum/size', True,
                       bool(matched) and all(comments(n).get('MD5', comments(n).get('FASTQ_MD5')) == raw['md5']
                           and comments(n).get('BYTES', comments(n).get('FILE SIZE', comments(n).get('FASTQ_BYTES'))) == raw['size'] for n in matched), [evidence],
                       reason='Matches are verified by sample/run and exact URI by the production resolver; names are used here to enumerate the independently frozen originals.')
    # ENA conventions apply only to explicit generated-file records for paired libraries.
    case = 'ena-SRP033351'
    data = read(candidate / 'cases' / case / 'off/saved-0.json')
    actual = Counter()
    for node in nodes(data, 'array_data_file'):
        facts = comments(node)
        actual[facts.get('READ_ROLE', 'unknown')] += 1
    record('ena-read-roles', 'off', case, 'paired mates versus unpaired reads', {'paired': 32, 'unpaired': 7}, dict(actual),
           [{'path': str(source / 'cases' / case / 'http/off.manifest.json'),
             'url': 'https://ena-docs.readthedocs.io/en/latest/retrieval/file-download/archive-generated-files.html'}],
           reason='16 paired-library runs contain 32 mate files plus seven explicitly unsuffixed generated read files; file count is not library layout.')
    for mode in modes:
        data = read(candidate / 'cases/sra-ERP185509' / mode / 'saved-0.json')
        values = {c.get('term_accession_number') for n in nodes(data, 'source') for c in n.get('characteristics', [])
                  if c['name'] == 'organism' and c.get('term_accession_number')}
        record('taxonomy-projection', mode, 'ERP185509', 'selected generated taxonomy identifiers', ['83332'], sorted(values),
               [{'path': str(source / 'cases/sra-ERP185509/source-snapshots.json')}],
               reason='Species-level authored assertions remain separately supplied annotations; no unsupported species/strain relation is invented.')
    pairs = [m['direct_vs_saved_json'] for c in summary['cases'] for m in c['modes'].values() if m.get('direct_vs_saved_json')]
    record('route-equivalence', 'all available', 'additional cohort', 'equivalent direct/fresh-JSON pairs', 18,
           sum(p['sdrf_verdict'] == 'equivalent' and not p['idf_discrepancies'] for p in pairs), [{'path': str(candidate / 'summary.json')}])
    result = {'checks': records, 'denominators': dict(Counter(r['issue'] for r in records)),
              'passed': sum(r['passed'] for r in records), 'failed': sum(not r['passed'] for r in records),
              'input_seals': seals, 'scope': 'Concrete semantic repair targets supplement the complete IDF/SDRF field comparisons, not a claim of complete scientific interpretation.'}
    assert all(verify_seal(p) == seals[str(p)] for p in (source, candidate))
    with output.open('x') as stream: json.dump(result, stream, indent=2, ensure_ascii=False)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('source', 'candidate', 'output'): parser.add_argument('--' + name, type=Path, required=True)
    result = check(**vars(parser.parse_args()))
    print(json.dumps({k: result[k] for k in ('denominators', 'passed', 'failed')}))
    if result['failed']: raise SystemExit(1)
