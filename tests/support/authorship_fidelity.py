"""Independent full-IDF and paired-source audit oracles (no MSC imports).

Differences are evidence statements, not a truth vote against a reference.
Repeated occurrences and group ownership survive normalization. In particular,
reordering a publication column is legal; moving its DOI to another PMID is not.
"""
from collections import Counter, defaultdict
import csv
import io
import json
from pathlib import Path
import re

from tests.support.deposited_fidelity import _scoped_table, HEADER_ALIASES

EQUIVALENT = 'equivalent-or-documented-transformation'
LOSS = 'source-supported-loss-alteration-or-misassociation'
UNSUPPORTED = 'unsupported-introduction'
CURATION = 'reference-only-annotation-requiring-assessment'
ENRICHMENT = 'evidence-supported-correction-or-enrichment'
UNRESOLVED = 'source-version-reference-defect-or-unresolved'

IDF_ALIASES = {'Person Roles Term Source Ref': 'Person Roles Term Source REF',
               'Status Term Source Ref': 'Publication Status Term Source REF',
               'Status Term Accession Number': 'Publication Status Term Accession Number',
               **{f'Comment [{name}]': f'Comment[{name}]' for name in (
                   'AEExperimentType', 'AdditionalFile:fa', 'SecondaryAccession', 'SequenceDataURI', 'Submitted Name', 'RelatedExperiment')}}
GROUPS = {
    'publication': ('PubMed ID', 'Publication '),
    'contact': ('Person ',), 'protocol': ('Protocol ',),
    'factor': ('Experimental Factor ',), 'design': ('Experimental Design',),
    'ontology': ('Term Source ',), 'quality': ('Quality Control ',),
    'replicate': ('Replicate ',), 'normalization': ('Normalization ',),
}


def present(value):
    return value is not None and value != ''


def table(path):
    return list(csv.reader(io.StringIO(Path(path).read_bytes().decode('utf-8-sig')), delimiter='\t'))


def idf_cells(rows, *, imported_scalars=False):
    """Populated cell inventory keyed by semantic group and repeated identity.

    Trailing empty padding and physical row/column order are immaterial. Internal
    blanks remain in groups before identity assignment. No text case, whitespace,
    authors, protocol or ontology values are normalized here.
    """
    fields = defaultdict(list)
    for row in rows:
        if not row:
            continue
        label = IDF_ALIASES.get(row[0], row[0])
        if label == 'SDRF File':
            continue  # independently verified companion relationship, never embedded data
        values = list(row[1:])
        while values and not present(values[-1]):
            values.pop()
        fields[label].extend(values)
    grouped = defaultdict(dict)
    for label, values in fields.items():
        group = next((key for key, prefixes in GROUPS.items() if label.startswith(prefixes)), 'study')
        grouped[group][label] = values
    cells = {}
    for group, values in grouped.items():
        if group == 'study':
            for label, items in values.items():
                occurrence = Counter()
                for value in items:
                    if present(value):
                        # Related studies are independent links, not positional
                        # attributes of another IDF group. Keep multiplicity
                        # without requiring physical repeated-row order.
                        identity = label + ':' + str(value) if label == 'Comment[RelatedExperiment]' else label
                        occurrence[identity] += 1
                        cells[(group, identity, occurrence[identity], label)] = value
            continue
        identities = Counter()
        for index in range(max(map(len, values.values()), default=0)):
            record = {key: items[index] for key, items in values.items()
                      if index < len(items) and present(items[index])}
            if not record:
                continue
            if group == 'publication':
                identity = ('pmid', record['PubMed ID']) if record.get('PubMed ID') else (
                    ('doi', record['Publication DOI']) if record.get('Publication DOI') else ('unidentified', index))
            elif group == 'contact':
                identity = tuple(record.get('Person ' + key, '').strip() if imported_scalars else record.get('Person ' + key, '')
                                 for key in ('First Name', 'Mid Initials', 'Last Name', 'Email'))
            else:
                keys = {'protocol': 'Protocol Name', 'factor': 'Experimental Factor Name',
                        'ontology': 'Term Source Name', 'design': 'Experimental Design',
                        'quality': 'Quality Control Type', 'replicate': 'Replicate Type',
                        'normalization': 'Normalization Type'}
                identity = (record.get(keys[group], ''),)
            identity = json.dumps(identity, ensure_ascii=False)
            identities[identity] += 1
            for label, value in record.items():
                cells[(group, identity, identities[identity], label)] = value
    return cells


def classify(source, generated, reference=None, *, supplied=True, source_complete=False,
             transformation=None, evidence=None):
    """Conservative triad classification; absent evidence is not evidence of absence."""
    if transformation:
        return EQUIVALENT, transformation
    if source == generated and present(source):
        if present(reference) and reference != generated:
            return UNRESOLVED, 'Source preserved; deposited reference differs and needs version/curation assessment.'
        return EQUIVALENT, 'Exact source value preserved.'
    if supplied and present(source) and source != generated:
        return LOSS, 'Value supplied to the converter is missing, changed or bound to a different identity.'
    if present(reference) and reference == generated and not present(source) and not evidence:
        return EQUIVALENT, 'Reference agreement; independence depends on the request trace.'
    if present(generated) and not present(source):
        if evidence:
            return ENRICHMENT, evidence
        if source_complete:
            return UNSUPPORTED, 'No support in the completely assessed source scope.'
        return UNRESOLVED, 'Generated value requires additional source/projection evidence.'
    if present(reference) and not present(source) and not present(generated):
        return CURATION, 'Reference supplies this annotation; no counterpart was established in the assessed mapping. Human curation is not proven.'
    if present(reference) and reference == generated:
        return EQUIVALENT, 'Reference agreement; independence depends on the request trace.'
    return UNRESOLVED, 'Available evidence does not resolve the discrepancy.'


def idf_comparison(source_rows, actual_rows, reference_rows=None, *, same_input=True, documented=False, source_sdrf=None):
    source = idf_cells(source_rows, imported_scalars=documented)
    actual = idf_cells(actual_rows, imported_scalars=documented)
    reference = idf_cells(reference_rows or [], imported_scalars=documented)
    result = []
    for key in sorted(source.keys() | actual.keys() | reference.keys(), key=str):
        left, right, ref = source.get(key), actual.get(key), reference.get(key)
        classification, reason = classify(left, right, ref, supplied=same_input)
        if documented and same_input and left != right:
            group, identity, occurrence, field = key
            if isinstance(left, str) and isinstance(right, str) and left.strip() == right and (
                    group == 'contact' or field == 'Experiment Description'):
                classification, reason = EQUIVALENT, 'AEParser strips boundary whitespace in contact scalars and study description; internal text is retained.'
            elif group == 'ontology' and occurrence > 1 and right is None:
                previous = {k[-1]: v for k, v in source.items() if k[:3] == (group, identity, 1)}
                current = {k[-1]: v for k, v in source.items() if k[:3] == (group, identity, occurrence)}
                retained = {k[-1]: v for k, v in actual.items() if k[:3] == (group, identity, 1)}
                if previous == current == retained:
                    classification, reason = EQUIVALENT, 'Identical named ontology declaration is consolidated; every declaration value is retained.'
            elif left is None and field == 'MAGE-TAB Version' and right == '1.1':
                classification, reason = EQUIVALENT, 'Exporter declares its MAGE-TAB output version.'
            elif left is None and field in {'Investigation Accession', 'Comment[ArrayExpressAccession]'} and right in source.values():
                classification, reason = EQUIVALENT, 'Existing explicit accession is also emitted in the canonical accession field.'
            elif left is None and field in {'Term Source Name', 'Comment[SecondaryAccessionTermSourceRef]'} and right in {'ArrayExpress', 'GEO', 'ENA', 'SRA'}:
                patterns = {'ArrayExpress': r'E-[A-Z]+-\d+', 'GEO': r'GSE\d+', 'ENA': r'ERP\d+', 'SRA': r'SRP\d+'}
                if any(re.fullmatch(patterns[right], str(v)) for v in source.values()):
                    classification, reason = EQUIVALENT, 'Repository vocabulary label is declared for an explicit source accession.'
            elif left is None and field == 'Comment[SequenceDataURI]' and isinstance(right, str):
                match = re.fullmatch(r'http://www\.ebi\.ac\.uk/ena/data/view/([SED]RR)(\d+)-\1(\d+)', right)
                if match and source_sdrf:
                    prefix, start, end = match.groups()
                    runs = {v for row in source_sdrf[1:] for v in row if re.fullmatch(r'[SED]RR\d+', str(v))}
                    if 0 <= int(end) - int(start) < 10000:
                        declared = {prefix + str(i).zfill(len(start)) for i in range(int(start), int(end) + 1)}
                        if declared == runs:
                            classification, reason = ENRICHMENT, 'Sequence browser range denotes exactly the source SDRF run identities.'
        result.append({'scope': 'IDF', 'identity': list(key[:-1]), 'field': key[-1],
                       'source': left, 'generated': right, 'reference': ref,
                       'classification': classification, 'reason': reason})
    return result


def assert_idf_preserved(source, actual, *, transformations=None, documented=False):
    """Source preservation, with narrowly declared per-cell transformations."""
    differences = idf_comparison(source, actual, documented=documented)
    transformations = transformations or {}
    failures = []
    for item in differences:
        if item['classification'] != LOSS:
            continue
        key = (*item['identity'], item['field'])
        if key in transformations and (item['source'], item['generated']) == transformations[key]:
            continue
        failures.append(item)
    assert not failures, json.dumps(failures, ensure_ascii=False, indent=2)
    return {'source_populated_cells': len(idf_cells(source)), 'compared_cells': len(differences)}


def scoped_sdrf_cells(rows, accession):
    """One record per populated occurrence, including node/qualifier ownership."""
    graph, descriptors, values, headers = _scoped_table(rows, accession)
    result = []
    for index, row in enumerate(values, 1):
        for descriptor, header, (_, value) in zip(descriptors, headers, row):
            if present(value):
                result.append({'row': index, 'owner': descriptor, 'field': header, 'value': value})
    return graph, result


ACCESSION = re.compile(r'(?<![A-Z0-9])(?:GSM\d+|[SED]R[SRX]\d+|SAM[END][A-Z]?\d+)(?!\d)')


def row_identity(headers, row):
    """Public sample/experiment/run identities plus explicit label and file keys."""
    ids = {'sample': set(), 'assay': set(), 'run': set(), 'file': set(), 'channel': set(), 'channel_source': set()}
    for header, value in zip(headers, row):
        text = str(value or '')
        header = HEADER_ALIASES.get(header, header)
        identity_column = (header in {'Source Name', 'Sample Name', 'Extract Name', 'Labeled Extract Name',
            'Assay Name', 'Hybridization Name', 'Scan Name', 'Array Data File', 'Derived Array Data File',
            'Array Data Matrix File', 'Derived Array Data Matrix File', 'Comment[FASTQ_URI]', 'Comment[File URI]'}
            or bool(re.fullmatch(r'Comment\[(?:ENA|SRA|BioSD)_(?:SAMPLE|EXPERIMENT|RUN)\]', header)))
        for match in ACCESSION.findall(text) if identity_column else []:
            kind = 'run' if re.match(r'[SED]RR', match) else 'assay' if re.match(r'[SED]RX', match) else 'sample'
            ids[kind].add(match)
        if header == 'Label' and text:
            ids['channel'].add(text.casefold())
        if header == 'Comment[Sample_source_name]' and text:
            ids['channel_source'].add(text)
        if ('Data File' in header or 'FASTQ_URI' in header) and text:
            ids['file'].add(text.rsplit('/', 1)[-1])
    return {key: sorted(value) for key, value in ids.items()}


def align_sdrf(reference, actual):
    """Conservative one-to-one joins; unresolved multiplicity stays visible.

    Match by file then run then assay then sample, with channel agreement where
    both declare it. Never infer alignment from row order or accession names.
    Exact duplicate identities retain separate occurrences, but ambiguous
    nonidentical rows are not arbitrarily assigned.
    """
    left = [row_identity(reference[0], row) for row in reference[1:]]
    right = [row_identity(actual[0], row) for row in actual[1:]]
    available = set(range(len(right)))
    matches, unresolved = [], []
    for index, identity in enumerate(left):
        selected, method = [], None
        for kind in ('file', 'run', 'assay', 'sample'):
            candidates = [j for j in available if set(identity[kind]) & set(right[j][kind])
                          and (not identity['channel'] or not right[j]['channel'] or identity['channel'] == right[j]['channel'])
                          and (not identity['channel_source'] or not right[j]['channel_source'] or identity['channel_source'] == right[j]['channel_source'])
                          and (not identity['sample'] or not right[j]['sample'] or set(identity['sample']) & set(right[j]['sample']))]
            if candidates:
                selected, method = candidates, kind
                break
        if len(selected) > 1 and len({tuple(actual[j + 1]) for j in selected}) == 1:
            selected = [min(selected)]
        if len(selected) == 1:
            j = selected[0]
            available.remove(j)
            matches.append({'reference_row': index + 1, 'generated_row': j + 1, 'method': method,
                            'reference_identity': identity, 'generated_identity': right[j]})
        else:
            unresolved.append({'reference_row': index + 1, 'identity': identity,
                               'candidate_rows': [j + 1 for j in selected],
                               'reason': 'ambiguous multiplicity' if selected else 'no shared identity'})
    return {'matches': matches, 'unresolved': unresolved, 'generated_unmatched_rows': [j + 1 for j in sorted(available)],
            'reference_rows': len(left), 'generated_rows': len(right)}


def pair_sdrf(reference, actual, accession):
    """Inventory every populated cell; classify only comparisons actually joined.

    Graph ownership descriptors must also agree to compare a field automatically.
    Different graphs or unresolved row mappings are a verification limitation,
    not proof of metadata loss or of curator intervention.
    """
    alignment = align_sdrf(reference, actual)
    lgraph, lcells = scoped_sdrf_cells(reference, accession)
    rgraph, rcells = scoped_sdrf_cells(actual, accession)
    by_left, by_right = defaultdict(list), defaultdict(list)
    for cell in lcells: by_left[cell['row']].append(cell)
    for cell in rcells: by_right[cell['row']].append(cell)
    compared, used_left, used_right = [], set(), set()
    lcounts, rcounts = Counter(lgraph), Counter(rgraph)
    def key(cell, graph):
        node, owner, qualifier = cell['owner']
        if node < 0:
            return ('row', owner, qualifier)
        name = graph[node]
        if lgraph != rgraph and (lcounts[name] != 1 or rcounts[name] != 1):
            return None  # no guessed correspondence between repeated processing steps
        return (name, graph[:node].count(name), owner, qualifier)
    for match in alignment['matches']:
        expected = by_left[match['reference_row']]
        generated = by_right[match['generated_row']]
        fields = {str(key(c, rgraph)): c for c in generated if key(c, rgraph) is not None}
        for cell in expected:
            signature = key(cell, lgraph)
            other = fields.get(str(signature)) if signature is not None else None
            classification = EQUIVALENT if other and cell['value'] == other['value'] else UNRESOLVED
            if signature is not None and other is None:
                classification = CURATION
            used_left.add((cell['row'], str(cell['owner'])))
            if other: used_right.add((other['row'], str(other['owner'])))
            compared.append({'scope': 'SDRF', 'identity': match, 'field': cell['field'], 'owner': cell['owner'],
                'source': None, 'generated': other['value'] if other else None, 'reference': cell['value'],
                'classification': classification, 'reason': 'Exact value and owner agreement.' if classification == EQUIVALENT
                else ('Reference supplies this annotation at a corresponding owner; human curation is not established.' if classification == CURATION
                      else 'Source provenance or graph correspondence requires assessment; not counted as converter failure.')})
    for inventory, used, side in [(lcells, used_left, 'reference'), (rcells, used_right, 'generated')]:
        for cell in inventory:
            if (cell['row'], str(cell['owner'])) in used: continue
            compared.append({'scope': 'SDRF', 'identity': {side + '_row': cell['row']},
                'field': cell['field'], 'owner': cell['owner'], 'source': None,
                'generated': cell['value'] if side == 'generated' else None,
                'reference': cell['value'] if side == 'reference' else None,
                'classification': UNRESOLVED, 'reason': 'Unaligned row or owner: no positional inference; inspect source and both inventories.'})
    return {'alignment': alignment, 'same_graph': lgraph == rgraph,
            'reference_populated_cells': len(lcells), 'generated_populated_cells': len(rcells),
            'reference_inventory': lcells, 'generated_inventory': rcells, 'comparisons': compared}
