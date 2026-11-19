# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Independent deposited SDRF oracle; no converter or renderer imports.

Rows are multisets of exact typed cells keyed by node, attribute occurrence and
qualifier owner. Factors are row variables (MAGE-TAB 1.1 section 2.3.2(11)),
bound to that entire scientific row rather than their incidental node carrier.
Only the named, recorded source compatibility rules below apply.
"""
from collections import Counter
import csv
import hashlib
import io
import json
from pathlib import Path


DEPOSITED_PANEL = (
    'E-GEOD-100', 'E-GEOD-48034', 'E-GEOD-48035', 'E-GEOD-52564', 'E-GEOD-60450',
    'E-MTAB-13131', 'E-MTAB-13662', 'E-MTAB-14416', 'E-MTAB-14560', 'E-MTAB-14566',
    'E-MTAB-14960', 'E-MTAB-16192', 'E-MTAB-16253', 'E-MTAB-16847', 'E-MTAB-16856',
    'E-MTAB-308', 'E-MTAB-5481', 'E-MTAB-5485',
)

# Enumerated deposited spellings; no generic whitespace/case normalization.
HEADER_ALIASES = {
    **{f'Comment [{name}]': f'Comment[{name}]' for name in (
        'ENA_SAMPLE', 'Sample_description', 'Sample_source_name', 'Sample_title',
        'INSTRUMENT_MODEL', 'LIBRARY_LAYOUT', 'LIBRARY_SELECTION', 'LIBRARY_SOURCE',
        'LIBRARY_STRATEGY', 'ENA_EXPERIMENT', 'Platform_title', 'ENA_RUN', 'FASTQ_URI',
        'Derived ArrayExpress FTP file', 'Derived Array Data FTP File')},
    **{f'Characteristics [{name}]': f'Characteristics[{name}]' for name in (
        'developmental stage', 'genotype', 'immunophenotype', 'organism', 'organism part',
        'sex', 'strain/background', 'reference rna', 'reference rna 2',
        'reference rna 2 cat. #', 'reference rna 2 lot #', 'reference rna cat. #',
        'reference rna lot #', 'strain', 'cell type', 'Organism')},
    **{f'FactorValue [{name}]': f'Factor Value[{name}]' for name in (
        'developmental stage', 'immunophenotype', 'ORGANISM PART', 'REFERENCE RNA',
        'REFERENCE RNA CAT. #', 'cell type')},
    'Factor Value [single cell identifier]': 'Factor Value[single cell identifier]',
}

BLANK_LITERALS = {
    'E-GEOD-48034': {'Term Source REF': ('  ',), 'Term Accession Number': ('  ',),
        **{f'Characteristics[{name}]': ('  ',) for name in (
            'reference rna lot #', 'reference rna 2', 'reference rna 2 cat. #', 'reference rna 2 lot #')}},
    'E-GEOD-48035': {'Term Source REF': ('  ',), 'Term Accession Number': ('  ',),
        **{f'Characteristics[{name}]': ('  ',) for name in (
            'reference rna lot #', 'reference rna 2', 'reference rna 2 cat. #', 'reference rna 2 lot #')}},
    'E-GEOD-52564': {'Term Source REF': ('  ',), 'Term Accession Number': ('  ', '        ')},
    'E-GEOD-60450': {'Term Source REF': ('  ',), 'Term Accession Number': ('  ',)},
    'E-MTAB-308': {'Characteristics[Age]': ('  ',), 'Unit[TimeUnit]': ('  ',), 'Factor Value[AGE]': ('  ',)},
}
OPTIONAL_EMPTY_COLUMNS = {'E-MTAB-16192': frozenset({'Derived Array Data File'})}
NODE_HEADERS = frozenset({'Source Name', 'Sample Name', 'Extract Name', 'Labeled Extract Name',
    'Hybridization Name', 'Assay Name', 'Scan Name', 'Normalization Name',
    'Array Data File', 'Derived Array Data File', 'Array Data Matrix File',
    'Derived Array Data Matrix File', 'Image File', 'Protocol REF'})
QUALIFIER_HEADERS = frozenset({'Term Source REF', 'Term Accession Number', 'Performer', 'Date'})


def without_derived_files(table):
    """Independent expected export projection; original deposited bytes stay sealed.

    Remove file nodes and their owned metadata only. Explicit scientific fields
    (and all their qualifiers) survive. No values or row occurrences normalize.
    """
    indices, excluded = [], False
    for index, original in enumerate(table[0]):
        header = HEADER_ALIASES.get(original, original)
        if header in NODE_HEADERS:
            excluded = header in {'Derived Array Data File', 'Derived Array Data Matrix File'}
        elif header.startswith(('Factor Value[', 'Characteristics[', 'Parameter Value[')):
            excluded = False
        if not excluded and not header.startswith('Comment[SAMPLE_FILE_'):
            indices.append(index)
    return [[row[i] for i in indices] for row in table]


def load_deposited_reference(root, accession):
    """Read one source IDF/SDRF pair only after its recorded hashes match."""
    root = Path(root)
    if accession not in DEPOSITED_PANEL:
        raise ValueError(f'undeclared deposited accession: {accession}')
    provenance = json.loads((root/'provenance.json').read_text())
    entries = {item['file']: item for item in provenance['files']}
    result = {'accession': accession, 'sha256': {}}
    for kind in ('idf','sdrf'):
        name = f'{accession}.{kind}.txt'
        expected = entries[name]
        raw = (root/name).read_bytes()
        digest = hashlib.sha256(raw).hexdigest()
        if digest != expected['sha256'] or len(raw) != expected['bytes']:
            raise ValueError(f'deposited source hash/size mismatch: {name}')
        result[kind] = list(csv.reader(io.StringIO(raw.decode('utf-8')), delimiter='\t'))
        result['sha256'][kind] = digest
    return result


def _scoped_table(table, accession, *, normalize_blanks=True):
    assert isinstance(table, (tuple,list)) and table, 'missing SDRF table'
    headers = [HEADER_ALIASES.get(value,value) for value in table[0]]
    assert headers and all(isinstance(h,str) for h in headers), 'invalid SDRF headers'
    assert all(len(row)==len(headers) for row in table[1:]), 'nonrectangular SDRF'
    omitted = {i for i,h in enumerate(headers) if h in OPTIONAL_EMPTY_COLUMNS.get(accession,())
               and all(row[i]=='' for row in table[1:])
               and (i+1==len(headers) or headers[i+1] in NODE_HEADERS
                    or headers[i+1].startswith('Factor Value['))}
    indices = [i for i in range(len(headers)) if i not in omitted]
    descriptors, graph = [], []
    attributes, qualifiers = Counter(), Counter()
    node_index, owner_node, owner = -1, -1, None
    factors = Counter()
    for index in indices:
        header = headers[index]
        if header in NODE_HEADERS:
            graph.append(header)
            node_index += 1
            owner_node = node_index
            attributes, qualifiers = Counter(), Counter()
            owner = ('node',)
            descriptor = (node_index, owner, ('value',))
        elif header.startswith('Factor Value['):
            factors[header] += 1
            owner_node = -1  # Experimental variable of the complete SDRF row.
            owner = ('factor',header,factors[header])
            descriptor = (owner_node, owner, ('value',))
        elif header.startswith('Comment['):
            assert owner is not None, f'orphan comment {header}'
            qualifiers[(owner,header)] += 1
            descriptor = (owner_node, owner, (header,qualifiers[(owner,header)]))
        elif header in QUALIFIER_HEADERS:
            assert owner is not None, f'orphan qualifier {header}'
            qualifiers[(owner,header)] += 1
            descriptor = (owner_node, owner, (header,qualifiers[(owner,header)]))
        elif header == 'Unit' or header.startswith('Unit['):
            assert owner is not None, f'orphan unit {header}'
            qualifiers[(owner,header)] += 1
            owner = (*owner,'unit',header,qualifiers[(owner,header)])
            descriptor = (owner_node, owner, ('value',))
        else:
            assert node_index>=0, f'attribute before source node: {header}'
            attributes[header] += 1
            owner_node = node_index
            owner = ('attribute',header,attributes[header])
            descriptor = (node_index, owner, ('value',))
        descriptors.append(descriptor)
    assert len(set(descriptors))==len(descriptors), 'ambiguous column ownership'
    normalized = []
    for row in table[1:]:
        values = []
        for index in indices:
            value = row[index]
            if normalize_blanks and type(value) is str and value in BLANK_LITERALS.get(accession,{}).get(headers[index],()):
                value = ''
            # Tag types: numeric zero must not compare equal to False or '0'.
            values.append((type(value).__name__,value))
        normalized.append(tuple(values))
    return graph, descriptors, normalized, [headers[i] for i in indices]


def assert_deposited_sdrf(expected, actual, *, accession):
    """Assert exact source row multiplicity and every node/qualifier-owned cell.

    Factor blocks may move within the same complete row; this scientific check
    does not replace the separate MAGE-TAB terminal-factor layout validator.
    Attribute order within one owner may change. Node/protocol order, duplicate
    attributes, qualifier ownership and exact nonblank literals may not change.
    Unexpected additions are reported explicitly and require a reviewed rule.
    """
    egraph, ecolumns, erows, _ = _scoped_table(expected,accession)
    agraph, acolumns, arows, _ = _scoped_table(actual,accession)
    assert egraph==agraph, f'{accession}: node/protocol sequence changed: {egraph!r} != {agraph!r}'
    missing, added = set(ecolumns)-set(acolumns), set(acolumns)-set(ecolumns)
    assert not missing and not added, f'{accession}: missing scoped columns {missing!r}; added scoped columns {added!r}'
    positions = [acolumns.index(column) for column in ecolumns]
    expected_rows = Counter(erows)
    actual_rows = Counter(tuple(row[i] for i in positions) for row in arows)
    assert expected_rows==actual_rows, (f'{accession}: row/cell occurrences differ; '
        f'lost={sum((expected_rows-actual_rows).values())}, added={sum((actual_rows-expected_rows).values())}; '
        f'first lost={next(iter(expected_rows-actual_rows),None)!r}; '
        f'first added={next(iter(actual_rows-expected_rows),None)!r}')
    # Blank compatibility belongs to known source cells, not every output
    # occurrence of an otherwise permitted header/literal. Bind it to the full
    # scientific row, allowing a deposited blank to remain or become empty.
    _, _, eraw, eheaders = _scoped_table(expected,accession,normalize_blanks=False)
    _, _, araw, _ = _scoped_table(actual,accession,normalize_blanks=False)
    def retained_blanks(rows, raw_rows):
        return Counter((row,ecolumns[i],raw[i]) for row,raw in zip(rows,raw_rows)
            for i,header in enumerate(eheaders)
            if raw[i][0]=='str' and raw[i][1] in BLANK_LITERALS.get(accession,{}).get(header,()))
    expected_blanks = retained_blanks(erows,eraw)
    aligned_rows = [tuple(row[i] for i in positions) for row in arows]
    aligned_raw = [tuple(row[i] for i in positions) for row in araw]
    extra_blanks = retained_blanks(aligned_rows,aligned_raw)-expected_blanks
    assert not extra_blanks, f'{accession}: blank literals introduced outside their source row/column'
    protocol_columns = [i for i,(node,owner,leaf) in enumerate(ecolumns)
        if node>=0 and egraph[node]=='Protocol REF' and owner==('node',) and leaf[0] in {'Term Source REF','Term Accession Number'}]
    return {'accession':accession, 'row_count':len(erows), 'distinct_rows':len(expected_rows),
            'node_sequence':egraph, 'column_count':len(ecolumns),
            'protocol_qualifier_cells':sum(row[i] != ('str','') for row in erows for i in protocol_columns)}


def assert_deposited_protocol_qualifiers(expected, actual, *, accession):
    """Check protocol qualifiers against exact complete graph-row identity.

    This narrower check does not certify unrelated attribute/comment ownership.
    It remains useful when that separate strict comparison exposes a defect.
    """
    def occurrences(table):
        graph, columns, rows, _ = _scoped_table(table,accession)
        selected = [i for i,(node,owner,leaf) in enumerate(columns)
            if node>=0 and owner==('node',) and (leaf==('value',) or
               (graph[node]=='Protocol REF' and leaf[0] in {'Term Source REF','Term Accession Number'}))]
        qualifiers = [i for i in selected if columns[i][2]!=('value',)]
        return (graph, [columns[i] for i in selected],
                Counter(tuple(row[i] for i in selected) for row in rows),
                sum(row[i]!=('str','') for row in rows for i in qualifiers))
    wanted, observed = occurrences(expected), occurrences(actual)
    assert wanted==observed, f'{accession}: protocol qualifier/graph row occurrences differ'
    return {'accession':accession,'protocol_qualifier_cells':wanted[3]}


ORDERED_FIELD_PANEL = ["E-MTAB-13662", "E-MTAB-14560", "E-MTAB-14566", "E-MTAB-16192",
         "E-MTAB-16253", "E-MTAB-13131", "E-MTAB-16856", "E-MTAB-14416",
         "E-MTAB-16847", "E-MTAB-14960"]

PROTOCOL_QUALIFIER_COUNTS = {
    'E-GEOD-48034': 20, 'E-GEOD-48035': 844,
    'E-GEOD-52564': 34, 'E-GEOD-60450': 48,
}


def assert_ordered_archive_fields(expected, actual):
    """Preserve repeated columns, normalized row multiplicity and graph order."""
    normalize = lambda row: tuple("" if v is None else str(v).strip() for v in row)
    def columns(headers):
        counts = Counter()
        keys = []
        for h in headers:
            counts[h] += 1
            keys.append((h, counts[h]))
        return keys
    akeys, ekeys = columns(actual[0]), columns(expected[0])
    assert set(akeys) == set(ekeys)
    positions = [akeys.index(key) for key in ekeys]
    assert Counter(normalize([r[i] for i in positions]) for r in actual[1:]) == Counter(map(normalize, expected[1:]))
    # Attribute ordering within a node may normalize; graph node/edge order may not.
    graph = lambda headers: [h for h in headers if h.endswith(" Name") or h.endswith(" File") or h == "Protocol REF"]
    assert graph(actual[0]) == graph(expected[0])


def protocol_occurrences(table):
    headers = table[0]
    columns, owner, occurrence = [], None, 0
    for i, header in enumerate(headers):
        if header == 'Protocol REF':
            occurrence += 1
            owner = ('protocol', occurrence, i)
        elif header not in {'Term Source REF', 'Term Accession Number'}:
            owner = None
        if owner and header in {'Term Source REF', 'Term Accession Number'}:
            columns.append((owner, header, i))
    return Counter((owner[1], header, row[owner[2]], row[i])
                   for row in table[1:] for owner, header, i in columns if row[i])


def assert_protocol_occurrence_counts(source, output, expected):
    """Keep the independent application-occurrence census alongside row identity."""
    original = protocol_occurrences(source)
    assert sum(original.values()) == expected
    assert protocol_occurrences(output) == original
    assert len(output) == len(source)
