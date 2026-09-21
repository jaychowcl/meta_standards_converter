# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
import copy

import pytest

from meta_standards_converter.magetab.constructor import AEConstructor
from meta_standards_converter.magetab.parser import AEParser
from tests.converters.test_ae2json import resolved_input


VALUE = 'Comment[SecondaryAccession]'
SOURCE = 'Comment[SecondaryAccessionTermSourceRef]'


@pytest.mark.parametrize('comments', [
    f'{VALUE}\tPRJNA229389\tSRP033200\tGSE52564\n{SOURCE}\tBioProject\tSRA\tGEO\n',
    f'{VALUE}\tPRJNA229389\n{SOURCE}\tBioProject\n{VALUE}\tSRP033200\n{SOURCE}\tSRA\n{VALUE}\tGSE52564\n{SOURCE}\tGEO\n',
])
def test_secondary_accessions_survive_semantic_overlay_and_roundtrip(comments):
    package = AEParser().parse(resolved_input(idf='Investigation Accession\tE-GEOD-52564\nExperimental Factor Name\tdisease\nProtocol Name\tP-extract\n' + comments))
    before = copy.deepcopy(package.to_mapping())
    rows = AEConstructor().miniml2magetab(package)
    pairs = []
    for i, row in enumerate(rows):
        if row[0] == VALUE:
            assert len(row) == 2
            assert rows[i + 1][0] == SOURCE
            assert len(rows[i + 1]) == 2
            pairs.append((row[1], rows[i + 1][1]))
    # Core accession order is the order emitted, including the parser's GEO-first identity rule.
    assert pairs == [('GSE52564', 'GEO'), ('PRJNA229389', 'BioProject'), ('SRP033200', 'SRA')]
    assert not any(value.startswith('E-') for value, _ in pairs)
    comments_started = False
    for row in rows:
        if row[0].startswith('Comment['):
            comments_started = True
        else:
            assert not comments_started
    text = '\n'.join('\t'.join(str(v) for v in r) for r in rows if r[0] != 'SDRF File')
    recovered = AEParser().parse(resolved_input(idf=text)).to_mapping()
    assert recovered['series']['accession'] == before['series']['accession']
    assert package.to_mapping() == before


@pytest.mark.parametrize('comments', [
    f'{VALUE}\tCUSTOM-A\t\tCUSTOM-B\tCUSTOM-C\n{SOURCE}\tSourceA\tIGNORE\t\tSourceC\n',
    f'{VALUE}\tCUSTOM-A\n{SOURCE}\tSourceA\n{VALUE}\tCUSTOM-B\n{VALUE}\tCUSTOM-C\n{SOURCE}\tSourceC\n',
    f'{VALUE}\tCUSTOM-A\n{SOURCE}\tSourceA\n{VALUE}\t\n{SOURCE}\tIGNORE\n{VALUE}\tCUSTOM-B\n{SOURCE}\t\n{VALUE}\tCUSTOM-C\n{SOURCE}\tSourceC\n',
])
def test_blank_accessions_and_missing_sources_never_shift_other_associations(comments):
    data = AEParser().parse(resolved_input(idf='Investigation Accession\tE-MTAB-1\n' + comments)).to_mapping()
    found = {a['value']: a.get('database') for a in data['series']['accession']}
    assert found['CUSTOM-A'] == 'SourceA'
    assert found['CUSTOM-C'] == 'SourceC'
    assert found['CUSTOM-B'] not in {'SourceA', 'SourceC', 'IGNORE'}


def test_final_accession_partition_is_idempotent_and_preserves_conflicting_labels():
    from meta_standards_converter.magetab.accession_rows import expand_secondary_accession_rows
    rows = [[VALUE, 'CUSTOM', 'CUSTOM'], [SOURCE, 'A', 'B']]
    expanded = [[VALUE, 'CUSTOM'], [SOURCE, 'A'], [VALUE, 'CUSTOM'], [SOURCE, 'B']]
    assert expand_secondary_accession_rows(rows) == expanded
    assert expand_secondary_accession_rows(expanded) == expanded


@pytest.mark.parametrize('sources', [('SourceA', 'SourceB'), ('SourceB', 'SourceA')])
@pytest.mark.parametrize('accession', ['CUSTOM', 'SRP123'])
def test_conflicting_sources_survive_complete_parse_export(sources, accession):
    comments = ''.join(f'{VALUE}\t{accession}\n{SOURCE}\t{source}\n' for source in sources)
    data = AEParser().parse(resolved_input(idf='Investigation Accession\tE-MTAB-1\nExperimental Factor Name\tdisease\nProtocol Name\tP-extract\n' + comments))
    pairs = [(a['value'], a.get('database')) for a in data.to_mapping()['series']['accession'] if a['value'] == accession]
    assert pairs == [(accession, source) for source in sources]
    from meta_standards_converter.magetab.accession_rows import secondary_accession_pairs
    assert secondary_accession_pairs(AEConstructor().miniml2magetab(data)) == pairs
