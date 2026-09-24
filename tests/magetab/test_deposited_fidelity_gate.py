# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Independent deposited-source acceptance with row and qualifier ownership."""
from copy import deepcopy
from pathlib import Path

import pytest

from tests.support.deposited_fidelity import (
    DEPOSITED_PANEL, assert_deposited_sdrf, assert_deposited_protocol_qualifiers,
    load_deposited_reference,
)
from tests.support.recorded_replay import ReplayTransport


FIXTURES = Path(__file__).parents[1] / 'fixtures/fidelity'
def scoped_table():
    return [
        ['Source Name', 'Characteristics[dose]', 'Unit', 'Term Source REF',
         'Protocol REF', 'Term Source REF', 'Extract Name', 'Comment[bytes]',
         'Protocol REF', 'Term Source REF', 'Assay Name', 'Factor Value[dose]'],
        ['s1', '0', 'hour', 'UO', 'P1', 'ArrayExpress', 'e1', 0, 'P2', 'EFO', 'a1', ' 0 '],
        ['s2', '1', 'hour', 'UO', 'P1', 'EFO', 'e2', 0, 'P2', 'ArrayExpress', 'a2', ' 1 '],
    ]


@pytest.mark.parametrize('mutation', ['drop-row', 'duplicate-row', 'qualifier-row-swap',
    'qualifier-protocol-swap', 'unit-owner', 'move-owner', 'zero-blank', 'zero-type', 'literal-whitespace',
    'drop-column', 'extra-column', 'unknown-alias'])
def test_comparator_rejects_scientific_and_scope_changes(mutation):
    source = scoped_table()
    actual = deepcopy(source)
    if mutation == 'drop-row': actual.pop()
    elif mutation == 'duplicate-row': actual.append(actual[1][:])
    elif mutation == 'qualifier-row-swap': actual[1][5], actual[2][5] = actual[2][5], actual[1][5]
    elif mutation == 'qualifier-protocol-swap': actual[1][5], actual[1][9] = actual[1][9], actual[1][5]
    elif mutation == 'unit-owner':
        actual = [[row[i] for i in (0,1,2,4,3,5,6,7,8,9,10,11)] for row in actual]
    elif mutation == 'move-owner':
        for row in actual: row[6], row[7] = row[7], row[6]
    elif mutation == 'zero-blank': actual[1][7] = ''
    elif mutation == 'zero-type': actual[1][7] = '0'
    elif mutation == 'literal-whitespace': actual[1][-1] = '0'
    elif mutation == 'drop-column':
        for row in actual: row.pop(7)
    elif mutation == 'extra-column':
        actual[0].append('Comment[extra]')
        for row in actual[1:]: row.append('')
    elif mutation == 'unknown-alias': actual[0][1] = 'Characteristics [dose]'
    with pytest.raises(AssertionError): assert_deposited_sdrf(source, actual, accession='synthetic')


def test_exact_duplicate_occurrences_and_attribute_reordering_are_retained():
    source = [['Source Name','Comment[a]','Comment[b]','Comment[a]'],
              ['s1',' 0 ',0,'second'],['s1',' 0 ',0,'second']]
    actual = [[row[i] for i in (0,2,1,3)] for row in source]
    report = assert_deposited_sdrf(source, actual, accession='synthetic')
    assert report['row_count'] == 2
    assert report['distinct_rows'] == 1


def row_factor_table():
    return [
        ['Source Name', 'Assay Name', 'Array Data File', 'Comment[file]',
         'Factor Value[dose]', 'Term Source REF', 'Unit[time]', 'Term Source REF',
         'Comment[unit note]', 'Factor Value[dose]', 'Term Source REF'],
        ['s1', 'a1', 'f1', 'file1', '0', 'EFO', 'hour', 'UO', 'note1', 'second1', 'EFO2'],
        ['s2', 'a2', 'f2', 'file2', '1', 'EFO2', 'day', 'UO2', 'note2', 'second2', 'EFO'],
    ]


def test_factor_occurrences_and_qualifiers_have_whole_row_scope():
    source = row_factor_table()
    # Factor placement does not change its experimental row association.
    actual = [[row[i] for i in (0,1,4,5,6,7,8,9,10,2,3)] for row in source]
    assert_deposited_sdrf(source, actual, accession='synthetic')


def test_comment_cannot_move_between_a_characteristic_and_its_node():
    source = [['Source Name','Comment[detail]','Characteristics[dose]'],['s1','zero dose','0']]
    actual = [[row[i] for i in (0,2,1)] for row in source]
    with pytest.raises(AssertionError):
        assert_deposited_sdrf(source,actual,accession='synthetic')


@pytest.mark.parametrize('header',['Array Data Matrix File','Derived Array Data Matrix File'])
def test_matrix_files_are_graph_nodes_not_reorderable_attributes(header):
    source=[['Source Name',header,'Comment[file]','Assay Name'],['s1','matrix','note','a1']]
    actual=[[row[i] for i in (0,3,1,2)] for row in source]
    with pytest.raises(AssertionError,match='sequence'):
        assert_deposited_sdrf(source,actual,accession='synthetic')


@pytest.mark.parametrize('mutation', ['factor-row-swap', 'assay-swap', 'file-swap',
    'duplicate-factor-swap', 'factor-qualifier-swap', 'factor-unit-owner', 'factor-comment-owner'])
def test_row_factors_cannot_cross_experiments_or_qualifier_owners(mutation):
    source = row_factor_table()
    actual = deepcopy(source)
    if mutation == 'factor-row-swap':
        actual[1][4:], actual[2][4:] = actual[2][4:], actual[1][4:]
    elif mutation in {'assay-swap', 'file-swap'}:
        index = 1 if mutation == 'assay-swap' else 2
        actual[1][index], actual[2][index] = actual[2][index], actual[1][index]
    elif mutation == 'duplicate-factor-swap':
        actual[1][4], actual[1][9] = actual[1][9], actual[1][4]
    elif mutation == 'factor-qualifier-swap':
        actual[1][5], actual[1][10] = actual[1][10], actual[1][5]
    elif mutation == 'factor-unit-owner':
        actual = [[row[i] for i in (0,1,2,3,4,5,9,10,6,7,8)] for row in actual]
    elif mutation == 'factor-comment-owner':
        actual = [[row[i] for i in (0,1,2,3,4,5,6,7,9,10,8)] for row in actual]
    with pytest.raises(AssertionError):
        assert_deposited_sdrf(source, actual, accession='synthetic')


def test_known_blank_compatibility_is_accession_header_and_literal_scoped():
    source = [['Source Name','Characteristics [reference rna 2]'],['s1','  ']]
    actual = [['Source Name','Characteristics[reference rna 2]'],['s1','']]
    assert_deposited_sdrf(source, actual, accession='E-GEOD-48035')
    with pytest.raises(AssertionError): assert_deposited_sdrf(source, actual, accession='E-MTAB-14416')
    source[1][1] = ' '  # Only the observed two-space literal is permitted.
    with pytest.raises(AssertionError): assert_deposited_sdrf(source, actual, accession='E-GEOD-48035')



def test_blank_exceptions_cannot_introduce_whitespace_into_another_source_row():
    source = [['Source Name','Characteristics [reference rna 2]'],['s1','  '],['s2','']]
    actual = [['Source Name','Characteristics[reference rna 2]'],['s1',''],['s2','  ']]
    with pytest.raises(AssertionError): assert_deposited_sdrf(source,actual,accession='E-GEOD-48035')


def test_only_entirely_empty_optional_file_column_may_disappear():
    source = [['Source Name','Derived Array Data File','Factor Value[dose]'],['s1','','0']]
    actual = [['Source Name','Factor Value[dose]'],['s1','0']]
    assert_deposited_sdrf(source, actual, accession='E-MTAB-16192')
    source[1][1] = 'file.txt'
    with pytest.raises(AssertionError): assert_deposited_sdrf(source, actual, accession='E-MTAB-16192')


def test_empty_file_with_owned_annotation_cannot_disappear():
    source=[['Source Name','Derived Array Data File','Comment[file]'],['s1','','retained annotation']]
    actual=[['Source Name','Comment[file]'],['s1','retained annotation']]
    with pytest.raises(AssertionError): assert_deposited_sdrf(source,actual,accession='E-MTAB-16192')


def test_reference_loader_rejects_bytes_not_matching_recorded_source(tmp_path):
    import shutil
    for name in ('provenance.json','E-MTAB-16847.idf.txt','E-MTAB-16847.sdrf.txt'):
        shutil.copyfile(FIXTURES/name,tmp_path/name)
    load_deposited_reference(tmp_path,'E-MTAB-16847')
    with (tmp_path/'E-MTAB-16847.sdrf.txt').open('a') as stream: stream.write('\n')
    with pytest.raises(ValueError, match='hash'): load_deposited_reference(tmp_path,'E-MTAB-16847')


@pytest.mark.parametrize('accession', DEPOSITED_PANEL)
def test_full_deposited_panel_preserves_scoped_source_occurrences(accession):
    reference = load_deposited_reference(FIXTURES, accession)
    isolation = ReplayTransport(Path(__file__).parents[1]/'fixtures/recorded_provider/success.manifest.json')
    with isolation.offline():
        from meta_standards_converter.converters import Converter
        result = Converter().convert(FIXTURES/(accession+'.idf.txt'), out_type='magetab',
            enrichment='off', options={'expand_studies':False})
    isolation.assert_complete()
    assert len(result.items) == 1
    item = result.items[0]
    assert item.execution == 'succeeded' and item.validation == 'valid', result.to_dict()
    assert item.completeness == ('partial' if accession=='E-GEOD-48035' else 'complete'), result.to_dict()
    assert all('source-text encoding' in d.message for d in item.diagnostics if d.severity in {'warning','error'})
    assert len(item.payload) == 1
    tables = item.payload[0]
    outputs = [row[1] for row in tables if row[0]=='SDRF File']
    assert len(outputs) == 1
    from tests.support.deposited_fidelity import NODE_HEADERS
    header = outputs[0][0]
    last_node = max(i for i,value in enumerate(header) if value in NODE_HEADERS)
    assert all(i > last_node for i,value in enumerate(header) if value.startswith('Factor Value['))
    expected = {'E-GEOD-48034':20,'E-GEOD-48035':844,'E-GEOD-52564':34,'E-GEOD-60450':48}
    if accession in expected:
        report=assert_deposited_protocol_qualifiers(reference['sdrf'],outputs[0],accession=accession)
        assert report['protocol_qualifier_cells'] == expected[accession]
    assert_deposited_sdrf(reference['sdrf'], outputs[0], accession=accession)


@pytest.mark.parametrize('index', [5,9])
def test_protocol_only_gate_retains_row_and_application_scope(index):
    source=scoped_table()
    actual=deepcopy(source)
    actual[1][index],actual[2][index]=actual[2][index],actual[1][index]
    with pytest.raises(AssertionError):
        assert_deposited_protocol_qualifiers(source,actual,accession='synthetic')


@pytest.mark.parametrize('accession,count', [('E-GEOD-48034',20),('E-GEOD-48035',844),
                                           ('E-GEOD-52564',34),('E-GEOD-60450',48)])
def test_actual_exports_retain_946_protocol_qualifiers_with_graph_row_identity(accession,count):
    reference=load_deposited_reference(FIXTURES,accession)
    isolation=ReplayTransport(Path(__file__).parents[1]/'fixtures/recorded_provider/success.manifest.json')
    with isolation.offline():
        from meta_standards_converter.converters import Converter
        result=Converter().convert(FIXTURES/(accession+'.idf.txt'),out_type='magetab',
            enrichment='off',options={'expand_studies':False})
    isolation.assert_complete()
    assert len(result.items)==1 and len(result.items[0].payload)==1
    tables=[row[1] for row in result.items[0].payload[0] if row[0]=='SDRF File']
    assert len(tables)==1
    report=assert_deposited_protocol_qualifiers(reference['sdrf'],tables[0],accession=accession)
    assert report['protocol_qualifier_cells']==count


def test_factor_comment_cannot_become_unit_comment():
    source=[['Source Name','Factor Value[dose]','Comment[note]','Unit[time]'],['s1','0','factor note','hour']]
    actual=[[row[i] for i in (0,1,3,2)] for row in source]
    with pytest.raises(AssertionError): assert_deposited_sdrf(source,actual,accession='synthetic')


def test_deposited_cohort_and_protocol_totals_are_explicit():
    assert len(DEPOSITED_PANEL) == len(set(DEPOSITED_PANEL)) == 18
    reports = []
    for accession in ('E-GEOD-48034','E-GEOD-48035','E-GEOD-52564','E-GEOD-60450'):
        source = load_deposited_reference(FIXTURES,accession)['sdrf']
        reports.append(assert_deposited_sdrf(source,source,accession=accession))
    assert sum(report['protocol_qualifier_cells'] for report in reports) == 946
