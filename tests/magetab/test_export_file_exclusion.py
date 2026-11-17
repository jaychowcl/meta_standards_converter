"""Processed-file exclusion is a projection, never a JSON mutation."""
from copy import deepcopy
import pytest
from tests.test_geo_source_fidelity import geo_package, rows, row


@pytest.mark.parametrize('handler', ['sequencing', 'bulk_sequencing', 'array', 'generic'])
def test_fallback_derived_files_are_excluded_and_json_is_unchanged(handler):
    data = geo_package(); before = deepcopy(data)
    output = rows(data, handler)
    assert data == before
    assert 'processed1.txt' not in str(output)
    assert 'Derived Array Data' not in str(output)


@pytest.mark.parametrize('kind', ['derived_array_data_file', 'derived_array_data_matrix_file'])
def test_canonical_exclusion_preserves_rows_documents_and_independent_factors(kind):
    data = geo_package()
    data['series']['assay_paths'] = [
        {'document': f'doc{i}.sdrf.txt', 'steps': [
            {'kind': 'source', 'name': 'GSM1', 'sample_ref': 'GSM1'},
            {'kind': 'assay', 'name': 'assay1'},
            {'kind': 'array_data_file', 'name': 'raw.fastq', 'comments': [{'name': 'MD5', 'value': 'raw-hash'}]},
            {'kind': 'image_file', 'name': 'slide.png'},
            {'kind': kind, 'name': f'processed{i}.txt', 'comments': [{'name': 'MD5', 'value': 'derived-hash'}],
             'factor_values': [{'name': 'dose', 'value': '1', 'unit': {'value': 'mg'},
                                'comments': [{'name': 'note', 'value': 'authored factor'}]}]}
        ]} for i in (1, 2)]
    before = deepcopy(data); output = rows(data); table = row(output, 'SDRF File')[1]
    assert data == before and len(table) == 3
    assert 'Derived Array Data' not in str(output) and 'derived-hash' not in str(output)
    assert table[1][table[0].index('Array Data File')] == 'raw.fastq'
    assert table[1][table[0].index('Image File')] == 'slide.png'
    assert table[1][table[0].index('Comment[MD5]')] == 'raw-hash'
    assert table[1][table[0].index('Factor Value[dose]')] == '1'
    assert 'authored factor' in table[1] and 'mg' in table[1]


def test_native_sample_results_do_not_become_file_comments():
    data = geo_package(); data['source']['format'] = 'SRA'
    data['series']['assay_paths'] = [{'steps': [
        {'kind': 'source', 'name': 'GSM1', 'sample_ref': 'GSM1'},
        {'kind': 'derived_array_data_file', 'name': 'analysis.txt',
         'factor_values': [{'name': 'dose', 'value': '2'}]}]}]
    output = rows(data)
    assert 'analysis.txt' not in str(output)
    assert 'SAMPLE_FILE_' not in str(output)
    assert row(output, 'SDRF File')[1][1][-1] == '2'


def test_raw_inventory_and_source_prose_are_preserved():
    from meta_standards_converter.magetab.constructor import AEConstructor
    from meta_standards_converter.miniml import MINiMLCodec
    data = geo_package()
    data['sample'][0]['raw_data'] = [{'value': 'https://example.test/raw.sra'}]
    data['series']['summary'] = 'The authors processed processed1.txt.'
    plan = AEConstructor().prepare(MINiMLCodec().decode(data).package)
    assert plan.data['sample'][0]['raw_data'] == data['sample'][0]['raw_data']
    assert 'processed1.txt' in row(AEConstructor().render(plan), 'Experiment Description')[1]
