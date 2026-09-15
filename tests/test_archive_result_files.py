# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from copy import deepcopy

import pytest

from meta_standards_converter.magetab.native_files import _record, _checksums
from meta_standards_converter.magetab.semantics import _miniml_path_columns
from meta_standards_converter.miniml.insdc_support import result_file
from meta_standards_converter.miniml.archive_results import normalize_result_bundles


@pytest.mark.parametrize('algorithm,key,digest', [('MD5', 'checksum', 'a'*32), ('SHA256', 'file_checksum', 'b'*64)])
def test_structured_checksums_are_exported_and_available_for_file_comparison(algorithm, key, digest):
    node = {'kind': 'derived_array_data_file', 'name': 'aligned.cram',
            'link': {'value': 'https://example.org/aligned.cram', 'type': 'CRAM',
                     key: digest, 'checksum_method': algorithm, 'bytes': '123', 'role': 'SUBMISSION_FILE'}}
    record = _record(node)
    assert _checksums(record) == {algorithm: digest}
    assert record['BYTES'] == '123'
    columns = _miniml_path_columns([node])
    assert ('Comment[CHECKSUM]', digest) in columns
    assert ('Comment[CHECKSUM_METHOD]', algorithm) in columns
    assert ('Comment[File URI]', node['link']['value']) in columns


def acquisition():
    sample = {'iid': 'SRS1', 'channel': [], 'sra_run': [{'run': 'SRR1', 'experiment': 'SRX1'}]}
    path = {'steps': [{'kind': 'source', 'name': 'SRS1', 'sample_ref': 'SRS1'},
                      {'kind': 'assay', 'name': 'SRX1', 'sample_ref': 'SRS1'},
                      {'kind': 'scan', 'name': 'SRR1', 'sample_ref': 'SRS1'},
                      {'kind': 'array_data_file', 'name': 'read.fastq'}]}
    return sample, [path]


def test_projection_of_later_files_cannot_inherit_prior_analysis_applications():
    sample, paths = acquisition()
    for name, method in [('a.cram','ERZ1:analysis'), ('a.cram.crai','ERZ1:analysis'), ('b.cram','ERZ2:analysis')]:
        result_file(sample, {'uri': 'https://example.org/'+name, 'format': 'CRAI' if name.endswith('crai') else 'CRAM'},
                    paths, ['SRR1'], {'name': method})
    assert [[s['protocol_ref'] for s in p['steps'] if s.get('protocol_ref')] for p in paths[1:]] == [
        ['ERZ1:analysis'], ['ERZ1:analysis'], ['ERZ2:analysis']]


def companion_data():
    sample, paths = acquisition()
    for name, fmt, role in [('a.cram','CRAM','SUBMISSION_FILE'), ('a.cram.crai','CRAI','INDEX_FILE')]:
        result_file(sample, {'uri': 'https://example.org/'+name, 'format': fmt, 'role': role,
                            'analysis_accession': 'ERZ1', 'md5': 'a'*32}, paths, ['SRR1'])
    return {'source': {'format': 'ENA'}, 'series': {'iid': 'ERP1', 'protocols': [], 'assay_paths': paths}, 'sample': [sample]}


def test_verified_index_is_a_companion_not_a_second_result():
    data = companion_data()
    normalize_result_bundles(data)
    results = [s for p in data['series']['assay_paths'] for s in p['steps'] if s['kind'].startswith('derived_')]
    assert len(results) == 1
    assert results[0]['name'] == 'a.cram'
    member = results[0]['link']['companion_files'][0]
    assert member['role'] == 'index' and member['node']['name'] == 'a.cram.crai'
    assert ('Comment[INDEX_URI]', 'https://example.org/a.cram.crai') in _miniml_path_columns(results)
    before = deepcopy(data)
    normalize_result_bundles(data)
    assert data == before


@pytest.mark.parametrize('change', ['analysis', 'run', 'missing analysis', 'multiple crams'])
def test_ambiguous_or_unrelated_indexes_are_not_paired(change):
    data = companion_data(); paths = data['series']['assay_paths']
    if change == 'analysis': paths[-1]['steps'][-1]['link']['analysis_accession'] = 'ERZ2'
    elif change == 'run': next(s for s in paths[-1]['steps'] if s['kind'] == 'scan')['name'] = 'SRR2'
    elif change == 'missing analysis':
        for path in paths:
            path['steps'][-1].get('link', {}).pop('analysis_accession', None)
    else:
        other = deepcopy(paths[1]); other['steps'][-1]['name'] = 'b.cram'; paths.append(other)
    normalize_result_bundles(data)
    assert not any(s.get('link', {}).get('companion_files') for p in paths for s in p['steps'])


def test_structured_and_comment_checksum_conflict_cannot_establish_file_identity():
    node = {'kind':'derived_array_data_file', 'name':'x', 'link':{'checksum':'a'*32},
            'comments':[{'name':'MD5','value':'b'*32}]}
    assert not _checksums(_record(node))
    columns = _miniml_path_columns([node])
    assert ('Comment[CHECKSUM]', 'a'*32) in columns and ('Comment[MD5]', 'b'*32) in columns


def test_two_explicit_checksum_algorithms_remain_associated_and_exported():
    node = {'kind':'derived_array_data_file', 'name':'x',
            'link':{'checksum':'a'*32, 'file_checksum':'b'*64, 'checksum_method':'SHA256'}}
    columns = _miniml_path_columns([node])
    assert ('Comment[MD5]', 'a'*32) in columns
    assert ('Comment[CHECKSUM]', 'b'*64) in columns
    assert ('Comment[CHECKSUM_METHOD]', 'SHA256') in columns
    node['link'].update(file_checksum='b'*32, checksum_method='MD5')
    columns = _miniml_path_columns([node])
    assert ('Comment[MD5]', 'a'*32) in columns
    assert ('Comment[CHECKSUM]', 'b'*32) in columns
    assert not _checksums(_record(node))


def test_index_comparison_restores_parallel_branches_for_occurrence_pruning():
    from meta_standards_converter.miniml.archive_results import comparison_view
    data = companion_data(); original = deepcopy(data)
    normalize_result_bundles(data)
    view = comparison_view(data)
    assert view['series']['assay_paths'] == original['series']['assay_paths']
    assert len(data['series']['assay_paths']) == 2


def test_ambiguous_companion_reports_without_modifying_paths(caplog):
    from meta_standards_converter.miniml.archive_results import normalize_index_companions
    data = companion_data()
    other = deepcopy(data['series']['assay_paths'][1]); other['steps'][-1]['name'] = 'b.cram'
    data['series']['assay_paths'].append(other)
    before = deepcopy(data); issues = []
    normalize_index_companions(data, issues=issues, check_only=True)
    assert issues and 'ambiguous' in issues[0] and 'ERZ1' in issues[0]
    assert data == before


def test_saved_explicit_filename_hash_is_encoded_in_export_without_mutation():
    node = {'kind':'array_data_file','name':'18858_6#91_1.fastq.gz',
            'link':{'value':'https://files.example.org/18858_6#91_1.fastq.gz','type':'fastq'}}
    original = deepcopy(node)
    assert _record(node)['URI'] == 'https://files.example.org/18858_6%2391_1.fastq.gz'
    assert ('Comment[File URI]', 'https://files.example.org/18858_6%2391_1.fastq.gz') in _miniml_path_columns([node])
    assert node == original
