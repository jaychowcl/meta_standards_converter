# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Source-derived acceptance gaps found by the complete cohort replay."""
from copy import deepcopy
import pytest

from tests.magetab.test_metadata_integrity import SOURCE
from meta_standards_converter.miniml import MINiMLCodec
from meta_standards_converter.metadata.workflows import prepare_workflows


def package():
    source = SOURCE['fastq']
    file = {'uri': source['URI'], 'filename': source['NAME'], 'format': 'FASTQ',
            'md5': source['MD5'], 'bytes': source['BYTES']}
    return {'miniml_schema_version': '3.1', 'source': {'format': 'SRA'},
            'sample': [{'iid': 'ERS15958371', 'channel': [], 'sra_run': [
                {'run': 'ERR11634701', 'experiment': 'ERX11036790', 'files': [file]}]}],
            'series': {'iid': 'ERP148790', 'assay_paths': [{'steps': [
                {'kind': 'source', 'name': 'ERS15958371', 'sample_ref': 'ERS15958371'},
                {'kind': 'assay', 'name': 'ERX11036790'},
                {'kind': 'scan', 'name': 'ERR11634701'},
                {'kind': 'array_data_file', 'name': source['NAME'],
                 'link': {'value': source['URI'], 'type': 'fastq'},
                 'comments': [{'name': 'READ_TYPE', 'value': 'sample_barcode'}]}]}]}}


def file_comments(data):
    result = prepare_workflows(MINiMLCodec().decode(data).package).to_mapping()
    return result, {c['name']: c['value'] for c in result['series']['assay_paths'][0]['steps'][-1].get('comments', [])}


def test_exact_run_file_completes_authored_checksum_size_and_keeps_owner():
    data = package()
    before = deepcopy(data)
    result, comments = file_comments(data)
    assert comments.get('MD5') == SOURCE['fastq']['MD5']
    assert comments.get('BYTES') == SOURCE['fastq']['BYTES']
    assert comments['READ_TYPE'] == 'sample_barcode'
    assert data == before
    assert prepare_workflows(result).to_mapping() == result


@pytest.mark.parametrize('mutation', ['checksum', 'size', 'uri', 'run', 'name', 'ambiguous'])
def test_incompatible_or_ambiguous_run_file_cannot_complete(mutation):
    data = package()
    supplied = data['sample'][0]['sra_run'][0]['files'][0]
    node = data['series']['assay_paths'][0]['steps'][-1]
    if mutation == 'checksum':
        node['comments'].append({'name': 'MD5', 'value': '0' * 32})
    elif mutation == 'size':
        node['comments'].append({'name': 'BYTES', 'value': '1'})
    elif mutation == 'uri':
        node['link']['value'] = 'https://example.org/' + node['name']
    elif mutation == 'run':
        data['sample'][0]['sra_run'][0]['run'] = 'ERR99999999'
    elif mutation == 'name':
        node['name'] = 'different.fastq.gz'
    else:
        data['sample'][0]['sra_run'][0]['files'].append(dict(supplied, md5='0' * 32))
    _, comments = file_comments(data)
    assert comments == {c['name']: c['value'] for c in node['comments']}


def test_provider_read_roles_reach_the_canonical_file_occurrence():
    data = package()
    file = data['sample'][0]['sra_run'][0]['files'][0]
    file.update(read_role='unpaired')
    node = data['series']['assay_paths'][0]['steps'][-1]
    node['comments'] = []
    _, comments = file_comments(data)
    assert comments.get('READ_ROLE') == 'unpaired'
    assert not comments.get('READ_INDEX')
