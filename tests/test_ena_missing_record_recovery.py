# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Recover exact ENA identities without confusing membership with field coverage."""
from copy import deepcopy
from types import SimpleNamespace
from unittest.mock import Mock
import xml.etree.ElementTree as ET

import pytest
import requests

from meta_standards_converter.helpers.request_helper import RateLimitedRequester, RequestSettings
from meta_standards_converter.sources.archive_support import StudyRecords, StudySeed, ArchiveHTTP
from meta_standards_converter.sources.ena import ENASource
from meta_standards_converter.miniml.ena_parser import ENAParser
from tests.test_native_archive_sources import HTTP


def records():
    return StudyRecords(StudySeed('ERP1', 'PRJEB1'))


def xml_samples(ids):
    root = ET.Element('SAMPLE_SET')
    for acc in ids:
        node = ET.SubElement(root, 'SAMPLE', accession=acc)
        ET.SubElement(ET.SubElement(node, 'IDENTIFIERS'), 'PRIMARY_ID').text = acc
    return root


def test_truncated_batches_split_and_successful_children_are_not_downloaded_again():
    calls = []
    def fetch(url, params, fmt):
        ids = url.rsplit('/', 1)[-1].split(','); calls.append(ids)
        if len(ids) > 25:
            raise requests.exceptions.ChunkedEncodingError('truncated')
        return xml_samples(ids)
    target = records()
    root = ENASource(http=HTTP(fetch)).verified_xml([f'SAMEA{i}' for i in range(100)], target, 'sample')
    assert len(root) == 100 and not target.issues
    assert [len(c) for c in calls] == [100, 50, 25, 25, 50, 25, 25]


def test_incomplete_well_formed_batch_retries_only_missing_ids():
    calls = []
    def fetch(url, params, fmt):
        ids = url.rsplit('/', 1)[-1].split(','); calls.append(ids)
        return xml_samples(ids[:1] if len(calls) == 1 else ids)
    target = records()
    root = ENASource(http=HTTP(fetch)).verified_xml(['SAMEA1', 'SAMEA2'], target, 'sample')
    assert [n.get('accession') for n in root] == ['SAMEA1', 'SAMEA2']
    assert calls == [['SAMEA1', 'SAMEA2'], ['SAMEA2']]
    assert not target.issues


def test_recovery_is_bounded_and_rejects_unrequested_records():
    calls = []
    def fetch(url, params, fmt):
        calls.append(url)
        return xml_samples(['SAMEA999'])
    target = records()
    assert ENASource(http=HTTP(fetch)).verified_xml(['SAMEA1'], target, 'sample') is None
    assert len(calls) <= 8
    assert any('SAMEA1' in issue for issue in target.issues)


def test_verified_batches_preserve_duplicate_source_occurrences():
    value=xml_samples(['SAMEA1','SAMEA1'])
    value[0].set('alias','first');value[1].set('alias','second')
    root=ENASource(http=HTTP(lambda *args:value)).verified_xml(['SAMEA1'],records(),'sample')
    assert [node.get('alias') for node in root]==['first','second']


def test_study_accepts_explicit_bioproject_alias():
    value=ET.fromstring('<STUDY_SET><STUDY accession="ERP1"><IDENTIFIERS><PRIMARY_ID>ERP1</PRIMARY_ID><EXTERNAL_ID namespace="BioProject">PRJEB1</EXTERNAL_ID></IDENTIFIERS></STUDY></STUDY_SET>')
    root=ENASource(http=HTTP(lambda *args:value)).verified_xml(['PRJEB1'],records(),'project')
    assert root is not None and root[0].get('accession')=='ERP1'


def test_matching_identity_in_wrong_entity_type_is_rejected():
    value=ET.fromstring('<EXPERIMENT_SET><EXPERIMENT accession="SAMEA1"/></EXPERIMENT_SET>')
    target=records()
    assert ENASource(http=HTTP(lambda *args:value)).verified_xml(['SAMEA1'],target,'sample') is None
    assert any('missing XML record SAMEA1' in issue for issue in target.issues)


def test_physical_batch_limit_includes_connection_retries_and_redirects():
    def get(url, **kwargs):
        if url.endswith('/redirect'):
            raise requests.ReadTimeout('stalled')
        return Mock(status_code=302, headers={'Location': 'https://www.ebi.ac.uk/redirect'})
    transport = Mock(side_effect=get)
    requester = RateLimitedRequester('ena_portal', RequestSettings(request_delay=0),
        get=transport, host_gate=Mock(wait=Mock(return_value=0)), random_value=lambda: 0)
    target = records()
    assert ENASource(http=ArchiveHTTP('ena_portal', requester=requester)).verified_xml(
        ['SAMEA1'], target, 'sample') is None
    assert transport.call_count <= 8
    assert transport.call_count == requester.provider_attempts


def test_unchanged_small_batch_does_not_multiply_transport_retries():
    from collections import Counter
    from tests.test_metadata_stream_recovery import response,client
    calls=[]
    def get(url,**kwargs):
        calls.append(url)
        return response([requests.exceptions.ChunkedEncodingError('truncated')])
    target=records()
    source=ENASource(http=ArchiveHTTP('ena_portal',requester=client(get)))
    source.verified_xml([f'SAMEA{i}' for i in range(100)],target,'sample')
    assert len(calls)<=8
    assert max(Counter(calls).values())<=4


def test_pre_header_timeout_does_not_fan_out_batches():
    from tests.test_metadata_stream_recovery import client
    get=Mock(side_effect=requests.ReadTimeout('headers unavailable'))
    source=ENASource(http=ArchiveHTTP('ena_portal',requester=client(get)))
    source.verified_xml([f'SAMEA{i}' for i in range(100)],records(),'sample')
    assert get.call_count==4
    assert len({call.args[0] for call in get.call_args_list})==1


def missing_sample_records():
    value = records()
    value.xml = [ET.fromstring('''<STUDY_SET><STUDY accession="ERP1"><DESCRIPTOR><STUDY_TITLE>Recovery fixture</STUDY_TITLE></DESCRIPTOR></STUDY></STUDY_SET>'''),
        ET.fromstring('''<EXPERIMENT_SET><EXPERIMENT accession="ERX1"><STUDY_REF accession="ERP1"/><DESIGN><SAMPLE_DESCRIPTOR accession="ERS1"/><LIBRARY_DESCRIPTOR><LIBRARY_STRATEGY>RNA-Seq</LIBRARY_STRATEGY><LIBRARY_SOURCE>TRANSCRIPTOMIC</LIBRARY_SOURCE><LIBRARY_SELECTION>cDNA</LIBRARY_SELECTION><LIBRARY_LAYOUT><SINGLE/></LIBRARY_LAYOUT></LIBRARY_DESCRIPTOR></DESIGN><PLATFORM><ILLUMINA><INSTRUMENT_MODEL>Illumina HiSeq 2000</INSTRUMENT_MODEL></ILLUMINA></PLATFORM></EXPERIMENT></EXPERIMENT_SET>'''),
        ET.fromstring('<RUN_SET><RUN accession="ERR1"><EXPERIMENT_REF accession="ERX1"/></RUN></RUN_SET>')]
    row = {'study_accession': 'PRJEB1', 'secondary_study_accession': 'ERP1',
        'run_accession': 'ERR1', 'experiment_accession': 'ERX1', 'sample_accession': 'SAMEA1',
        'secondary_sample_accession': 'ERS1', 'sample_title': 'sample one',
        'scientific_name': 'Homo sapiens', 'tax_id': '9606',
        'fastq_ftp': 'ftp.sra.ebi.ac.uk/one.fastq.gz', 'fastq_md5': 'a'*32, 'fastq_bytes': '42'}
    value.indexed = {'read_run': [row], 'read_experiment': [deepcopy(row)]}
    value.linked = [{'provider': 'biosamples', 'kind': 'sample', 'accession': 'SAMEA1',
        'metadata': {'accession': 'SAMEA1', 'characteristics': {'tissue': [{'text': 'pancreas'}]}}},
        {'provider': 'biosamples', 'kind': 'sample', 'accession': 'SAMEA999',
        'metadata': {'accession': 'SAMEA999', 'characteristics': {'tissue': [{'text': 'unrelated'}]}}}]
    return value


def test_run_bound_fallback_restores_sample_and_files_with_partial_field_coverage():
    value = missing_sample_records()
    data = ENAParser().parse(value).to_mapping()
    assert [s['iid'] for s in data['sample']] == ['SAMEA1']
    sample = data['sample'][0]
    assert sample['sra_run'][0]['run'] == 'ERR1'
    assert sample['sra_run'][0]['files'][0]['md5'] == 'a'*32
    assert {'name': 'tissue', 'value': 'pancreas'} in sample['channel'][0]['characteristics']
    coverage = next(r['metadata'] for r in data['extensions']['insdc']['records'] if r['kind'] == 'retrieval_coverage')
    assert coverage['complete'] is False
    assert coverage['membership_complete'] is True
    assert coverage['entities']['sample']['expected'] == ['SAMEA1']
    assert coverage['entities']['sample']['fallback'] == ['SAMEA1']
    assert any('full Browser' in issue for issue in value.issues)


def test_wrong_study_rows_cannot_create_fallback_samples():
    value = missing_sample_records()
    for rows in value.indexed.values():
        rows[0]['study_accession'] = 'PRJEB999'
        rows[0]['secondary_study_accession'] = 'ERP999'
    assert all(s.iid != 'SAMEA1' for s in ENAParser().parse(value).samples)


def test_explicit_pool_members_remain_in_membership_when_sample_xml_is_missing():
    value = missing_sample_records()
    descriptor = value.xml[1].find('EXPERIMENT/DESIGN/SAMPLE_DESCRIPTOR')
    pool = ET.SubElement(descriptor, 'POOL')
    ET.SubElement(pool, 'MEMBER', accession='ERS1')
    ET.SubElement(pool, 'MEMBER', accession='ERS2')
    data = ENAParser().parse(value).to_mapping()
    assert {s['iid'] for s in data['sample']} == {'SAMEA1', 'ERS2'}
    assert all(s['sra_run'][0]['run']=='ERR1' for s in data['sample'])
    assert next(r['metadata'] for r in data['extensions']['insdc']['records']
                if r['kind']=='retrieval_coverage')['membership_complete'] is True


def test_ambiguous_alias_cannot_attach_a_run_to_an_arbitrary_sample():
    value = missing_sample_records()
    for rows in value.indexed.values():
        second = deepcopy(rows[0]); second['sample_accession'] = 'SAMEA2'; rows.append(second)
    data = ENAParser().parse(value).to_mapping()
    assert not any(s['sra_run'] for s in data['sample'])
    assert any('ambiguous' in issue for issue in value.issues)
