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
from meta_standards_converter.metadata.archive_enrichment import merge_archive_metadata
from meta_standards_converter.miniml import MINiMLCodec
from tests.test_archive_fidelity_enrichment import workflow, enriched_native


def merge(extra, prefer=True):
    return merge_archive_metadata(enriched_native(), MINiMLCodec().decode(extra).package, prefer=prefer)


def test_paired_file_scan_aliases_share_one_verified_workflow():
    extra=workflow().to_mapping()
    path=extra['series']['assay_paths'][0]
    scan=next(s for s in path['steps'] if s['kind']=='scan')
    scan['name']='SRR11192680_1.fastq.gz'
    scan['comments'].append({'name':'FASTQ_URI','value':'ftp.sra.ebi.ac.uk/vol1/SRR11192680_1.fastq.gz'})
    second=deepcopy(path)
    scan2=next(s for s in second['steps'] if s['kind']=='scan')
    scan2['name']='SRR11192680_2.fastq.gz'
    scan2['comments'][-1]['value']='ftp.sra.ebi.ac.uk/vol1/SRR11192680_2.fastq.gz'
    for scan, checksum in [(scan, 'a' * 32), (scan2, 'b' * 32)]:
        scan['comments'].extend([{'name':'FASTQ_MD5','value':checksum}, {'name':'FASTQ_FILE_NAME','value':scan['name']}, {'name':'FASTQ_BYTES','value':'1234'}])
    extra['series']['assay_paths'].append(second)
    result,issues=merge(extra)
    assert not issues
    assert all(any(s['kind']=='extract' for s in p['steps']) for p in result.to_mapping()['series']['assay_paths'])
    uris=[s['link']['value'] for p in result.to_mapping()['series']['assay_paths'] for s in p['steps'] if s.get('link')]
    assert 'ftp://ftp.sra.ebi.ac.uk/vol1/SRR11192680_1.fastq.gz' in uris


@pytest.mark.parametrize('conflict',['name','comment'])
def test_conflicting_run_identifiers_are_rejected(conflict):
    extra=workflow().to_mapping();scan=next(s for s in extra['series']['assay_paths'][0]['steps'] if s['kind']=='scan')
    if conflict=='name':scan['name']='SRR999999'
    else:scan['comments'].append({'name':'SRA_RUN','value':'SRR999999'})
    result,issues=merge(extra)
    assert any('conflicting' in i for i in issues)
    assert not any(s['kind']=='extract' for p in result.to_mapping()['series']['assay_paths'] for s in p['steps'])


def test_raw_processing_derived_order_is_preserved():
    extra=workflow().to_mapping();steps=extra['series']['assay_paths'][0]['steps']
    steps.extend([{'kind':'array_data_file','name':'one.fastq.gz','link':{'value':'https://example.org/one.fastq.gz','type':'fastq'}},
                  {'kind':'normalization','name':'normalized counts'},
                  {'kind':'derived_array_data_file','name':'counts.tsv','link':{'value':'https://example.org/counts.tsv'}}])
    result,issues=merge(extra);assert not issues
    paths=result.to_mapping()['series']['assay_paths']
    derived=[p for p in paths if any(s.get('name')=='counts.tsv' for s in p['steps'])]
    assert derived
    assert [s['kind'] for s in derived[0]['steps']][-3:]==['array_data_file','normalization','derived_array_data_file']
    for path in paths:
        if not any(s['kind']=='derived_array_data_file' for s in path['steps']):
            assert not any(s['kind']=='normalization' for s in path['steps'])


def test_consistent_peer_only_run_keeps_workflow():
    extra=enriched_native().to_mapping()
    extra['sample'][0]['sra_run'][0]['run']='SRR99'
    for p in extra['series']['assay_paths']:
        for s in p['steps']:
            if s['kind']=='scan':
                s['name']='SRR99'
                for c in s.get('comments',[]):
                    if c['name'] in ('ENA_RUN','SRA_RUN'):c['value']='SRR99'
    result,issues=merge(extra,prefer=False)
    assert not issues
    assert any(s.get('name')=='SRR99' for p in result.to_mapping()['series']['assay_paths'] for s in p['steps'])


def test_different_paired_workflow_methods_remain_ambiguous():
    extra=workflow().to_mapping();other=deepcopy(extra['series']['assay_paths'][0])
    next(s for s in other['steps'] if s['kind']=='extract')['material_type']={'value':'DNA'}
    extra['series']['assay_paths'].append(other)
    result,issues=merge(extra)
    assert any('ambiguous' in i for i in issues)
    assert not any(s['kind']=='extract' for p in result.to_mapping()['series']['assay_paths'] for s in p['steps'])


def test_unmatched_sample_file_cannot_suppress_matched_sample_result():
    extra=workflow().to_mapping();sample=extra['sample'][0]
    uri='https://example.org/shared.tsv'
    sample['supplementary_data']=[{'value':uri}]
    other=deepcopy(extra['series']['assay_paths'][0])
    for step in other['steps']:
        if step.get('sample_ref'):step['sample_ref']='UNMATCHED'
        if step['kind'] in ('scan','assay'):step['comments']=[];step['name']='unknown'
    other['steps'].append({'kind':'derived_array_data_file','name':'shared.tsv','link':{'value':uri}})
    extra['series']['assay_paths'].append(other)
    result,issues=merge(extra)
    paths=[p for p in result.to_mapping()['series']['assay_paths'] if any(s.get('link',{}).get('value')==uri for s in p['steps'])]
    assert paths
    assert not any(s['kind']=='scan' for p in paths for s in p['steps'])
    assert any(c['name']=='organism' for c in paths[0]['steps'][0]['characteristics'])


def test_ftp_location_reserved_filename_is_preserved():
    from meta_standards_converter.metadata.archive_workflows import file_node
    node=file_node({'uri':'ftp.sra.ebi.ac.uk/vol1/read#1.fastq.gz'})
    assert node['name']=='read#1.fastq.gz'
    assert node['link']['value']=='ftp://ftp.sra.ebi.ac.uk/vol1/read%231.fastq.gz'


def test_explicit_unmatched_sample_ref_is_not_recovered_by_run():
    extra=workflow().to_mapping()
    for step in extra['series']['assay_paths'][0]['steps']:
        if step.get('sample_ref'):step['sample_ref']='unmatched-sample'
    result,issues=merge(extra)
    assert any('sample binding' in i for i in issues)
    assert not any(s['kind']=='extract' for p in result.to_mapping()['series']['assay_paths'] for s in p['steps'])


def test_fileless_processing_is_not_applied_to_native_raw_files():
    extra=workflow().to_mapping()
    extra['series']['assay_paths'][0]['steps'].append({'kind':'normalization','name':'counts'})
    result,issues=merge(extra)
    assert not issues
    assert not any(s['kind']=='normalization' for p in result.to_mapping()['series']['assay_paths'] for s in p['steps'])


@pytest.mark.parametrize('scope',['sample','experiment'])
def test_sparse_linked_paths_cannot_remove_native_assay_or_run(scope):
    extra=workflow().to_mapping()
    if scope=='experiment':
        extra['series']['assay_paths'][0]['steps']=[s for s in extra['series']['assay_paths'][0]['steps'] if s['kind']!='scan']
    else:
        first=extra['series']['assay_paths'][0]['steps'][0]
        extra['series']['assay_paths']=[{'steps':[first,{'kind':'derived_array_data_file','name':'counts.tsv','link':{'value':'https://example.org/counts.tsv'}}]}]
    result,issues=merge(extra);assert not issues
    paths=result.to_mapping()['series']['assay_paths']
    raw=[p for p in paths if any(s['kind']=='array_data_file' for s in p['steps'])]
    assert raw
    for p in raw:
        assert any(s['kind']=='scan' and s['name']=='SRR11192680' for s in p['steps'])
        assert any(s['kind']=='assay' and s['name']=='SRX7812918' for s in p['steps'])
