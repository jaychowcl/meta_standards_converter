# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Only fully covered generated acquisitions may be absorbed into continuations."""
from copy import deepcopy
import pytest
from meta_standards_converter.metadata.archive_workflows import merge_workflows


def inputs():
    steps = [{'kind': 'source', 'name': 'SRS1', 'sample_ref': 'SRS1'},
             {'kind': 'extract', 'name': 'aliquot-A'},
             {'kind': 'assay', 'name': 'SRX1'}, {'kind': 'scan', 'name': 'SRR1'},
             {'kind': 'array_data_file', 'name': 'reads.fastq', 'link': {'value': 'https://example.org/reads.fastq'},
              'comments': [{'name': 'MD5', 'value': '123'}]}]
    sample = {'iid': 'SRS1', 'channel': [{}], 'sra_run': [{'run': 'SRR1', 'experiment': 'SRX1'}]}
    data = {'source': {'format': 'SRA'}, 'sample': [sample], 'series': {'iid': 'SRP1', 'assay_paths': [{'steps': steps}]}}
    extra = deepcopy(data);extra['source']['format'] = 'MAGE-TAB'
    extra['series']['assay_paths'][0]['steps'][-1].pop('comments')
    extra['series']['assay_paths'][0]['steps'].append({'kind': 'derived_array_data_file', 'name': 'counts.tsv',
        'link': {'value': 'https://example.org/counts.tsv'},
        'factor_values': [{'name': 'compound', 'value': 'drug', 'unit': {'value': 'mg'}}]})
    return data, extra


def merge(data, extra):
    issues=[]
    merge_workflows(data, extra, {'SRS1': 'SRS1'}, {}, True, issues)
    return data['series']['assay_paths'], issues


def test_generated_prefix_is_covered_and_complementary_file_evidence_retained():
    paths, issues = merge(*inputs())
    assert not issues
    assert len(paths) == 1
    assert paths[0]['steps'][-1]['factor_values'][0]['value'] == 'drug'
    assert {'name': 'MD5', 'value': '123'} in paths[0]['steps'][-2]['comments']


@pytest.mark.parametrize('difference', ['file', 'processing', 'authored', 'conflict'])
def test_distinct_or_authored_paths_are_not_pruned(difference):
    data, extra=inputs()
    if difference=='file':data['series']['assay_paths'][0]['steps'][-1]['link']['value']='https://example.org/other.fastq'
    elif difference=='processing':data['series']['assay_paths'][0]['steps'].append({'kind':'derived_array_data_file','name':'other.tsv','link':{'value':'https://example.org/other.tsv'}})
    elif difference=='authored':data['source']['format']='MAGE-TAB'
    else:extra['series']['assay_paths'][0]['steps'][-2]['comments']=[{'name':'MD5','value':'conflicting'}]
    paths, issues=merge(data,extra)
    assert len(paths)>=2
    if difference=='conflict': assert any('annotation scope' in x for x in issues)


def test_distinct_processing_continuations_remain_and_authored_occurrences_survive():
    data, extra=inputs()
    other=deepcopy(extra['series']['assay_paths'][0]);other['steps'][-1]['name']='other.tsv';other['steps'][-1]['link']['value']='https://example.org/other.tsv'
    extra['series']['assay_paths'] += [other,deepcopy(other)]
    paths,issues=merge(data,extra)
    assert len(paths)==3
    assert all(p['steps'][-1]['factor_values'] for p in paths)


def test_acquisition_annotations_are_shared_before_file_alternatives_are_rendered():
    data, extra=inputs()
    alternative=deepcopy(data['series']['assay_paths'][0])
    alternative['steps'][-1]['name']='submitted.fastq'
    alternative['steps'][-1]['link']['value']='https://example.org/submitted.fastq'
    alternative['steps'][0]['comments']=[{'name':'native alias','value':'SAMEA1'}]
    data['series']['assay_paths'].append(alternative)
    paths,issues=merge(data,extra)
    assert len(paths)==2
    assert paths[0]['steps'][:-1]==paths[1]['steps'][:-2]
    assert {'name':'native alias','value':'SAMEA1'} in paths[1]['steps'][0]['comments']


def test_missing_library_run_binding_does_not_authorize_reconciliation():
    data,extra=inputs()
    for document in (data,extra):
        for path in document['series']['assay_paths']:
            path['steps']=[s for s in path['steps'] if s['kind'] not in {'assay','scan'}]
    paths,issues=merge(data,extra)
    assert len(paths)==2
