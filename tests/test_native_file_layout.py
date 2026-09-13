# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from copy import deepcopy
import json
import pytest
from tests.test_archive_export_cleanup import legacy
from tests.test_protocol_export import render


def file(uri, format='fastq', **values):
    return {'kind':'array_data_file','name':uri.rsplit('/',1)[-1], 'link':{'value':uri,'type':format},
            'comments':[{'name':k,'value':str(v)} for k,v in {'FORMAT':format,**values}.items()]}


def package():
    data=legacy();base=deepcopy(data['series']['assay_paths'][0]['steps']);base=[s for s in base if s['kind']!='array_data_file']
    paths=[]
    for f in [file('ftp://reads/one_1.fastq.gz',MD5='abc',BYTES=12),file('ftp://reads/one_2.fastq.gz'),
              file('ftp://reads/one_I1.fastq.gz'),file('ftp://reads/original.srf','SRF',ROLE='SUBMISSION_FILE'),
              file('ftp://reads/one_1.fastq.gz')]:
        paths.append({'steps':deepcopy(base)+[f]})
    data['series']['assay_paths']=paths
    return data


def projected(data):
    from meta_standards_converter.magetab.native_files import project_native_files
    project_native_files(data)
    return data['series']['assay_paths']


def comments(path):
    return [(c['name'],c['value']) for s in path['steps'] for c in s.get('comments',[])]


def test_read_relationships_keep_index_reads_and_submitted_columns():
    data=package();paths=projected(data)
    assert len(paths)==3
    assert all(('SUBMITTED_FILE_NAME','original.srf') in comments(p) for p in paths)
    assert all(('SUBMITTED_FILE_URI','ftp://reads/original.srf') in comments(p) for p in paths)
    assert {v for p in paths for k,v in comments(p) if k=='FASTQ_URI'}=={'ftp://reads/one_1.fastq.gz','ftp://reads/one_2.fastq.gz','ftp://reads/one_I1.fastq.gz'}
    assert ('FASTQ_MD5','abc') in comments(paths[0]) and ('FASTQ_BYTES','12') in comments(paths[0])
    once=deepcopy(data);projected(data);assert data==once


def test_no_fastq_is_one_run_with_grouped_archive_fields():
    data=package();base=data['series']['assay_paths'][0]['steps'][:-1]
    data['series']['assay_paths']=[{'steps':base+[file('https://reads/one.sra','SRA Normalized',MD5='abc')]}, {'steps':base+[file('https://reads/one.lite','SRA Lite')]}]
    paths=projected(data);assert len(paths)==1
    c=comments(paths[0]);assert ('FASTQ_URI','') in c
    assert [v for k,v in c if k=='ARCHIVE_FILE_URI']==['https://reads/one.sra','https://reads/one.lite']
    assert [v for k,v in c if k=='ARCHIVE_FILE_MD5']==['abc','']


def test_conflicting_files_and_different_workflows_are_not_collapsed():
    data=package();first=deepcopy(data['series']['assay_paths'][0]);first['steps'][-1]['comments'][1]['value']='different-checksum'
    data['series']['assay_paths'].append(first)
    other=deepcopy(first);other['steps'][0]['description']='different workflow';data['series']['assay_paths'].append(other)
    assert len(projected(data))==5


def test_processed_paths_and_multiple_runs_remain_scoped():
    data=package();processed=deepcopy(data['series']['assay_paths'][0]);processed['steps'].append({'kind':'derived_array_data_file','name':'counts.h5','link':{'value':'https://results/counts.h5'}})
    data['series']['assay_paths'].append(processed)
    sample_result={'steps':[deepcopy(processed['steps'][0]),{'kind':'derived_array_data_file','name':'sample.tsv','link':{'value':'https://results/sample.tsv'}}]}
    data['series']['assay_paths'].append(sample_result)
    another=deepcopy(data['series']['assay_paths'][0]);next(s for s in another['steps'] if s['kind']=='scan')['name']='SRR999999'
    next(s for s in another['steps'] if s['kind']=='scan')['comments']=[];data['series']['assay_paths'].append(another)
    paths=projected(data)
    assert len(paths)==6
    assert sum('counts.h5' in str(p) for p in paths)==1 and sample_result in paths
    other=next(p for p in paths if 'SRR999999' in str(p));assert not any(k.startswith('SUBMITTED_FILE') for k,v in comments(other))


@pytest.mark.parametrize('format',['GEO','MAGE-TAB'])
def test_non_native_layout_unchanged(format):
    data=package();data['source']['format']=format;before=deepcopy(data);projected(data);assert data==before


def test_ordinary_export_changes_only_private_copy():
    data=package();before=deepcopy(data);rows=render(data);table=next(r[1] for r in rows if r[0]=='SDRF File')
    assert len(table)==4 and 'Comment[FASTQ_URI]' in table[0]
    assert 'Array Data File' not in table[0]
    assert data==before


def test_saved_json_no_enrich_has_no_requests_or_input_changes(tmp_path):
    from meta_standards_converter.converters.json2ae import JSON2AEConverter
    from meta_standards_converter.magetab.constructor import AEConstructor
    class Forbidden:
        def __getattr__(self, name):
            raise AssertionError('Unexpected retrieval: ' + name)
    path=tmp_path/'native.json';path.write_text(json.dumps(package()));before=path.read_bytes()
    rows=JSON2AEConverter(enricher=Forbidden(),ae_constructor=AEConstructor(pubmed_client=Forbidden(),insdc_client=Forbidden())).convert(str(path),enrich=False)[0]
    table=next(r[1] for r in rows if r[0]=='SDRF File')
    assert len(table)==4 and path.read_bytes()==before


def test_archive_groups_and_versioned_identities_remain_aligned():
    data=package();base=data['series']['assay_paths'][0]['steps'][:-1]
    data['series']['assay_paths']=[{'steps':deepcopy(base)+[file('https://r/v1.sra','SRA',MD5='abc')]},
        {'steps':deepcopy(base)+[file('https://r/v2.sra','SRA',CHECKSUM='def',CHECKSUM_METHOD='SHA256')]}]
    rows=render(data);table=next(r[1] for r in rows if r[0]=='SDRF File');pairs=list(zip(table[0],table[1]))
    assert [v for k,v in pairs if k=='Comment[ARCHIVE_FILE_URI]']==['https://r/v1.sra','https://r/v2.sra']
    assert [v for k,v in pairs if k=='Comment[ARCHIVE_FILE_MD5]']==['abc','']
    assert [v for k,v in pairs if k=='Comment[ARCHIVE_FILE_CHECKSUM_METHOD]']==['','SHA256']


def test_enrichment_file_aliases_complete_the_same_read_relationship():
    data=package();base=data['series']['assay_paths'][0]['steps'][:-1]
    uri='ftp://reads/one_1.fastq.gz'
    enriched={'kind':'array_data_file','name':'one_1.fastq.gz','link':{'value':uri},
              'comments':[{'name':'FASTQ_URI','value':uri}]}
    untyped={'kind':'array_data_file','name':'one_1.fastq.gz','link':{'value':uri}}
    data['series']['assay_paths'].extend([{'steps':deepcopy(base)+[enriched]}, {'steps':deepcopy(base)+[untyped]}])
    paths=projected(data)
    assert len(paths)==3
    assert not any(k=='FASTQ_FASTQ_URI' for p in paths for k,v in comments(p))
    assert not any(k=='ARCHIVE_FILE_URI' and v==uri for p in paths for k,v in comments(p))
    # An explicit FASTQ_URI label is sufficient even without a format field.
    data=package();data['series']['assay_paths']=[{'steps':deepcopy(base)+[enriched]}]
    paths=projected(data);assert ('FASTQ_URI',uri) in comments(paths[0])


def test_repeated_sample_processed_relationship_is_not_an_extra_workflow():
    data=package();source=deepcopy(data['series']['assay_paths'][0]['steps'][0])
    result={'steps':[source,{'kind':'derived_array_data_file','name':'sample.tsv','link':{'value':'https://results/sample.tsv'}}]}
    different=deepcopy(result);different['steps'].insert(1,{'kind':'protocol_application','protocol_ref':'different-processing'})
    data['series']['assay_paths'].extend([result,deepcopy(result),different])
    paths=projected(data)
    assert len(paths)==5 and paths.count(result)==1 and different in paths
