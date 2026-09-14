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
    assert len(paths)==5  # The complete counts branch also represents its raw prefix.
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


def test_checksum_verified_mirrors_keep_alternate_metadata_without_extra_rows():
    data=package();base=data['series']['assay_paths'][0]['steps'][:-1]
    one=file('https://ebi/reads.fastq.gz',MD5='a'*32,BYTES=100,ROLE='GENERATED_FILE')
    mirror=file('s3://bucket/reads.fastq.gz.1',CHECKSUM='a'*32,CHECKSUM_METHOD='MD5',BYTES=100,ROLE='GENERATED_FILE')
    mirror['name']=one['name']
    alias=deepcopy(one);alias['comments']=[];alias['link'].pop('type')
    data['series']['assay_paths']=[{'steps':deepcopy(base)+[f]} for f in [one,mirror,alias]]
    before=deepcopy(data);table=next(r[1] for r in render(data) if r[0]=='SDRF File')
    assert len(table)==2
    row=list(zip(table[0],table[1]))
    assert ('Comment[FASTQ_URI]','https://ebi/reads.fastq.gz') in row
    assert ('Comment[FASTQ_ALTERNATIVE_URI]','s3://bucket/reads.fastq.gz.1') in row
    assert ('Comment[FASTQ_ALTERNATIVE_CHECKSUM_METHOD]','MD5') in row
    assert data==before


@pytest.mark.parametrize('change', ['checksum','algorithm','size','name','workflow','run','invalid'])
def test_mirror_matching_rejects_conflicts_and_different_scopes(change):
    data=package();first=deepcopy(data['series']['assay_paths'][0]);first['steps'][-1]=file('https://ebi/one.fastq.gz',MD5='a'*32,BYTES=100)
    second=deepcopy(first);second['steps'][-1]['link']['value']='s3://bucket/one.fastq.gz'
    node=second['steps'][-1]
    if change=='checksum':node['comments'][1]['value']='b'*32
    if change=='algorithm':node['comments'][1]={'name':'CHECKSUM','value':'a'*32};node['comments'].append({'name':'CHECKSUM_METHOD','value':'unknown'})
    if change=='size':node['comments'][2]['value']='101'
    if change=='name':node['name']='one.v2.fastq.gz'
    if change=='workflow':second['steps'][0]['description']='different'
    if change=='run':next(s for s in second['steps'] if s['kind']=='scan')['name']='SRR999999'
    if change=='invalid':
        first['steps'][-1]['comments'][1]['value']='abc';node['comments'][1]['value']='abc'
    data['series']['assay_paths']=[first,second]
    assert len(projected(data))==2


def test_sparse_raw_path_is_subsumed_by_complete_result_branches_with_metadata():
    data=package();raw=deepcopy(data['series']['assay_paths'][0])
    complete=[]
    for name in ('counts.tsv','abundance.tsv'):
        p=deepcopy(raw);p['steps'][-1]['comments']=[{'name':'FORMAT','value':'fastq'}]
        p['steps'] += [{'kind':'protocol_application','protocol_ref':'processing'},
                       {'kind':'derived_array_data_file','name':name,'link':{'value':'https://results/'+name}}]
        complete.append(p)
    data['series']['assay_paths']=[raw,*complete]
    paths=projected(data)
    assert len(paths)==2
    assert all(('FASTQ_MD5','abc') in comments(p) and ('FASTQ_BYTES','12') in comments(p) for p in paths)
    once=deepcopy(data);projected(data);assert data==once


def test_sparse_sample_result_requires_unique_file_identity_and_same_processing():
    data=package();raw=deepcopy(data['series']['assay_paths'][0]);source=deepcopy(raw['steps'][0])
    proc={'kind':'protocol_application','protocol_ref':'processing'}
    result={'kind':'derived_array_data_file','name':'counts.tsv','link':{'value':'https://results/v1/counts.tsv'}}
    full={'steps':raw['steps']+[proc,result]}
    sparse={'steps':[source,deepcopy(proc),{'kind':'derived_array_data_file','name':'counts.tsv','link':{'value':'counts.tsv'}}]}
    compressed={'steps':[source,deepcopy(proc),{'kind':'derived_array_data_file','name':'counts.tsv.gz','link':{'value':'https://results/counts.tsv.gz'}}]}
    data['series']['assay_paths']=[deepcopy(sparse),deepcopy(full),compressed]
    paths=projected(data);assert len(paths)==2 and compressed in paths
    other=deepcopy(full);other['steps'][-1]['link']['value']='https://results/v2/counts.tsv'
    for sequence in ([sparse,full,other],[other,full,sparse]):
        data['series']['assay_paths']=deepcopy(sequence)
        assert len(projected(data))==3


@pytest.mark.parametrize('change',['checksum','processing','sample','workflow'])
def test_sparse_consolidation_preserves_conflicts_and_boundaries(change):
    data=package();raw=deepcopy(data['series']['assay_paths'][0]);full=deepcopy(raw)
    full['steps'].append({'kind':'derived_array_data_file','name':'counts.tsv'})
    if change=='checksum':full['steps'][-2]['comments'][1]['value']='different'
    if change=='processing':raw['steps'].insert(-1,{'kind':'protocol_application','protocol_ref':'different'})
    if change=='sample':full['steps'][0]['sample_ref']='other-sample'
    if change=='workflow':full['steps'][0]['description']='different'
    data['series']['assay_paths']=[raw,full]
    assert len(projected(data))==2


@pytest.mark.parametrize('conflicting',[False,True])
def test_sparse_file_completion_does_not_drop_empty_or_repeated_metadata(conflicting):
    data=package();raw=deepcopy(data['series']['assay_paths'][0]);full=deepcopy(raw)
    full['steps'].append({'kind':'derived_array_data_file','name':'counts.tsv'})
    if conflicting:raw['steps'][-1]['comments'].append({'name':'MD5','value':'contradiction'})
    else:full['steps'][-2]['comments'][1]['value']=''
    data['series']['assay_paths']=[raw,full]
    paths=projected(data)
    if conflicting:
        assert len(paths)==2 and 'contradiction' in str(paths)
    else:
        assert len(paths)==1 and ('FASTQ_MD5','abc') in comments(paths[0])


def test_sparse_explicit_uri_completes_filename_only_longer_branch():
    data=package();raw=deepcopy(data['series']['assay_paths'][0]);full=deepcopy(raw)
    full['steps'][-1]['link']['value']=full['steps'][-1]['name']
    full['steps'].append({'kind':'derived_array_data_file','name':'counts.tsv'})
    for sequence in ([raw,full],[full,raw]):
        data['series']['assay_paths']=deepcopy(sequence)
        paths=projected(data)
        assert len(paths)==1 and ('FASTQ_URI','ftp://reads/one_1.fastq.gz') in comments(paths[0])


@pytest.mark.parametrize('sparse_first',[False,True])
def test_native_sdrf_columns_preserve_each_acquisition_and_result_path(sparse_first):
    data=package();base=deepcopy(data['series']['assay_paths'][0]['steps']);source=deepcopy(base[0])
    data['series']['protocols']=[{'name':'prep','description':'Prepare explicit material','type':{'value':'nucleic acid extraction protocol'}},
                               {'name':'process','description':'Process explicit result','type':{'value':'normalization data transformation protocol'}}]
    base=[s for s in base if s['kind']!='protocol_application']
    base[1:1]=[{'kind':'protocol_application','protocol_ref':'prep'},{'kind':'extract','name':'explicit extract','material_type':{'value':'RNA'}}]
    result={'kind':'derived_array_data_file','name':'counts.tsv','link':{'value':'https://results/counts.tsv'}}
    processing={'kind':'protocol_application','protocol_ref':'process'}
    full={'steps':base+[deepcopy(processing),deepcopy(result)]}
    sparse={'steps':[source,processing,{**result,'name':'sample.tsv','link':{'value':'https://results/sample.tsv'}}]}
    data['series']['assay_paths']=[sparse,full] if sparse_first else [full,sparse]
    rows={r[0]:r[1:] for r in render(data)};table=rows['SDRF File'][0];header=table[0]
    assert header.index('Extract Name')<header.index('Assay Name')<header.index('Scan Name')<header.index('Derived Array Data File')
    descriptions=dict(zip(rows['Protocol Name'],rows['Protocol Description']))
    for row in table[1:]:
        methods=[descriptions[v] for k,v in zip(header,row) if k=='Protocol REF' and v]
        if row[header.index('Assay Name')]:assert methods==['Prepare explicit material','Process explicit result']
        else:assert methods==['Process explicit result']


def test_ordered_native_renderer_keeps_protocol_blocks_repeated_materials_and_end():
    from meta_standards_converter.magetab.semantics import render_miniml_assay_documents
    paths=[{'steps':[{'kind':'source','name':'s'},{'kind':'protocol_application','protocol_ref':'a'},
                     {'kind':'protocol_application','protocol_ref':'b'},{'kind':'extract','name':'e1'},
                     {'kind':'protocol_application','protocol_ref':'c'},{'kind':'extract','name':'e2'},
                     {'kind':'protocol_application','protocol_ref':'end'}]},
           {'steps':[{'kind':'source','name':'s2'},{'kind':'protocol_application','protocol_ref':'only'},
                     {'kind':'extract','name':'e3'}]}]
    table=next(iter(render_miniml_assay_documents(paths,preserve_order=True).values()))
    assert [v for v in table[1] if v]==['s','a','b','e1','c','e2','end']
    assert [v for v in table[2] if v]==['s2','only','e3']
    crossing=[{'steps':[{'kind':'source','name':'s'},{'kind':'assay','name':'a'},{'kind':'scan','name':'r'}]},
              {'steps':[{'kind':'source','name':'s2'},{'kind':'scan','name':'r2'},{'kind':'assay','name':'a2'}]}]
    with pytest.raises(ValueError,match='Incompatible native SDRF path ordering'):
        render_miniml_assay_documents(crossing,preserve_order=True)
