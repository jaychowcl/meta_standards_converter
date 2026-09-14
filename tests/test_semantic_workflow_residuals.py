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
from meta_standards_converter.miniml.archive_residuals import finalize
from tests.test_archive_fidelity_enrichment import workflow, enriched_native


def aliased():
    extra=workflow().to_mapping();p=extra['series']['assay_paths'][0]
    scan=next(s for s in p['steps'] if s['kind']=='scan');scan['name']='submitted.srf'
    scan['comments'].extend([{'name':'FASTQ_URI','value':'ftp.sra.ebi.ac.uk/vol1/read_1.fastq.gz'},
                             {'name':'FASTQ_MD5','value':'a'*32},
                             {'name':'SUBMITTED_FILE_NAME','value':'submitted.srf'}])
    second=deepcopy(p)
    for c in second['steps'][-1]['comments']:
        if c['name']=='FASTQ_URI':c['value']='ftp.sra.ebi.ac.uk/vol1/read_2.fastq.gz'
        if c['name']=='FASTQ_MD5':c['value']='b'*32
    extra['series']['assay_paths'].append(second)
    data,issues=merge_archive_metadata(enriched_native(),MINiMLCodec().decode(extra).package,prefer=True)
    assert not issues
    return data.to_mapping()


def paths(data):
    return [p for r in data['extensions']['insdc']['records'] if r['kind']=='MINiML' and r['accession']=='E-MTAB-1'
            for p in r['metadata'].get('series',{}).get('assay_paths',[])]


def test_verified_aliased_workflows_retain_names_not_mapped_payload():
    data=aliased();residuals=paths(data)
    assert len(residuals)==2
    for p in residuals:
        assert p['document']
        assert p['steps'][0]['name']=='embryo'
        assert next(s['name'] for s in p['steps'] if s['kind']=='scan')=='submitted.srf'
        assert next(s['name'] for s in p['steps'] if s['kind']=='assay')=='assay'
        assert not any(k in s for s in p['steps'] for k in ('characteristics','material_type','factor_values'))
        assert not any(c['name'] in ('FASTQ_MD5','SUBMITTED_FILE_NAME') for s in p['steps'] for c in s.get('comments',[]))
    assert {c['value'] for p in residuals for s in p['steps'] for c in s.get('comments',[]) if c['name']=='FASTQ_URI'}=={
        'ftp.sra.ebi.ac.uk/vol1/read_1.fastq.gz','ftp.sra.ebi.ac.uk/vol1/read_2.fastq.gz'}
    again=finalize(deepcopy(data)).to_mapping();assert again==data


def legacy_source(data):
    # Recreate the complete original incoming snapshot independently of the
    # residual pruner, preserving the parser's document and descriptive names.
    source=workflow().to_mapping();p=source['series']['assay_paths'][0]
    sid=data['sample'][0]['iid']
    for s in p['steps']:
        if s.get('sample_ref'):s['sample_ref']=sid
        if s.get('protocol_ref'):s['protocol_ref']='E-MTAB-1:'+s['protocol_ref']
    scan=next(s for s in p['steps'] if s['kind']=='scan');scan['name']='submitted.srf'
    scan['comments'].extend([{'name':'FASTQ_URI','value':'ftp.sra.ebi.ac.uk/vol1/read_1.fastq.gz'},
                            {'name':'FASTQ_MD5','value':'a'*32},
                            {'name':'SUBMITTED_FILE_NAME','value':'submitted.srf'}])
    return {'provider':'MAGE-TAB','kind':'MINiML','accession':'E-MTAB-1',
            'metadata':{'series':{'iid':data['series']['iid'],'assay_paths':[p]}}}


def test_alias_projection_preserves_unmapped_siblings_and_coupled_values():
    data=aliased();record=legacy_source(data);p=record['metadata']['series']['assay_paths'][0]
    p['steps'][0]['comments']=[{'name':'unknown detail','value':'retain me'}]
    for path in data['series']['assay_paths']:
        for step in path['steps']:
            for value in step.get('factor_values',[]):value['unit']={'value':'hours'}
    core={k:deepcopy(v) for k,v in data.items() if k!='extensions'}
    out=finalize(deepcopy(data),[record]).to_mapping();remaining=paths(out)[0]
    assert remaining['steps'][0]['comments']==p['steps'][0]['comments']
    assert 'characteristics' not in remaining['steps'][0]
    assert next(s['factor_values'] for s in remaining['steps'] if s.get('factor_values'))==next(s['factor_values'] for s in p['steps'] if s.get('factor_values'))
    assert {k:v for k,v in out.items() if k!='extensions'}==core


@pytest.mark.parametrize('case',['source_duplicate','sparse_source_duplicate','target_duplicate','wrong_run','file_conflict','procedure_order','unbound_sample','unknown_file_fact','qualified_file_fact','qualified_filename','qualified_uri'])
def test_alias_projection_declines_ambiguous_or_unverified_correspondence(case):
    data=aliased();record=legacy_source(data);source=record['metadata']['series']['assay_paths'][0]
    matching=[p for p in data['series']['assay_paths'] if any(s.get('link',{}).get('value')=='ftp://ftp.sra.ebi.ac.uk/vol1/read_1.fastq.gz' for s in p['steps'])]
    assert matching
    if case=='source_duplicate':record['metadata']['series']['assay_paths'].append(deepcopy(source))
    elif case=='sparse_source_duplicate':
        sparse=deepcopy(source)
        sparse['steps'][-1]['comments']=[c for c in sparse['steps'][-1]['comments'] if c['name']!='FASTQ_MD5']
        record['metadata']['series']['assay_paths'].append(sparse)
    elif case=='target_duplicate':data['series']['assay_paths'].append(deepcopy(matching[0]))
    elif case=='wrong_run':data['sample'][0]['sra_run'][0]['experiment']='SRX99999'
    elif case=='unbound_sample':
        for s in source['steps']:
            if s.get('sample_ref'):s['sample_ref']='unbound'
    elif case=='unknown_file_fact':
        source['steps'][-1]['comments'].append({'name':'FASTQ_BYTES','value':'unknown'})
    elif case=='qualified_file_fact':
        next(c for c in source['steps'][-1]['comments'] if c['name']=='FASTQ_MD5')['qualification']='source-only'
    elif case=='qualified_filename':
        source['steps'][-1]['comments'].append({'name':'FASTQ_FILE_NAME','value':'read_1.fastq.gz','assertion_context':'source-only-qualifier'})
    elif case=='qualified_uri':
        next(c for c in source['steps'][-1]['comments'] if c['name']=='FASTQ_URI')['assertion_context']='source-only-qualifier'
    elif case=='file_conflict':
        for p in matching:
            for s in p['steps']:
                for c in s.get('comments',[]):
                    if c['name']=='MD5':c['value']='c'*32
    else:
        for p in matching:
            indices=[i for i,s in enumerate(p['steps']) if s.get('protocol_ref')]
            a,b=indices[:2];p['steps'][a],p['steps'][b]=p['steps'][b],p['steps'][a]
    out=finalize(data,[record]).to_mapping()
    assert paths(out)[0]==source


def test_plain_filename_projection_and_repeated_application_positions():
    data=aliased();record=legacy_source(data);source=record['metadata']['series']['assay_paths'][0]
    source['steps'][-1]['comments'].append({'name':'FASTQ_FILE_NAME','value':'read_1.fastq.gz'})
    first=next(s for s in source['steps'] if s.get('protocol_ref'))
    second=deepcopy(first);first['comments']=[{'name':'duration','value':'10'}];second['comments']=[{'name':'duration','value':'20'}]
    source['steps'].insert(source['steps'].index(first)+1,second)
    for p in data['series']['assay_paths']:
        i=next(i for i,s in enumerate(p['steps']) if s.get('protocol_ref'))
        p['steps'][i:i+1]=[deepcopy(second),deepcopy(first)]
    out=finalize(data,[record]).to_mapping();remaining=paths(out)[0]
    assert 'characteristics' not in remaining['steps'][0]
    assert [s for s in remaining['steps'] if s.get('protocol_ref')==first['protocol_ref']]==[first,second]
    assert not any(c['name']=='FASTQ_FILE_NAME' for s in remaining['steps'] for c in s.get('comments',[]))


def mex_alias():
    data=aliased();record=legacy_source(data);source=record['metadata']['series']['assay_paths'][0]
    tail=[]
    for role in ('barcodes.tsv.gz','features.tsv.gz','matrix.mtx.gz'):
        tail.extend([{'kind':'protocol_application','protocol_ref':'processing'},
                     {'kind':'derived_array_data_file','name':'pool_'+role,'comments':[{'name':'MD5','value':'d'*32}]}])
    source['steps'].extend(deepcopy(tail))
    data['series']['protocols'].append({'name':'processing','type':{'value':'alignment protocol'},'description':'Cell Ranger generates the matrix.'})
    for path in data['series']['assay_paths']:
        if any(s.get('link',{}).get('value','').endswith('read_1.fastq.gz') for s in path['steps']):
            path['steps'].extend(deepcopy(tail))
            for s in path['steps']:
                if s['kind']=='derived_array_data_file':s['link']={'value':s['name'],'repository':'ArrayExpress','source_accession':'E-MTAB-1'}
    return data,record


def test_filename_only_mex_members_prune_only_against_verified_same_path_bundle():
    data,record=mex_alias()
    out=finalize(data,[record]).to_mapping()
    left=paths(out)[0]
    assert not any('characteristics' in s or 'material_type' in s or 'comments' in s and any(c['name']=='MD5' for c in s['comments']) for s in left['steps'])
    assert any(r['kind']=='result_layout' for r in out['extensions']['insdc']['records'])
    assert sum(s['kind']=='derived_array_data_file' for p in out['series']['assay_paths'] for s in p['steps'])==1
    assert finalize(out).to_mapping()==out


@pytest.mark.parametrize('case',['method','run','directory','prefix','checksum','qualifier','opaque','wrong_role'])
def test_filename_only_mex_cannot_borrow_an_unproved_companion(case):
    from meta_standards_converter.miniml.archive_results import normalize_result_bundles
    data,record=mex_alias();original=deepcopy(record['metadata']['series']['assay_paths'][0])
    if case=='method':data['series']['protocols'][-1]['description']='unrelated'
    normalize_result_bundles(data)
    matrix=next(s for p in data['series']['assay_paths'] for s in p['steps'] if s.get('link',{}).get('companion_files')) if case!='method' else None
    if case=='run':data['sample'][0]['sra_run'][0]['experiment']='SRX999'
    if case=='directory':matrix['link']['companion_files'][0]['node']['link']['value']='https://other/pool_barcodes.tsv.gz'
    if case=='prefix':matrix['link']['companion_files'][0]['node']['name']='other_barcodes.tsv.gz'
    if case=='checksum':matrix['link']['companion_files'][0]['node']['comments'][0]['value']='conflict'
    if case=='qualifier':record['metadata']['series']['assay_paths'][0]['steps'][-5]['comments'][0]['qualification']='unmapped'
    if case=='opaque':matrix['link']['companion_files'][0]['node']['link']='opaque'
    if case=='wrong_role':matrix['link']['companion_files'][0]['role']='unrelated'
    original=deepcopy(record['metadata']['series']['assay_paths'][0])
    out=finalize(data,[record]).to_mapping()
    assert paths(out)[0]==original


def test_organism_label_projection_requires_complete_unique_coupled_value():
    data=aliased();record=legacy_source(data);source=record['metadata']['series']['assay_paths'][0]
    source['steps'][0]['characteristics']=[{'name':'Organism','value':'Danio rerio','term_source_ref':'NCBITaxon','term_accession_number':'7955'}]
    for path in data['series']['assay_paths']:
        path['steps'][0]['characteristics']=[{'name':'organism','value':'Danio rerio','term_source_ref':'NCBITaxon','term_accession_number':'7955'}]
    out=finalize(deepcopy(data),[record]).to_mapping()
    assert not paths(out)[0]['steps'][0].get('characteristics')
    source['steps'][0]['characteristics'][0]['term_accession_number']='conflict'
    out=finalize(data,[record]).to_mapping()
    assert paths(out)[0]['steps'][0]['characteristics']==source['steps'][0]['characteristics']
