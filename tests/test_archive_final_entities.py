# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from copy import deepcopy
from tests.test_native_archive_enrichment import native
from tests.test_protocol_export import render
from meta_standards_converter.miniml.archive_residuals import finalize


def test_native_platforms_follow_explicit_instruments_and_mixed_scope():
    data=native().to_mapping();sample=data['sample'][0]
    assert data['platform']
    platform={p['iid']:p for p in data['platform']}[sample['platform_ref']['ref']]
    assert platform['title']==sample['sra_run'][0]['instrument_model']
    assert not platform.get('accession')
    run=deepcopy(sample['sra_run'][0]);run.update(run='SRR99',experiment='SRX99',instrument_model='Sequel II')
    sample['sra_run'].append(run)
    data=finalize(data,[]).to_mapping();sample=data['sample'][0]
    assert len(data['platform'])==2
    assert not sample.get('platform_ref')
    assert len({r['platform_ref']['ref'] for r in sample['sra_run']})==2


def test_native_database_is_repository_only_and_used_ontologies_export():
    data=native().to_mapping()
    assert not {'INSDC','NCBITaxon','EFO','UO'} & {d['iid'] for d in data['database']}
    data['database'].extend([{'iid':'UBERON','name':'UBERON','url':'https://example.org/uberon.owl','version':'source-version'}, {'iid':'UNUSED','name':'UNUSED','version':'do not export'}])
    data['series']['assay_paths'][0]['steps'][0]['characteristics'].append({'name':'tissue','value':'cortex','term_source_ref':'UBERON','term_accession_number':'UBERON:0000956'})
    data=finalize(data,[]).to_mapping()
    assert 'UBERON' not in {d['iid'] for d in data['database']}
    rows={r[0]:r[1:] for r in render(data)}
    assert 'UBERON' in rows['Term Source Name']
    assert 'UNUSED' not in rows['Term Source Name']
    assert 'SRA' not in rows['Term Source Name']
    i=rows['Term Source Name'].index('UBERON')
    assert rows['Term Source Version'][i]=='source-version'


def test_identical_organizations_coalesce_without_merging_contacts_or_scopes():
    data=native().to_mapping();sid=data['sample'][0]['iid']
    data['organization']=[{'iid':'o1','name':'Lab','role':'owner','sample_accession':sid},
                          {'iid':'o2','name':'Lab','role':'owner','sample_accession':'SAMN2'},
                          {'iid':'o3','name':'Lab','role':'owner','address':{'city':'Elsewhere'}}]
    data['contributor']=[{'iid':'c1','person':{'last':'Contact'},'organization_ref':{'ref':'o1'}},
                         {'iid':'c2','person':{'last':'Contact'},'organization_ref':{'ref':'o2'}}]
    data['series']['contributor_ref']=[{'ref':'c2'}];data['sample'][0]['contact_ref']=[{'ref':'c1'}]
    data=finalize(data,[]).to_mapping()
    assert len(data['organization'])==2
    assert len(data['contributor'])==2
    assert data['contributor'][0]['organization_ref']==data['contributor'][1]['organization_ref']
    assert data['sample'][0]['contact_ref']==[{'ref':'c1'}]
    assert data['series']['contributor_ref']==[{'ref':'c2'}]
    org=next(o for o in data['organization'] if o['iid']=='o1')
    assert {o['iid'] for o in org['source_occurrences']}=={'o1','o2'}
    assert finalize(data,[]).to_mapping()==data


def test_local_vocabulary_name_survives_repository_separation():
    data=native().to_mapping();data['database'].append({'iid':'local-vocabulary','name':'Supplied vocabulary label'})
    data['sample'][0]['channel'][0]['characteristics'].append({'name':'condition','value':'x','term_source_ref':'local-vocabulary'})
    data=finalize(data,[]).to_mapping()
    assert any(r['kind']=='term_source_declaration' and r['metadata']['name']=='Supplied vocabulary label' for r in data['extensions']['insdc']['records'])


def test_missing_instrument_markers_do_not_create_platforms():
    from meta_standards_converter.miniml.archive_entities import local_platforms
    for literal in ('not provided', 'unspecified', 'missing: not reported'):
        data={'source':{'format':'ENA'},'sample':[{'iid':'SAM1','sra_run':[{'run':'ERR1','instrument_model':literal}]}]}
        local_platforms(data)
        assert not data['platform']
        assert not data['sample'][0].get('platform_ref')
        assert data['sample'][0]['sra_run'][0]['instrument_model']==literal


def test_saved_native_vocabulary_conflict_is_namespaced_during_enrichment():
    from tests.test_native_archive_enrichment import linked
    from meta_standards_converter.miniml import MINiMLCodec
    from meta_standards_converter.metadata.archive_enrichment import merge_archive_metadata
    data=native().to_mapping();data['database'].append({'iid':'LOCAL','name':'Native vocabulary','url':'https://native/vocab'})
    data['sample'][0]['channel'][0]['characteristics'].append({'name':'native trait','value':'a','term_source_ref':'LOCAL'})
    package=finalize(data,[]);extra=linked(package,'E-MTAB-1','title').to_mapping()
    extra['source']['format']='MAGE-TAB';extra['database']=[{'iid':'LOCAL','name':'Incoming vocabulary','url':'https://incoming/vocab'}]
    extra['sample'][0]['channel'][0]['characteristics']=[{'name':'incoming trait','value':'b','term_source_ref':'LOCAL'}]
    result,issues=merge_archive_metadata(package,MINiMLCodec().decode(extra).package,prefer=True);assert not issues
    data=result.to_mapping();traits={c['name']:c for c in data['sample'][0]['channel'][0]['characteristics']}
    assert traits['native trait']['term_source_ref']=='LOCAL'
    assert traits['incoming trait']['term_source_ref']=='E-MTAB-1:LOCAL'
    rows={r[0]:r[1:] for r in render(data)}
    assert len(rows['Term Source Name'])==len(set(rows['Term Source Name']))


def test_compatible_persisted_vocabulary_details_complete_once():
    from tests.test_native_archive_enrichment import linked
    from meta_standards_converter.miniml import MINiMLCodec
    from meta_standards_converter.metadata.archive_enrichment import merge_archive_metadata
    data=native().to_mapping();data['database'].append({'iid':'LOCAL','name':'Supplied vocabulary'})
    data['sample'][0]['channel'][0]['characteristics'].append({'name':'trait','value':'a','term_source_ref':'LOCAL'})
    package=finalize(data,[]);extra=linked(package,'E-MTAB-1','title').to_mapping()
    extra['database']=[{'iid':'LOCAL','name':'Supplied vocabulary','url':'https://supplied/vocab'}]
    result,_=merge_archive_metadata(package,MINiMLCodec().decode(extra).package,prefer=True)
    records=[r for r in result.to_mapping()['extensions']['insdc']['records'] if r['kind']=='term_source_declaration' and r['accession']=='LOCAL']
    assert len(records)==1 and records[0]['metadata']['url']=='https://supplied/vocab'


def test_saved_peer_vocabulary_conflicts_remap_both_core_and_retained_declarations():
    from meta_standards_converter.metadata.archive_enrichment import merge_archive_metadata
    packages=[]
    for label in ('Native','Peer'):
        data=native().to_mapping();data['database'].append({'iid':'LOCAL','name':label+' vocabulary','url':'https://'+label.lower()+'/vocab'})
        data['sample'][0]['channel'][0]['characteristics'].append({'name':label+' trait','value':label,'term_source_ref':'LOCAL'})
        packages.append(finalize(data,[]))
    result,issues=merge_archive_metadata(*packages,prefer=False);assert not issues
    data=result.to_mapping();traits={c['name']:c for c in data['sample'][0]['channel'][0]['characteristics']}
    assert traits['Native trait']['term_source_ref']=='LOCAL'
    assert traits['Peer trait']['term_source_ref']==data['series']['iid']+':LOCAL'
    rows={r[0]:r[1:] for r in render(data)}
    assert len(rows['Term Source Name'])==len(set(rows['Term Source Name']))
    peer=rows['Term Source Name'].index(data['series']['iid']+':LOCAL')
    assert rows['Term Source File'][peer]=='https://peer/vocab'


def test_repeated_enrichment_namespaces_changed_vocabulary_versions_without_collision():
    from tests.test_native_archive_enrichment import linked
    from meta_standards_converter.miniml import MINiMLCodec
    from meta_standards_converter.metadata.archive_enrichment import merge_archive_metadata
    data=native().to_mapping();data['database'].append({'iid':'LOCAL','name':'Native vocabulary'})
    package=finalize(data,[])
    refs=[]
    for version in ('1','2','2'):
        extra=linked(native(),'E-MTAB-1','title').to_mapping()
        extra['database']=[{'iid':'LOCAL','name':'Incoming vocabulary','url':'https://incoming/vocab','version':version}]
        extra['sample'][0]['channel'][0]['characteristics']=[{'name':'trait version '+version,'value':version,'term_source_ref':'LOCAL'}]
        package,issues=merge_archive_metadata(package,MINiMLCodec().decode(extra).package,prefer=True)
        assert not issues
        data=package.to_mapping();traits={c['name']:c for c in data['sample'][0]['channel'][0]['characteristics']}
        refs.append(traits['trait version '+version]['term_source_ref'])
    assert refs[0]!=refs[1] and refs[1]==refs[2]
    records=[r['metadata'] for r in data['extensions']['insdc']['records'] if r['kind']=='term_source_declaration']
    assert len({r['iid'] for r in records})==len(records)
    rows={r[0]:r[1:] for r in render(data)}
    for ref,version in zip(refs[:2],('1','2')):
        i=rows['Term Source Name'].index(ref)
        assert rows['Term Source Version'][i]==version
        assert rows['Term Source File'][i]=='https://incoming/vocab'


def test_native_blank_term_references_do_not_declare_a_vocabulary():
    data=native().to_mapping();data['database'].append({'iid':'  ','name':'  '})
    node=data['series']['assay_paths'][0]['steps'][0]
    node['factor_values']=[{'name':'mixture','value':'mix','term_source_ref':'  ','term_accession_number':'\t'}]
    original=deepcopy(data)
    rows={r[0]:r[1:] for r in render(data)}
    assert all(str(v).strip() for v in rows['Term Source Name'])
    assert data==original
    result=finalize(data,[]).to_mapping()
    assert 'term_source_ref' not in result['series']['assay_paths'][0]['steps'][0]['factor_values'][0]


def test_native_non_rna_single_cell_libraries_do_not_assert_coding_rna_type():
    data=native().to_mapping();data['sample'][0]['title']='single-cell RNA sequencing and ADT libraries'
    data['sample'][0]['library_strategy']='OTHER'
    for run in data['sample'][0]['sra_run']:run['library_strategy']='OTHER'
    rows={r[0]:r[1:] for r in render(data)}
    assert not rows.get('Comment[AEExperimentType]')


def test_native_sequence_data_links_do_not_span_unrelated_run_accessions():
    data=native().to_mapping();template=data['sample'][0]['sra_run'][0]
    accessions=['ERR6054545','ERR6054546','ERR6286716','ERR6548408','ERR10123638']
    data['sample'][0]['sra_run']=[{**deepcopy(template),'run':acc} for acc in accessions]
    rows={r[0]:r[1:] for r in render(data)}
    assert rows['Comment[SequenceDataURI]']==['http://www.ebi.ac.uk/ena/data/view/'+s for s in
        ['ERR6054545-ERR6054546','ERR6286716','ERR6548408','ERR10123638']]
