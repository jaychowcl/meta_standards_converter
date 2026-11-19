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
from meta_standards_converter.miniml.archive_entities import declare_ontologies
from tests.support.archive import native
from tests.support.magetab import render


from meta_standards_converter.miniml.sra_parser import SRAParser


from tests.support.archive import fixture_records

from meta_standards_converter.miniml.archive_residuals import finalize


@pytest.mark.parametrize('uri,expected',[
 ('http://purl.org/obo/owl/NCBITaxon#NCBITaxon_4932','NCBITaxon'),
 ('https://purl.org/obo/owl/UBERON#UBERON_0000955','UBERON'),
 ('http://purlXorg/obo/owl/NCBITaxon#NCBITaxon_4932','EFO'),
 ('http://purl.org/obo/owl/NCBITaxon#UBERON_4932','EFO'),
 ('http://purlXobolibraryYorg/obo/NCBITaxon_4932','EFO')])
def test_only_explicit_matching_ontology_namespaces_are_corrected(uri,expected):
    data=native().to_mapping();group={'name':'organism','value':'yeast','term_source_ref':'EFO','term_accession_number':uri}
    data['sample'][0]['channel'][0]['characteristics']=[deepcopy(group)]
    declare_ontologies(data)
    assert data['sample'][0]['channel'][0]['characteristics'][0]['term_source_ref']==expected
    assert data['sample'][0]['channel'][0]['characteristics'][0]['term_accession_number']==uri
    if expected!='EFO':assert any(r['kind']=='annotation' and r['metadata']['annotation']==group for r in data['extensions']['insdc']['records'])
    before=deepcopy(data);declare_ontologies(data);assert data==before


def test_saved_native_export_repairs_legacy_term_source_without_input_mutation():
    data=native().to_mapping();uri='http://purl.org/obo/owl/NCBITaxon#NCBITaxon_4932'
    group={'name':'organism','value':'yeast','term_source_ref':'EFO','term_accession_number':uri}
    for path in data['series']['assay_paths']:path['steps'][0]['characteristics']=[deepcopy(group)]
    before=deepcopy(data);output=render(data);assert data==before
    table=next(r[1] for r in output if r[0]=='SDRF File')
    for row in table[1:]:
        i=row.index(uri);assert table[0][i]=='Term Accession Number';assert row[i-1]=='NCBITaxon'


def test_native_ontology_references_have_declarations():
    data = SRAParser().parse(fixture_records('sra')).to_mapping()
    assert 'NCBITaxon' not in {d['iid'] for d in data['database']}
    from tests.support.magetab import render
    assert 'NCBITaxon' in {r[0]: r[1:] for r in render(data)}['Term Source Name']


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


def test_local_vocabulary_name_survives_repository_separation():
    data=native().to_mapping();data['database'].append({'iid':'local-vocabulary','name':'Supplied vocabulary label'})
    data['sample'][0]['channel'][0]['characteristics'].append({'name':'condition','value':'x','term_source_ref':'local-vocabulary'})
    data=finalize(data,[]).to_mapping()
    assert any(r['kind']=='term_source_declaration' and r['metadata']['name']=='Supplied vocabulary label' for r in data['extensions']['insdc']['records'])


def test_saved_native_vocabulary_conflict_is_namespaced_during_enrichment():
    from tests.support.archive import linked
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
    from tests.support.archive import linked
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
    from tests.support.archive import linked
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
