# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from copy import deepcopy
from meta_standards_converter.magetab.parser import AEParser
from meta_standards_converter.miniml import MINiMLCodec
from meta_standards_converter.metadata.archive_enrichment import merge_archive_metadata
from tests.converters.test_ae2json import resolved_input
from tests.test_native_archive_enrichment import native


def workflow():
    header=['Source Name','Characteristics[organism]','Characteristics[Age]','Unit[TimeUnit]','Term Source REF','Term Accession Number','Material Type','Protocol REF','Extract Name','Material Type','Protocol REF','Extract Name','Material Type','Assay Name','Comment[ENA_EXPERIMENT]','Scan Name','Comment[ENA_RUN]','Factor Value[AGE]','Unit[TimeUnit]']
    row=['embryo','Danio rerio','1','days','UO','UO:0000033','organism part','P1','RNA','total RNA','P2','mRNA','3 prime end of mRNA','assay','SRX7812918','scan','SRR11192680','1','d']
    idf='Investigation Accession\tE-MTAB-1\nInvestigation Title\tStudy\nComment[SecondaryAccession]\tSRP250911\nExperimental Factor Name\tAGE\nExperimental Factor Type\tage\nProtocol Name\tP1\tP2\nProtocol Type\textraction\tlibrary construction\nProtocol Description\textract RNA\tselect mRNA\nPerson First Name\tJane\nPerson Last Name\tDoe\nPerson Email\tjane@example.org\n'
    return AEParser().parse(resolved_input(idf=idf,sdrfs=['\t'.join(header)+'\n'+'\t'.join(row)+'\n']))


def test_source_only_rows_bind_samples_and_keep_scoped_units_materials():
    data=workflow().to_mapping(); sample=data['sample'][0]; steps=data['series']['assay_paths'][0]['steps']
    assert {s.get('sample_ref') for s in steps if s.get('sample_ref')}=={sample['iid']}
    age=next(a for a in sample['channel'][0]['characteristics'] if a['name']=='Age')
    assert age['unit']=={'value':'days','term_source_ref':'UO','term_accession_number':'UO:0000033'}
    assert age==next(a for a in steps[0]['characteristics'] if a['name']=='Age')
    assert not any(a['value']=='3 prime end of mRNA' for a in sample['channel'][0]['characteristics'])
    assert [s['material_type']['value'] for s in steps if s['kind']=='extract']==['total RNA','3 prime end of mRNA']


def enriched_native():
    data=native().to_mapping()
    data['sample'][0]['channel'][0]['organism']=[{'value':'Danio rerio','taxid':'7955'}]
    return MINiMLCodec().decode(data).package


def test_complete_workflows_factors_people_taxonomy_and_native_files():
    original=enriched_native(); linked=workflow()
    data,issues=merge_archive_metadata(original,linked,prefer=True)
    data=data.to_mapping()
    assert not issues
    assert data['contributor'] and all(p['iid'].startswith('E-MTAB-1:') for p in data['contributor'])
    assert data['series']['contributor_ref']
    assert data['sample'][0]['channel'][0]['organism']==[{'value':'Danio rerio','taxid':'7955'}]
    assert data['sample'][0]['iid']==original.samples[0].iid
    original_links=[s['link'] for p in original.to_mapping()['series']['assay_paths'] for s in p['steps'] if s.get('link')]
    new_links=[s['link'] for p in data['series']['assay_paths'] for s in p['steps'] if s.get('link')]
    assert all(link in new_links for link in original_links)
    for path in data['series']['assay_paths']:
        steps=path['steps']
        assert len([s for s in steps if s['kind']=='extract'])==2
        assert any(s.get('factor_values') for s in steps)
        assert {s['protocol_ref'] for s in steps if s.get('protocol_ref')}=={'E-MTAB-1:P1','E-MTAB-1:P2'}
        assert next(s['name'] for s in steps if s['kind']=='scan')=='SRR11192680'


def test_missing_refs_recover_only_by_verified_run_and_database_ids_merge():
    package=enriched_native().to_mapping(); package['database'].append({'iid':'GEO','name':'GEO'})
    extra=workflow().to_mapping(); extra['database'].append({'iid':'GEO','name':'Gene Expression Omnibus'})
    for path in extra['series']['assay_paths']:
        for step in path['steps']: step.pop('sample_ref',None)
    data,issues=merge_archive_metadata(MINiMLCodec().decode(package).package,MINiMLCodec().decode(extra).package,prefer=True)
    assert not issues
    assert len([d for d in data.to_mapping()['database'] if d['iid']=='GEO'])==1
    assert any(s['kind']=='extract' for p in data.to_mapping()['series']['assay_paths'] for s in p['steps'])


def test_sample_processed_file_does_not_acquire_a_run():
    extra=workflow().to_mapping(); extra['sample'][0]['supplementary_data']=[{'value':'https://example.org/counts.tsv'}]
    data,_=merge_archive_metadata(enriched_native(),MINiMLCodec().decode(extra).package,prefer=True)
    paths=[p for p in data.to_mapping()['series']['assay_paths'] if any(s.get('link',{}).get('value')=='https://example.org/counts.tsv' for s in p['steps'])]
    assert paths and all(not any(s['kind'] in ('scan','assay') for s in p['steps']) for p in paths)


def test_archive_broker_header_spellings_are_semantic():
    source=resolved_input(idf='Investigation Accession\tE-MTAB-1\n',sdrfs=['Source Name\tAssay Name\tComment [ENA_EXPERIMENT]\tScan Name\tComment [ENA_RUN]\tFactorValue [age]\tUnit [TimeUnit]\na\tassay\tSRX1\tscan\tSRR1\t2\tdays\n'])
    data=AEParser().parse(source).to_mapping()
    assert data['sample'][0]['sra_run'][0]['experiment']=='SRX1'
    scan=next(s for s in data['series']['assay_paths'][0]['steps'] if s['kind']=='scan')
    assert scan['factor_values'][0]['unit']['value']=='days'


def test_fastq_alternatives_do_not_make_identical_workflows_ambiguous():
    extra=workflow().to_mapping(); path=extra['series']['assay_paths'][0]
    path['steps'][-1].setdefault('comments',[]).append({'name':'FASTQ_URI','value':'ftp://example.org/one.fastq'})
    second=deepcopy(path);second['steps'][-1]['comments'][-1]['value']='ftp://example.org/two.fastq'
    extra['series']['assay_paths'].append(second)
    data,issues=merge_archive_metadata(enriched_native(),MINiMLCodec().decode(extra).package,prefer=True)
    assert not issues
    assert any(s['kind']=='extract' for p in data.to_mapping()['series']['assay_paths'] for s in p['steps'])


def test_conflicting_organism_does_not_keep_old_taxid_and_local_refs_are_scoped():
    extra=workflow().to_mapping()
    extra['sample'][0]['channel'][0]['organism']=[{'value':'Mus musculus'}]
    extra['database'].append({'iid':'LOCAL','name':'incoming vocabulary'})
    extra['sample'][0]['channel'][0]['characteristics'].append({'name':'condition','value':'x','term_source_ref':'LOCAL'})
    extra['organization']=[{'iid':'org1','name':'Lab'}]
    extra['contributor'][0]['organization_ref']={'ref':'org1'}
    package=enriched_native().to_mapping();package['database'].append({'iid':'LOCAL','name':'different vocabulary'})
    data,_=merge_archive_metadata(MINiMLCodec().decode(package).package,MINiMLCodec().decode(extra).package,prefer=True)
    data=data.to_mapping()
    assert not data['sample'][0]['channel'][0]['organism'][0].get('taxid')
    assert next(c for c in data['sample'][0]['channel'][0]['characteristics'] if c['name']=='condition')['term_source_ref']=='E-MTAB-1:LOCAL'
    assert data['organization'][0]['iid']=='E-MTAB-1:org1'
    assert data['contributor'][0]['organization_ref']=={'ref':'E-MTAB-1:org1'}


def test_mismatched_run_workflow_is_reported_and_retained():
    extra=workflow().to_mapping()
    for step in extra['series']['assay_paths'][0]['steps']:
        for comment in step.get('comments',[]):
            if comment['name']=='ENA_RUN': comment['value']='SRR999'
    data,issues=merge_archive_metadata(enriched_native(),MINiMLCodec().decode(extra).package,prefer=True)
    assert any('scope' in i for i in issues)
    assert not any(s['kind']=='extract' for p in data.to_mapping()['series']['assay_paths'] for s in p['steps'])
    assert any(r['kind']=='MINiML' for r in data.to_mapping()['extensions']['insdc']['records'])
