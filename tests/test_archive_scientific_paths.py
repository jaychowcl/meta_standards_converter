# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from copy import deepcopy
from meta_standards_converter.metadata.archive_enrichment import merge_archive_metadata
from meta_standards_converter.miniml import MINiMLCodec
from tests.test_native_archive_enrichment import native, linked
from tests.test_protocol_export import render


def rich():
    package=native();extra=linked(package,'GSE1','explicit title').to_mapping()
    extra['series']['assay_paths']=[]
    extra['sample'][0]['channel'][0].update(source={'value':'cortex'},molecule={'value':'polyA RNA'},extract_protocol='Purified RNA and prepared libraries.')
    extra['sample'][0]['data_processing']='Mapped and counted reads.'
    extra['sample'][0]['supplementary_data']=[{'value':'https://example.org/counts.tsv'}]
    return package,MINiMLCodec().decode(extra).package


def test_geo_fields_complete_scoped_acquisition_and_result_paths():
    package,extra=rich();result,issues=merge_archive_metadata(package,extra,prefer=True);assert not issues
    data=result.to_mapping();protocols={p['name']:p for p in data['series']['protocols']}
    for path in data['series']['assay_paths']:
        steps=path['steps'];source=steps[0]
        assert {'name':'Sample_title','value':'explicit title'} in source['comments']
        assert {'name':'Sample_source_name','value':'cortex'} in source['comments']
        applied=[protocols[s['protocol_ref']]['description'] for s in steps if s.get('protocol_ref') in protocols]
        if any(s['kind']=='assay' for s in steps):
            extract=next(s for s in steps if s['kind']=='extract')
            assert extract['material_type']=={'value':'polyA RNA'}
            assert steps.index(extract)<next(i for i,s in enumerate(steps) if s['kind']=='assay')
            assert 'Purified RNA and prepared libraries.' in applied
            assert 'Mapped and counted reads.' not in applied
        if any(s['kind']=='derived_array_data_file' for s in steps):
            assert 'Mapped and counted reads.' in applied
            assert not any(s['kind']=='scan' for s in steps)


def test_saved_native_projection_is_idempotent_and_study_files_are_scoped():
    package,extra=rich();result,_=merge_archive_metadata(package,extra,prefer=True)
    data=result.to_mapping();data['series']['supplementary_data']=[{'value':'https://example.org/study-matrix.txt','type':'TXT'}]
    before=deepcopy(data);rows=render(data)
    values={r[0]:r[1:] for r in rows}
    assert values['Comment[Study supplementary file]']==['https://example.org/study-matrix.txt']
    assert values['Comment[Study supplementary file type]']==['TXT']
    assert 'https://example.org/study-matrix.txt' not in str(values['SDRF File'])
    assert rows==render(data)
    assert data==before


def test_later_geo_updates_generated_material_and_protocol_together():
    package,extra=rich();first,_=merge_archive_metadata(package,extra,prefer=True)
    update=extra.to_mapping();update['sample'][0]['channel'][0].update(molecule={'value':'total RNA'},extract_protocol='New extraction.')
    result,issues=merge_archive_metadata(first,MINiMLCodec().decode(update).package,prefer=True);assert not issues
    data=result.to_mapping();protocols={p['name']:p for p in data['series']['protocols']}
    for p in data['series']['assay_paths']:
        if any(s['kind']=='assay' for s in p['steps']):
            assert next(s for s in p['steps'] if s['kind']=='extract')['material_type']['value']=='total RNA'
            assert any(protocols.get(s.get('protocol_ref'),{}).get('description')=='New extraction.' for s in p['steps'])


def test_partial_ae_workflow_does_not_supply_material_to_other_experiment():
    from tests.test_archive_fidelity_enrichment import workflow, enriched_native
    package=enriched_native().to_mapping();second=deepcopy(package['series']['assay_paths'][0])
    for s in second['steps']:
        if s['kind']=='assay':s['name']='SRX99'
        if s['kind']=='scan':s['name']='SRR99';s['comments']=[]
    package['series']['assay_paths'].append(second)
    run=deepcopy(package['sample'][0]['sra_run'][0]);run.update(run='SRR99',experiment='SRX99');package['sample'][0]['sra_run'].append(run)
    extra=workflow().to_mapping();extra['sample'][0]['channel'][0].update(molecule={'value':'total RNA'},extract_protocol='First experiment only')
    result,issues=merge_archive_metadata(MINiMLCodec().decode(package).package,MINiMLCodec().decode(extra).package,prefer=True);assert not issues
    for p in result.to_mapping()['series']['assay_paths']:
        if any(s.get('name')=='SRX99' for s in p['steps']):assert not any(s['kind']=='extract' for s in p['steps'])


def test_explicit_scan_protocol_is_retained_at_acquisition_scope():
    package,extra=rich();extra=extra.to_mapping();extra['sample'][0]['scan_protocol']='Sequenced on supplied instrument.'
    result,issues=merge_archive_metadata(package,MINiMLCodec().decode(extra).package,prefer=True);assert not issues
    data=result.to_mapping();protocols={p['name']:p for p in data['series']['protocols']}
    for p in data['series']['assay_paths']:
        if any(s['kind']=='scan' for s in p['steps']):
            assert any(protocols.get(s.get('protocol_ref'),{}).get('description')=='Sequenced on supplied instrument.' for s in p['steps'])


def test_partial_ae_workflow_does_not_project_sample_protocols_onto_other_experiment():
    from tests.test_archive_fidelity_enrichment import workflow, enriched_native
    package=enriched_native().to_mapping();second=deepcopy(package['series']['assay_paths'][0])
    for s in second['steps']:
        if s['kind']=='assay':s['name']='SRX99'
        if s['kind']=='scan':s['name']='SRR99';s['comments']=[]
    second['steps'].append({'kind':'derived_array_data_file','name':'second.tsv'})
    package['series']['assay_paths'].append(second)
    run=deepcopy(package['sample'][0]['sra_run'][0]);run.update(run='SRR99',experiment='SRX99');package['sample'][0]['sra_run'].append(run)
    extra=workflow().to_mapping();extra['sample'][0].update(scan_protocol='First experiment sequencing only',data_processing='First experiment processing only')
    result,issues=merge_archive_metadata(MINiMLCodec().decode(package).package,MINiMLCodec().decode(extra).package,prefer=True);assert not issues
    data=result.to_mapping();protocols={p['name']:p for p in data['series']['protocols']}
    for path in data['series']['assay_paths']:
        if any(s.get('name')=='SRX99' for s in path['steps']):
            assert not any('First experiment' in protocols.get(s.get('protocol_ref'),{}).get('description','') for s in path['steps'])


def test_authored_extraction_and_scanning_are_not_projected_as_duplicate_methods():
    from meta_standards_converter.miniml.archive_paths import complete_native_paths
    data=native().to_mapping();sid=data['sample'][0]['iid'];steps=data['series']['assay_paths'][0]['steps']
    definitions=[{'name':'P-MTAB-1','type':{'value':'nucleic_acid_extraction'},'description':'Extract RNA. '},
                 {'name':'P-MTAB-2','type':{'value':'scanning'},'description':'Call bases.'}]
    data['series']['assay_paths']=[data['series']['assay_paths'][0]]
    data['series'].setdefault('protocols',[]).extend(definitions)
    steps[1:1]=[{'kind':'protocol_application','protocol_ref':'P-MTAB-1'},{'kind':'extract','name':'authored','material_type':{'value':'RNA'},'sample_ref':sid},
                {'kind':'protocol_application','protocol_ref':'P-MTAB-2'}]
    data['sample'][0]['channel'][0]['extract_protocol']='Extract RNA.'
    data['sample'][0]['scan_protocol']='Call bases.'
    before=deepcopy(data['series']['protocols'])
    complete_native_paths(data)
    assert data['series']['protocols']==before
    assert sum(s.get('protocol_ref')=='P-MTAB-2' for s in steps)==1
    # Previously saved generated copies are safely removed only when authored
    # applications already represent them on all relevant sample acquisitions.
    duplicate=data['series']['iid']+':'+sid+':scan_protocol'
    data['series']['protocols'].append({'name':duplicate,'type':{'value':'nucleic acid sequencing protocol'},'description':'Call bases.'})
    steps.insert(next(i for i,s in enumerate(steps) if s['kind']=='scan'),{'kind':'protocol_application','protocol_ref':duplicate})
    complete_native_paths(data)
    assert data['series']['protocols']==before
    assert not any(s.get('protocol_ref')==duplicate for s in steps)


def test_authored_processing_does_not_hide_uncovered_derived_matrix_branch():
    from meta_standards_converter.miniml.archive_paths import complete_native_paths
    data=native().to_mapping();sid=data['sample'][0]['iid']
    data['sample'][0]['data_processing']='Count reads.'
    data['series'].setdefault('protocols',[]).append({'name':'P-MTAB-1','type':{'value':'normalization'},'description':'Count reads.'})
    data['series']['assay_paths']=[
        {'steps':[{'kind':'source','name':sid,'sample_ref':sid},
                  {'kind':'protocol_application','protocol_ref':'P-MTAB-1'},
                  {'kind':'derived_array_data_file','name':'counts.tsv'}]},
        {'steps':[{'kind':'source','name':sid,'sample_ref':sid},
                  {'kind':'derived_array_data_matrix_file','name':'counts.mtx'}]}]
    complete_native_paths(data)
    methods={p['name']:p for p in data['series']['protocols']}
    for path in data['series']['assay_paths']:
        applications=[s['protocol_ref'] for s in path['steps'] if s.get('protocol_ref')]
        assert len(applications)==1
        assert methods[applications[0]]['description']=='Count reads.'
    assert data['series']['assay_paths'][0]['steps'][1]['protocol_ref']=='P-MTAB-1'
    before=deepcopy(data);complete_native_paths(data);assert data==before


def test_native_paths_export_all_explicit_biosample_aliases_without_changing_identity():
    from meta_standards_converter.miniml.archive_paths import complete_native_paths
    data=native().to_mapping();sample=data['sample'][0];sid=sample['iid']
    sample['accession']=[{'value':sid,'database':'SRA'},
                         {'value':'SAMN123','database':'BioSample'},
                         {'value':'SAMEA456','database':'BioSample'}]
    paths=data['series']['assay_paths'];source=paths[0]['steps'][0]
    source.setdefault('comments',[]).append({'name':'BioSD_SAMPLE','value':'SAMN789'})
    expected={c['value'] for c in source['comments'] if c['name']=='BioSD_SAMPLE'} | {'SAMN123','SAMEA456'}
    pooled=deepcopy(paths[0]);pooled['steps'].insert(1,{'kind':'sample','name':'different','sample_ref':'different'})
    paths.append(pooled);original_pool=deepcopy(pooled)
    complete_native_paths(data)
    assert {c['value'] for c in source['comments'] if c['name']=='BioSD_SAMPLE'}==expected
    assert source['name']==sid and source['sample_ref']==sid
    assert pooled==original_pool
    paths.pop();before=deepcopy(data);complete_native_paths(data);assert data==before
    # Rendering an existing native package projects aliases on its private copy.
    for path in paths:
        path['steps'][0]['comments']=[c for c in path['steps'][0].get('comments',[]) if c['name']!='BioSD_SAMPLE']
    before=deepcopy(data)
    rows=next(r[1] for r in render(data) if r[0]=='SDRF File')
    values={r[i] for r in rows[1:] for i,h in enumerate(rows[0]) if h=='Comment[BioSD_SAMPLE]'}
    assert {'SAMN123','SAMEA456'}<=values
    assert data==before
