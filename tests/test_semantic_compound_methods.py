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
from tests.test_native_archive_enrichment import native
from tests.test_protocol_export import render
from meta_standards_converter.miniml.archive_paths import complete_native_paths


TEXT = ('Cord blood mononuclear cells (CBMCs) were isolated from cord blood as described in Breton et al. 2015. '
        'Library was constructed as described in the protocol of the manufacturer (10x Genomics Single Cell) '
        'or for Dropseq (Macosko et al., 2015) with the modifications described in Stoeckius et al. 2017')


def compound():
    data=native().to_mapping(); sample=data['sample'][0]
    experiment=sample['sra_run'][0]['experiment'];ref=experiment+':library'
    data['series']['protocols']=[{'name':ref,'type':{'value':'library construction protocol'},'description':TEXT}]
    sample['channel'][0].update(extract_protocol=TEXT,molecule={'value':'polyA RNA'})
    for path in data['series']['assay_paths']:
        path['steps'].insert(1,{'kind':'protocol_application','protocol_ref':ref})
    return data,ref


def test_compound_method_is_one_application_before_generated_material():
    data,ref=compound();complete_native_paths(data)
    assert len(data['series']['protocols'])==2  # both supplied type-specific definitions
    for path in data['series']['assay_paths']:
        steps=path['steps'];apps=[s for s in steps if s['kind']=='protocol_application']
        assert apps==[{'kind':'protocol_application','protocol_ref':ref}]
        assert steps[steps.index(apps[0])+1]['kind']=='extract'
    before=deepcopy(data);complete_native_paths(data);assert data==before


def test_saved_compound_projection_is_repaired_on_export_copy():
    data,ref=compound();sid=data['sample'][0]['iid'];generated=data['series']['iid']+':'+sid+':extract_protocol'
    data['series']['protocols'].append({'name':generated,'type':{'value':'nucleic acid extraction protocol'},'description':TEXT})
    for p in data['series']['assay_paths']:
        p['steps'][1:1]=[{'kind':'protocol_application','protocol_ref':generated},
                        {'kind':'extract','name':sid+':extract','sample_ref':sid,'material_type':{'value':'polyA RNA'}}]
    before=deepcopy(data);idf=render(data);sdrf=next(r[1] for r in idf if r[0]=='SDRF File')
    for row in sdrf[1:]:
        assert sum(bool(row[i]) for i,h in enumerate(sdrf[0]) if h=='Protocol REF')==1
    assert data==before and idf==render(data)


@pytest.mark.parametrize('case',['authored_material','decorated','ambiguous_assay','ambiguous_application','wrong_run','different_text'])
def test_compound_method_requires_native_identity_and_generated_boundary(case):
    data,ref=compound();path=data['series']['assay_paths'][0];data['series']['assay_paths']=[path];steps=path['steps']
    if case=='authored_material':
        steps[1:1]=[{'kind':'protocol_application','protocol_ref':ref},
                    {'kind':'extract','name':'authored-extract','sample_ref':data['sample'][0]['iid']}]
    elif case=='decorated':steps[1]['comments']=[{'name':'performed by','value':'A'}]
    elif case=='ambiguous_assay':steps.insert(2,deepcopy(next(s for s in steps if s['kind']=='assay')))
    elif case=='ambiguous_application':steps.insert(2,deepcopy(steps[1]))
    elif case=='wrong_run':next(s for s in steps if s['kind']=='scan')['name']='SRR999'
    else:data['sample'][0]['channel'][0]['extract_protocol']='Extract RNA independently.'
    original=deepcopy([s for s in steps if s.get('protocol_ref')==ref]);complete_native_paths(data)
    assert [s for s in steps if s.get('protocol_ref')==ref]==original
    assert len([s for s in steps if s.get('protocol_ref')])>=2


def test_native_method_does_not_cross_an_authored_treatment():
    data,ref=compound();data['series']['protocols'].append({'name':'P-MTAB-99','type':{'value':'treatment protocol'},'description':'Stimulate cells.'})
    for p in data['series']['assay_paths']:p['steps'].insert(1,{'kind':'protocol_application','protocol_ref':'P-MTAB-99'})
    complete_native_paths(data)
    for p in data['series']['assay_paths']:
        refs=[s['protocol_ref'] for s in p['steps'] if s.get('protocol_ref')]
        assert refs.index('P-MTAB-99')<refs.index(ref)


@pytest.mark.parametrize('incoming', ['Use kit,  then extract 20 ug RNA.', 'Use kit, then extract 20 µg RNA.', 'Use kit, then extract 20 μg RNA.'])
def test_compound_projection_uses_permitted_description_normalization(incoming):
    data,ref=compound()
    data['series']['protocols'][0]['description']='Use kit, then extract 20 ug RNA.'
    data['sample'][0]['channel'][0]['extract_protocol']=incoming
    complete_native_paths(data)
    assert all([s['protocol_ref'] for s in p['steps'] if s.get('protocol_ref')]==[ref] for p in data['series']['assay_paths'])
    assert data['series']['protocols'][0]['description']=='Use kit, then extract 20 ug RNA.'
