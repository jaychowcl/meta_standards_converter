# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Source-position and typed-owner contracts for field and unit comments."""
from copy import deepcopy

import pytest

from meta_standards_converter.miniml.model import OntologyValue, NamedComment, MINiMLModelError
from meta_standards_converter.miniml import MINiMLCodec
from meta_standards_converter.magetab.semantics import render_miniml_assay_documents
from tests.test_magetab_blank_nodes import parse
from tests.support.deposited_fidelity import assert_deposited_sdrf


def test_ontology_comments_roundtrip_and_old_positional_constructor():
    original={'value':'RNA','term_source_ref':'EFO','comments':[{'name':'note','value':' 0 '}]}
    value=OntologyValue.from_value(original)
    assert value.comments==(NamedComment('note',' 0 '),)
    assert value.to_mapping()==original
    assert OntologyValue('RNA','EFO',None,()).to_mapping()=={'value':'RNA','term_source_ref':'EFO'}
    original['comments'][0]['value']='changed'
    assert value.comments[0].value==' 0 '
    with pytest.raises(MINiMLModelError): OntologyValue('RNA',comments=('wrong',))
    with pytest.raises(MINiMLModelError): OntologyValue.from_value({'value':'RNA','comments':[{'name':'','value':'bad'}]})


@pytest.mark.parametrize('field', ['Material Type','Technology Type','Label'])
def test_scalar_field_comment_does_not_rebind_to_previous_named_attribute(field):
    header=['Source Name','Characteristics[dose]',field,'Comment[detail]','Assay Name']
    row=['s1','0','RNA','field note','a1']
    data=parse(header,[row])
    paths=data.to_mapping()['series']['assay_paths']
    source=paths[0]['steps'][0]
    key={'Material Type':'material_type','Technology Type':'technology_type','Label':'label'}[field]
    assert source[key]['comments']==[{'name':'detail','value':'field note'}]
    assert 'comments' not in source['characteristics'][0]
    output=render_miniml_assay_documents(paths)['x.sdrf.txt']
    assert_deposited_sdrf([header,row],output,accession='synthetic')
    assert MINiMLCodec().decode(data.to_mapping(),strict=True).package.to_mapping()['series']['assay_paths']==paths


@pytest.mark.parametrize('prefix,container', [('Characteristics','characteristics'),('Factor Value','factor_values'),('Parameter Value','parameter_values')])
def test_named_and_unit_comments_keep_repeated_occurrence_owners(prefix,container):
    leading=['Source Name','Protocol REF'] if prefix=='Parameter Value' else ['Source Name']
    header=[*leading,f'{prefix}[dose]','Comment[note]','Unit[time]','Comment[note]',f'{prefix}[dose]','Unit[time]','Comment[note]']
    row=['s1',*(['P1'] if len(leading)==2 else []),'0','factor note','hour','unit one','1','day','unit two']
    data=parse(header,[row])
    paths=data.to_mapping()['series']['assay_paths']
    values=next(s[container] for s in paths[0]['steps'] if container in s)
    assert values[0]['comments']==[{'name':'note','value':'factor note'}]
    assert values[0]['unit']['comments']==[{'name':'note','value':'unit one'}]
    assert values[1]['unit']['comments']==[{'name':'note','value':'unit two'}]
    output=render_miniml_assay_documents(paths)['x.sdrf.txt']
    assert_deposited_sdrf([header,row],output,accession='synthetic')


def test_empty_ontology_field_preserves_meaningful_comment():
    data=parse(['Source Name','Material Type','Comment[note]'],[['s1','','unknown material']])
    source=data.to_mapping()['series']['assay_paths'][0]['steps'][0]
    assert source['material_type']=={'value':'','comments':[{'name':'note','value':'unknown material'}]}
    assert not source.get('comments')


def test_blank_unit_comment_is_explicit_unbound_evidence_not_named_comment():
    data=parse(['Source Name','Characteristics[dose]','Unit[time]','Comment[note]'],[['s1','0','','unknown unit']])
    source=data.to_mapping()['series']['assay_paths'][0]['steps'][0]
    assert not source['characteristics'][0].get('comments')
    residual=data.to_mapping()['extensions']['magetab']['unbound_annotations']
    assert any(v.get('value')=='unknown unit' and v.get('column_index')==3 for v in residual)


def test_duplicate_scalar_field_and_its_comment_are_explicit_residuals():
    data=parse(['Source Name','Material Type','Comment[note]','Material Type','Comment[note]'],
               [['s1','RNA','first','DNA','second']])
    source=data.to_mapping()['series']['assay_paths'][0]['steps'][0]
    assert source['material_type']=={'value':'RNA','comments':[{'name':'note','value':'first'}]}
    residual=data.to_mapping()['extensions']['magetab']['unbound_annotations']
    assert {(v['column_index'],v['value']) for v in residual}>={(3,'DNA'),(4,'second')}


def test_duplicate_unit_and_comment_are_not_overwritten_or_rebound():
    data=parse(['Source Name','Characteristics[dose]','Unit[time]','Unit[time]','Comment[note]'],
               [['s1','0','hour','day','second unit']])
    source=data.to_mapping()['series']['assay_paths'][0]['steps'][0]
    assert source['characteristics'][0]['unit']=={'value':'hour'}
    assert not source.get('comments')
    residual=data.to_mapping()['extensions']['magetab']['unbound_annotations']
    assert {(v['column_index'],v['value']) for v in residual}>={(3,'day'),(4,'second unit')}


def test_merge_retains_typed_field_comments_without_mutating_input():
    from meta_standards_converter.metadata.archive_enrichment import merge_archive_metadata
    data=parse(['Source Name','Material Type','Comment[note]','Assay Name'],[['s1','RNA','literal','a1']])
    package=MINiMLCodec().decode(data.to_mapping(),strict=True).package
    before=deepcopy(package.to_mapping())
    merged,_=merge_archive_metadata(package,package,prefer=True)
    source=merged.to_mapping()['series']['assay_paths'][0]['steps'][0]
    assert source['material_type']['comments']==[{'name':'note','value':'literal'}]
    assert package.to_mapping()==before


@pytest.mark.parametrize('existing', [[],[{'name':'note','value':'native'}]])
def test_merge_completes_comments_only_on_same_identity_field(existing):
    from meta_standards_converter.metadata.archive_enrichment import merge_archive_metadata
    base=parse(['Source Name','Extract Name','Material Type','Assay Name'],[['GSM1','e1','RNA','ERX1']]).to_mapping()
    base['database'].append({'iid':'GEO','name':'GEO'})
    base['series']['assay_paths'][0]['steps'][1]['material_type']['comments']=existing
    peer=deepcopy(base)
    peer['series']['assay_paths'][0]['steps'][1]['material_type']['comments']=[{'name':'note','value':'linked'}]
    codec=MINiMLCodec()
    left,right=codec.decode(base,strict=True).package,codec.decode(peer,strict=True).package
    before=(left.to_mapping(),right.to_mapping())
    merged,_=merge_archive_metadata(left,right,prefer=True)
    fields=[s['material_type'] for p in merged.to_mapping()['series']['assay_paths'] for s in p['steps'] if s['kind']=='extract']
    assert fields and all(v['comments']==[{'name':'note','value':'linked'},*existing] for v in fields)
    assert (left.to_mapping(),right.to_mapping())==before


def test_verified_join_preserves_max_comment_occurrences_and_input_isolation():
    from meta_standards_converter.metadata.archive_workflows import _bind_native
    note={'name':'note','value':''}
    preferred=[{'kind':'extract','name':'e1','material_type':{'value':'RNA','comments':[note]}}]
    native=[{'kind':'extract','name':'e1','material_type':{'value':'RNA','comments':[note,note,{'name':'note','value':'native'}]}}]
    before=deepcopy((preferred,native))
    result=_bind_native(deepcopy(preferred),'GSM1',native)
    assert result[0]['material_type']['comments']==[note,note,{'name':'note','value':'native'}]
    result[0]['material_type']['comments'][-1]['value']='changed'
    assert (preferred,native)==before


@pytest.mark.parametrize('mismatch',['value','term_source_ref','node-name','ambiguous-native','ambiguous-incoming'])
def test_verified_join_does_not_complete_unmatched_field_comments(mismatch):
    from meta_standards_converter.metadata.archive_workflows import _bind_native
    preferred=[{'kind':'extract','name':'e1','material_type':{'value':'RNA','term_source_ref':'EFO'}}]
    native=deepcopy(preferred)
    native[0]['material_type']['comments']=[{'name':'note','value':'native'}]
    if mismatch=='value': native[0]['material_type']['value']='DNA'
    elif mismatch=='term_source_ref': native[0]['material_type']['term_source_ref']='OTHER'
    elif mismatch=='node-name': native[0]['name']='e2'
    elif mismatch=='ambiguous-native': native.append(deepcopy(native[0]))
    elif mismatch=='ambiguous-incoming': preferred.append(deepcopy(preferred[0]))
    result=_bind_native(preferred,'GSM1',native)
    assert all('comments' not in step['material_type'] for step in result)


def test_scalar_field_comments_remain_visible_to_identity_and_library_consumers():
    from meta_standards_converter.metadata.archive_workflows import path_ids
    from meta_standards_converter.miniml.archive_libraries import path_facts,apply_library_facts
    steps=[{'kind':'extract','name':'e','material_type':{'value':'RNA','comments':[{'name':'LIBRARY_LAYOUT','value':'PAIRED'}]}},
           {'kind':'assay','name':'a','technology_type':{'value':'sequencing assay','comments':[{'name':'ENA_EXPERIMENT','value':'ERX1'}]}}]
    assert path_ids({'steps':steps})=={'ERX1'}
    assert path_facts(steps)==[{'library_layout':'PAIRED'}]
    apply_library_facts({},steps,{'library_layout':'SINGLE'})
    assert steps[0]['material_type']['comments']==[{'name':'LIBRARY_LAYOUT','value':'SINGLE'}]
    assert not steps[0].get('comments')
    assert path_ids({'steps':steps})=={'ERX1'}


@pytest.mark.parametrize('nested', ['characteristics','factor_values'])
def test_named_and_unit_descendant_comments_do_not_supply_node_identity_or_library(nested):
    from meta_standards_converter.metadata.archive_workflows import path_ids
    from meta_standards_converter.miniml.archive_libraries import path_facts
    comments=[{'name':'ENA_EXPERIMENT','value':'ERX999'},{'name':'LIBRARY_LAYOUT','value':'SINGLE'}]
    step={'kind':'assay','name':'a',nested:[{'name':'dose','value':'0','comments':comments,'unit':{'value':'hour','comments':comments}}]}
    assert path_ids({'steps':[step]})==set()
    assert path_facts([step])==[]


def test_channel_unit_matches_path_with_intervening_named_comment():
    data=parse(['Source Name','Characteristics[dose]','Comment[note]','Unit[time]','Term Source REF'],
               [['s1','0','named note','hour','UO']]).to_mapping()
    channel=data['sample'][0]['channel'][0]['characteristics'][0]
    named=data['series']['assay_paths'][0]['steps'][0]['characteristics'][0]
    assert channel['unit']==named['unit']=={'value':'hour','term_source_ref':'UO'}
    assert channel['unit_type']==named['unit_type']=='time'


@pytest.mark.parametrize('boundary',['Characteristics[other]','Factor Value[other]','Parameter Value[other]','Protocol REF','Assay Name'])
def test_channel_companion_scan_does_not_cross_owner_boundaries(boundary):
    data=parse(['Source Name','Characteristics[dose]','Comment[note]',boundary,'Unit'],
               [['s1','0','named note','other','hour']]).to_mapping()
    dose=next(v for v in data['sample'][0]['channel'][0]['characteristics'] if v['name']=='dose')
    assert 'unit' not in dose


def test_channel_repeated_unit_does_not_overwrite_first_occurrence():
    data=parse(['Source Name','Characteristics[dose]','Unit','Term Source REF','Unit','Term Source REF'],
               [['s1','0','hour','UO','day','OTHER']]).to_mapping()
    assert data['sample'][0]['channel'][0]['characteristics'][0]['unit']=={'value':'hour','term_source_ref':'UO'}
