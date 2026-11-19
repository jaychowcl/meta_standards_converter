# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""MAGE-TAB 1.1 section 2.3.3 places complete factor blocks last."""
from copy import deepcopy

import pytest

from meta_standards_converter.magetab.sdrf.model import SDRFAttr, SDRFNode, SDRFPath
from meta_standards_converter.magetab.sdrf.renderer import SDRFRenderer
from meta_standards_converter.magetab.semantics import render_miniml_assay_documents
from meta_standards_converter.magetab.validation import MAGETabValidationError, validate_magetab


@pytest.mark.parametrize('preserve_order', [False, True])
def test_terminal_projection_preserves_factor_subtrees_keys_and_input(preserve_order):
    paths = [SDRFPath([
        SDRFNode('Source Name', 's', 's1'),
        SDRFNode('Assay Name', 'a', 'a1', [
            SDRFAttr('Comment[node]', 'assay'),
            SDRFAttr('Factor Value[dose]', 0, [
                SDRFAttr('Term Source REF', 'EFO'),
                SDRFAttr('Comment[factor]', 'literal  '),
                SDRFAttr('Unit[time]', 'hour', [SDRFAttr('Term Source REF', 'UO')])]),
            SDRFAttr('Factor Value[dose]', '', required=True)]),
        SDRFNode('Array Data File', 'f', 'f1', [SDRFAttr('Comment[node]', 'file')]),
    ])]
    renderer = SDRFRenderer(preserve_order=preserve_order)
    columns = renderer.plan_columns(paths)
    before = deepcopy((paths, columns))
    assert renderer.render_paths(columns, paths) == [
        ['Source Name', 'Assay Name', 'Comment[node]', 'Array Data File', 'Comment[node]',
         'Factor Value[dose]', 'Term Source REF', 'Comment[factor]', 'Unit[time]', 'Term Source REF',
         'Factor Value[dose]'],
        ['s1', 'a1', 'assay', 'f1', 'file', 0, 'EFO', 'literal  ', 'hour', 'UO', ''],
    ]
    assert (paths, columns) == before


def test_typed_factor_comments_after_units_do_not_take_node_comments():
    paths = [{'document':'x.sdrf.txt', 'steps':[
        {'kind':'source', 'name':'s1'},
        {'kind':'assay', 'name':'a1', 'comments':[{'name':'detail','value':'assay'}],
         'factor_values':[{'name':'dose','value':'0', 'unit':{'value':'hour','term_source_ref':'UO'},
                           'comments':[{'name':'detail','value':'factor'}]}]},
        {'kind':'array_data_file','name':'f1','comments':[{'name':'detail','value':'file'}]},
    ]}]
    before = deepcopy(paths)
    table = render_miniml_assay_documents(paths)['x.sdrf.txt']
    assert table == [
        ['Source Name','Assay Name','Comment[detail]','Array Data File','Comment[detail]',
         'Factor Value[dose]','Comment[detail]','Unit','Term Source REF'],
        ['s1','a1','assay','f1','file','0','factor','hour','UO'],
    ]
    assert paths == before
    from tests.support.magetab import parse
    parsed = parse(table[0],table[1:]).to_mapping()['series']['assay_paths'][0]['steps']
    factor = next(v for step in parsed for v in step.get('factor_values', []))
    assert factor['comments'] == [{'name':'detail','value':'factor'}]
    assert next(s for s in parsed if s['kind']=='assay')['comments'] == [{'name':'detail','value':'assay'}]
    assert next(s for s in parsed if s['kind']=='array_data_file')['comments'] == [{'name':'detail','value':'file'}]


def test_no_factor_projection_is_unchanged():
    paths = [SDRFPath([SDRFNode('Source Name','s','s1',[
        SDRFAttr('Characteristics[dose]',0,[SDRFAttr('Unit','hour')])]),
        SDRFNode('Assay Name','a','a1'),SDRFNode('Array Data File','f','f1')])]
    renderer=SDRFRenderer()
    assert renderer.render_paths(renderer.plan_columns(paths),paths)==[
        ['Source Name','Characteristics[dose]','Unit','Assay Name','Array Data File'],
        ['s1',0,'hour','a1','f1']]


@pytest.mark.parametrize('following', ['Array Data File','Protocol REF','Characteristics[dose]','Technology Type'])
def test_validator_rejects_nonfactor_columns_after_terminal_factor_block(following):
    rows=[['Experimental Factor Name','dose'],['SDRF File',[
        ['Source Name','Assay Name','Factor Value[dose]',following],['s1','a1','0','']]]]
    with pytest.raises(MAGETabValidationError,match='Factor Value.*after'):
        validate_magetab(rows)


def test_validator_accepts_terminal_factor_qualifiers_and_comments():
    rows=[['Experimental Factor Name','dose'],['Term Source Name','UO'],['SDRF File',[
        ['Source Name','Assay Name','Factor Value[dose]','Unit[time]','Term Source REF','Comment[detail]'],
        ['s1','a1','0','hour','UO','zero']]]]
    assert validate_magetab(rows) is rows


@pytest.mark.parametrize('boundary', ['Factor Value[other]','Protocol REF','Assay Name','Parameter Value[time]'])
def test_factor_companion_scan_stops_at_next_owner(boundary):
    from meta_standards_converter.magetab.semantics import _assay_path
    headers=['Source Name','Factor Value[dose]','Comment[note]',boundary,'Unit']
    path=_assay_path('x',1,headers,['s1','0','factor comment','next','hour'])
    factor=next(s for s in path['steps'] if s.get('header')=='Factor Value[dose]')
    assert 'unit' not in factor
    assert next(s for s in path['steps'] if s.get('header')=='Comment[note]')['value']=='factor comment'


def test_factor_companion_scan_keeps_post_unit_comment_position_explicit():
    from meta_standards_converter.magetab.semantics import _assay_path
    headers=['Source Name','Factor Value[dose]','Comment[factor]','Unit','Comment[unit]','Factor Value[other]']
    path=_assay_path('x',1,headers,['s1','0','factor comment','hour','unit comment','other'])
    factor=next(s for s in path['steps'] if s.get('header')=='Factor Value[dose]')
    assert factor['unit']=='hour' and factor['companion_columns']['unit']==3
    # The source model retains both positions for typed-owner migration.
    comments=[s for s in path['steps'] if s['kind']=='comment']
    assert [(s['column_index'],s['value']) for s in comments]==[(2,'factor comment'),(4,'unit comment')]
