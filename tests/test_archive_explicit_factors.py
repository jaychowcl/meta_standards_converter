# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from copy import deepcopy
import xml.etree.ElementTree as ET
import pytest
from meta_standards_converter.miniml import MINiMLCodec
from meta_standards_converter.miniml.ena_parser import ENAParser
from meta_standards_converter.miniml.sra_parser import SRAParser
from tests.test_native_archive_parsers import fixture_records
from tests.test_protocol_export import render
from meta_standards_converter.miniml.archive_residuals import source_records


def records(provider):
    result=fixture_records(provider)
    experiment=next(n for root in result.xml for n in root.iter('EXPERIMENT'))
    attrs=ET.SubElement(experiment,'EXPERIMENT_ATTRIBUTES')
    for tag,value,unit in [('Experimental Factor: AGE','2 d',None),
                            ('Experimental Factor: Custom dose','5','mg'),
                            ('Experimental Factor: Custom dose','5','mg'),
                            ('tissue','brain',None)]:
        attr=ET.SubElement(attrs,'EXPERIMENT_ATTRIBUTE')
        ET.SubElement(attr,'TAG').text=tag;ET.SubElement(attr,'VALUE').text=value
        if unit:ET.SubElement(attr,'UNITS').text=unit
        if tag.endswith('AGE'):ET.SubElement(attr,'UNKNOWN').text='retain sibling'
    return result


@pytest.mark.parametrize('provider,parser',[('ena',ENAParser),('sra',SRAParser)])
def test_explicit_experiment_factors_reach_assays_idf_and_residuals(provider,parser):
    package=parser().parse(records(provider));data=package.to_mapping()
    declared={v['name']:v for v in data['series'].get('variable',[])}
    assert set(declared)=={'AGE','Custom dose'}
    assert declared['AGE']['factor']=='age' and declared['Custom dose']['factor']=='other'
    for path in data['series']['assay_paths']:
        assay=next(n for n in path['steps'] if n['kind']=='assay')
        assert assay['factor_values']==[{'name':'AGE','value':'2 d'},
            {'name':'Custom dose','value':'5','unit':{'value':'mg'}},
            {'name':'Custom dose','value':'5','unit':{'value':'mg'}}]
        assert all(not s.get('factor_values') for s in path['steps'] if s['kind']!='assay')
    residual=str(data['extensions']['insdc'])
    assert 'retain sibling' in residual and 'brain' in residual
    assert 'Experimental Factor: Custom dose' not in residual
    output=render(data);assert set(next(r[1:] for r in output if r[0]=='Experimental Factor Name'))==set(declared)
    table=next(r[1] for r in output if r[0]=='SDRF File')
    assert 'Factor Value[AGE]' in table[0] and 'Factor Value[Custom dose]' in table[0]
    assert any('2 d' in row for row in table[1:])


def test_saved_factor_recovery_is_idempotent_scoped_and_keeps_preferred_groups():
    from meta_standards_converter.miniml.archive_factors import restore_native_factors
    package=SRAParser().parse(records('sra'));data=package.to_mapping()
    retained=source_records(package)
    data['series'].pop('variable',None)
    for path in data['series']['assay_paths']:
        for node in path['steps']:node.pop('factor_values',None)
    data['extensions']['insdc']['records']=retained+deepcopy(retained)
    unrelated=deepcopy(data['series']['assay_paths'][0])
    next(n for n in unrelated['steps'] if n['kind']=='assay')['name']='SRX999999'
    data['series']['assay_paths'].append(unrelated)
    restore_native_factors(data)
    assert not any(n.get('factor_values') for n in unrelated['steps'])
    before=deepcopy(data);restore_native_factors(data);assert data==before
    for path in data['series']['assay_paths'][:-1]:
        assay=next(n for n in path['steps'] if n['kind']=='assay')
        assay['factor_values']=[{'name':'AGE','value':'3','unit':{'value':'days'}},
                               {'name':'Custom dose','value':'not provided'}]
    restore_native_factors(data)
    for path in data['series']['assay_paths'][:-1]:
        values=next(n for n in path['steps'] if n['kind']=='assay')['factor_values']
        assert [v for v in values if v['name']=='AGE']==[{'name':'AGE','value':'3','unit':{'value':'days'}}]
        assert [v for v in values if v['name']=='Custom dose']==[{'name':'Custom dose','value':'5','unit':{'value':'mg'}}]*2


def test_factor_residual_multiplicity_does_not_count_file_projections():
    from meta_standards_converter.miniml.archive_residuals import finalize
    package=SRAParser().parse(records('sra'));data=package.to_mapping()
    # One represented dose copied into several file/run paths cannot consume two source occurrences.
    for path in data['series']['assay_paths']:
        assay=next(n for n in path['steps'] if n['kind']=='assay')
        assay['factor_values']=[{'name':'Custom dose','value':'5','unit':{'value':'mg'}}]
    data['series']['assay_paths']*=3
    result=finalize(data,source_records(package)).to_mapping()
    assert 'Experimental Factor: Custom dose' in str(result['extensions']['insdc'])


@pytest.mark.parametrize('child', ['TAG', 'VALUE', 'UNITS'])
def test_mapped_factor_children_keep_unknown_attributes_and_nested_content(child):
    source = records('sra')
    attribute = source.xml[0].find('.//EXPERIMENT_ATTRIBUTE[UNITS]')
    target = attribute.find(child)
    target.set('reviewed_unit_context', 'retained detail')
    ET.SubElement(target, 'NOTE').text = 'unmapped nested detail'
    data = SRAParser().parse(source).to_mapping()
    residual = str(data['extensions']['insdc'])
    assert 'retained detail' in residual and 'unmapped nested detail' in residual
    assert 'Experimental Factor: Custom dose' in residual


@pytest.mark.parametrize('bad_scope', ['run', 'sample', 'unknown_scan'])
def test_factor_recovery_and_residuals_reject_incompatible_bindings(bad_scope):
    from meta_standards_converter.miniml.archive_factors import restore_native_factors
    from meta_standards_converter.miniml.archive_residuals import Projection
    package = SRAParser().parse(records('sra'))
    data = package.to_mapping()
    path = data['series']['assay_paths'][0]
    data['series']['assay_paths'] = [path]
    for step in path['steps']:
        step.pop('factor_values', None)
        if bad_scope == 'sample' and step.get('sample_ref'):
            step['sample_ref'] = 'SRS99999999'
        if step['kind'] == 'scan' and bad_scope != 'sample':
            step['name'] = 'SRR99999999' if bad_scope == 'run' else 'unbound scan'
            if bad_scope == 'unknown_scan':
                step.pop('comments', None)
    restore_native_factors(data, source_records(package))
    assert not any(s.get('factor_values') for s in path['steps'])
    next(s for s in path['steps'] if s['kind'] == 'assay')['factor_values'] = [{'name': 'AGE', 'value': '2 d'}]
    assert not Projection(data).factor('SRX7812918', {'name': 'AGE', 'value': '2 d'})


def test_true_assay_only_path_and_verified_scan_alias_can_recover_factors():
    from meta_standards_converter.miniml.archive_factors import restore_native_factors
    package = SRAParser().parse(records('sra'))
    for assay_only in (True, False):
        data = package.to_mapping()
        for path in data['series']['assay_paths']:
            for step in path['steps']:
                step.pop('factor_values', None)
                if step['kind'] == 'scan':
                    step['name'] = 'supplied scan alias'
            if assay_only:
                path['steps'] = [s for s in path['steps'] if s['kind'] != 'scan']
        restore_native_factors(data, source_records(package))
        assert all(any(s.get('factor_values') for s in p['steps']) for p in data['series']['assay_paths'])


@pytest.mark.parametrize('second', [{'name': 'AGE', 'value': '2 d'},
                                   {'name': 'AGE', 'value': '2', 'unit': {'value': 'd'}}])
def test_equivalent_factor_projections_share_occurrence_capacity(second):
    from meta_standards_converter.miniml.archive_residuals import Projection
    data = SRAParser().parse(records('sra')).to_mapping()
    for path in data['series']['assay_paths']:
        for step in path['steps']:
            step.pop('factor_values', None)
            if step['kind'] == 'assay':
                step['factor_values'] = [{'name': 'AGE', 'value': '2 d'}]
            if step['kind'] == 'scan':
                step['factor_values'] = [{'name': 'AGE', 'value': '2', 'unit': {'value': 'd'}}]
    projection = Projection(data)
    assert projection.factor('SRX7812918', {'name': 'AGE', 'value': '2 d'})
    assert not projection.factor('SRX7812918', second)

    assert not projection.factor('SRX7812918', {'name': 'AGE', 'value': '2', 'unit': {'value': 'h'}})
    assert not projection.factor('SRX999999', second)
    # Two supplied occurrences on one node do support two source occurrences.
    assay = next(n for n in data['series']['assay_paths'][0]['steps'] if n['kind'] == 'assay')
    assay['factor_values'] = [deepcopy(second), deepcopy(second)]
    projection = Projection(data)
    assert projection.factor('SRX7812918', second)
    assert projection.factor('SRX7812918', second)
    assert not projection.factor('SRX7812918', second)


@pytest.mark.parametrize('reverse', [False, True])
def test_mixed_factor_forms_together_are_consumed_independently_of_source_order(reverse):
    from meta_standards_converter.miniml.archive_residuals import Projection
    data = SRAParser().parse(records('sra')).to_mapping()
    values = [{'name': 'AGE', 'value': '2 d'}, {'name': 'AGE', 'value': '2', 'unit': {'value': 'd'}}]
    for path in data['series']['assay_paths']:
        for node in path['steps']:
            node.pop('factor_values', None)
            if node['kind'] == 'assay':
                node['factor_values'] = deepcopy(values)
    projection = Projection(data)
    ordered = values[::-1] if reverse else values
    assert [projection.factor('SRX7812918', v) for v in ordered] == [True, True]
    assert not projection.factor('SRX7812918', values[0])
