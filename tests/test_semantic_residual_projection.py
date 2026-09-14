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
from tests.test_native_file_layout import package
from meta_standards_converter.miniml.archive_dates import normalize_archive_dates
from meta_standards_converter.miniml.archive_residuals import finalize


def source_record(source, workflows=False):
    return {'provider':'MAGE-TAB','kind':'MINiML_workflows' if workflows else 'MINiML',
            'accession':source['series']['iid'], 'metadata':deepcopy(source['series'] if workflows else source)}


def retained(data, record):
    return [r for r in finalize(data,[record]).to_mapping()['extensions']['insdc']['records'] if r['kind']==record['kind']]


@pytest.mark.parametrize('workflows',[False,True])
def test_mapped_dates_do_not_retain_whole_source_workflows(workflows):
    source=package();target=deepcopy(source);normalize_archive_dates(target)
    result=retained(target,source_record(source,workflows))
    assert all(set(r['metadata']) <= {'database'} for r in result)


@pytest.mark.parametrize('status',['missing','conflict','one_occurrence','precision'])
def test_unrepresented_date_occurrences_retain_original_source_shapes(status):
    source=package();target=deepcopy(source);normalize_archive_dates(target)
    ss=target['sample'][0]['status']
    dates=[s for s in ss if s.get('database')=='ENA' and s.get('release_date')=='2020-01-01']
    assert len(dates)==2
    if status=='missing':target['sample'][0]['status']=[s for s in ss if s not in dates]
    if status=='conflict':
        for s in dates:s['release_date']='2021-01-01'
    if status=='precision':
        for s in dates:s['release_date']='2020-01-01T00:00:00Z'
    if status=='one_occurrence':ss.remove(dates[0])
    result=retained(target,source_record(source))
    chars=result[0]['metadata']['sample'][0]['channel'][0]['characteristics']
    assert sum(c.get('name')=='ENA-FIRST-PUBLIC' for c in chars)==(1 if status=='one_occurrence' else 2)


def test_scoped_path_residual_keeps_unknown_sibling_and_event_position_only():
    source=package();normalize_archive_dates(source)
    path=source['series']['assay_paths'][0];source['series']['assay_paths']=[path]
    path['steps'].insert(1,{'kind':'protocol_application','protocol_ref':'unknown_method'})
    path['steps'].insert(2,{'kind':'protocol_application','protocol_ref':'unknown_method','comments':[{'name':'operator note','value':'source-only detail'}]})
    source['series'].setdefault('protocols',[]).append({'name':'unknown_method','type':{'value':'treatment protocol'}})
    target=deepcopy(source);target['series']['assay_paths'][0]['steps'][2].pop('comments')
    result=retained(target,source_record(source))
    steps=result[0]['metadata']['series']['assay_paths'][0]['steps']
    assert len(steps)==len(path['steps'])
    assert steps[2]['comments']==path['steps'][2]['comments']
    assert set(steps[0]) <= {'kind','name','sample_ref'}
    assert steps[1]=={'kind':'protocol_application','protocol_ref':'unknown_method'}
    assert 'link' not in steps[-1]


def test_ambiguous_path_signatures_do_not_pair_unknown_occurrences():
    source=package();normalize_archive_dates(source)
    a=source['series']['assay_paths'][0];b=deepcopy(a)
    a['unknown']='first';b['unknown']='second';source['series']['assay_paths']=[a,b]
    target=deepcopy(source)
    for p in target['series']['assay_paths']:p.pop('unknown')
    result=retained(target,source_record(source))
    assert result[0]['metadata']['series']['assay_paths']==source['series']['assay_paths']


def test_date_annotations_and_unknown_sample_bindings_are_not_lost():
    source=package();target=deepcopy(source);normalize_archive_dates(target)
    source['sample'][0]['channel'][0]['characteristics'][0].update(unit='unknown unit',extra_note='keep')
    source['series']['assay_paths'][0]['steps'][0]['sample_ref']='not-in-output'
    result=retained(target,source_record(source))
    assert 'unknown unit' in str(result) and 'not-in-output' in str(result)


@pytest.mark.parametrize('same_reference',[False,True])
def test_residual_containment_never_erases_changed_procedure_order(same_reference):
    source=package();normalize_archive_dates(source)
    path=source['series']['assay_paths'][0];source['series']['assay_paths']=[path]
    first={'kind':'protocol_application','protocol_ref':'treat','comments':[{'name':'duration','value':'10'}]}
    second={'kind':'protocol_application','protocol_ref':'treat' if same_reference else 'wash',
            'comments':[{'name':'duration','value':'20'}]}
    path['steps'][1:1]=[first,second]
    target=deepcopy(source);target['series']['assay_paths'][0]['steps'][1:3]=[deepcopy(second),deepcopy(first)]
    result=retained(target,source_record(source))
    assert result[0]['metadata']['series']['assay_paths'][0]['steps'][1:3]==[first,second]


def test_linked_biosamples_does_not_recopied_exact_underscore_fields():
    import xml.etree.ElementTree as E
    from tests.test_native_archive_parsers import fixture_records
    from meta_standards_converter.miniml.ena_parser import ENAParser
    records=fixture_records('ena');sample=records.xml[1].find('SAMPLE')
    attrs=sample.find('SAMPLE_ATTRIBUTES')
    for value in ('native','native'):
        attr=E.SubElement(attrs,'SAMPLE_ATTRIBUTE');E.SubElement(attr,'TAG').text='GAL_sample_id';E.SubElement(attr,'VALUE').text=value
    records.linked.append({'provider':'biosamples','kind':'sample','accession':sample.get('accession'), 'metadata':{
        'accession':sample.get('accession'),'characteristics':{'GAL_sample_id':[{'text':'linked-conflict'}], 'new_field':[{'text':'extra'},{'text':'second'}]}}})
    data=ENAParser().parse(records).to_mapping();cs=data['sample'][0]['channel'][0]['characteristics']
    assert [c['value'] for c in cs if c['name']=='GAL_sample_id']==['native','native']
    assert [c['value'] for c in cs if c['name']=='new_field']==['extra','second']
    assert 'linked-conflict' in str(data['extensions'])
