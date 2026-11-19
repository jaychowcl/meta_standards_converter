# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from copy import deepcopy
from meta_standards_converter.miniml.archive_administration import normalize_administration
from meta_standards_converter.miniml.archive_residuals import Projection


def test_explicit_administration_is_scoped_and_idempotent():
    attrs = [
        {'name':'External Id','value':'SAMN123'}, {'name':'SRA accession','value':'SRS456'},
        {'name':'Submitter Id','value':'local-sample'}, {'name':'ENA-CHECKLIST','value':'ERC000011'},
        {'name':'BioSampleModel','value':'Generic'}, {'name':'broker name','value':'ArrayExpress'},
        {'name':'organism','value':'mouse'}, {'name':'organism','value':'mouse'},
    ]
    data = {'source':{'format':'ENA'}, 'series':{'iid':'ERP1','assay_paths':[{'steps':[
        {'kind':'source','name':'s','sample_ref':'s','characteristics':deepcopy(attrs)}]}]},
        'sample':[{'iid':'s','channel':[{'characteristics':deepcopy(attrs)}]}]}
    normalize_administration(data)
    sample = data['sample'][0]
    assert sample['channel'][0]['characteristics'] == attrs[-2:]
    assert data['series']['assay_paths'][0]['steps'][0]['characteristics'] == attrs[-2:]
    assert {a['value'] for a in sample['accession']} == {'SAMN123','SRS456'}
    assert any(r['target'] == 'local-sample' and r['namespace'] == 'ENA submitter' for r in sample['relation'])
    assert any(o['name'] == 'ArrayExpress' and o['role'] == 'broker name' for o in data['organization'])
    assert {'name':'ENA-CHECKLIST','value':'ERC000011'} in sample['comments']
    before = deepcopy(data); normalize_administration(data); assert data == before
    projection = Projection(data)
    assert all(projection.character(sample, a['name'], a['value']) for a in attrs[:-2])


def test_invalid_external_identifier_is_not_declared_as_verified_accession():
    data = {'source':{'format':'SRA'},'series':{},'sample':[{'iid':'s','channel':[{'characteristics':[
        {'name':'External Id','value':'unknown,identifier'}, {'name':'growth model','value':'biological'}]}]}]}
    normalize_administration(data)
    sample = data['sample'][0]
    assert not sample.get('accession')
    assert sample['comments'][0]['value'] == 'unknown,identifier'
    assert sample['channel'][0]['characteristics'] == [{'name':'growth model','value':'biological'}]


def test_coalesced_organization_role_occurrences_are_residual_matches():
    from meta_standards_converter.miniml.archive_entities import coalesce_organizations
    data = {'source':{'format':'ENA'}, 'series':{},'organization':[
        {'iid':'a','name':'Lab','web_link':'https://lab.org','role':'INSDC center name','sample_accession':'s'},
        {'iid':'b','name':'Lab','web_link':'https://lab.org','role':'broker name','sample_accession':'s'}],
        'sample':[{'iid':'s','relation':[{'type':'archive center','target':'a'},{'type':'archive broker','target':'b'}]}]}
    coalesce_organizations(data)
    projection = Projection(data)
    assert projection.character(data['sample'][0], 'INSDC center name', 'Lab')
    assert not projection.character(data['sample'][0], 'INSDC center name', 'Lab')
    assert projection.character(data['sample'][0], 'broker name', 'Lab')


def test_snapshot_residual_does_not_repeat_mapped_administrative_characteristics():
    from meta_standards_converter.miniml.archive_residuals import _mapped_date_view
    original = {'source':{'format':'ENA'}, 'series':{},'sample':[{'iid':'s','channel':[{'characteristics':[
        {'name':'ENA-CHECKLIST','value':'ERC000011'}, {'name':'unknown','value':'keep'}]}]}]}
    target = deepcopy(original); normalize_administration(target)
    view = _mapped_date_view(original, target)
    assert view['sample'][0]['channel'][0]['characteristics'] == [{'name':'unknown','value':'keep'}]
    assert original['sample'][0]['channel'][0]['characteristics'][0]['name'] == 'ENA-CHECKLIST'


def identifier_example():
    attrs = [{'name':'Sample Name','value':'ERS000087','note':'first occurrence'},
             {'name':'Sample Name','value':'ERS000087','note':'second occurrence'},
             {'name':'sample name','value':'zebrafish embryo'},
             {'name':'Sample Name','value':'ERS999999'},
             {'name':'Alias','value':'E-MTAB-308:Zebrafish embryo 1 dpf 2'},
             {'name':'Alias','value':'ZF_male_sample1'}]
    return {'source':{'format':'ENA'}, 'series':{'iid':'ERP1','assay_paths':[{'steps':[
        {'kind':'source','name':'s','sample_ref':'s','characteristics':deepcopy(attrs)},
        {'kind':'extract','name':'RNA','characteristics':[deepcopy(attrs[0]),deepcopy(attrs[-1])]}]}]},
        'sample':[{'iid':'s','accession':[{'database':'ENA','value':'ERS000087'}],
                   'channel':[{'characteristics':deepcopy(attrs)}]}]}


def test_verified_sample_names_and_local_aliases_keep_occurrences_and_scope():
    from meta_standards_converter.miniml.archive_residuals import _mapped_date_view
    original = identifier_example()
    data = deepcopy(original); normalize_administration(data)
    sample = data['sample'][0]
    kept = original['sample'][0]['channel'][0]['characteristics'][2:4]
    assert sample['channel'][0]['characteristics'] == kept
    steps = data['series']['assay_paths'][0]['steps']
    assert steps[0]['characteristics'] == kept
    assert steps[1] == original['series']['assay_paths'][0]['steps'][1]
    assert [a['note'] for a in sample['accession'] if a.get('label') == 'Sample Name'] == ['first occurrence','second occurrence']
    assert {'type':'Alias','target':'E-MTAB-308:Zebrafish embryo 1 dpf 2','namespace':'E-MTAB-308'} in sample['relation']
    assert {'name':'Alias','value':'ZF_male_sample1'} in sample['comments']
    assert not data.get('organization')
    assert _mapped_date_view(original, data)['series']['assay_paths'][0]['steps'][1] == steps[1]
    view = _mapped_date_view(original, data)
    assert view['sample'][0]['channel'][0]['characteristics'] == kept
    assert view['series']['assay_paths'][0]['steps'][0]['characteristics'] == kept
    before = deepcopy(data); normalize_administration(data); assert data == before


def test_unverified_or_free_text_sample_name_is_never_an_organization_or_identifier():
    data = identifier_example()
    data['sample'][0]['accession'] = []
    normalize_administration(data)
    sample = data['sample'][0]
    assert len(sample['channel'][0]['characteristics']) == 4
    assert not sample.get('accession') and not data.get('organization')


def test_alias_without_explicit_namespace_stays_a_labelled_comment():
    data = identifier_example()
    data['sample'][0]['channel'][0]['characteristics'].append({'name':'Alias','value':'arbitrary:local'})
    normalize_administration(data)
    assert {'name':'Alias','value':'arbitrary:local'} in data['sample'][0]['comments']


def test_displaced_administrative_unit_identifier_remains_residual():
    from meta_standards_converter.miniml.archive_residuals import _mapped_date_view
    original = identifier_example()
    original['series']['assay_paths'] = []
    original['sample'][0]['channel'][0]['characteristics'][0]['unit'] = {
        'value':'unknown-unit','term_source_ref':'UO','term_accession_number':'UO:1'}
    target = deepcopy(original); normalize_administration(target)
    assert not any(c.get('note') == 'first occurrence' for c in _mapped_date_view(original,target)['sample'][0]['channel'][0]['characteristics'])
    next(a for a in target['sample'][0]['accession'] if a.get('note') == 'first occurrence')['unit']['term_accession_number'] = 'UO:2'
    view = _mapped_date_view(original,target)
    assert view['sample'][0]['channel'][0]['characteristics'][0] == original['sample'][0]['channel'][0]['characteristics'][0]
