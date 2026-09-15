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
