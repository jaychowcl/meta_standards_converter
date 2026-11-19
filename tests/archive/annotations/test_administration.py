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
from tests.support.archive import legacy
from tests.support.archive import fixture_records
from meta_standards_converter.miniml.ena_parser import ENAParser

from meta_standards_converter.miniml import MINiMLCodec

from meta_standards_converter.miniml.sra_parser import SRAParser


from meta_standards_converter.magetab.idf import IDFConstructor

from tests.support.archive import native


from meta_standards_converter.miniml.archive_residuals import finalize


def test_administration_preserves_occurrences_and_cleans_copies():
    from meta_standards_converter.miniml.archive_administration import normalize_administration
    data=legacy();sample=data['sample'][0]
    items=[{'name':'INSDC status','value':'public'},{'name':'INSDC status','value':'public'},
           {'name':'INSDC center name','value':'SC'},{'name':'INSDC center alias','value':'SC'}]
    sample['channel'][0]['characteristics']+=deepcopy(items)
    for path in data['series']['assay_paths']:path['steps'][0]['characteristics']+=deepcopy(items)
    normalize_administration(data)
    statuses=[s for s in sample['status'] if s.get('comment')==[items[0]]]
    assert len(statuses)==2
    orgs={o['iid']:o for o in data['organization']}
    links=[r for r in sample['relation'] if r['type']=='archive center' and orgs[r['target']].get('role','').startswith('INSDC')]
    assert len(links)==2
    assert {orgs[r['target']]['role'] for r in links}=={'INSDC center name','INSDC center alias'}
    assert all(c not in items for c in sample['channel'][0]['characteristics'])
    assert all(c not in items for p in data['series']['assay_paths'] for st in p['steps'] for c in st.get('characteristics',[]))
    before=deepcopy(data);normalize_administration(data);assert data==before


def test_fresh_independent_center_and_status_sources_remain_scoped():
    records=fixture_records('ena');node=records.xml[1].find('SAMPLE');acc=node.get('accession');node.set('center_name','SC')
    attrs=node.find('SAMPLE_ATTRIBUTES');a=ET.SubElement(attrs,'SAMPLE_ATTRIBUTE')
    ET.SubElement(a,'TAG').text='INSDC status';ET.SubElement(a,'VALUE').text='public';ET.SubElement(a,'CUSTOM').text='keep'
    records.linked.append({'provider':'biosamples','kind':'sample','accession':acc,'metadata':{'accession':acc,'status':'PUBLIC','characteristics':{'INSDC center name':[{'text':'SC'}]}}})
    data=ENAParser().parse(records).to_mapping();sample=data['sample'][0]
    assert len([r for r in sample['relation'] if r['type']=='archive center'])==2
    assert any(s.get('database')=='INSDC' and s.get('comment')==[{'name':'INSDC status','value':'public'}] for s in sample['status'])
    assert any(s.get('database')=='BioSamples' and s.get('comment')==[{'name':'status','value':'PUBLIC'}] for s in sample['status'])
    assert 'keep' in str(data['extensions'])
    assert not any(c['name']=='INSDC status' for c in sample['channel'][0]['characteristics'])


def test_owner_and_center_relations_follow_scoped_source_roles_on_saved_json():
    from meta_standards_converter.miniml.archive_administration import normalize_administration
    data=legacy(); sample=data['sample'][0]; sid=sample['iid']
    data['organization']=[
        {'iid':'owner','name':'Owner','role':'owner','sample_accession':sid},
        {'iid':'center','name':'Center','role':'center_name','sample_accession':sid},
        {'iid':'shared','name':'Shared','roles':['owner','center'], 'source_occurrences':[
            {'iid':'shared-owner','role':'owner','sample_accession':sid},
            {'iid':'other-center','role':'center','sample_accession':'SAMN999'}]}]
    sample['relation']=[{'type':'archive center','target':'owner'},
                        {'type':'archive center','target':'shared'}]
    normalize_administration(data)
    assert {'type':'sample owner','target':'owner'} in sample['relation']
    assert {'type':'sample owner','target':'shared'} in sample['relation']
    assert {'type':'archive center','target':'center'} in sample['relation']
    assert not any(r['type']=='archive center' and r['target'] in {'owner','shared'} for r in sample['relation'])
    before=deepcopy(data); normalize_administration(data); assert data==before


def test_legacy_roleless_occurrences_recover_only_a_sole_aggregate_role():
    from meta_standards_converter.miniml.archive_administration import normalize_administration
    data=legacy();sample=data['sample'][0];sid=sample['iid']
    data['organization']=[
        {'iid':'old','role':'owner','source_occurrences':[{'iid':'a','sample_accession':sid}]},
        {'iid':'mixed','roles':['owner','center'],'source_occurrences':[{'iid':'b','sample_accession':sid}]}]
    normalize_administration(data)
    assert {'type':'sample owner','target':'old'} in sample['relation']
    assert {'type':'organization','target':'mixed'} in sample['relation']
    assert not any(r['target']=='mixed' and r['type']!='organization' for r in sample['relation'])


def test_native_organizations_contacts_and_invalid_email():
    records = fixture_records('sra')
    records.xml.append(ET.fromstring('''<RecordSet><DocumentSummary uid="9"><Project><ProjectID><ArchiveID accession="PRJNA609050"/></ProjectID></Project><Organization role="owner" type="institute" url="https://lab.test"><Name>Lab</Name><Contact email="a@lab.test"><Name><First>A</First><Last>B</Last></Name></Contact><Contact email="https://lab.test/contact"/></Organization><Organization role="owner"><Name>Unstaffed</Name></Organization></DocumentSummary></RecordSet>'''))
    data = SRAParser().parse(records).to_mapping()
    org = next(o for o in data['organization'] if o['name'] == 'Lab')
    assert org['web_link'] == 'https://lab.test'
    assert any(o['name'] == 'Unstaffed' for o in data['organization'])
    person = next(c for c in data['contributor'] if c.get('email') == 'a@lab.test')
    assert person['organization_ref'] == {'ref': org['iid']}
    assert not any(c.get('email', '').startswith('http') for c in data['contributor'])
    assert 'https://lab.test/contact' in str(data['extensions'])


def test_idf_inline_contacts_affiliations_and_roles():
    data = {'organization': [{'iid': 'O', 'name': 'Lab'}], 'contributor': [], 'series': {'contact': [{'person': {'first': 'Jane', 'last': 'Doe'}, 'organization_ref': {'ref': 'O'}, 'roles': [{'value': 'submitter', 'term_source_ref': 'EFO', 'term_accession_number': 'EFO:1'}]}]}}
    rows = dict((r[0], r[1:]) for r in IDFConstructor()._idf_persons(data))
    assert rows['Person Last Name'] == ['Doe']
    assert rows['Person Affiliation'] == ['Lab']
    assert rows['Person Roles'] == ['submitter']
    assert rows['Person Roles Term Source Ref'] == ['EFO']


def test_declared_archive_centers_are_organizations_without_invented_people():
    from meta_standards_converter.miniml.ena_parser import ENAParser
    records=fixture_records('ena')
    records.xml[0].find('STUDY').set('center_name','Declared sequencing centre')
    data=ENAParser().parse(records).to_mapping()
    assert any(o['name']=='Declared sequencing centre' and o.get('role')=='center_name' for o in data['organization'])
    assert not data['contributor']


def test_known_address_fields_are_projected_as_address_fields():
    records=fixture_records('sra')
    data=SRAParser().parse(records).to_mapping()
    person=next(c for c in data['contributor'] if c.get('email')=='bdaisley@uwo.ca')
    assert person['address']['postal_code']=='N6A3K7'
    assert person['address']['city']=='London'
    assert person['address']['country']=='Canada'


def test_biosample_owner_contacts_are_scoped_and_residuals_are_pruned():
    from tests.support.archive import nodes
    from meta_standards_converter.miniml.archive_residuals import finalize
    records = fixture_records('sra')
    owner = records.xml[1].find('.//BioSample/Owner')
    contact = owner.find('Contacts/Contact')
    contact.set('sec_email', 'second@lab.test')
    contact.set('role', 'submitter')
    ET.SubElement(contact, 'Unknown').text = 'retain this detail'
    data = SRAParser().parse(records).to_mapping()
    sample = data['sample'][0]
    refs = {r['ref'] for r in sample['contact_ref']}
    people = [c for c in data['contributor'] if c['iid'] in refs]
    assert len(people) == 1
    person = people[0]
    assert person['person'] == {'first': 'Brendan', 'last': 'Daisley'}
    assert person['extensions']['secondary_email'] == 'second@lab.test'
    assert person['roles'] == [{'value': 'submitter'}]
    assert person['iid'] not in {r['ref'] for r in data['series']['contributor_ref']}
    org = next(o for o in data['organization'] if o['iid'] == person['organization_ref']['ref'])
    assert org['name'] == 'University of Western Ontario'
    assert org['web_link'] == 'http://www.uwo.ca'
    residual = [r for r in data['extensions']['insdc']['records'] if r['kind'] == 'BioSample']
    assert 'retain this detail' in str(residual)
    assert 'Brendan' not in str(residual) and 'second@lab.test' not in str(residual)
    assert finalize(data).to_mapping() == data
    rows = dict((r[0], r[1:]) for r in IDFConstructor()._idf_persons(data))
    assert rows['Person First Name'].count('Brendan') == 2  # independent source occurrences
    MINiMLCodec().decode(data, strict=True)


def test_owner_without_contact_and_malformed_secondary_email_are_preserved():
    records = fixture_records('sra')
    owner = records.xml[1].find('.//BioSample/Owner')
    owner.find('Contacts/Contact').set('sec_email', 'https://not-an-email.test')
    ET.SubElement(owner, 'Unknown').text = 'owner detail'
    data = SRAParser().parse(records).to_mapping()
    assert 'https://not-an-email.test' in str(data['extensions'])
    owner.remove(owner.find('Contacts'))
    data = SRAParser().parse(records).to_mapping()
    assert any(o.get('sample_accession') == 'SAMN14218700' for o in data['organization'])
    assert not data['sample'][0].get('contact_ref')
    assert 'owner detail' in str(data['extensions'])


def test_identical_organizations_coalesce_without_merging_contacts_or_scopes():
    data=native().to_mapping();sid=data['sample'][0]['iid']
    data['organization']=[{'iid':'o1','name':'Lab','web_link':'https://lab.example.org','role':'owner','sample_accession':sid},
                          {'iid':'o2','name':'Lab','web_link':'https://lab.example.org','role':'owner','sample_accession':'SAMN2'},
                          {'iid':'o3','name':'Lab','role':'owner','address':{'city':'Elsewhere'}}]
    data['contributor']=[{'iid':'c1','person':{'last':'Contact'},'organization_ref':{'ref':'o1'}},
                         {'iid':'c2','person':{'last':'Contact'},'organization_ref':{'ref':'o2'}}]
    data['series']['contributor_ref']=[{'ref':'c2'}];data['sample'][0]['contact_ref']=[{'ref':'c1'}]
    data=finalize(data,[]).to_mapping()
    assert len(data['organization'])==2
    assert len(data['contributor'])==2
    assert data['contributor'][0]['organization_ref']==data['contributor'][1]['organization_ref']
    assert data['sample'][0]['contact_ref']==[{'ref':'c1'}]
    assert data['series']['contributor_ref']==[{'ref':'c2'}]
    org=next(o for o in data['organization'] if o['iid']=='o1')
    assert {o['iid'] for o in org['source_occurrences']}=={'o1','o2'}
    assert finalize(data,[]).to_mapping()==data


def test_organization_roles_do_not_define_identity_and_names_alone_do_not_merge():
    from meta_standards_converter.miniml.archive_entities import coalesce_organizations
    data = {'source':{'format':'ENA'}, 'organization':[
        {'iid':'a','name':'Lab','web_link':'https://lab.org','role':'owner'},
        {'iid':'b','name':'Lab','web_link':'https://lab.org','role':'center'},
        {'iid':'c','name':'Unknown lab'}, {'iid':'d','name':'Unknown lab'}]}
    coalesce_organizations(data)
    assert len(data['organization']) == 3
    shared = data['organization'][0]
    assert {v['role'] for v in shared['source_occurrences']} == {'owner','center'}
    assert set(shared['roles']) == {'owner','center'}
    before = deepcopy(data); coalesce_organizations(data); assert data == before


def test_list_organization_roles_survive_consolidation():
    from meta_standards_converter.miniml.archive_entities import coalesce_organizations
    data = {'source':{'format':'ENA'}, 'organization':[
        {'iid':'a','name':'Lab','web_link':'https://lab.org','roles':['owner']},
        {'iid':'b','name':'Lab','web_link':'https://lab.org','roles':['centre']}]}
    coalesce_organizations(data)
    assert set(data['organization'][0]['roles']) == {'owner','centre'}
