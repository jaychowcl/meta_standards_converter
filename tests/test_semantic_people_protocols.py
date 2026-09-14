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

from meta_standards_converter.magetab.idf import IDFConstructor
from meta_standards_converter.magetab.protocol_export import prepare_protocols
from tests.test_native_file_layout import package
from tests.test_protocol_export import definition


def people(data):
    rows = IDFConstructor()._idf_persons(data)
    return [dict(zip([r[0] for r in rows], column)) for column in zip(*(r[1:] for r in rows))]


def test_people_scope_exact_facts_and_original_occurrences():
    data = package()
    one = {'iid': 'one', 'person': {'first': 'Jane', 'last': 'Doe'}, 'organization': 'Lab', 'email': 'jane@lab.org'}
    data['contributor'] = [one, {**deepcopy(one), 'iid': 'two'},
        {'iid': 'platform', 'person': {'first': 'Platform', 'last': 'Owner'}},
        {'iid': 'admin', 'person': {'first': 'Geo', 'last': 'Curators'}, 'organization': 'NCBI', 'email': 'geo-group@ncbi.nlm.nih.gov'}]
    data['series']['contributor_ref'] = [{'ref': 'admin'}]
    data['sample'][0]['contact_ref'] = [{'ref': 'one'}, {'ref': 'two'}]
    data['platform'] = [{'iid': 'P', 'contact_ref': [{'ref': 'platform'}]}]
    before = deepcopy(data)
    assert [p['Person First Name'] for p in people(data)] == ['Jane']
    assert data == before
    data['series']['contact_ref'] = [{'ref': 'platform'}]
    assert [p['Person First Name'] for p in people(data)] == ['Jane', 'Platform']


def test_people_complement_only_with_compatible_identity_and_keep_conflicts():
    data = {'series': {}, 'contributor': [
        {'iid': 'a', 'person': {'first': 'Jane', 'last': 'Doe'}, 'email': 'jane@lab.org'},
        {'iid': 'b', 'person': {'first': 'Jane', 'last': 'Doe'}, 'email': 'jane@lab.org', 'organization': 'Lab'},
        {'iid': 'c', 'person': {'first': 'Jane', 'last': 'Doe'}, 'email': 'other@lab.org', 'organization': 'Lab'},
        {'iid': 'd', 'person': {'first': 'Jane', 'last': 'Doe'}, 'email': 'jane@lab.org', 'organization': 'Different lab'},
        {'iid': 'e', 'person': {'first': 'Other', 'last': 'Person'}, 'email': 'jane@lab.org'}]}
    # The incomplete first record cannot choose between two conflicting affiliations.
    assert len(people(data)) == 5
    del data['contributor'][3]
    result = people(data)
    assert len(result) == 3
    assert result[0]['Person Affiliation'] == 'Lab'


def test_same_name_is_not_identity_and_protocol_contacts_are_eligible():
    same = {'person': {'first': 'Alex', 'last': 'Kim'}}
    data = {'series': {'protocols': [{'name': 'P', 'contacts': ['a']}]},
        'platform': [{'iid': 'platform', 'contact_ref': [{'ref': 'a'}]}],
        'contributor': [{**same, 'iid': 'a', 'email': 'a@lab.org'}, {**same, 'iid': 'b', 'email': 'b@lab.org'}]}
    assert len(people(data)) == 2


def test_protocol_alias_definitions_do_not_erase_real_applications():
    data = package()
    a, b, c = definition('a'), definition('b'), definition('c')
    b['type'] = {'value': 'nucleic acid library construction protocol'}
    c['type'] = {'value': 'nucleic acid extraction protocol'}
    data['series']['protocols'] = [a, b, c]
    path = data['series']['assay_paths'][0]
    path['steps'].insert(1, {'kind': 'protocol_application', 'protocol_ref': 'a', 'performer': 'first'})
    path['steps'].insert(2, {'kind': 'protocol_application', 'protocol_ref': 'b', 'performer': 'second'})
    prepare_protocols(data)
    assert len(data['series']['protocols']) == 2
    applications = [s for s in path['steps'] if s.get('performer')]
    assert len(applications) == 2 and applications[0]['protocol_ref'] == applications[1]['protocol_ref']
    assert {s['performer'] for s in applications} == {'first', 'second'}


def test_protocol_aliases_keep_ontology_conflicts_and_registered_accessions():
    data = package()
    a, b = definition('a'), definition('b')
    a['type'] = {'value': 'treatment protocol', 'term_accession_number': 'EFO:1'}
    b['type'] = {'value': 'sample treatment protocol', 'term_accession_number': 'EFO:2'}
    data['series']['protocols'] = [a, b]
    prepare_protocols(data)
    assert len(data['series']['protocols']) == 2
    a['type'].pop('term_accession_number'); b['type'].pop('term_accession_number')
    data['series']['protocols'] = [a, b]
    prepare_protocols(data)
    assert len(data['series']['protocols']) == 1
    a['name'], b['name'] = 'P-MTAB-1', 'P-MTAB-2'
    data['series']['protocols'] = [a, b]
    prepare_protocols(data)
    assert len(data['series']['protocols']) == 2


def test_nested_contacts_same_id_conflicts_and_extra_affiliations():
    jane = {'iid': 'same', 'person': {'first': 'Jane', 'last': 'Doe'}, 'email': 'jane@lab.org'}
    data = {'series': {}, 'sample': [{'contact': [jane, {**jane, 'department': 'Biology'},
        {**jane, 'department': 'Chemistry'}]}], 'contributor': []}
    rows = people(data)
    assert len(rows) == 3
    assert {r['Person Affiliation'] for r in rows} == {None, 'Biology', 'Chemistry'}


def test_opaque_contacts_and_unscoped_saved_native_people_survive():
    a = {'iid': 'a', 'person': {'first': 'Jane', 'last': 'Doe'}, 'web_link': 'https://lab.org/Jane'}
    data = {'source': {'format': 'SRA'}, 'series': {}, 'contributor': [a, {**a, 'iid': 'b', 'web_link': 'https://lab.org/jane'}]}
    assert len(people(data)) == 2


@pytest.mark.parametrize('reverse', [False, True])
def test_local_protocol_cannot_choose_between_registered_identities(reverse):
    data = package()
    registered = [definition('P-MTAB-1'), definition('P-MTAB-2')]
    data['series']['protocols'] = [definition('local'), *(reversed(registered) if reverse else registered)]
    prepare_protocols(data)
    assert len(data['series']['protocols']) == 3


@pytest.mark.parametrize('reverse', [False, True])
def test_sparse_protocol_cannot_choose_between_conflicting_instruments(reverse):
    data = package()
    sparse, a, b = [definition(n) for n in ('sparse', 'a', 'b')]
    sparse.pop('hardware'); a['hardware'] = ['A']; b['hardware'] = ['B']
    data['series']['protocols'] = [sparse, *([b, a] if reverse else [a, b])]
    data['series']['assay_paths'][0]['steps'].insert(1, {'kind': 'protocol_application', 'protocol_ref': 'sparse'})
    prepare_protocols(data)
    assert len(data['series']['protocols']) == 3
    ref = data['series']['assay_paths'][0]['steps'][1]['protocol_ref']
    assert not next(p for p in data['series']['protocols'] if p['name'] == ref).get('hardware')
