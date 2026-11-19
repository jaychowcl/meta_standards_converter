# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Contact presentation must preserve evidence without inventing identity joins."""
from copy import deepcopy
import json
from pathlib import Path

import pytest

from meta_standards_converter.magetab.idf import IDFConstructor


def frozen_contacts():
    return json.loads((Path(__file__).parents[1] /
        'fixtures/edge_cases/contact-consolidation/contact-consolidation.json').read_text())


def people(data):
    before = deepcopy(data)
    rows = IDFConstructor()._idf_persons(data)
    assert data == before
    return [dict(zip([r[0] for r in rows], values)) for values in zip(*(r[1:] for r in rows))]


def matched_sources(data):
    return [c for c in data['contributor'] if c.get('email') == 'smyth@wehi.edu.au']


@pytest.mark.parametrize('reverse', [False, True])
def test_frozen_contacts_consolidate_only_supported_identity(reverse):
    data = frozen_contacts()
    if reverse:
        data['contributor'].reverse()
    result = people(data)
    assert len(result) == 3
    matched = [p for p in result if p['Person Email'] == 'smyth@wehi.edu.au']
    assert len(matched) == 1
    assert '3052' in matched[0]['Person Address']
    assert matched[0]['Person Roles'] == 'submitter'
    assert matched[0]['Person Affiliation'] == 'Smyth, Bioinformatics, Walter and Eliza Hall Institute of Medical Research'
    assert len([p for p in result if not p['Person Email']]) == 2


@pytest.mark.parametrize('change', ['street', 'postcode', 'organization', 'department', 'opaque_suffix',
                                  'opaque_substring', 'middle_name', 'first_name', 'email', 'persistent_id'])
def test_conflicting_contact_facts_and_identities_remain_separate(change):
    data = frozen_contacts()
    geo, ae = matched_sources(data)
    if change == 'street': geo['address']['line'] = ['9 Different Road']
    elif change == 'postcode': ae['address'] = ae['address'].replace('Vic, Australia', 'Vic, 9999, Australia')
    elif change == 'organization': ae['organization'] = 'Other Institute'
    elif change == 'department': ae['department'] = 'Chemistry'
    elif change == 'opaque_suffix': ae['address'] += ', Unknown building'
    elif change == 'opaque_substring': ae['address'] = ae['address'].replace('Royal Pde', 'Royal Pde Annex')
    elif change == 'middle_name': ae['person']['middle'] = 'J'
    elif change == 'first_name': ae['person']['first'] = 'Other'
    elif change == 'email': ae['email'] = 'other@wehi.edu.au'
    elif change == 'persistent_id':
        geo['accession'] = [{'database': 'ORCID', 'value': '0000-0001-2345-6789'}]
        ae['accession'] = [{'database': 'ORCID', 'value': '0000-0002-2345-6789'}]
    assert len(people(data)) == 4


def test_generic_mailbox_cannot_supply_identity():
    data = frozen_contacts()
    for c in matched_sources(data): c['email'] = 'info@wehi.edu.au'
    assert len(people(data)) == 4


@pytest.mark.parametrize('reverse', [False, True])
def test_explicit_identity_can_complete_structured_address_and_qualified_role(reverse):
    data = frozen_contacts()
    geo, ae = matched_sources(data)
    ae['address'] = {k: deepcopy(v) for k, v in geo['address'].items() if k != 'postal_code'}
    ae['roles'] = [{'value': 'submitter', 'term_source_ref': 'EFO', 'term_accession_number': 'EFO:0000001'}]
    for c in (geo, ae):
        c.pop('email')
        c['accession'] = [{'database': 'ORCID', 'value': '0000-0001-2345-6789'}]
    if reverse: data['contributor'].reverse()
    result = people(data)
    assert len(result) == 3
    contact = next(p for p in result if p['Person Roles'])
    assert '3052' in contact['Person Address']
    assert (contact['Person Roles'], contact['Person Roles Term Source Ref'],
            contact['Person Roles Term Accession Number']) == ('submitter', 'EFO', 'EFO:0000001')


@pytest.mark.parametrize('reverse', [False, True])
def test_sparse_contact_does_not_bridge_two_incompatible_addresses(reverse):
    original = matched_sources(frozen_contacts())[0]
    sparse = {'iid': 'sparse', 'person': deepcopy(original['person']), 'email': original['email']}
    conflict = deepcopy(original)
    conflict['iid'] = 'conflict'
    conflict['address']['line'] = ['9 Different Road']
    values = [original, sparse, conflict]
    if reverse: values.reverse()
    assert len(people({'series': {}, 'contributor': values})) == 3


def test_conflicting_persistent_ids_survive_even_identical_display_facts():
    original = matched_sources(frozen_contacts())[0]
    other = deepcopy(original)
    original['accession'] = [{'database': 'ORCID', 'value': '0000-0001-2345-6789'}]
    other['iid'] = 'other'
    other['accession'] = [{'database': 'ORCID', 'value': '0000-0002-2345-6789'}]
    assert len(people({'series': {}, 'contributor': [original, other]})) == 2


def test_conflicting_structured_components_survive_identical_rendered_addresses():
    original = matched_sources(frozen_contacts())[0]
    other = deepcopy(original)
    original['address'] = {'line': ['1 Street'], 'city': 'London', 'country': 'UK'}
    other['address'] = {'line': ['1 Street', 'London'], 'city': 'UK'}
    other['iid'] = 'other'
    result = people({'series': {}, 'contributor': [original, other]})
    assert len(result) == 2
    assert result[0]['Person Address'] == result[1]['Person Address']
