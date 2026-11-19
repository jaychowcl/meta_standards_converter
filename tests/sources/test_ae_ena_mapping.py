# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Exact alias discovery is only a lead; public record links establish identity."""
from xml.etree import ElementTree as ET
from unittest.mock import Mock
import pytest
from meta_standards_converter.sources.ena import ENASource

from tests.support.ae import MAPPINGS


def project(ae, study, primary, *, link=True, public=True):
    node = ET.fromstring(f'''<PROJECT_SET><PROJECT accession="{primary}" alias="{ae}">
      <IDENTIFIERS><PRIMARY_ID>{primary}</PRIMARY_ID><SECONDARY_ID>{study}</SECONDARY_ID></IDENTIFIERS>
      <PROJECT_ATTRIBUTES><PROJECT_ATTRIBUTE><TAG>ENA-STATUS</TAG><VALUE>{'public' if public else 'private'}</VALUE></PROJECT_ATTRIBUTE></PROJECT_ATTRIBUTES>
      <PROJECT_LINKS /></PROJECT></PROJECT_SET>''')
    if link:
        links = node.find('.//PROJECT_LINKS')
        links.append(ET.fromstring(f'<PROJECT_LINK><URL_LINK><URL>https://www.ebi.ac.uk/biostudies/arrayexpress/studies/{ae}</URL></URL_LINK></PROJECT_LINK>'))
    return node


def source(ae, study, primary, **options):
    service = ENASource(http=Mock())
    service.search = Mock(return_value=[{'study_alias': ae, 'study_accession': primary, 'secondary_study_accession': study}])
    service.xml = Mock(return_value=project(ae, study, primary, **options))
    return service


@pytest.mark.parametrize('ae,study,primary', MAPPINGS)
def test_reported_studies_require_explicit_public_mapping(ae, study, primary):
    service = source(ae, study, primary)
    mapping = service.resolve_arrayexpress(ae)
    assert (mapping.study, mapping.project) == (study, primary)
    assert mapping.evidence == (f'https://www.ebi.ac.uk/biostudies/arrayexpress/studies/{ae}',)
    service.search.assert_called_once_with('study', f'study_alias="{ae}"',
        'study_accession,secondary_study_accession,study_alias')


@pytest.mark.parametrize('options', [{'link': False}, {'public': False}])
def test_alias_and_nonpublic_evidence_are_insufficient(options):
    service = source(*MAPPINGS[0], **options)
    with pytest.raises(ValueError) as caught:
        service.resolve_arrayexpress(MAPPINGS[0][0])
    assert caught.value.code == 'fallback_unverified'


def test_mismatched_record_or_link_is_rejected():
    ae, study, primary = MAPPINGS[0]
    for root in (project('E-MTAB-999', study, primary), project(ae, 'ERP999', 'PRJEB999')):
        service = source(ae, study, primary)
        service.xml.return_value = root
        with pytest.raises(ValueError) as caught:
            service.resolve_arrayexpress(ae)
        assert caught.value.code == 'fallback_unverified'


def test_project_and_read_study_aliases_deduplicate_using_verified_xml():
    service = source(*MAPPINGS[0])
    service.search.return_value.append({'study_alias': MAPPINGS[0][0], 'study_accession': MAPPINGS[0][1]})
    assert service.resolve_arrayexpress(MAPPINGS[0][0]).study == MAPPINGS[0][1]


def test_two_distinct_verified_studies_are_ambiguous():
    ae, study, primary = MAPPINGS[0]
    service = source(ae, study, primary)
    service.search.return_value.append({'study_alias': ae, 'study_accession': 'PRJEB999', 'secondary_study_accession': 'ERP999'})
    service.xml.side_effect = lambda ids: project(ae, study, primary) if ids == [primary] else project(ae, 'ERP999', 'PRJEB999')
    with pytest.raises(ValueError) as caught:
        service.resolve_arrayexpress(ae)
    assert caught.value.code == 'fallback_ambiguous'


@pytest.mark.parametrize('rows', [{}, [None], [{'study_alias': MAPPINGS[0][0]}]])
def test_malformed_discovery_is_not_absence(rows):
    service = source(*MAPPINGS[0])
    service.search.return_value = rows
    with pytest.raises(ValueError) as caught:
        service.resolve_arrayexpress(MAPPINGS[0][0])
    assert caught.value.code == 'provider_response_invalid'


def test_link_from_nested_record_does_not_verify_parent():
    ae, study, primary = MAPPINGS[0]
    service = source(ae, study, primary, link=False)
    nested = ET.SubElement(service.xml.return_value.find('PROJECT'), 'SAMPLE')
    nested.append(project(ae, 'ERP999', 'PRJEB999').find('.//PROJECT_LINKS'))
    with pytest.raises(ValueError) as caught:
        service.resolve_arrayexpress(ae)
    assert caught.value.code == 'fallback_unverified'


def test_invalid_xml_response_is_classified():
    service = source(*MAPPINGS[0])
    service.xml.side_effect = ET.ParseError('PRIVATE XML')
    with pytest.raises(ValueError) as caught:
        service.resolve_arrayexpress(MAPPINGS[0][0])
    assert caught.value.code == 'provider_response_invalid'
    assert 'PRIVATE' not in str(caught.value)
