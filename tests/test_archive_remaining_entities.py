# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
import xml.etree.ElementTree as ET

from meta_standards_converter.miniml import MINiMLCodec
from meta_standards_converter.miniml.sra_parser import SRAParser
from meta_standards_converter.metadata.archive_enrichment import merge_archive_metadata
from meta_standards_converter.magetab.idf import IDFConstructor
from tests.test_native_archive_parsers import fixture_records
from tests.test_archive_fidelity_enrichment import enriched_native, workflow


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


def test_platform_declarations_and_all_sample_references_are_imported():
    extra = workflow().to_mapping()
    extra['platform'] = [{'iid': 'GPL13112', 'accession': [{'value': 'GPL13112', 'database': 'GEO'}], 'title': 'Illumina', 'technology': 'high-throughput sequencing', 'contact_ref': [{'ref': extra['contributor'][0]['iid']}]}]
    extra['database'].append({'iid': 'GEO', 'name': 'GEO'})
    extra['sample'][0]['platform_ref'] = {'ref': 'GPL13112'}
    data, issues = merge_archive_metadata(enriched_native(), MINiMLCodec().decode(extra).package, prefer=True)
    data = data.to_mapping()
    assert not issues
    platform = data['platform'][0]
    assert data['sample'][0]['platform_ref'] == {'ref': platform['iid']}
    assert platform['accession'][0]['value'] == 'GPL13112'
    assert platform['contact_ref'][0]['ref'] in {c['iid'] for c in data['contributor']}


def test_idf_inline_contacts_affiliations_and_roles():
    data = {'organization': [{'iid': 'O', 'name': 'Lab'}], 'contributor': [], 'series': {'contact': [{'person': {'first': 'Jane', 'last': 'Doe'}, 'organization_ref': {'ref': 'O'}, 'roles': [{'value': 'submitter', 'term_source_ref': 'EFO', 'term_accession_number': 'EFO:1'}]}]}}
    rows = dict((r[0], r[1:]) for r in IDFConstructor()._idf_persons(data))
    assert rows['Person Last Name'] == ['Doe']
    assert rows['Person Affiliation'] == ['Lab']
    assert rows['Person Roles'] == ['submitter']
    assert rows['Person Roles Term Source Ref'] == ['EFO']


def test_native_ontology_references_have_declarations():
    data = SRAParser().parse(fixture_records('sra')).to_mapping()
    assert 'NCBITaxon' in {d['iid'] for d in data['database']}


def test_declared_archive_centers_are_organizations_without_invented_people():
    from meta_standards_converter.miniml.ena_parser import ENAParser
    records=fixture_records('ena')
    records.xml[0].find('STUDY').set('center_name','Declared sequencing centre')
    data=ENAParser().parse(records).to_mapping()
    assert any(o['name']=='Declared sequencing centre' and o.get('role')=='center_name' for o in data['organization'])
    assert not data['contributor']
