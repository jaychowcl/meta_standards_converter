"""Source-derived reproductions; no production repair is part of this audit."""
import xml.etree.ElementTree as ET

import pytest

from meta_standards_converter.converters import Converter, PreparationPolicy
from meta_standards_converter.miniml.ena_parser import ENAParser
from meta_standards_converter.miniml.sra_parser import SRAParser
from tests.support.archive import fixture_records


# Constructed PubMed boundary, deliberately independent of MSC output. NLM's
# Author model permits Initials without ForeName, Suffix and CollectiveName.
ARTICLE = """<PubmedArticleSet><PubmedArticle><MedlineCitation>
<PMID>123456</PMID><Article><ArticleTitle>Author boundary audit</ArticleTitle>
<AuthorList>
<Author><LastName>Smith</LastName><Initials>AB</Initials></Author>
<Author><LastName>Jones</LastName><ForeName>Alice</ForeName><Initials>A</Initials><Suffix>Jr</Suffix></Author>
<Author><CollectiveName>Example Consortium</CollectiveName></Author>
<Author><LastName>Börsch</LastName><ForeName>Anna</ForeName><Initials>A</Initials></Author>
</AuthorList></Article></MedlineCitation><PubmedData>
<PublicationStatus>ppublish</PublicationStatus><ArticleIdList>
<ArticleId IdType="doi">10.1234/audit</ArticleId></ArticleIdList>
</PubmedData></PubmedArticle></PubmedArticleSet>"""


def source_package(provider):
    records = fixture_records(provider)
    study = next(node for root in records.xml for node in root.iter('STUDY'))
    links = ET.SubElement(study, 'STUDY_LINKS')
    link = ET.SubElement(ET.SubElement(links, 'STUDY_LINK'), 'XREF_LINK')
    ET.SubElement(link, 'DB').text = 'PubMed'
    ET.SubElement(link, 'ID').text = '123456'
    records.xml.append(ET.fromstring(ARTICLE))
    return (SRAParser if provider == 'sra' else ENAParser)().parse(records)


@pytest.mark.parametrize('provider', ['sra', 'ena'])
@pytest.mark.parametrize('stage', ['json', 'magetab'])
@pytest.mark.parametrize('component,expected', [('initials', 'AB Smith'), ('suffix', 'Alice Jones Jr')])
@pytest.mark.xfail(strict=True, reason='MSC-AUTH-001: native PubMed projection drops Initials-only names and Suffix; audit excludes fixes')
def test_native_publication_retains_supplied_author_components(provider, stage, component, expected):
    package = source_package(provider)
    if stage == 'json':
        authors = package.to_mapping()['series']['pubmed_publication'][0]['author_list']
    else:
        result = Converter().convert(package, out_type='magetab', preparation=PreparationPolicy(
            enrichment='off', expand_studies=False, linked_metadata=False, peer_provider=False))
        assert result.items[0].execution == 'succeeded', result.to_dict()
        authors = next(row[1] for row in result.items[0].payload[0] if row[0] == 'Publication Author List')
    assert expected in authors.split(', '), f'{provider}/{stage} lost {component}: {authors}'


@pytest.mark.parametrize('provider', ['sra', 'ena'])
def test_native_publication_retains_group_unicode_order_and_identifier_alignment(provider):
    pub = source_package(provider).to_mapping()['series']['pubmed_publication'][0]
    assert pub['author_list'].split(', ')[2:] == ['Example Consortium', 'Anna Börsch']
    assert (pub['pubmed_id'], pub['doi'], pub['title']) == ('123456', '10.1234/audit', 'Author boundary audit')


@pytest.mark.parametrize('stage', ['json', 'magetab'])
@pytest.mark.xfail(strict=True, reason='MSC-AUTH-002: GEO scalar characteristics disappear during v3 migration; audit excludes fixes')
def test_geo_retains_untagged_characteristic_on_its_channel(stage):
    from meta_standards_converter.miniml.geo_parser import GEOParser
    xml = '''<MINiML xmlns="http://www.ncbi.nlm.nih.gov/geo/info/MINiML">
    <Sample iid="GSM1"><Accession database="GEO">GSM1</Accession><Channel-Count>2</Channel-Count>
      <Channel position="1"><Source>treated</Source><Molecule>total RNA</Molecule>
        <Characteristics>exposure for 6 hours</Characteristics></Channel>
      <Channel position="2"><Source>control</Source><Molecule>total RNA</Molecule>
        <Characteristics>untreated reference</Characteristics></Channel></Sample>
    <Series iid="GSE1"><Accession database="GEO">GSE1</Accession><Title>Untyped annotation audit</Title>
      <Sample-Ref ref="GSM1"/></Series></MINiML>'''
    package = GEOParser().parse(xml)[0]
    if stage == 'json':
        channels = package.to_mapping()['sample'][0]['channel']
        actual = [[item['value'] for item in channel.get('characteristics', [])] for channel in channels]
        assert actual == [['exposure for 6 hours'], ['untreated reference']]
    else:
        result = Converter().convert(package, out_type='magetab', preparation=PreparationPolicy(
            enrichment='off', expand_studies=False, linked_metadata=False, peer_provider=False))
        assert result.items[0].execution == 'succeeded'
        sdrf = next(row[1] for row in result.items[0].payload[0] if row[0] == 'SDRF File')
        assert 'Characteristics[unspecified]' in sdrf[0]
        column = sdrf[0].index('Characteristics[unspecified]')
        assert [row[column] for row in sdrf[1:]] == ['exposure for 6 hours', 'untreated reference']


@pytest.mark.parametrize('stage', ['json', 'magetab'])
@pytest.mark.xfail(strict=True, reason='MSC-AUTH-003: repeated recognized IDF RelatedExperiment rows retain only the first row; audit excludes fixes')
def test_all_repeated_related_experiment_rows_survive_import_and_export(stage, tmp_path):
    source = tmp_path / 'E-MTAB-1.idf.txt'
    source.write_text('Investigation Accession\tE-MTAB-1\nInvestigation Title\tRepeated relationship audit\n'
        'Comment [RelatedExperiment]\tE-MTAB-2\nComment [RelatedExperiment]\tE-MTAB-3\n'
        'Comment [RelatedExperiment]\tE-MTAB-4\nSDRF File\tE-MTAB-1.sdrf.txt\n')
    (tmp_path / 'E-MTAB-1.sdrf.txt').write_text('Source Name\nSample 1\n')
    result = Converter().convert(source, out_type=stage, preparation=PreparationPolicy(
        enrichment='off', expand_studies=False, linked_metadata=False, peer_provider=False))
    assert result.items[0].execution == 'succeeded'
    if stage == 'json':
        package = result.items[0].payload[0].to_mapping()
        actual = [r['target'] for r in package['series'].get('relation', []) if r['type'] == 'related experiment']
    else:
        actual = [value for row in result.items[0].payload[0] if row[0] == 'Comment[RelatedExperiment]' for value in row[1:]]
    assert actual == ['E-MTAB-2', 'E-MTAB-3', 'E-MTAB-4']
