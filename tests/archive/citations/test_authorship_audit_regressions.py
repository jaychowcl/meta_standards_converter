# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Source-derived regressions for the defects reproduced by the authorship audit."""
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


def source_package(provider, article=ARTICLE, pmids=('123456',)):
    records = fixture_records(provider)
    study = next(node for root in records.xml for node in root.iter('STUDY'))
    links = ET.SubElement(study, 'STUDY_LINKS')
    for pmid in pmids:
        link = ET.SubElement(ET.SubElement(links, 'STUDY_LINK'), 'XREF_LINK')
        ET.SubElement(link, 'DB').text = 'PubMed'
        ET.SubElement(link, 'ID').text = pmid
    records.xml.append(ET.fromstring(article))
    return (SRAParser if provider == 'sra' else ENAParser)().parse(records)


@pytest.mark.parametrize('provider', ['sra', 'ena'])
@pytest.mark.parametrize('stage', ['json', 'magetab'])
@pytest.mark.parametrize('component,expected', [('initials', 'AB Smith'), ('suffix', 'Alice Jones Jr')])
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


@pytest.mark.parametrize('provider', ['sra', 'ena'])
@pytest.mark.parametrize('author,expected', [
    ('<Author><ForeName>Alice Quinn</ForeName><Initials>AQ</Initials><LastName>Jones</LastName><Suffix>III</Suffix></Author>', 'Alice Quinn Jones III'),
    ('<Author><ForeName> </ForeName><Initials>AJ</Initials><LastName>Börsch</LastName></Author>', 'AJ Börsch'),
    ('<Author><CollectiveName>Étude Consortium</CollectiveName><LastName>Unused</LastName></Author>', 'Étude Consortium'),
])
def test_author_projection_prefers_supplied_full_or_collective_name(provider, author, expected):
    root = ET.fromstring(ARTICLE)
    authors = root.find('.//AuthorList')
    authors.clear()
    authors.append(ET.fromstring(author))
    pub = source_package(provider, ET.tostring(root, encoding='unicode')).to_mapping()['series']['pubmed_publication'][0]
    assert pub['author_list'] == expected


@pytest.mark.parametrize('provider', ['sra', 'ena'])
def test_multiple_publications_keep_authors_on_their_paper_and_allow_blank_list(provider):
    root = ET.fromstring(ARTICLE)
    second = ET.fromstring(ARTICLE.replace('123456', '654321').replace('10.1234/audit', '10.1234/second'))[0]
    second.find('.//AuthorList').clear()
    root.append(second)
    data = source_package(provider, ET.tostring(root, encoding='unicode'), ('123456', '654321')).to_mapping()
    pubs = {p['pubmed_id']: p for p in data['series']['pubmed_publication']}
    assert pubs['123456']['author_list'] == 'AB Smith, Alice Jones Jr, Example Consortium, Anna Börsch'
    assert not pubs['654321'].get('author_list')
    assert [pubs[p]['doi'] for p in ('123456', '654321')] == ['10.1234/audit', '10.1234/second']


def test_mixed_scalar_characteristics_keep_occurrences_and_qualifiers():
    from meta_standards_converter.miniml.migration import MINiMLV1Migrator
    from copy import deepcopy
    tagged = {'name': 'duration', 'value': '6', 'unit': {'value': 'hour', 'term_accession_number': 'UO:0000032'},
              'term_source_ref': 'EFO', 'comments': [{'name': 'evidence', 'value': 'source'}]}
    data = {'source': {'format': 'GEO'}, 'series': {'iid': 'GSE1'}, 'sample': [{'iid': 'GSM1', 'channel': [
        {'characteristics': ['none', tagged, 'none', 0, None]}, {'characteristics': ['untreated reference']}]}]}
    before = deepcopy(data)
    migrated = MINiMLV1Migrator().migrate(data).package.to_mapping()
    assert data == before
    channels = migrated['sample'][0]['channel']
    assert channels[0]['characteristics'] == [{'name': 'unspecified', 'value': 'none'}, tagged,
        {'name': 'unspecified', 'value': 'none'}, {'name': 'unspecified', 'value': '0'}]
    assert channels[1]['characteristics'] == [{'name': 'unspecified', 'value': 'untreated reference'}]


def test_related_rows_keep_json_occurrences_and_export_unique_targets(tmp_path):
    source = tmp_path / 'E-MTAB-1.idf.txt'
    source.write_text('Investigation Accession\tE-MTAB-1\nInvestigation Title\tRepeated links\n'
        'Comment [RelatedExperiment]\tE-MTAB-2\t\tE-MTAB-3\n'
        'Comment[RelatedExperiment]\t E-MTAB-2 \tE-MTAB-4\n'
        'PubMed ID\t123\t456\nPublication DOI\t10.1234/first\t10.1234/second\n'
        'Publication Author List\tFirst Author\tSecond Author\nSDRF File\tE-MTAB-1.sdrf.txt\n')
    source.with_name('E-MTAB-1.sdrf.txt').write_text('Source Name\nSample 1\n')
    policy = PreparationPolicy(enrichment='off', expand_studies=False, linked_metadata=False, peer_provider=False)
    imported = Converter().convert(source, out_type='json', preparation=policy).items[0].payload[0]
    series = imported.to_mapping()['series']
    assert [r['target'] for r in series['relation'] if r['type'] == 'related experiment'] == [
        'E-MTAB-2', 'E-MTAB-3', 'E-MTAB-2', 'E-MTAB-4']
    assert [(p['pubmed_id'], p['doi'], p['author_list']) for p in series['pubmed_publication']] == [
        ('123', '10.1234/first', 'First Author'), ('456', '10.1234/second', 'Second Author')]
    exported = Converter().convert(imported, out_type='magetab', preparation=policy).items[0].payload[0]
    assert next(r[1:] for r in exported if r[0] == 'Comment[RelatedExperiment]') == [
        'E-MTAB-2', 'E-MTAB-3', 'E-MTAB-4']
