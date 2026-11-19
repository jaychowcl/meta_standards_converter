# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Source oracles must reject lost author components and mistaken run scope."""
import xml.etree.ElementTree as ET
import pytest
from tests.support.authorship_sources import citation_facts, author_representation_equivalent, native_coverage, native_contact_checks


def document(xml=None, data=None, url='https://example.org/source'):
    return {'kind': 'xml' if xml else 'json', 'xml': ET.fromstring(xml) if xml else None,
            'data': data, 'url': url, 'sha256': 'a' * 64, 'body': 'source.body',
            'member': None, 'reference_assisted': False}


def article():
    return document('''<PubmedArticle><MedlineCitation><PMID>1</PMID><Article>
    <ArticleTitle>A <i>scoped</i> title</ArticleTitle><AuthorList>
    <Author><LastName>Smith</LastName><Initials>AB</Initials></Author>
    <Author><LastName>Jones</LastName><ForeName>Alice</ForeName><Initials>A</Initials><Suffix>Jr</Suffix></Author>
    <Author><CollectiveName>Börsch Consortium</CollectiveName></Author>
    </AuthorList></Article></MedlineCitation><PubmedData><ArticleIdList>
    <ArticleId IdType="doi">10.1/a</ArticleId></ArticleIdList></PubmedData></PubmedArticle>''')


def test_raw_pubmed_oracle_preserves_initials_suffix_group_unicode_and_inline_title():
    pub = citation_facts([article()])['1'][0]
    assert pub['author_list'] == 'AB Smith, Alice Jones Jr, Börsch Consortium'
    assert pub['title'] == 'A scoped title' and pub['doi'] == '10.1/a'
    assert author_representation_equivalent(pub, 'Smith AB, Jones A Jr, Börsch Consortium.')


@pytest.mark.parametrize('actual', ['Smith, Alice Jones Jr, Börsch Consortium',
    'AB Smith, Alice Jones, Börsch Consortium', 'Alice Jones Jr, AB Smith, Börsch Consortium',
    'AB Smith, Alice Jones Jr, Borsch Consortium'])
def test_author_equivalence_does_not_hide_identity_order_suffix_or_unicode_loss(actual):
    assert not author_representation_equivalent(citation_facts([article()])['1'][0], actual)


def test_assembly_run_links_do_not_expand_the_read_study_denominator():
    docs = [document(data=[{'run_accession': 'ERR1', 'secondary_sample_accession': 'ERS1'}],
                     url='https://www.ebi.ac.uk/ena/portal/api/filereport?accession=ERP1&result=read_run'),
            document(data=[{'run_accession': 'ERR1;ERR2'}],
                     url='https://www.ebi.ac.uk/ena/portal/api/search?result=assembly'),
            document('<RUN_SET><RUN accession="ERR1"/><RUN accession="ERR2"/></RUN_SET>')]
    package = {'sample': [{'accession': [{'value': 'ERS1'}], 'sra_run': [{'run': 'ERR1'}]}]}
    result = native_coverage(docs, package)
    assert result['source_runs'] == ['ERR1']
    assert result['additional_retrieved_xml_runs_not_in_study_listing'] == ['ERR2']
    assert result['runs_equal'] and result['samples_equal']


def test_native_owner_contact_cannot_move_to_another_sample():
    docs = [document('''<BioSampleSet><BioSample accession="SAMN1"><Owner><Name>Lab A</Name>
      <Contacts><Contact email="a@example.org"><Name><First>Alice</First><Last>A</Last></Name></Contact></Contacts>
      </Owner></BioSample><BioSample accession="SAMN2"><Owner><Name>Lab B</Name>
      <Contacts><Contact><Name><First>Bob</First><Last>B</Last></Name></Contact></Contacts></Owner></BioSample></BioSampleSet>''')]
    package = {'sample': [
        {'accession': [{'value': 'SAMN1'}], 'contact_ref': [{'ref': 'b'}]},
        {'accession': [{'value': 'SAMN2'}], 'contact_ref': [{'ref': 'a'}]}],
        'contributor': [{'iid': 'a', 'person': {'first': 'Alice', 'last': 'A'}, 'email': 'a@example.org'},
                        {'iid': 'b', 'person': {'first': 'Bob', 'last': 'B'}}]}
    result = native_contact_checks(docs, package, [])
    assert len(result['checks']) == 5
    assert all(item['classification'] == 'source-supported-loss-alteration-or-misassociation' for item in result['checks'])


def test_organization_only_owner_is_not_a_person_or_paper_author():
    docs = [document('<BioSample accession="SAMN1"><Owner><Name>Research Institute</Name></Owner></BioSample>')]
    result = native_contact_checks(docs, {'sample': []}, [])
    assert not result['checks']
    assert result['organization_only_evidence'][0]['organization'] == 'Research Institute'
