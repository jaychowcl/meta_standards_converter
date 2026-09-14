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
import pytest
import requests
from meta_standards_converter.miniml import MINiMLCodec
from meta_standards_converter.miniml.sra_parser import SRAParser
from meta_standards_converter.miniml.ena_parser import ENAParser
from meta_standards_converter.metadata.enrichment import MINiMLEnricher
from tests.test_native_archive_parsers import fixture_records
from tests.test_native_archive_converters import Source
from tests.test_archive_residuals import nodes


class PubMed:
    def __init__(self, fail=False): self.calls=[]; self.fail=fail
    def pubmed_summary(self, pubmed_id):
        self.calls.append(pubmed_id)
        if self.fail: raise requests.ConnectionError('unavailable')
        return ('10.1/fetched','Fetched authors','Fetched title','published','EFO','EFO_0001796')


class NoRuns:
    def extract_sra_accessions(self, *args): raise AssertionError('native sequencing lookup')
    def fetch_sra_runs(self, *args, **kwargs): raise AssertionError('native sequencing lookup')


def package():
    d=SRAParser().parse(fixture_records('sra')).to_mapping()
    d['series']['pubmed_id']=['123','456']
    d['series']['pubmed_publication']=[{'pubmed_id':'123','title':'Native title','status':'published'},
                                      {'pubmed_id':'456','title':'Second title','status':'preprint','status_term_source_ref':'EFO','status_term_accession_number':'EFO_0010558'}]
    return MINiMLCodec().decode(d).package


def test_native_hydration_fills_only_missing_fields_without_sequencing_changes():
    original=package();client=PubMed();service=MINiMLEnricher(pubmed_fetcher=client,insdc_fetcher=NoRuns())
    data=service.enrich(original).to_mapping();pubs=data['series']['pubmed_publication']
    assert client.calls==['123','456']
    assert pubs[0]['title']=='Native title' and pubs[0]['doi']=='10.1/fetched'
    assert pubs[0]['status_term_accession_number']=='EFO_0001796'
    assert pubs[1]['status']=='preprint' and pubs[1]['status_term_accession_number']=='EFO_0010558'
    assert data['sample']==original.to_mapping()['sample']
    assert data['series']['assay_paths']==original.to_mapping()['series']['assay_paths']
    assert original.to_mapping()['series']['pubmed_publication'][0].get('doi') is None
    assert 'EFO' not in {d['iid'] for d in data['database']}
    assert pubs[0]['status_term_source_ref']=='EFO'
    service.enrich(MINiMLCodec().decode(data).package)
    assert client.calls==['123','456']


def test_missing_pubmed_ids_do_not_call_remote_and_failure_retains_native(caplog):
    client=PubMed(fail=True);service=MINiMLEnricher(pubmed_fetcher=client,insdc_fetcher=NoRuns())
    service.enrich(SRAParser().parse(fixture_records('sra')))
    assert not client.calls
    data=service.enrich(package()).to_mapping()
    assert data['series']['pubmed_publication'][0]['title']=='Native title'
    assert service.publication_issues==['PubMed 123: ConnectionError','PubMed 456: ConnectionError']
    assert 'PubMed 123' in caplog.text


@pytest.mark.parametrize('provider,parser',[('sra',SRAParser),('ena',ENAParser)])
def test_fetched_status_is_mapped_and_removed_from_residual(provider,parser):
    records=fixture_records(provider)
    study=next(n for root in records.xml for n in root.iter('STUDY'))
    ET.SubElement(ET.SubElement(study,'XREF_LINK'),'DB').text='PubMed'
    ET.SubElement(study.find('XREF_LINK'),'ID').text='123'
    records.xml.append(ET.fromstring('<PubmedArticleSet><PubmedArticle><MedlineCitation><PMID>123</PMID><Article><ArticleTitle>Title</ArticleTitle></Article></MedlineCitation><PubmedData><PublicationStatus>ppublish</PublicationStatus><Unmapped>keep</Unmapped></PubmedData></PubmedArticle></PubmedArticleSet>'))
    data=parser().parse(records).to_mapping()
    pub=data['series']['pubmed_publication'][0]
    assert pub['status']=='published' and pub['status_term_accession_number']=='EFO_0001796'
    assert not any(n.get('tag')=='PublicationStatus' for n in nodes(data['extensions']))
    assert 'keep' in str(data['extensions'])


@pytest.mark.parametrize('provider',['sra','ena'])
def test_converter_hydrates_after_linked_enrichment_and_reports_failure(provider):
    from meta_standards_converter.converters.sra2json import SRA2JSONConverter
    from meta_standards_converter.converters.ena2json import ENA2JSONConverter
    class Linked:
        def enrich(self, native):
            data=native.to_mapping();data['series']['pubmed_id']=['123']
            return MINiMLCodec().decode(data).package,[]
    service=MINiMLEnricher(pubmed_fetcher=PubMed(fail=True),insdc_fetcher=NoRuns())
    cls=SRA2JSONConverter if provider=='sra' else ENA2JSONConverter
    result=cls(source=Source(provider),linked_enricher=Linked(),publication_enricher=service).convert('SRP250911',enrich_from_geo_ae=True)
    assert result.packages and result.studies[0].status=='partial'
    assert 'PubMed 123: ConnectionError' in result.studies[0].issues


def test_publication_discovery_uses_explicit_identifiers_and_urls_only():
    from meta_standards_converter.sources.archive_support import publication_ids
    records=fixture_records('ena')
    records.xml.append(ET.fromstring('<PROJECT_SET><PROJECT><Publication id="123"><DbType>ePubmed</DbType></Publication><URL_LINK><URL>https://pubmed.ncbi.nlm.nih.gov/456/</URL></URL_LINK><TITLE>Paper about 789</TITLE></PROJECT></PROJECT_SET>'))
    records.linked.append({'provider':'ena','kind':'cross_references','accession':records.seed.study,'metadata':[{'url':'https://pubmed.ncbi.nlm.nih.gov/987/','description':'paper 999'}]})
    assert publication_ids(records)=={'123','456','987'}


def test_empty_pubmed_response_is_reported_as_unavailable():
    class Empty(PubMed):
        def pubmed_summary(self, pubmed_id): return (None,)*6
    service=MINiMLEnricher(pubmed_fetcher=Empty(),insdc_fetcher=NoRuns())
    service.enrich(package())
    assert service.publication_issues==['PubMed 123: ValueError','PubMed 456: ValueError']


def test_json2ae_native_publications_and_no_enrich_are_separate(tmp_path):
    import json
    from meta_standards_converter.converters.json2ae import JSON2AEConverter
    from meta_standards_converter.magetab.constructor import AEConstructor
    path=tmp_path/'native.json';path.write_text(json.dumps(package().to_mapping()));original=path.read_bytes()
    client=PubMed();service=MINiMLEnricher(pubmed_fetcher=client,insdc_fetcher=NoRuns())
    converter=JSON2AEConverter(enricher=service,ae_constructor=AEConstructor(pubmed_client=client,insdc_client=NoRuns()))
    rows=converter.convert(str(path),enrich=False)[0]
    assert client.calls==[] and path.read_bytes()==original
    rows=converter.convert(str(path))[0]
    assert client.calls==['123','456'] and path.read_bytes()==original
    assert next(r for r in rows if r[0]=='Publication DOI')[1:]==['10.1/fetched','10.1/fetched']
    assert next(r for r in rows if r[0]=='Publication Title')[1:]==['Native title','Second title']
