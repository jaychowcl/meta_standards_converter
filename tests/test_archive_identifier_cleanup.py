# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from copy import deepcopy

from meta_standards_converter.miniml import MINiMLCodec
from meta_standards_converter.miniml.model import PubMedPublication
from meta_standards_converter.metadata.enrichment import MINiMLEnricher
from tests.test_archive_publications import PubMed, NoRuns
from tests.test_native_archive_enrichment import native


def test_null_pmid_does_not_become_literal_none():
    assert PubMedPublication.from_mapping({'pubmed_id':None,'title':'Known paper'}).pubmed_id == ''


def test_invalid_saved_pmids_never_reach_lookup_and_remain_residual():
    data = native().to_mapping()
    data['series']['pubmed_id'] = ['None', 'https://doi.org/10.1234/paper']
    data['series']['pubmed_publication'] = [{'pubmed_id':'None','title':'Known paper'}]
    client = PubMed()
    service = MINiMLEnricher(pubmed_fetcher=client, insdc_fetcher=NoRuns())
    result = service.enrich(MINiMLCodec().decode(data).package).to_mapping()
    assert client.calls == []
    assert result['series'].get('pubmed_id',[]) == []
    assert result['series']['pubmed_publication'][0]['title'] == 'Known paper'
    assert result['series']['pubmed_publication'][0]['pubmed_id'] == ''
    assert 'None' in str(result['extensions']['insdc'])


def test_doi_resolution_is_cached_and_sample_citations_stay_sample_scoped():
    data = native().to_mapping()
    citation = {'pubmed_id':'', 'doi':'10.1234/source', 'title':'Supplied title'}
    data['sample'][0]['pubmed_publication'] = [deepcopy(citation)]
    data['sample'][0]['relation'].append({'type':'publication', 'target':'10.1234/source', 'publication':deepcopy(citation)})
    calls = []
    def resolve(kind, value):
        calls.append((kind,value)); return '123'
    service = MINiMLEnricher(pubmed_fetcher=PubMed(), insdc_fetcher=NoRuns(), publication_identifier_resolver=resolve)
    result = service.enrich(MINiMLCodec().decode(data).package).to_mapping()
    assert calls == [('doi','10.1234/source')]
    assert not result['series'].get('pubmed_id')
    assert result['sample'][0]['pubmed_publication'][0]['pubmed_id'] == '123'
    assert result['sample'][0]['pubmed_publication'][0]['title'] == 'Supplied title'


def test_status_only_publication_does_not_trigger_identifier_resolution():
    data = native().to_mapping()
    data['series']['pubmed_publication'] = [{'pubmed_id':'','status':'published'}]
    def forbidden(*args): raise AssertionError('No citation identifier was supplied')
    client = PubMed()
    service = MINiMLEnricher(pubmed_fetcher=client, insdc_fetcher=NoRuns(), publication_identifier_resolver=forbidden)
    service.enrich(MINiMLCodec().decode(data).package)
    assert not client.calls


def test_invalid_relation_keeps_scoped_citation_and_unknown_siblings():
    from meta_standards_converter.miniml.publication_identifiers import clean_publication_identifiers
    data = native().to_mapping()
    data['sample'][0]['relation'] = [{'type':'PubMed', 'target':'None', 'unknown':'retained',
        'publication':{'pubmed_id':'None','title':'Useful citation','doi':'10.1234/known'}}]
    residuals = clean_publication_identifiers(data)
    assert data['sample'][0]['pubmed_publication'][0]['title'] == 'Useful citation'
    assert data['sample'][0]['pubmed_publication'][0]['doi'] == '10.1234/known'
    assert 'retained' in str(residuals)
    assert not data['series'].get('pubmed_publication')


def test_conflicting_doi_and_pmcid_do_not_select_an_arbitrary_pmid():
    data = native().to_mapping()
    data['series']['pubmed_publication'] = [{'pubmed_id':'', 'pmcid':'PMC123', 'doi':'10.1234/other', 'title':'Source title'}]
    calls = []
    def resolve(kind, value):
        calls.append(kind)
        return {'pmcid':'111', 'doi':'222'}[kind]
    client = PubMed()
    service = MINiMLEnricher(pubmed_fetcher=client, insdc_fetcher=NoRuns(), publication_identifier_resolver=resolve)
    output = service.enrich(MINiMLCodec().decode(data).package).to_mapping()
    assert set(calls) == {'pmcid','doi'}
    assert not client.calls
    assert output['series']['pubmed_publication'][0]['pubmed_id'] == ''
    assert service.publication_issues


def test_local_identifier_cleanup_preserves_valid_repeated_source_occurrences():
    from meta_standards_converter.miniml.publication_identifiers import clean_publication_identifiers
    data={'source':{'format':'GEO'},'series':{'iid':'GSE1','pubmed_id':['123','None','123',' 456 ']}}
    records=clean_publication_identifiers(data)
    assert data['series']['pubmed_id']==['123','123','456']
    assert records[0]['metadata']['pubmed_id']=='None'
    before=deepcopy(data);assert not clean_publication_identifiers(data);assert data==before
