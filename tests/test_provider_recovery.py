# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from meta_standards_converter.sources.metadata import NativeMetadataSource
"""Provider recovery facts cross the public conversion boundary without prose parsing."""
from meta_standards_converter.converters import PreparationPolicy
import json
from types import SimpleNamespace
import xml.etree.ElementTree as ET

import pytest
import requests

from meta_standards_converter.converters import Converter
from meta_standards_converter.converters.unified.contracts import Diagnostic, ConversionItemResult
from meta_standards_converter.sources.archive_support import ArchiveHTTP
from meta_standards_converter.sources.ena import ENASource
from meta_standards_converter.converters.ena2json import ENA2JSONConverter


class Requester:
    def __init__(self, body=b'[broken', error=None):
        self.body, self.error = body, error
    def get(self, *args, **kwargs):
        if self.error:
            raise self.error
        response = requests.Response()
        response.status_code = 200
        response._content = self.body
        response.headers['Content-Length'] = str(len(self.body))
        response.iter_content = lambda **kw: iter([self.body])
        response.close = lambda: None
        return response


def issues(item):
    return [d['provider_issue'] for d in item.to_dict()['diagnostics'] if 'provider_issue' in d]


@pytest.mark.parametrize('error,reason', [
    (None, 'provider_response_invalid'),
    (requests.Timeout('secret URL'), 'provider_transport_failure'),
    (requests.HTTPError(response=SimpleNamespace(status_code=429)), 'provider_http_failure'),
])
def test_native_resolution_preserves_failure_cause(error, reason):
    source = ENASource(http=ArchiveHTTP('ena_portal', requester=Requester(error=error)))
    item = Converter(services={'ena_source': NativeMetadataSource('ena', source=source)}).convert('SRP250911', in_type='ena_accession', out_type='json', preparation=PreparationPolicy(enrichment='off')).items[0]
    found = issues(item)
    assert len(found) == 1
    assert found[0]['reason'] == reason
    assert found[0]['operation'] == 'study_resolution'
    assert found[0]['error']['provider'] == 'ena'
    assert found[0]['error']['retry_category'] == 'retryable'
    assert item.error.correlation_id == found[0]['error']['correlation_id']
    assert item.error.error_type != 'InputError'
    assert 'secret URL' not in json.dumps(item.to_dict())


def test_legacy_diagnostic_serialization_is_unchanged():
    item = ConversionItemResult('one', 'local', diagnostics=[Diagnostic('warning', 'text')])
    assert item.to_dict()['diagnostics'] == [{'code':'warning', 'message':'text', 'stage':'input', 'severity':'error'}]


def test_empty_provider_result_is_not_a_provider_failure():
    source = ENASource(http=ArchiveHTTP('ena_portal', requester=Requester(body=b'[]')))
    item = Converter(services={'ena_source': NativeMetadataSource('ena', source=source)}).convert('SRP250911', in_type='ena_accession', out_type='json', preparation=PreparationPolicy(enrichment='off')).items[0]
    assert issues(item) == []


@pytest.mark.parametrize('fmt,body', [('json', b'{bad'), ('xml', b'<bad>')])
def test_remote_parse_failure_but_not_local_parse_failure(fmt, body):
    from meta_standards_converter.provider_issues import collect_provider_issues
    item = ConversionItemResult('one', 'remote')
    with collect_provider_issues(item):
        with pytest.raises((ValueError, ET.ParseError)):
            ArchiveHTTP('ena_portal', requester=Requester(body)).get('https://www.ebi.ac.uk/ena/test', fmt=fmt)
    assert issues(item)[0]['reason'] == 'provider_response_invalid'
    local = Converter().convert('<bad>', in_type='geo_xml', out_type='json', preparation=PreparationPolicy(enrichment='off')).items[0]
    assert issues(local) == []


def test_handled_and_recovered_issues_do_not_make_success_partial():
    from meta_standards_converter.provider_issues import collect_provider_issues, provider_boundary, recover
    item = ConversionItemResult('one', 'remote', execution='succeeded', completeness='complete', status='complete')
    @provider_boundary(provider='ena', operation='run_inventory')
    def failed():
        raise requests.Timeout('ignored')
    with collect_provider_issues(item):
        try:
            failed()
        except requests.Timeout as error:
            recover(error)
        try:
            {}['optional']
        except KeyError:
            pass
    assert item.completeness == 'complete'
    assert [d.code for d in item.diagnostics] == ['provider_recovered']


def test_cached_partial_result_replays_its_structured_issues():
    from meta_standards_converter.provider_issues import collect_provider_issues, provider_boundary
    from meta_standards_converter.metadata.session import PreparationSession
    @provider_boundary(provider='ena', operation='study_resolution')
    def fail():
        raise requests.Timeout('ignored')
    def partial():
        try:
            fail()
        except requests.Timeout:
            return []
    results = [ConversionItemResult(str(i), 'remote') for i in range(2)]
    with PreparationSession() as session:
        for result in results:
            with collect_provider_issues(result):
                session.retrieve('shared', partial)
    assert len(issues(results[0])) == len(issues(results[1])) == 1
    assert issues(results[0]) == issues(results[1])


def test_optional_ena_file_failure_retains_native_fallback_and_reason():
    from meta_standards_converter.sources.insdc import INSDCWebfetcher
    from meta_standards_converter.provider_issues import collect_provider_issues
    item = ConversionItemResult('one', 'remote', execution='succeeded', status='complete', completeness='complete')
    with collect_provider_issues(item):
        assert INSDCWebfetcher(ena_requester=Requester()).fetch_ena_fastq_files('SRP1') == {}
    assert item.completeness == 'partial'
    assert issues(item)[0]['operation'] == 'run_inventory'


def test_biostudies_wrapped_parse_failure_keeps_cause(monkeypatch):
    from meta_standards_converter.sources.magetab import AEWebFetcher
    from meta_standards_converter.provider_issues import collect_provider_issues
    from meta_standards_converter.converters.unified.outcomes import failure
    source = AEWebFetcher()
    monkeypatch.setattr(source, '_request', lambda *a, **k: (Requester().get(), 'https://www.ebi.ac.uk/biostudies/test'))
    item = ConversionItemResult('one', 'remote')
    with collect_provider_issues(item):
        with pytest.raises(ValueError) as caught:
            source._fetch_json('https://www.ebi.ac.uk/biostudies/test')
    envelope, _ = failure(caught.value, stage='load', provider=None, item_id='one')
    assert envelope.error_type == 'JSONDecodeError'
    assert envelope.retry_category.value == 'retryable'
    assert len(issues(item)) == 1
    assert envelope.correlation_id == issues(item)[0]['error']['correlation_id']


@pytest.mark.parametrize('status,retryable', [(400,False),(429,True),(503,True)])
def test_http_classification(status, retryable):
    from meta_standards_converter.provider_issues import collect_provider_issues
    item = ConversionItemResult('one', 'remote')
    error = requests.HTTPError(response=SimpleNamespace(status_code=status))
    with collect_provider_issues(item):
        with pytest.raises(requests.HTTPError):
            ArchiveHTTP('ena_portal', requester=Requester(error=error)).get('https://www.ebi.ac.uk/ena/test')
    assert (issues(item)[0]['error']['retry_category'] == 'retryable') == retryable


def test_recovered_429_has_no_unresolved_issue(tmp_path):
    from meta_standards_converter.helpers.request_helper import RateLimitedRequester, RequestSettings, HostRequestGate
    from meta_standards_converter.provider_issues import collect_provider_issues
    responses = [requests.Response(), requests.Response()]
    responses[0].status_code, responses[1].status_code = 429, 200
    for response in responses:
        response.close = lambda: None
    replies = iter(responses)
    requester = RateLimitedRequester(service='ena_portal', get=lambda *a, **k: next(replies),
        settings=RequestSettings(max_retries=1, request_delay=0), sleep=lambda _:None,
        host_gate=HostRequestGate(directory=tmp_path))
    item = ConversionItemResult('one', 'remote')
    with collect_provider_issues(item):
        assert requester.get('https://www.ebi.ac.uk/ena/test').status_code == 200
    assert issues(item) == []
    assert requester.retry_count == 1


def test_optional_publication_provider_failure_has_operation():
    from meta_standards_converter.provider_issues import collect_provider_issues
    from meta_standards_converter.metadata.enrichment import MINiMLEnricher
    def unavailable(**kwargs):
        raise requests.Timeout('private request context')
    item = ConversionItemResult('one','remote')
    with collect_provider_issues(item):
        pub = MINiMLEnricher(pubmed_fetcher=SimpleNamespace(pubmed_summary=unavailable))._pubmed_publication('26221022')
    assert pub['pubmed_id'] == '26221022'
    assert issues(item)[0]['operation'] == 'publication_retrieval'
    assert issues(item)[0]['error']['provider'] == 'pubmed'


@pytest.mark.parametrize('url,provider', [
    ('https://www.ebi.ac.uk/biosamples/samples/SAMEA1','biosamples'),
    ('https://pmc.ncbi.nlm.nih.gov/tools/idconv/api/v1/articles/','pubmed'),
])
def test_auxiliary_provider_identity(url, provider):
    from meta_standards_converter.provider_issues import collect_provider_issues
    item=ConversionItemResult('one','remote')
    with collect_provider_issues(item):
        with pytest.raises(ValueError):
            ArchiveHTTP('ena_portal',requester=Requester()).get(url,fmt='json')
    assert issues(item)[0]['error']['provider']==provider


@pytest.mark.parametrize('params', [
    {'db': 'pubmed', 'api_key': 'test-secret'},
    [('db', 'pubmed'), ('api_key', 'test-secret')],
    'db=pubmed&api_key=test-secret', b'db=pubmed&api_key=test-secret',
])
def test_provider_failure_accepts_requests_parameter_forms_without_leaking_them(params):
    from meta_standards_converter.provider_issues import collect_provider_issues, provider_boundary
    class Client:
        @provider_boundary()
        def get(self, url, params=None):
            raise requests.Timeout('private transport detail')
    item = ConversionItemResult('one', 'remote')
    with collect_provider_issues(item), pytest.raises(requests.Timeout):
        Client().get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi', params=params)
    assert issues(item)[0]['error']['provider'] == 'pubmed'
    assert 'test-secret' not in str(issues(item))
    assert 'private transport detail' not in str(issues(item))
