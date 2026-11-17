# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from meta_standards_converter.sources.metadata import NativeMetadataSource
"""Optional top-level budgets count network attempts and retain valid candidates."""
from meta_standards_converter.converters import PreparationPolicy, RuntimeOptions
from unittest.mock import Mock
import pytest
import requests
from meta_standards_converter.converters import Converter, InputSpec
from meta_standards_converter.helpers.request_helper import RateLimitedRequester, RequestSettings, HostRequestGate
from tests.converters.test_audit_source_precedence import package, NoPublication


@pytest.fixture
def requester(tmp_path):
    RateLimitedRequester.reset_service_state()
    response = Mock(status_code=200, headers={})
    call = Mock(return_value=response)
    return RateLimitedRequester('test', get=call, settings=RequestSettings(request_delay=0, max_retries=2), host_gate=HostRequestGate(tmp_path))


@pytest.mark.parametrize('name,value', [('metadata_max_requests',0), ('metadata_max_requests',True), ('metadata_max_requests',1.5), ('metadata_timeout_seconds',-1), ('metadata_timeout_seconds',float('nan')), ('metadata_timeout_seconds',True)])
def test_options_must_be_positive(name, value):
    with pytest.raises((TypeError, ValueError), match=name):
        RuntimeOptions(**{name:value})


def test_retries_and_redirects_consume_actual_attempts(requester):
    from meta_standards_converter.metadata.retrieval_budget import _budget_scope, _BudgetExhausted
    requester._get.side_effect = [requests.Timeout(), Mock(status_code=302, headers={'Location':'https://example.test/next'}), Mock(status_code=200, headers={})]
    with _budget_scope(max_requests=2) as budget:
        with pytest.raises(_BudgetExhausted): requester.get('https://example.test/start')
        assert budget.requests == 2
    assert requester._get.call_count == 2


def test_nested_scopes_share_allowance_and_new_inputs_are_isolated(requester):
    from meta_standards_converter.metadata.retrieval_budget import _budget_scope, _BudgetExhausted
    with _budget_scope(max_requests=1) as outer:
        requester.get('https://example.test/one')
        with _budget_scope(max_requests=100) as inner:
            assert inner is outer
            with pytest.raises(_BudgetExhausted): requester.get('https://example.test/two')
    with _budget_scope(max_requests=1) as other:
        requester.get('https://example.test/three')
        assert other.requests == 1


def test_deadline_clamps_request_timeout_and_provider_wait(requester):
    from meta_standards_converter.metadata.retrieval_budget import _budget_scope, _BudgetExhausted
    now = [100.0]
    requester._wall_clock = lambda: now[0]
    requester._sleep = lambda delay: now.__setitem__(0, now[0] + delay)
    with _budget_scope(timeout_seconds=2, clock=lambda: now[0]):
        requester.get('https://example.test/one', timeout=(10, 20))
        timeout = requester._get.call_args.kwargs['timeout']
        assert timeout.total <= 2
        requester._host_gate.defer('example.test', delay_seconds=10, clock=requester._wall_clock)
        with pytest.raises(_BudgetExhausted): requester.get('https://example.test/two')
        assert requester._get.call_count == 1
        assert now[0] <= 102


def test_enrichment_cutoff_retains_valid_input_and_persists_coverage(requester):
    class Enricher:
        def enrich_selected(self, data, **kwargs):
            requester.get('https://example.test/one')
            requester.get('https://example.test/two')
            raise AssertionError('unreachable')
    result = Converter(services={'enricher': Enricher()}).convert(package('MAGE-TAB'), out_type='json', runtime=RuntimeOptions(**{'metadata_max_requests': 1}), preparation=PreparationPolicy(expand_studies=False, enrichment='curators'))
    item = result.items[0]
    assert item.status == 'partial' and item.validation == 'valid', result.to_dict()
    assert item.payload[0]['sample'][0]['iid'] == 'SRS1'
    assert any('budget' in d.message and 'standard' in d.message for d in item.diagnostics)
    saved = Converter().convert(item.payload[0], out_type='json', preparation=PreparationPolicy(enrichment='off'))
    assert saved.status == 'partial'


def test_core_cutoff_fails_without_publication_but_supplementary_cutoff_retains_core(requester, tmp_path):
    from meta_standards_converter.converters.sra2json import SRA2JSONConverter
    from meta_standards_converter.sources.archive_support import Resolution, StudySeed, StudyRecords
    from meta_standards_converter.miniml import MINiMLCodec
    class Source:
        def resolve(self, accession): return Resolution([StudySeed(accession,accession)])
        def fetch(self, seed, *, _on_core=None, publications=True):
            requester.get('https://example.test/core')
            requester.get('https://example.test/core2')
            records = StudyRecords(seed); records.indexed['read_run'] = [{'run_accession':'SRR1'}]
            if _on_core: _on_core(records)
            requester.get('https://example.test/supplement')
            return records
    class Parser:
        def parse(self, records): return MINiMLCodec().decode(package(),strict=True).package
    facade = Converter(services={'sra_source': NativeMetadataSource('sra', source=Source(),parser=Parser())})
    failed = facade.convert('SRP1', in_type='sra', out_type='json', outdir=tmp_path / 'failed', runtime=RuntimeOptions(metadata_max_requests=1), preparation=PreparationPolicy(enrichment='off', expand_studies=False))
    assert failed.status == 'failed' and not failed.items[0].artifacts
    result = facade.convert('SRP1', in_type='sra', out_type='json', runtime=RuntimeOptions(metadata_max_requests=2), preparation=PreparationPolicy(enrichment='off', expand_studies=False))
    assert result.status == 'partial' and result.items[0].validation == 'valid', result.to_dict()
    assert any('supplement' in d.message and 'budget' in d.message for d in result.items[0].diagnostics)


def test_batch_input_budget_isolation(requester):
    class Enricher:
        def enrich_selected(self, data, **kwargs):
            requester.get('https://example.test/publication')
            return data
    result = Converter(services={'enricher': Enricher()}).convert([InputSpec(package('MAGE-TAB'), id='a'), InputSpec(package('GEO'), id='b')], out_type='json', runtime=RuntimeOptions(metadata_max_requests=1), preparation=PreparationPolicy(enrichment='curators', expand_studies=False))
    assert [i.status for i in result.items] == ['complete','complete'], result.to_dict()
    assert requester._get.call_count == 2


def test_cached_metadata_consumes_no_requests_and_redirects_do(tmp_path):
    from meta_standards_converter.metadata.retrieval_budget import _budget_scope, _BudgetExhausted
    from meta_standards_converter.retrieval import RetrievalService
    from tests.test_retrieval import _policy, _Session, _Response
    session = _Session([_Response(status_code=302, headers={'Location':'/final'}), _Response()])
    service = RetrievalService(tmp_path, policy=_policy(), session=session)
    with _budget_scope(max_requests=2) as budget:
        path = service.localize('https://data.example.org/metadata')
        assert budget.requests == 2
        assert service.localize('https://data.example.org/metadata') == path
        assert budget.requests == 2
        with pytest.raises(_BudgetExhausted): service.localize('https://data.example.org/other')


def test_budget_partial_memo_does_not_poison_a_later_input(requester):
    from meta_standards_converter.metadata.session import PreparationSession
    from meta_standards_converter.metadata.retrieval_budget import _budget_scope, _BudgetExhausted
    def optional():
        try:
            requester.get('https://example.test/one'); requester.get('https://example.test/two')
            return 'complete'
        except _BudgetExhausted: return 'partial'
    with PreparationSession() as session:
        with _budget_scope(max_requests=1): assert session.retrieve('same', optional) == 'partial'
        with _budget_scope(max_requests=2): assert session.retrieve('same', optional) == 'complete'


def test_deadline_during_body_cannot_cache_partial_bytes(tmp_path):
    from meta_standards_converter.metadata.retrieval_budget import _budget_scope, _BudgetExhausted
    from meta_standards_converter.retrieval import RetrievalService
    from tests.test_retrieval import _policy, _Session, _Response
    now = [0.0]
    class Slow(_Response):
        def iter_content(self, **kwargs):
            yield b'a'; now[0] = 2; yield b'bc'
    service = RetrievalService(tmp_path, policy=_policy(), session=_Session([Slow()]))
    with _budget_scope(timeout_seconds=1, clock=lambda:now[0]):
        with pytest.raises(_BudgetExhausted): service.localize('https://data.example.org/metadata')
    assert not service._destination('https://data.example.org/metadata').exists()


def test_ena_partial_xref_is_not_cached_as_complete():
    from meta_standards_converter.metadata.retrieval_budget import _BudgetExhausted
    from meta_standards_converter.sources.ena import ENASource
    from meta_standards_converter.sources.archive_support import StudyRecords, StudySeed
    http = Mock(); source = ENASource(http=http); source.xref_page_size = 1
    row = {'Target Primary Accession':'ERP1','Source Primary Accession':'E-MTAB-1'}
    http.get.side_effect = [[row], _BudgetExhausted('cutoff'), [row], []]
    records = StudyRecords(StudySeed('ERP1','ERP1'))
    with pytest.raises(_BudgetExhausted): source.cross_references('ERP1',records)
    assert 'ERP1' not in source._xrefs
    source.cross_references('ERP1', records)
    assert source._xrefs['ERP1'] == ([row], [], [])


def test_unbounded_options_keep_existing_fetch_signature():
    from tests.test_native_archive_converters import Source as RecordedSource
    class Source(RecordedSource):
        def fetch(self, seed, *, publications=True):
            assert publications is False
            return super().fetch(seed)
    from meta_standards_converter.converters.sra2json import SRA2JSONConverter
    result = Converter(services={'sra_source': NativeMetadataSource('sra', source=Source())}).convert('SRP250911', in_type='sra', out_type='json', runtime=RuntimeOptions(metadata_max_requests=None, metadata_timeout_seconds=None), preparation=PreparationPolicy(enrichment='off', expand_studies=False))
    assert result.items[0].payload and result.items[0].validation == 'valid'


def test_native_sources_offer_core_before_supplementary_discovery():
    import xml.etree.ElementTree as ET
    from meta_standards_converter.sources.ena import ENASource
    from meta_standards_converter.sources.sra import SRASource
    from meta_standards_converter.sources.archive_support import StudySeed
    from tests.test_native_archive_sources import HTTP
    class CoreSeen(Exception): pass
    def ena_http(url, params, fmt):
        assert params.get('result') not in {'analysis','assembly'}
        assert 'xref' not in url and 'biosamples' not in url
        if '/xml/' in url: return ET.fromstring('<ROOT/>')
        if 'returnFields' in url: return 'columnId\nrun_accession\n'
        return [{'run_accession':'ERR1'}] if params.get('result') == 'read_run' else []
    def ready(records):
        assert records.indexed.get('read_run') or records.xml
        raise CoreSeen()
    with pytest.raises(CoreSeen): ENASource(http=HTTP(ena_http)).fetch(StudySeed('ERP1','PRJEB1'),_on_core=ready)
    def sra_http(url, params, fmt):
        assert params.get('db') == 'sra'
        if 'esearch' in url: return {'esearchresult':{'count':'1','idlist':['1']}}
        return ET.fromstring('<EXPERIMENT_PACKAGE_SET><EXPERIMENT_PACKAGE><STUDY accession="SRP1"/><EXPERIMENT accession="SRX1"><STUDY_REF accession="SRP1"/></EXPERIMENT></EXPERIMENT_PACKAGE></EXPERIMENT_PACKAGE_SET>')
    with pytest.raises(CoreSeen): SRASource(http=HTTP(sra_http)).fetch(StudySeed('SRP1','SRP1'),_on_core=ready)
