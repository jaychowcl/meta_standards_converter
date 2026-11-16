"""Bounded publication retrieval preserves source ownership and export fidelity."""
from concurrent.futures import ThreadPoolExecutor
from copy import deepcopy
from dataclasses import replace
from datetime import date
import json
from pathlib import Path
from threading import Barrier, Event, Lock
import time
from types import SimpleNamespace
from unittest.mock import Mock
import xml.etree.ElementTree as ET

import pytest

from meta_standards_converter.conversion_logging import current, operation
from meta_standards_converter.helpers.request_helper import RateLimitedRequester, RequestSettings
from meta_standards_converter.metadata.retrieval_budget import (
    _active, _budget_scope, _BudgetExhausted, _cancelled, _RetrievalCancelled,
)
from meta_standards_converter.runtime_contracts import get_resource_profile
from meta_standards_converter.sources.archive_support import ArchiveHTTP, StudyRecords, StudySeed
from meta_standards_converter.sources.sra import SRASource
from tests.test_native_archive_parsers import fixture_records


def link_response(params, *, citations=True):
    root = ET.Element('eLinkResult')
    for uid in params['id']:
        item = ET.SubElement(root, 'LinkSet')
        ET.SubElement(item, 'DbFrom').text = params['dbfrom']
        ET.SubElement(ET.SubElement(item, 'IdList'), 'Id').text = uid
        if citations and params['db'] == 'pubmed':
            group = ET.SubElement(item, 'LinkSetDb')
            ET.SubElement(group, 'DbTo').text = params['db']
            ET.SubElement(group, 'LinkName').text = params['linkname']
            ET.SubElement(ET.SubElement(group, 'Link'), 'Id').text = {
                'sra': '123', 'bioproject': '124', 'biosample': '125',
            }[params['dbfrom']]
    return root


def fixture_http(workers, hook=None):
    calls = []
    def get(url, params=None, fmt='xml'):
        calls.append((url, deepcopy(params)))
        if hook:
            hook(url, params)
        if 'esummary' in url:
            return {'result': {'10': {'uid': '10', 'expxml': '<Experiment acc="SRX7812918"/>'}}}
        return link_response(params)
    return SimpleNamespace(profile=SimpleNamespace(network_workers=workers), get=get, calls=calls)


def fixture_output(workers, hook=None):
    from meta_standards_converter.miniml.sra_parser import SRAParser
    from meta_standards_converter.magetab.constructor import AEConstructor
    records = fixture_records('sra')
    source = SRASource(http=fixture_http(workers, hook))
    source.linked_publications(records, ['10'])
    package = SRAParser().parse(records)
    return {
        'linked': records.linked,
        'issues': records.issues,
        'json': package.to_mapping(),
        'magetab': AEConstructor().miniml2magetab(package, platform_handler='generic'),
    }


def test_two_workers_overlap_and_finish_retrieval_before_owner_lookup():
    barrier = Barrier(2)
    lock = Lock()
    active = maximum = completed = 0
    def hook(url, params):
        nonlocal active, maximum, completed
        if 'esummary' in url:
            assert active == 0
            assert completed == 5
            return
        with lock:
            active += 1
            maximum = max(maximum, active)
            first_pair = completed == 0
        try:
            if first_pair:
                barrier.wait(timeout=1)
            time.sleep(.005)
        finally:
            with lock:
                active -= 1
                completed += 1
    records = fixture_records('sra')
    SRASource(http=fixture_http(8, hook)).linked_publications(records, ['10'])
    assert maximum == 2
    assert completed == 5
    assert records.issues == []
    assert [r['metadata']['pubmed_id'] for r in records.linked if r['kind'] == 'publication_reference'] == ['123', '124', '125']


def test_batch_sizes_and_source_order_survive_reverse_completion():
    ids = [str(n) for n in range(201)]
    started, finished = [], []
    lock, release_first = Lock(), Event()
    def links(dbfrom, db, batch, name):
        key = (dbfrom, db, tuple(batch), name)
        with lock:
            started.append(key)
            position = len(started)
        if position == 1:
            assert release_first.wait(1)
        elif position == 2:
            with lock:
                finished.append(key)
            release_first.set()
        if position != 2:
            with lock:
                finished.append(key)
        return {uid: [] for uid in batch}
    source = SRASource(http=SimpleNamespace(profile=SimpleNamespace(network_workers=2)))
    source.publication_links = links
    records = StudyRecords(StudySeed('SRP1', 'SRP1'))
    source.linked_publications(records, ids)
    assert len(started) == 6
    expected = [(db, tuple(sorted(ids)[start:start + 100]))
                for db in ('pubmed', 'pmc') for start in range(0, 201, 100)]
    assert sorted((k[1], k[2]) for k in started) == sorted(expected)
    assert finished[0] == started[1]
    assert records.issues == []


@pytest.mark.parametrize('workers', [1, 2])
def test_outputs_match_frozen_serial_json_and_magetab_golden(workers, monkeypatch):
    class FrozenDate(date):
        @classmethod
        def today(cls):
            return cls(2026, 10, 20)
    monkeypatch.setattr('meta_standards_converter.magetab.idf.date', FrozenDate)
    def hook(url, params):
        if workers == 2 and 'elink' in url and params['db'] == 'pubmed':
            time.sleep(.01)
    expected = json.loads((Path(__file__).parent / 'fixtures/sra_publication_concurrency.json').read_text())
    actual = json.loads(json.dumps(fixture_output(workers, hook)))
    assert actual == expected


def test_one_worker_preserves_serial_request_order_and_avoids_executor(monkeypatch):
    monkeypatch.setattr('concurrent.futures.ThreadPoolExecutor', Mock(side_effect=AssertionError('executor')))
    http = fixture_http(1)
    SRASource(http=http).linked_publications(fixture_records('sra'), ['10'])
    assert [('owner' if 'esummary' in url else params['dbfrom'] + '_' + params['db'])
            for url, params in http.calls] == [
        'sra_pubmed', 'owner', 'sra_pmc', 'bioproject_pubmed', 'bioproject_pmc', 'biosample_pubmed',
    ]


def test_injected_http_without_profile_keeps_serial_compatibility(monkeypatch):
    monkeypatch.setattr('concurrent.futures.ThreadPoolExecutor', Mock(side_effect=AssertionError('executor')))
    http = fixture_http(1)
    del http.profile
    records = fixture_records('sra')
    SRASource(http=http).linked_publications(records, ['10'])
    assert not records.issues


def test_partial_failures_keep_successful_owners_and_original_issue_order():
    def links(dbfrom, db, ids, name):
        if db == 'pmc':
            if dbfrom == 'sra':
                time.sleep(.02)
            raise ValueError(name)
        return {uid: [dict(sra='123', bioproject='124', biosample='125')[dbfrom]] for uid in ids}
    source = SRASource(http=fixture_http(2))
    source.publication_links = links
    records = fixture_records('sra')
    source.linked_publications(records, ['10'])
    assert records.issues == ['sra publication links: ValueError', 'bioproject publication links: ValueError']
    assert [str(e) for e in records.failures] == ['sra_pmc', 'bioproject_pmc']
    assert len([r for r in records.linked if r['kind'] == 'publication_reference']) == 3


def test_parent_cancellation_stops_before_any_injected_request():
    stop = Event()
    stop.set()
    source = SRASource(http=fixture_http(2))
    source.publication_links = Mock(return_value={})
    token = _cancelled.set(stop)
    try:
        with pytest.raises(_RetrievalCancelled):
            source.linked_publications(fixture_records('sra'), ['10'])
    finally:
        _cancelled.reset(token)
    source.publication_links.assert_not_called()


def test_worker_fatal_error_cancels_sibling_without_scheduling_more():
    entered = Event()
    calls = []
    def links(dbfrom, db, ids, name):
        calls.append(name)
        if db == 'pubmed':
            assert entered.wait(1)
            raise KeyboardInterrupt('fatal source failure')
        entered.set()
        from meta_standards_converter.metadata.retrieval_budget import _check_cancelled
        deadline = time.monotonic() + 2
        while time.monotonic() < deadline:
            _check_cancelled()
            time.sleep(.005)
        pytest.fail('sibling did not receive cancellation')
    source = SRASource(http=fixture_http(2))
    source.publication_links = links
    with pytest.raises(KeyboardInterrupt, match='fatal source failure'):
        source.linked_publications(fixture_records('sra'), ['10'])
    assert calls == ['sra_pubmed', 'sra_pmc']
    assert _cancelled.get() is None


def test_existing_cancellation_reaches_workers_during_retrieval():
    parent = Event()
    calls = []
    barrier = Barrier(2)
    def links(dbfrom, db, ids, name):
        from meta_standards_converter.metadata.retrieval_budget import _check_cancelled
        calls.append(name)
        barrier.wait(timeout=1)
        parent.set()
        _check_cancelled()
        pytest.fail('parent cancellation lost')
    source = SRASource(http=fixture_http(2))
    source.publication_links = links
    token = _cancelled.set(parent)
    try:
        with pytest.raises(_RetrievalCancelled):
            source.linked_publications(fixture_records('sra'), ['10'])
        assert _cancelled.get() is parent
    finally:
        _cancelled.reset(token)
    assert len(calls) == 2


def test_first_fatal_failure_survives_secondary_worker_cleanup_error():
    entered = Barrier(2)
    origin = KeyboardInterrupt('original provider failure')
    def links(dbfrom, db, ids, name):
        entered.wait(timeout=1)
        if db == 'pmc':
            raise origin
        deadline = time.monotonic() + 2
        while not _cancelled.get().is_set() and time.monotonic() < deadline:
            time.sleep(.005)
        assert _cancelled.get().is_set()
        raise SystemExit('secondary cleanup failure')
    source = SRASource(http=fixture_http(2))
    source.publication_links = links
    with pytest.raises(KeyboardInterrupt) as caught:
        source.linked_publications(fixture_records('sra'), ['10'])
    assert caught.value is origin


def test_request_budget_and_operation_context_are_shared():
    seen = []
    barrier = Barrier(2)
    response = Mock(status_code=200, headers={})
    requester = RateLimitedRequester('ncbi_eutils', RequestSettings(max_in_flight=2),
        get=Mock(return_value=response), host_gate=Mock(wait=Mock(return_value=0)))
    def links(dbfrom, db, ids, name):
        seen.append((_active.get(), current()))
        barrier.wait(timeout=1)
        requester.get(SRASource.base + 'elink.fcgi')
        return {}
    source = SRASource(http=fixture_http(2))
    source.publication_links = links
    with operation('sra-publications-test', heartbeat=False) as op, _budget_scope(max_requests=1) as budget:
        with pytest.raises(_BudgetExhausted):
            source.linked_publications(StudyRecords(StudySeed('SRP1', 'SRP1')), ['10'])
        assert seen == [(budget, op), (budget, op)]
        assert budget.requests == requester.provider_attempts == 1
        assert op.counters['provider_attempts'] == 1


def test_expired_deadline_prevents_injected_requests():
    clock = Mock(return_value=0.)
    source = SRASource(http=fixture_http(2))
    source.publication_links = Mock(return_value={})
    with _budget_scope(timeout_seconds=1, clock=clock):
        clock.return_value = 2.
        with pytest.raises(_BudgetExhausted):
            source.linked_publications(fixture_records('sra'), ['10'])
    source.publication_links.assert_not_called()


def test_executor_start_failure_restores_parent_cancellation(monkeypatch):
    monkeypatch.setattr('concurrent.futures.ThreadPoolExecutor',
                        Mock(side_effect=RuntimeError('executor unavailable')))
    parent = Event()
    token = _cancelled.set(parent)
    try:
        with pytest.raises(RuntimeError, match='executor unavailable'):
            SRASource(http=fixture_http(2)).linked_publications(fixture_records('sra'), ['10'])
        assert _cancelled.get() is parent
        assert not parent.is_set()
    finally:
        _cancelled.reset(token)


def test_distinct_cache_misses_overlap_and_successful_entries_are_reused():
    barrier = Barrier(2)
    calls = []
    def get(url, params, fmt='xml'):
        calls.append(params['db'])
        barrier.wait(timeout=1)
        return link_response(params, citations=False)
    source = SRASource(http=SimpleNamespace(get=get))
    with ThreadPoolExecutor(max_workers=2) as pool:
        futures = [pool.submit(source.publication_links, 'sra', db, ['10'], 'sra_' + db)
                   for db in ('pubmed', 'pmc')]
        assert [f.result() for f in futures] == [{'10': []}, {'10': []}]
    assert source.publication_links('sra', 'pubmed', ['10'], 'sra_pubmed') == {'10': []}
    assert sorted(calls) == ['pmc', 'pubmed']


def test_duplicate_sources_are_rejected_without_poisoning_cache():
    params = {'dbfrom': 'sra', 'db': 'pubmed', 'id': ['10'], 'linkname': 'sra_pubmed'}
    bad = link_response(params)
    bad.append(deepcopy(bad[0]))
    get = Mock(side_effect=[bad, link_response(params)])
    source = SRASource(http=SimpleNamespace(get=get))
    with pytest.raises(ValueError, match='duplicate link source'):
        source.publication_links('sra', 'pubmed', ['10'], 'sra_pubmed')
    assert source.publication_links('sra', 'pubmed', ['10'], 'sra_pubmed') == {'10': ['123']}
    assert get.call_count == 2


def test_concurrent_streams_share_evidence_and_close_responses(tmp_path):
    barrier, lock = Barrier(2), Lock()
    responses, payloads = [], []
    def get(url, params=None, **kwargs):
        raw = ET.tostring(link_response(params, citations=False))
        def chunks(**kw):
            barrier.wait(timeout=1)
            yield raw
        response = Mock(status_code=200, headers={}, iter_content=chunks)
        with lock:
            responses.append(response)
            payloads.append(raw)
        return response
    requester = RateLimitedRequester('ncbi_eutils', RequestSettings(max_in_flight=2),
        get=get, host_gate=Mock(wait=Mock(return_value=0)))
    profile = get_resource_profile(overrides={'network_workers': 2})
    http = ArchiveHTTP('ncbi_eutils', requester=requester, resource_profile=profile, evidence_dir=tmp_path)
    records = StudyRecords(StudySeed('SRP1', 'SRP1'))
    with _budget_scope(max_requests=2) as budget:
        SRASource(http=http).linked_publications(records, ['10'])
        assert budget.requests == requester.provider_attempts == 2
    assert records.issues == []
    assert http.downloaded_bytes == sum(map(len, payloads))
    assert all(r.close.call_count == 1 for r in responses)
    assert {p.read_bytes() for p in tmp_path.iterdir()} == set(payloads)


def test_concurrent_stream_byte_exhaustion_is_partial_and_closes_responses():
    barrier = Barrier(2)
    responses = []
    def get(url, params=None, **kwargs):
        def chunks(**kw):
            barrier.wait(timeout=1)
            yield b'x' * 24
        response = Mock(status_code=200, headers={}, iter_content=chunks)
        responses.append(response)
        return response
    requester = RateLimitedRequester('ncbi_eutils', RequestSettings(max_in_flight=2),
        get=get, host_gate=Mock(wait=Mock(return_value=0)))
    profile = replace(get_resource_profile(overrides={'network_workers': 2}),
                      max_aggregate_download_bytes=18)
    http = ArchiveHTTP('ncbi_eutils', requester=requester, resource_profile=profile)
    records = StudyRecords(StudySeed('SRP1', 'SRP1'))
    SRASource(http=http).linked_publications(records, ['10'])
    assert requester.provider_attempts == 2
    assert http.downloaded_bytes == 18
    assert len(records.failures) == 2
    assert all('XMLSizeLimitError' in issue for issue in records.issues)
    assert all(r.close.call_count == 1 for r in responses)
    assert records.linked == []


def test_stage_timings_distinguish_retrieval_ownership_and_merge():
    emitter = Mock()
    with operation('sra-publications-test', heartbeat=False, emitter=emitter):
        SRASource(http=fixture_http(2)).linked_publications(fixture_records('sra'), ['10'])
    stages = [call.kwargs.get('stage') for call in emitter.emit.call_args_list
              if call.args[0] == 'stage_finished']
    assert 'publication_link_retrieval' in stages
    assert 'publication_owner_resolution' in stages
    assert 'publication_link_merge' in stages
