# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Concurrent evidence retrieval preserves ordering, limits and partial work."""
from concurrent.futures import ThreadPoolExecutor
from copy import deepcopy
from dataclasses import replace
import json
from threading import Barrier
from threading import Event
import time
from types import SimpleNamespace
from unittest.mock import Mock

import pytest

from meta_standards_converter.conversion_logging import context_callable, operation
from meta_standards_converter.helpers.request_helper import RateLimitedRequester, RequestSettings
from meta_standards_converter.metadata.retrieval_budget import _active, _budget_scope, _BudgetExhausted
from meta_standards_converter.runtime_contracts import get_resource_profile
from meta_standards_converter.sources.archive_support import ArchiveHTTP, StudyRecords, StudySeed
from meta_standards_converter.sources.ena import ENASource
from meta_standards_converter.xml_safety import XMLSizeLimitError


def http(workers, barrier=None):
    def get(url, params=None, fmt='json'):
        if barrier:
            barrier.wait(timeout=2)
        if '/xref/' in url:
            return [{'Target Primary Accession': params['accession'], 'Source': 'PubMed', 'Source Primary Accession': '123'}]
        acc = url.rsplit('/', 1)[-1]
        return {'accession': acc, 'characteristics': {'organism': [{'text': 'Homo sapiens'}]}}
    return SimpleNamespace(profile=SimpleNamespace(network_workers=workers), get=get)


def records():
    return StudyRecords(StudySeed('ERP1', 'PRJEB1'))


def test_two_lanes_overlap_and_merge_in_source_order():
    source = ENASource(http=http(2, Barrier(2)))
    result = records()
    source._fetch_independent_links(result, ['ERP1'], ['SAMEA1'])
    assert result.issues == []
    assert [r['provider'] for r in result.linked] == ['ena', 'biosamples']
    assert [r['accession'] for r in result.linked] == ['ERP1', 'SAMEA1']


def test_serial_and_parallel_records_and_exports_equal():
    from tests.test_native_archive_parsers import fixture_records
    from meta_standards_converter.miniml.ena_parser import ENAParser
    from meta_standards_converter.magetab.constructor import AEConstructor
    outputs = []
    for workers in (1, 2):
        result = deepcopy(fixture_records('ena'))
        ENASource(http=http(workers))._fetch_independent_links(result, ['SRP250911'], ['SAMN14218700'])
        package = ENAParser().parse(result)
        outputs.append((package.to_mapping(), AEConstructor().miniml2magetab(package, platform_handler='generic')))
    assert outputs[0] == outputs[1]
    assert outputs[0][0]['series']['assay_paths']


def test_one_worker_never_creates_executor(monkeypatch):
    monkeypatch.setattr('concurrent.futures.ThreadPoolExecutor', Mock(side_effect=AssertionError('parallel')))
    result = records()
    ENASource(http=http(1))._fetch_independent_links(result, ['ERP1'], ['SAMEA1'])
    assert len(result.linked) == 2


def test_budget_and_logging_context_propagate_to_workers():
    barrier = Barrier(2)
    get = Mock(return_value=Mock(status_code=200, headers={}))
    requester = RateLimitedRequester('test', RequestSettings(request_delay=0, max_in_flight=2),
                                    get=get, host_gate=Mock(wait=Mock(return_value=0)))
    seen = []
    def fetch(url, params=None, fmt=None):
        seen.append(_active.get())
        barrier.wait(timeout=2)
        requester.get(url)
        return [] if '/xref/' in url else {'accession': 'SAMEA1'}
    source = ENASource(http=SimpleNamespace(profile=SimpleNamespace(network_workers=2), get=fetch))
    with operation('test', heartbeat=False) as op, _budget_scope(max_requests=1) as budget:
        with pytest.raises(_BudgetExhausted):
            source._fetch_independent_links(records(), ['ERP1'], ['SAMEA1'])
        assert seen == [budget, budget]
        assert budget.requests == 1
        assert op.counters['provider_attempts'] == 1
    assert get.call_count == 1


def test_independent_failure_preserves_other_lane_and_issue_order():
    transport = http(2)
    original = transport.get
    def get(url, params=None, fmt=None):
        if '/biosamples/' in url:
            raise ValueError('unavailable')
        return original(url, params, fmt)
    transport.get = get
    result = records()
    ENASource(http=transport)._fetch_independent_links(result, ['ERP1'], ['SAMEA1'])
    assert [r['provider'] for r in result.linked] == ['ena']
    assert result.issues == ['linked SAMEA1: ValueError']
    assert len(result.failures) == 1 and isinstance(result.failures[0], ValueError)


def test_concurrent_streams_share_byte_budget_and_close(tmp_path):
    barrier = Barrier(2)
    responses = []
    def get(*args, **kwargs):
        def chunks(**kw):
            barrier.wait(timeout=2)
            yield b'123456'
            yield b'789012'
        response = Mock(headers={}, iter_content=chunks)
        responses.append(response)
        return response
    profile = replace(get_resource_profile(), max_aggregate_download_bytes=18)
    archive = ArchiveHTTP('ena', requester=SimpleNamespace(get=get), resource_profile=profile, evidence_dir=tmp_path)
    with ThreadPoolExecutor(max_workers=2) as executor:
        futures = [executor.submit(archive.get, 'https://example.org/', fmt='text') for _ in range(2)]
        failures = []
        values = []
        for future in futures:
            try: values.append(future.result(timeout=3))
            except XMLSizeLimitError as exc: failures.append(exc)
    assert failures
    assert archive.downloaded_bytes <= 18
    assert all(r.close.call_count == 1 for r in responses)
    assert all(p.read_text() in values for p in tmp_path.iterdir())


def test_fatal_worker_cancels_sibling_cooldown_without_an_extra_request():
    waiting = Event()
    def wait(*args, sleep, **kwargs):
        waiting.set()
        sleep(30)
        return 30
    get = Mock()
    requester = RateLimitedRequester('test', RequestSettings(max_in_flight=2),
                                    get=get, host_gate=SimpleNamespace(wait=wait))
    source = ENASource(requester=requester, resource_profile=get_resource_profile(overrides={'network_workers': 2}))
    def failed(*args):
        assert waiting.wait(2)
        raise KeyboardInterrupt('stop')
    source.cross_references = failed
    started = time.monotonic()
    with pytest.raises(KeyboardInterrupt, match='stop'):
        source._fetch_independent_links(records(), ['ERP1'], ['SAMEA1', 'SAMEA2'])
    assert time.monotonic() - started < 2
    get.assert_not_called()


def test_concurrent_metrics_and_request_budget_are_exact(tmp_path):
    response = Mock(status_code=200, headers={})
    requester = RateLimitedRequester('test', RequestSettings(request_delay=0, max_in_flight=2),
                                    get=Mock(return_value=response), host_gate=Mock(wait=Mock(return_value=.25)))
    with operation('test', heartbeat=False) as op, _budget_scope(max_requests=100) as budget:
        with ThreadPoolExecutor(max_workers=2) as executor:
            futures = [executor.submit(context_callable(requester.get), 'https://example.org/') for _ in range(100)]
            assert all(f.result() is response for f in futures)
        assert budget.requests == requester.provider_attempts == op.counters['provider_attempts'] == 100
        assert requester.rate_wait_seconds == op.counters['rate_wait_seconds'] == 25
