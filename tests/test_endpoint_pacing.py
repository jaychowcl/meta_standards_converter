# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Request behavior across service paths, credentials, aliases and workers."""
import multiprocessing
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from unittest.mock import Mock

import pytest
import requests

from meta_standards_converter.helpers.request_helper import (
    HostRequestGate, RateLimitedRequester, RequestSettings,
)
from meta_standards_converter.sources.archive_support import ArchiveHTTP


ENA = "https://www.ebi.ac.uk/ena/xref/rest/json/search"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
PUBMED = "https://www.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"


def response(status=200, headers=None):
    value = requests.Response()
    value.status_code = status
    value.headers.update(headers or {})
    value._content = b"{}"
    value._content_consumed = True
    value.url = EUTILS
    return value


class Clock:
    now = 100.0

    def time(self):
        return self.now

    def sleep(self, seconds):
        self.now += seconds


@pytest.fixture
def transport(tmp_path):
    RateLimitedRequester.reset_service_state()
    clock = Clock()
    gate = HostRequestGate(tmp_path)

    def create(service="ena_portal", settings=None, get=None):
        return RateLimitedRequester(service, settings=settings, get=get or Mock(return_value=response()),
                                    host_gate=gate, clock=clock.time, sleep=clock.sleep)

    return clock, create


@pytest.mark.parametrize("url,interval", [
    (ENA, .1),
    ("https://www.ebi.ac.uk/ena/browser/api/xml/ERP001058", .1),
    ("https://www.ebi.ac.uk/ena/portal/api/search", .1),
    ("https://www.ebi.ac.uk/biosamples/samples/SAMEA1094641", 1),
    ("https://www.ebi.ac.uk/biostudies/api/v1/studies/E-MTAB-1", 1),
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/file.tgz", 1),
    (EUTILS, .5), (PUBMED, .5),
    ("https://www.ebi.ac.uk/other-api/data", 1),
    ("https://unrelated.example.org/ena/data", 1),
])
def test_url_selects_policy_even_when_requester_label_is_ena(transport, url, interval):
    clock, create = transport
    requester = create()
    requester.get(url)
    requester.get(url)
    assert clock.now == pytest.approx(100 + interval)
    assert requester.rate_wait_seconds == pytest.approx(interval)


@pytest.mark.parametrize("suffix,params", [
    ("?api_key=test-secret", None), ("", {"api_key": "test-secret"}),
    ("", [("api_key", "test-secret")]), ("", "api_key=test-secret"),
])
def test_keyed_rate_requires_key_on_the_actual_request(transport, suffix, params):
    clock, create = transport
    requester = create("ncbi_eutils")
    for _ in range(2):
        requester.get(EUTILS + suffix, params=params)
    assert clock.now == pytest.approx(100.125)


def test_environment_key_alone_does_not_grant_keyed_rate(transport, monkeypatch):
    clock, create = transport
    monkeypatch.setenv("NCBI_API_KEY", "test-secret")
    requester = create("ncbi_eutils")
    requester.get(EUTILS)
    requester.get(EUTILS, params={"api_key": ""})
    assert clock.now == pytest.approx(100.5)


def test_ebi_services_do_not_share_pacing_or_cooldowns(transport):
    clock, create = transport
    requester = create()
    requester.get(ENA)
    requester.get("https://www.ebi.ac.uk/biosamples/samples/SAMEA1094641")
    requester.get("https://www.ebi.ac.uk/biostudies/api/v1/studies/E-MTAB-1")
    assert clock.now == 100
    requester.get(ENA)
    assert clock.now == pytest.approx(100.1)


def test_ncbi_aliases_share_pacing_across_requesters(transport):
    clock, create = transport
    create("sra").get(EUTILS)
    create("pubmed").get(PUBMED)
    assert clock.now == pytest.approx(100.5)


def test_last_attempt_still_persists_retry_after_for_other_alias(transport):
    clock, create = transport
    failed = create(settings=RequestSettings(max_retries=0),
                    get=Mock(return_value=response(429, {"Retry-After": "7"})))
    with pytest.raises(requests.HTTPError):
        failed.get(EUTILS)
    create().get(PUBMED)
    assert clock.now == 107
    assert failed.provider_attempts == 1 and failed.retry_count == 0


def test_ena_cooldown_does_not_delay_biosamples(transport):
    clock, create = transport
    failed = create(settings=RequestSettings(max_retries=0),
                    get=Mock(return_value=response(429, {"Retry-After": "7"})))
    with pytest.raises(requests.HTTPError):
        failed.get(ENA)
    create().get("https://www.ebi.ac.uk/biosamples/samples/SAMEA1094641")
    assert clock.now == 100
    create().get(ENA)
    assert clock.now == 107


def test_redirect_hops_are_counted_without_an_outer_budget(transport):
    clock, create = transport
    get = Mock(side_effect=[response(302, {"Location": PUBMED}), response()])
    requester = create(get=get)
    requester.get(EUTILS, params={"id": "1"}, headers={"Authorization": "secret"})
    assert get.call_count == 2
    assert requester.provider_attempts == 2
    assert clock.now == pytest.approx(100.5)
    assert get.call_args.kwargs["allow_redirects"] is False
    assert "params" not in get.call_args.kwargs
    assert "Authorization" not in get.call_args.kwargs["headers"]


def test_explicit_spacing_override_remains_supported(transport):
    clock, create = transport
    requester = create(settings=RequestSettings(request_delay=.25))
    requester.get(ENA)
    requester.get(ENA)
    assert clock.now == pytest.approx(100.25)


def test_native_archive_resource_settings_use_endpoint_defaults(transport):
    clock, create = transport
    requester = create(settings=ArchiveHTTP("ena_portal").requester.settings)
    requester.get(ENA)
    requester.get(ENA)
    assert clock.now == pytest.approx(100.1)


def _alias_process(directory, barrier, output, url):
    requester = RateLimitedRequester("worker", host_gate=HostRequestGate(directory),
                                     get=lambda *a, **k: (output.put(time.monotonic()), response())[1])
    barrier.wait(timeout=5)
    requester.get(url)


def test_ncbi_alias_budget_is_shared_across_processes(tmp_path):
    ctx = multiprocessing.get_context("fork")
    barrier, output = ctx.Barrier(2), ctx.Queue()
    workers = [ctx.Process(target=_alias_process, args=(tmp_path, barrier, output, url))
               for url in (EUTILS, PUBMED)]
    for worker in workers:
        worker.start()
    try:
        timestamps = sorted(output.get(timeout=5) for _ in workers)
        assert timestamps[1] - timestamps[0] >= .48
    finally:
        for worker in workers:
            worker.join(timeout=5)
            if worker.is_alive():
                worker.terminate()
            assert worker.exitcode == 0


def test_concurrent_callers_share_one_ena_start_budget(tmp_path):
    RateLimitedRequester.reset_service_state()
    starts, lock = [], threading.Lock()

    def get(*args, **kwargs):
        with lock:
            starts.append(time.monotonic())
        return response()

    clients = [RateLimitedRequester("ena_portal", get=get, host_gate=HostRequestGate(tmp_path))
               for _ in range(4)]
    with ThreadPoolExecutor(max_workers=4) as pool:
        list(pool.map(lambda client: client.get(ENA), clients))
    assert all(b - a >= .095 for a, b in zip(sorted(starts), sorted(starts)[1:]))
