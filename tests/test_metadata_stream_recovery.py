# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""A successful HTTP status does not establish a complete metadata response."""
from dataclasses import replace
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from threading import Thread, Event, Timer
from types import SimpleNamespace
from unittest.mock import Mock
import socket
import time

import pytest
import requests

from meta_standards_converter.helpers.request_helper import RateLimitedRequester, RequestSettings
from meta_standards_converter.metadata.retrieval_budget import _budget_scope, _BudgetExhausted, _cancelled, _RetrievalCancelled
from meta_standards_converter.runtime_contracts import get_resource_profile
from meta_standards_converter.sources.archive_support import ArchiveHTTP

_connect = socket.socket.connect


def response(chunks):
    def stream(**kwargs):
        for item in chunks:
            if isinstance(item, BaseException):
                raise item
            yield item
    return Mock(status_code=200, headers={}, iter_content=stream)


def client(get, retries=3):
    return RateLimitedRequester('ena_portal', RequestSettings(request_delay=0, max_retries=retries),
        get=get, host_gate=Mock(wait=Mock(return_value=0)), random_value=lambda: 0)


@pytest.mark.parametrize('error',[requests.exceptions.ChunkedEncodingError,requests.exceptions.ContentDecodingError])
def test_body_failure_retries_in_the_same_request_budget_and_retains_byte_charges(tmp_path,error):
    first = response([b'<SET>', error('truncated')])
    second = response([b'<SET/>'])
    get = Mock(side_effect=[first, second]); requester = client(get)
    archive = ArchiveHTTP('ena_portal', requester=requester, evidence_dir=tmp_path)
    assert archive.get('https://www.ebi.ac.uk/ena/browser/api/xml/ERP1').tag == 'SET'
    assert requester.provider_attempts == 2 and requester.retry_count == 1
    assert archive.downloaded_bytes == 11
    assert first.close.call_count == second.close.call_count == 1
    assert [p.read_bytes() for p in tmp_path.iterdir()] == [b'<SET/>']


def test_connection_and_body_failures_share_four_attempts():
    responses = [response([requests.exceptions.ChunkedEncodingError('truncated')]) for _ in range(3)]
    get = Mock(side_effect=[requests.ConnectionError('connect'), *responses])
    requester = client(get)
    with pytest.raises(requests.exceptions.ChunkedEncodingError):
        ArchiveHTTP('ena_portal', requester=requester).get('https://www.ebi.ac.uk/ena/browser/api/xml/ERP1')
    assert get.call_count == requester.provider_attempts == 4
    assert requester.retry_count == 3
    assert all(r.close.call_count == 1 for r in responses)


def test_body_retry_obeys_operation_request_allowance():
    first = response([requests.exceptions.ChunkedEncodingError('truncated')])
    get = Mock(return_value=first)
    with _budget_scope(max_requests=1), pytest.raises(_BudgetExhausted):
        ArchiveHTTP('ena_portal', requester=client(get)).get('https://www.ebi.ac.uk/ena/browser/api/xml/ERP1')
    assert get.call_count == 1 and first.close.call_count == 1


def test_injected_requester_without_consumption_extension_is_supported():
    value = response([b'<SET/>'])
    requester = SimpleNamespace(get=Mock(return_value=value))
    assert ArchiveHTTP('ena_portal', requester=requester).get('https://example.org').tag == 'SET'
    assert value.close.call_count == 1
    assert set(requester.get.call_args.kwargs) == {'params', 'stream'}


def test_oversized_chunk_exhausts_remaining_credit_before_any_later_request():
    profile=replace(get_resource_profile(),max_aggregate_download_bytes=4)
    get=Mock(return_value=response([b'<SET/>']))
    archive=ArchiveHTTP('ena',requester=client(get),resource_profile=profile)
    with pytest.raises(ValueError,match='aggregate'):
        archive.get('https://example.org/first')
    with pytest.raises(ValueError,match='budget exhausted'):
        archive.get('https://example.org/second')
    assert get.call_count==1


def test_metadata_response_deadline_defaults():
    assert get_resource_profile('standard').metadata_response_timeout_seconds == 300
    assert get_resource_profile('large').metadata_response_timeout_seconds == 900


def test_body_transfer_is_timed_inside_a_real_conversion_operation(caplog):
    import logging
    from meta_standards_converter.conversion_logging import operation
    from tests.test_conversion_logging import events
    with caplog.at_level(logging.DEBUG),operation('stream',heartbeat=False):
        value=ArchiveHTTP('ena',requester=client(Mock(return_value=response([b'<SET/>'])))).get('https://example.org')
    assert value.tag=='SET'
    assert any(e['stage']=='body_transfer' and e['event_type']=='stage_finished' for e in events(caplog))


@pytest.mark.parametrize('cancel',[False,True])
def test_slow_trickle_has_a_wall_clock_deadline_and_releases_connection(monkeypatch,cancel):
    profile = replace(get_resource_profile(), metadata_response_timeout_seconds=1.5 if cancel else .15)
    stopped=Event();token=_cancelled.set(stopped)
    canceller=Timer(.05,stopped.set) if cancel else None
    class Handler(BaseHTTPRequestHandler):
        def do_GET(self):
            self.send_response(200); self.send_header('Content-Length', '1000'); self.end_headers()
            try:
                for _ in range(1000):
                    self.wfile.write(b'x'); self.wfile.flush(); time.sleep(.02)
            except (BrokenPipeError, ConnectionResetError):
                pass
        def log_message(self, *args): pass
    server = ThreadingHTTPServer(('127.0.0.1', 0), Handler)
    server.daemon_threads = True
    def local_connect(sock, address):
        assert address == ('127.0.0.1', server.server_port)
        return _connect(sock, address)
    monkeypatch.setattr(socket.socket, 'connect', local_connect)
    thread = Thread(target=server.serve_forever, daemon=True); thread.start()
    requester = client(requests.get, retries=0)
    started = time.monotonic()
    try:
        if canceller:canceller.start()
        with pytest.raises(_RetrievalCancelled if cancel else requests.Timeout):
            ArchiveHTTP('ena_portal', requester=requester, resource_profile=profile).get(
                f'http://127.0.0.1:{server.server_port}/', fmt='text')
        assert time.monotonic() - started < 1.0
        assert requester.provider_attempts == 1
    finally:
        if canceller:canceller.cancel()
        _cancelled.reset(token)
        server.shutdown(); server.server_close(); thread.join()


def test_deadline_does_not_shutdown_a_socket_already_returned_to_a_pool():
    from types import SimpleNamespace
    from meta_standards_converter.helpers.metadata_stream import response_deadline
    owned, peer=socket.socketpair()
    value=SimpleNamespace(raw=SimpleNamespace(_connection=SimpleNamespace(sock=owned)))
    try:
        with pytest.raises(requests.Timeout),response_deadline(value,.02):
            # urllib3 has consumed the response and released its connection.
            # An injected shared session may now lend that socket elsewhere.
            value.raw._connection=None
            time.sleep(.06)
        peer.sendall(b'usable')
        assert owned.recv(6)==b'usable'
    finally:
        owned.close();peer.close()
