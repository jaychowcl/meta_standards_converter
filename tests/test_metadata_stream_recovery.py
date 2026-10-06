"""A successful HTTP status does not establish a complete metadata response."""
from dataclasses import replace
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from threading import Thread
from types import SimpleNamespace
from unittest.mock import Mock
import socket
import time

import pytest
import requests

from meta_standards_converter.helpers.request_helper import RateLimitedRequester, RequestSettings
from meta_standards_converter.metadata.retrieval_budget import _budget_scope, _BudgetExhausted
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


def test_body_failure_retries_in_the_same_request_budget_and_retains_byte_charges(tmp_path):
    first = response([b'<SET>', requests.exceptions.ChunkedEncodingError('truncated')])
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


def test_metadata_response_deadline_defaults():
    assert get_resource_profile('standard').metadata_response_timeout_seconds == 300
    assert get_resource_profile('large').metadata_response_timeout_seconds == 900


def test_slow_trickle_has_a_wall_clock_deadline_and_releases_connection(monkeypatch):
    profile = replace(get_resource_profile(), metadata_response_timeout_seconds=.15)
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
        with pytest.raises(requests.Timeout, match='metadata response deadline'):
            ArchiveHTTP('ena_portal', requester=requester, resource_profile=profile).get(
                f'http://127.0.0.1:{server.server_port}/', fmt='text')
        assert time.monotonic() - started < 1.0
        assert requester.provider_attempts == 1
    finally:
        server.shutdown(); server.server_close(); thread.join()
