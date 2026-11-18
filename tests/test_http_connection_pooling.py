# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Connection reuse must preserve the stateless metadata request contract."""
from concurrent.futures import ThreadPoolExecutor
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from threading import Barrier, Thread
from unittest.mock import Mock, patch
import socket

import pytest
import requests

from meta_standards_converter.conversion_logging import context_callable, logged_conversion
from meta_standards_converter.helpers.http_transport import transport_scope
from meta_standards_converter.helpers.request_helper import RateLimitedRequester, RequestSettings

_connect = socket.socket.connect


@pytest.fixture
def endpoint(monkeypatch):
    class Server(ThreadingHTTPServer):
        daemon_threads = True
        accepted = 0
        cookies = []

        def get_request(self):
            result = super().get_request()
            self.accepted += 1
            return result

    class Handler(BaseHTTPRequestHandler):
        protocol_version = 'HTTP/1.1'

        def do_GET(self):
            self.server.cookies.append(self.headers.get('Cookie'))
            body = b'ok'
            self.send_response(200)
            self.send_header('Content-Length', str(len(body)))
            self.send_header('Set-Cookie', 'server_state=must_not_persist; Path=/')
            self.end_headers()
            self.wfile.write(body)

        def log_message(self, *args):
            pass

    server = Server(('127.0.0.1', 0), Handler)
    # Permit only this test server; retain the repository's external-network ban.
    def local_connect(sock, address):
        assert address == ('127.0.0.1', server.server_port)
        return _connect(sock, address)
    monkeypatch.setattr(socket.socket, 'connect', local_connect)
    RateLimitedRequester.reset_service_state()
    thread = Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        yield server, f'http://127.0.0.1:{server.server_port}/'
    finally:
        server.shutdown()
        server.server_close()
        thread.join()


def client(**kwargs):
    return RateLimitedRequester('pool-test', RequestSettings(request_delay=0, max_in_flight=2),
                                host_gate=Mock(wait=Mock(return_value=0)), **kwargs)


def test_connections_reused_without_cookie_state(endpoint):
    server, url = endpoint
    with transport_scope():
        for _ in range(3):
            with client().get(url, stream=True) as response:
                assert response.content == b'ok'
    assert server.accepted == 1
    assert server.cookies == [None, None, None]


def test_scope_closes_sessions_on_exception_and_next_conversion_is_fresh(endpoint):
    server, url = endpoint
    sessions = []
    real = requests.Session

    def factory():
        session = real()
        session.close = Mock(wraps=session.close)
        sessions.append(session)
        return session

    @logged_conversion
    def convert(fail=False):
        with client().get(url) as response:
            assert response.content == b'ok'
        if fail:
            raise ValueError('test failure')

    with patch('meta_standards_converter.helpers.http_transport.requests.Session', factory):
        with pytest.raises(ValueError):
            convert(True)
        convert()
    assert len(sessions) == 2
    assert all(s.close.call_count == 1 for s in sessions)
    assert server.accepted == 2
    assert all(s.get_adapter(url).max_retries.total == 0 for s in sessions)


def test_context_copied_workers_have_separate_sessions(endpoint):
    server, url = endpoint
    barrier = Barrier(2)

    def worker():
        barrier.wait(timeout=5)
        for _ in range(2):
            with client().get(url) as response:
                assert response.content == b'ok'

    with transport_scope(), ThreadPoolExecutor(max_workers=2) as executor:
        futures = [executor.submit(context_callable(worker)) for _ in range(2)]
        for future in futures:
            future.result(timeout=5)
    assert server.accepted == 2


def test_injected_transport_unchanged():
    response = Mock(status_code=200, headers={})
    get = Mock(return_value=response)
    with transport_scope(), patch('requests.Session', side_effect=AssertionError('unexpected session')):
        assert client(get=get).get('https://example.org/') is response
    get.assert_called_once()


def test_terminal_throttle_closes_response():
    response = Mock(status_code=429, headers={'Retry-After': '1'})
    response.raise_for_status.side_effect = requests.HTTPError('429')
    requester = RateLimitedRequester('test', RequestSettings(max_retries=0),
                                    get=Mock(return_value=response), host_gate=Mock(wait=Mock(return_value=0)))
    with pytest.raises(requests.HTTPError):
        requester.get('https://example.org/')
    response.close.assert_called_once()


def test_pool_acquisition_is_bounded_by_connection_timeout(endpoint):
    from threading import Timer
    import time
    server, url = endpoint
    with transport_scope():
        first = client().get(url, stream=True)
        # Always release even on the buggy implementation; this reproducer cannot hang.
        rescue = Timer(.35, first.close)
        rescue.start()
        started = time.monotonic()
        try:
            with pytest.raises(requests.exceptions.ConnectionError):
                RateLimitedRequester("bounded-pool", RequestSettings(max_retries=0, request_delay=0),
                    host_gate=Mock(wait=Mock(return_value=0))).get(url, timeout=(.05, 2), stream=True)
            assert time.monotonic() - started < .3
        finally:
            first.close()
            rescue.cancel()
            rescue.join()


def test_pool_acquisition_honors_metadata_deadline_without_dispatch():
    from meta_standards_converter.helpers.http_transport import _HTTPPool
    from meta_standards_converter.metadata.retrieval_budget import _budget_scope, _BudgetExhausted
    pool = _HTTPPool('localhost', maxsize=1, block=True)
    pool.pool.get_nowait()
    with _budget_scope(timeout_seconds=.02) as budget:
        with pytest.raises(_BudgetExhausted):
            pool.urlopen("GET", "/")
        assert budget.requests == 0
    pool.close()


def test_pool_acquisition_honors_cancellation():
    from threading import Event, Timer
    from meta_standards_converter.helpers.http_transport import _HTTPPool
    from meta_standards_converter.metadata.retrieval_budget import _cancelled, _RetrievalCancelled
    pool = _HTTPPool('localhost', maxsize=1, block=True)
    pool.pool.get_nowait()
    event = Event()
    token = _cancelled.set(event)
    trigger = Timer(.02, event.set)
    trigger.start()
    try:
        with pytest.raises(_RetrievalCancelled):
            pool.urlopen("GET", "/", pool_timeout=.5)
    finally:
        trigger.join()
        _cancelled.reset(token)
        pool.close()
