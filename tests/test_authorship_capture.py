"""Capture bounds and immutable application-body replay controls."""
import io
import json
import requests
import pytest
from tests.support.authorship_capture import Capture, CaptureLimit
from tests.support.recorded_replay import ReplayTransport


def response(body=b'B\xc3\xb6rsch\n'):
    result = requests.Response()
    result.status_code = 200
    result.raw = io.BytesIO(body)
    result.headers['Content-Type'] = 'application/xml'
    result.headers['Set-Cookie'] = 'must-not-be-recorded'
    result.encoding = 'utf-8'
    return result


def test_capture_keeps_original_bytes_redacts_credentials_and_reuses_frozen_prior_phase(tmp_path):
    capture = Capture(tmp_path, 'synthetic', 'test')
    capture.original_send = lambda *a, **kw: response()
    request = requests.Request('GET', 'https://example.org/source?api_key=secret&id=1').prepare()
    capture.phase = 'off'
    first = capture.send(None, request)
    assert first.content == b'B\xc3\xb6rsch\n'
    capture.manifest('off')
    capture.phase = 'standard'
    assert capture.send(None, request).content == b'B\xc3\xb6rsch\n'
    capture.manifest('standard')
    capture.flush()
    assert capture.physical_requests == 1
    serialized = (tmp_path / 'trace.json').read_text()
    assert 'secret' not in serialized and 'must-not-be-recorded' not in serialized
    assert capture.events[1]['reused_from_phase'] == 'off'
    replay = ReplayTransport(tmp_path / 'standard.manifest.json')
    assert replay.send(None, request).content == b'B\xc3\xb6rsch\n'
    replay.assert_complete()


def test_capture_request_ceiling_is_global_to_phases(tmp_path):
    capture = Capture(tmp_path, 'synthetic', 'test', requests_limit=1)
    capture.original_send = lambda *a, **kw: response()
    request = requests.Request('GET', 'https://example.org/source').prepare()
    capture.phase = 'off'
    assert capture.send(None, request).content
    capture.phase = 'standard'
    with pytest.raises(CaptureLimit): capture.send(None, request)


def test_incomplete_stream_is_recorded_as_a_verification_limit_not_replayable_truth(tmp_path):
    capture = Capture(tmp_path, 'synthetic', 'test')
    capture.original_send = lambda *a, **kw: response()
    capture.phase = 'off'
    result = capture.send(None, requests.Request('GET', 'https://example.org/source').prepare())
    iterator = result.iter_content(2)
    next(iterator)
    result.close()
    capture.manifest('off')
    manifest = json.loads((tmp_path / 'off.manifest.json').read_text())
    assert manifest['incomplete_requests'] and not manifest['requests']
