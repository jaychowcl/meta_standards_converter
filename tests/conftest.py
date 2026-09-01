# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Fail-closed test guard for accidental network and child-process access."""

import os
from pathlib import Path
import socket
import subprocess
from urllib.parse import urlsplit

import pytest
import requests


_LIVE_ENV = "RUN_LIVE_API_TESTS"
_LIVE_DIR = (Path(__file__).parent / "live_api").resolve()
_LIVE_HOSTS = {"eutils.ncbi.nlm.nih.gov", "ftp.ebi.ac.uk", "www.ebi.ac.uk", "www.ncbi.nlm.nih.gov"}


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "live_api(max_requests=12): opt-in public provider contract",
    )


def pytest_collection_modifyitems(items):
    enabled = os.environ.get(_LIVE_ENV) == "1"
    for item in items:
        if item.get_closest_marker("live_api") is None:
            continue
        if not Path(str(item.path)).resolve().is_relative_to(_LIVE_DIR):
            raise pytest.UsageError("live_api tests are restricted to tests/live_api")
        if not enabled:
            item.add_marker(pytest.mark.skip(reason=f"set {_LIVE_ENV}=1 to run live API contracts"))


def _blocked(kind):
    def fail(*_args, **_kwargs):
        pytest.fail(f"Unreviewed {kind} access is blocked by tests/conftest.py")

    return fail


@pytest.fixture(autouse=True)
def block_external_side_effects(monkeypatch, request):
    live = request.node.get_closest_marker("live_api") is not None and os.environ.get(_LIVE_ENV) == "1"
    if live:
        live_marker = request.node.get_closest_marker("live_api")
        max_requests = int(live_marker.kwargs.get("max_requests", 12))
        original_request = requests.sessions.Session.request
        original_send = requests.sessions.Session.send
        sends = [0]

        def validate(url):
            parsed = urlsplit(str(url))
            if parsed.scheme not in {"http", "https"} or parsed.hostname not in _LIVE_HOSTS:
                pytest.fail(f"unapproved live API URL: {url}")

        def guarded_request(session, method, url, *args, **kwargs):
            validate(url)
            return original_request(session, method, url, *args, **kwargs)

        def guarded_send(session, prepared, *args, **kwargs):
            validate(prepared.url)
            sends[0] += 1
            if sends[0] > max_requests:
                pytest.fail(f"live API request budget exceeded ({max_requests})")
            return original_send(session, prepared, *args, **kwargs)

        monkeypatch.setattr(requests.sessions.Session, "request", guarded_request)
        monkeypatch.setattr(requests.sessions.Session, "send", guarded_send)
    else:
        monkeypatch.setattr(socket.socket, "connect", _blocked("network"))
        monkeypatch.setattr(socket, "create_connection", _blocked("network"))
    if request.node.get_closest_marker("fake_process") is None:
        for name in ("Popen", "run", "call", "check_call", "check_output"):
            monkeypatch.setattr(subprocess, name, _blocked("child-process"))
