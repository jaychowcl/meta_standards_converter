# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Fail-closed test guard for accidental network and child-process access."""

import socket
import subprocess

import pytest


def _blocked(kind):
    def fail(*_args, **_kwargs):
        pytest.fail(f"Unreviewed {kind} access is blocked by tests/conftest.py")

    return fail


@pytest.fixture(autouse=True)
def block_external_side_effects(monkeypatch, request):
    monkeypatch.setattr(socket.socket, "connect", _blocked("network"))
    monkeypatch.setattr(socket, "create_connection", _blocked("network"))
    if request.node.get_closest_marker("fake_process") is None:
        for name in ("Popen", "run", "call", "check_call", "check_output"):
            monkeypatch.setattr(subprocess, name, _blocked("child-process"))
