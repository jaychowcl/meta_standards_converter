# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
import os
import socket
import subprocess

import pytest


def test_network_access_is_blocked_by_default():
    with pytest.raises(pytest.fail.Exception, match="network access is blocked"):
        socket.create_connection(("127.0.0.1", 9))


def test_child_process_access_is_blocked_by_default():
    with pytest.raises(pytest.fail.Exception, match="child-process access is blocked"):
        subprocess.run(["should-not-run"], check=False)


@pytest.mark.fake_process
def test_reviewed_fake_process_uses_a_bounded_path(tmp_path, monkeypatch):
    executable = tmp_path / "fake-ok"
    executable.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    executable.chmod(0o755)
    monkeypatch.setenv("PATH", os.fspath(tmp_path))

    completed = subprocess.run(["fake-ok"], check=False)

    assert completed.returncode == 0
