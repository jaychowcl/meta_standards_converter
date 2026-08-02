# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from __future__ import annotations

from pathlib import Path
import subprocess
import sys

import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "benchmarks" / "provider_http_concurrency.py"


def test_provider_benchmark_is_bounded_and_model_free() -> None:
    source = SCRIPT.read_text(encoding="utf-8")

    assert '"warmups": 1' in source
    assert '"workers": [1, 2]' in source
    assert "speedup >= 0.20" in source
    assert "* 1.5" in source
    assert '"model_calls": 0' in source
    assert '"embedding_calls": 0' in source
    assert "openai" not in source.lower()


@pytest.mark.fake_process
def test_provider_benchmark_help_does_not_call_network(tmp_path: Path) -> None:
    completed = subprocess.run(
        [sys.executable, str(SCRIPT), "--help"],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )

    assert completed.returncode == 0
    assert "--provider" in completed.stdout
    assert not list(tmp_path.iterdir())


def test_provider_scenarios_use_fixed_public_identifiers() -> None:
    source = SCRIPT.read_text(encoding="utf-8")

    for identifier in ("27708338", "SRR390728", "E-MTAB-5061"):
        assert identifier in source
