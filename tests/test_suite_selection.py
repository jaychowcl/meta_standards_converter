# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from pathlib import Path

import pytest


@pytest.mark.parametrize(
    "arguments, passed, deselected",
    [
        ([], 1, 1),
        (["--full"], 2, 0),
        (["-m", "acceptance"], 0, 2),
        (["--full", "-m", "acceptance"], 1, 1),
    ],
    ids=["quick", "full", "acceptance-needs-full", "acceptance-only"],
)
def test_suite_selection(pytester, monkeypatch, arguments, passed, deselected):
    monkeypatch.delenv("RUN_LIVE_API_TESTS", raising=False)
    pytester.makeconftest(Path(__file__).with_name("conftest.py").read_text())
    pytester.makepyfile(
        test_examples="""
import pytest

def test_feature():
    pass

@pytest.mark.acceptance
def test_acceptance():
    pass
"""
    )
    live_dir = pytester.path / "live_api"
    live_dir.mkdir()
    (live_dir / "test_provider.py").write_text(
        "import pytest\n@pytest.mark.live_api\ndef test_provider():\n"
        "    pytest.fail('Provider must remain disabled even with --full')\n"
    )
    if "-m" in arguments:
        deselected += 1

    result = pytester.runpytest_inprocess("-q", "-p", "no:cacheprovider", *arguments)
    result.assert_outcomes(passed=passed, deselected=deselected, skipped=(0 if "-m" in arguments else 1))
    assert result.ret == (pytest.ExitCode.OK if passed or (0 if "-m" in arguments else 1) else pytest.ExitCode.NO_TESTS_COLLECTED)
