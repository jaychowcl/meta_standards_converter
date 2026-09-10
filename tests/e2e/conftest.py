# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
import pytest
from tests.support.contracts import install_replay

@pytest.fixture
def replay(monkeypatch):
    return install_replay(monkeypatch)

@pytest.fixture
def workspace(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    # Machine availability is an external input; keep the actual estimator/admission code.
    monkeypatch.setattr("meta_standards_converter.converters.json2h5ad._available_memory_bytes", lambda: 2**30)
    return tmp_path
