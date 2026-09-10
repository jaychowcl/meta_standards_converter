# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
import pytest
import json
import requests
from meta_standards_converter.converters import GEO2JSONConverter
from tests.support.contracts import FIXTURES, install_replay


def documents():
    case = FIXTURES / "edge_cases/geo-relations"
    return {"GSE1": (case / "child.xml").read_bytes(), "GSE2": (case / "parent.xml").read_bytes()}

@pytest.mark.parametrize("related,expected_ids", [(False, ["GSE1"]), (True, ["GSE1", "GSE2"])])
def test_reciprocal_parent_evidence_and_related_cardinality(workspace, monkeypatch, related, expected_ids):
    calls = install_replay(monkeypatch, documents=documents())
    packages = GEO2JSONConverter().convert("GSE1", related_series=related, out="out")
    assert [p.series.iid for p in packages] == expected_ids
    assert packages[0]["series"]["pubmed_id"] == ["42129775"]
    assert packages[0]["extensions"]["publication_inheritance"]["source_series"] == "GSE2"
    assert calls[:2] == ["geo:GSE1", "geo:GSE2"]
    assert calls.count("geo:GSE1") == calls.count("geo:GSE2") == 1
    assert calls.count("pubmed:42129775") == len(expected_ids)
    paths = list((workspace / "out").glob("*.json"))
    assert [p.name for p in paths] == ["GSE1.json"]
    assert [p["series"]["iid"] for p in json.loads(paths[0].read_text())] == expected_ids


def test_missing_publication_evidence_does_not_invent_inheritance(workspace, monkeypatch):
    docs = documents(); docs["GSE2"] = docs["GSE2"].replace(b"<Pubmed-ID>42129775</Pubmed-ID>", b"")
    calls = install_replay(monkeypatch, documents=docs)
    packages = GEO2JSONConverter().convert("GSE1")
    assert not packages[0]["series"].get("pubmed_id")
    assert "publication_inheritance" not in packages[0].get("extensions", {})
    assert calls == ["geo:GSE1", "geo:GSE2"]


def test_publication_failure_retains_identifier_and_package(workspace, monkeypatch):
    calls = install_replay(monkeypatch, documents=documents(), fail_on=("pubmed:42129775",))
    packages = GEO2JSONConverter().convert("GSE1", out="out")
    assert packages[0]["series"]["pubmed_publication"][0]["pubmed_id"] == "42129775"
    assert packages[0]["series"]["pubmed_publication"][0].get("doi") is None
    assert (workspace / "out/GSE1.json").is_file()
    assert calls[-1] == "pubmed:42129775"


def test_related_retrieval_failure_is_not_a_silent_sample_skip(workspace, monkeypatch):
    install_replay(monkeypatch, documents=documents(), fail_on=("geo:GSE2",))
    with pytest.raises(requests.RequestException):
        GEO2JSONConverter().convert("GSE1", related_series=True, out="out")
    assert not (workspace / "out").exists()
