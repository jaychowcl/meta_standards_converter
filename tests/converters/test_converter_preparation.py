"""Presets govern supplementary retrieval independently of the destination."""
from copy import deepcopy
from types import SimpleNamespace
import json

import pytest

from meta_standards_converter.converters import Converter
from meta_standards_converter.converters.archive_results import ArchiveImportResult, StudyImportOutcome
from meta_standards_converter.miniml import MINiMLCodec
from tests.converters.test_json2tsv import package


def native(provider="SRA"):
    value = package("SRP1", "SRS1")
    value["source"]["format"] = provider
    value["series"]["iid"] = "SRP1"
    value["database"][0]["iid"] = "GEO"
    value["series"]["relation"] = [{"type": "GEO", "target": "GSE2"}]
    return value


def services(calls, *, failure=False):
    class Enricher:
        publication_issues = []
        def enrich_selected(self, data, *, publications, run_metadata):
            calls.append(("standard", publications, run_metadata))
            if failure:
                raise OSError("unavailable")
            return data
    class Linked:
        def enrich(self, data):
            calls.append(("linked",))
            return data, []
    class Peer:
        def convert(self, accession, **options):
            calls.append(("peer", accession))
            p = MINiMLCodec().decode(native("ENA")).package
            return ArchiveImportResult(accession, [StudyImportOutcome("SRP1", "SRP1", "complete", p)])
    return {"enricher": Enricher(), "linked_enricher": Linked(), "ena2json": Peer()}


@pytest.mark.parametrize("preset,expected", [
    ("standard", ["peer", "linked", "standard"]),
    ("curators", ["standard"]),
    ("off", []),
])
@pytest.mark.parametrize("target", ["json", "csv"])
def test_preset_controls_all_enrichment(preset, expected, target):
    calls = []
    value = native()
    original = deepcopy(value)
    result = Converter(services=services(calls)).convert(value, out_type=target,
        enrichment=preset, options={"expand_studies": False})
    assert result.status == "complete", result.to_dict()
    assert [call[0] for call in calls] == expected
    assert value == original
    assert any(d.code.startswith("enrichment_") for d in result.items[0].diagnostics)


def test_file_and_object_share_preparation(tmp_path):
    source = tmp_path / "native.json"
    source.write_text(json.dumps(native()))
    a, b = [], []
    one = Converter(services=services(a)).convert(source, out_type="json", enrichment="curators")
    two = Converter(services=services(b)).convert(native(), out_type="json", enrichment="curators")
    assert one.status == two.status == "complete"
    assert a == b == [("standard", True, False)]
    assert one.items[0].payload == two.items[0].payload


def test_failed_enrichment_retains_valid_partial_payload():
    result = Converter(services=services([], failure=True)).convert(native(),
        out_type="json", enrichment="curators")
    assert result.status == "partial", result.to_dict()
    assert result.items[0].execution == "succeeded"
    assert result.items[0].payload[0].samples[0].iid == "SRS1"
    assert any(d.code == "source_partial" for d in result.items[0].diagnostics)


def test_off_magetab_does_not_resolve_hidden_evidence():
    from meta_standards_converter.converters.json2ae import JSON2AEConverter
    from meta_standards_converter.magetab.constructor import AEConstructor
    from meta_standards_converter.metadata.enrichment import MAGETabEvidenceResolver
    calls = []
    class Publications:
        def pubmed_summary(self, **kwargs):
            calls.append("unexpected")
            raise AssertionError("off must not retrieve publications")
    value = package()
    value["series"]["pubmed_id"] = ["123"]
    constructor = AEConstructor(evidence_resolver=MAGETabEvidenceResolver(pubmed_client=Publications()))
    result = Converter(services={"json2ae": JSON2AEConverter(ae_constructor=constructor)}).convert(
        value, out_type="magetab", enrichment="off", options={"expand_studies": False})
    assert result.status == "complete", result.to_dict()
    assert not calls


def test_curators_geo_excludes_sra_run_enrichment():
    calls = []
    value = package()
    value["source"]["format"] = "GEO"
    result = Converter(services=services(calls)).convert(value, out_type="json", enrichment="curators")
    assert result.status == "complete", result.to_dict()
    assert calls == [("standard", True, False)]
