"""Public-boundary regressions from the source-fidelity audit."""
from copy import deepcopy

import pytest

from meta_standards_converter.converters import Converter
from meta_standards_converter.converters.unified.families import recorded_neighbors
from meta_standards_converter.miniml import MINiMLCodec
from tests.converters.test_json2tsv import package


def ena_package():
    value = package("ERP000263", "ERS000001")
    value["source"]["format"] = "ENA"
    value["database"] = [{"iid": "ENA", "name": "ENA"}]
    value["series"]["accession"].append({"value": "PRJEB1"})
    value["extensions"] = {"insdc": {"version": "2.0", "records": [
        {"provider": "ena", "kind": "cross_references", "accession": "ERP000263",
         "metadata": [{"source": "ArrayExpress", "primary_id": "E-MTAB-308"}]},
        {"provider": "ena", "kind": "PROJECT", "accession": "PRJEB1", "metadata": {
            "tag": "PROJECT", "attributes": {"accession": "PRJEB1"}, "children": [
                {"tag": "CHILD_PROJECT", "attributes": {"accession": "PRJEB2"}}]}},
    ]}}
    return value


@pytest.mark.parametrize("target", ["json", "magetab", "tsv", "csv"])
def test_ena_list_records_are_retained_and_do_not_break_off(target):
    value = ena_package()
    original = deepcopy(value)
    result = Converter().convert(value, out_type=target, enrichment="off",
                                 options={"expand_studies": False})
    assert result.status == "complete", result.to_dict()
    assert value == original
    if target == "json":
        assert result.items[0].payload[0].to_mapping()["extensions"] == value["extensions"]


def test_family_traversal_uses_project_records_only():
    typed = MINiMLCodec().decode(ena_package(), strict=True).package
    assert recorded_neighbors(typed, "ena") == ["PRJEB2"]


def test_disabled_preparation_does_not_inspect_related_records(monkeypatch):
    def unexpected(*args):
        pytest.fail("disabled linked enrichment inspected related records")
    monkeypatch.setattr("meta_standards_converter.converters.unified.families.recorded_neighbors", unexpected)
    result = Converter().convert(ena_package(), out_type="json", enrichment="off",
                                 options={"expand_studies": False})
    assert result.status == "complete", result.to_dict()


@pytest.mark.parametrize("target", ["json", "magetab", "tsv", "csv", "h5ad", "obs"])
def test_invalid_references_fail_before_publication(target, tmp_path):
    value = package()
    value["series"]["sample_ref"] = [{"ref": "missing"}]
    result = Converter().convert(value, out_type=target, outdir=tmp_path,
                                 enrichment="off", options={"expand_studies": False})
    item = result.items[0]
    assert item.status == "failed", result.to_dict()
    assert item.validation == "invalid"
    assert any(d.code == "invalid_metadata" and "reference" in d.message for d in item.diagnostics)
    assert not item.artifacts and not list(tmp_path.iterdir())


def test_invalid_enrichment_is_rejected_with_diagnostic_and_valid_partial(tmp_path):
    class Enricher:
        def enrich_selected(self, data, **kwargs):
            value = data.to_mapping()
            value["series"]["sample_ref"] = [{"ref": "missing"}]
            value["sample"][0]["title"] = "should not survive"
            return value
    value = package()
    result = Converter(services={"enricher": Enricher()}).convert(
        value, out_type="json", outdir=tmp_path, enrichment="curators",
        options={"expand_studies": False})
    item = result.items[0]
    assert (item.status, item.execution, item.completeness, item.validation) == (
        "partial", "succeeded", "partial", "valid"), result.to_dict()
    MINiMLCodec().decode(item.payload[0], strict=True)
    assert item.payload[0].samples[0].title == "sample title"
    assert any(d.code == "source_partial" and "/series/sample_ref" in d.message
               and "standard" in d.message for d in item.diagnostics)
    assert item.artifacts


def test_strict_compatible_duplicate_titles_still_succeed():
    value = package()
    second = deepcopy(value["sample"][0])
    second["iid"] = "GSM2"
    value["sample"].append(second)
    result = Converter().convert(value, out_type="json", enrichment="off",
                                 options={"expand_studies": False})
    assert result.status == "complete", result.to_dict()
