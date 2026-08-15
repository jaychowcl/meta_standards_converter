from __future__ import annotations

import csv
import json

import anndata
import pandas
from scipy import sparse

from meta_standards_converter.converters.json2ae import json2ae
from meta_standards_converter.converters.json2h5ad import JSON2H5ADConverter
from meta_standards_converter.converters.json2tabular import JSON2TSVConverter
from meta_standards_converter.miniml import (
    MINiMLHarmonizationPatch,
    apply_miniml_harmonization_patch,
    miniml_source_fingerprint,
)


def _applied_package() -> dict:
    package = {
        "miniml_schema_version": "3.0",
        "source": {"format": "test"},
        "database": [],
        "organization": [],
        "contributor": [],
        "platform": [],
        "series": {"iid": "GSE1"},
        "sample": [
            {
                "iid": "GSM1",
                "channel": [
                    {
                        "source": {"value": "blood"},
                        "characteristics": [
                            {"name": "disease", "value": "raw case"}
                        ],
                    }
                ],
            }
        ],
    }
    path = "/sample/0/channel/0/characteristics/0"
    patch = MINiMLHarmonizationPatch(
        base_sha256=miniml_source_fingerprint(package),
        adds=(
            {
                "path": path,
                "harmonized_value": {
                    "field": "sample_disease_name",
                    "value": "disease",
                    "term_source_ref": "DOID",
                    "term_accession_number": "DOID:4",
                    "hierarchy_depth": 0,
                },
                "source_evidence": {
                    "target_id": "target-1",
                    "source_field": "disease",
                    "source_label": "raw case",
                    "source_field_path": f"{path}/name",
                    "source_value_path": f"{path}/value",
                    "derivation": "direct",
                    "match_kind": "exact_value",
                },
            },
        ),
    )
    return apply_miniml_harmonization_patch(package, patch)


def test_tsv_exposes_deterministic_patch_provenance_without_overrides(tmp_path) -> None:
    source = tmp_path / "package.json"
    destination = tmp_path / "manifest.tsv"
    source.write_text(json.dumps(_applied_package()), encoding="utf-8")

    JSON2TSVConverter().convert_source(source, destination)
    with destination.open(encoding="utf-8", newline="") as stream:
        row = next(csv.DictReader(stream, delimiter="\t"))

    prefix = "msc.harmonization.sample_disease_name"
    assert row[f"{prefix}.value"] == "disease"
    assert row[f"{prefix}.id"] == "DOID:4"
    assert row[f"{prefix}.ontology"] == "DOID"
    assert row[f"{prefix}.hierarchy_depth"] == "0"
    assert row[f"{prefix}.source_field"] == "disease"
    assert row[f"{prefix}.source_label"] == "raw case"
    assert row[f"{prefix}.source_path"].endswith("/characteristics/0/value")
    assert row[f"{prefix}.match_kind"] == "exact_value"
    assert len(row[f"{prefix}.patch_id"]) == 64


def test_magetab_keeps_raw_characteristic_and_adds_machine_comments(tmp_path) -> None:
    source = tmp_path / "package.json"
    source.write_text(json.dumps(_applied_package()), encoding="utf-8")
    magetab = json2ae().convert(str(source), enrich=False)[0]
    sdrf = next(row[1] for row in magetab if row and row[0] == "SDRF File")
    header, values = sdrf[0], sdrf[1]

    disease = header.index("Characteristics[disease]")
    assert values[disease] == "raw case"
    label = "Comment[msc_harmonization_sample_disease_name_source_label]"
    assert label in header
    assert values[header.index(label)] == "raw case"
    assert "Characteristics[sample_disease_name]" not in header


def test_h5ad_and_self_contained_miniml_retain_applied_fragments(tmp_path) -> None:
    source = tmp_path / "package.json"
    source.write_text(json.dumps(_applied_package()), encoding="utf-8")
    expression = tmp_path / "GSM1.h5ad"
    anndata.AnnData(
        X=sparse.csr_matrix([[1]]),
        obs=pandas.DataFrame(index=["cell-1"]),
        var=pandas.DataFrame(index=["gene-1"]),
    ).write_h5ad(expression)

    result = JSON2H5ADConverter().convert(
        str(source),
        out=str(tmp_path / "out"),
        asset_specs=[f"GSM1={expression}"],
    )
    converted = anndata.read_h5ad(result.sample_h5ads["GSM1"])
    prefix = "msc.harmonization.sample_disease_name"
    assert converted.obs[f"{prefix}.source_label"].iat[0] == "raw case"
    retained = converted.uns["msc_harmonization"]
    assert retained["schema_version"] == "1.0"
    assert json.loads(retained["patches_json"])[0]["status"] == "applied"
    miniml = json.loads(converted.uns["msc_miniml"]["packages_json"])
    assert miniml[0]["extensions"]["msc_harmonization"]["patches"]

