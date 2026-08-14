# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from __future__ import annotations

import pytest

from meta_standards_converter.miniml import (
    MINIML_SCHEMA_VERSION,
    MINiMLCodec,
    MINiMLModelError,
    HarmonizedValue,
    append_harmonized_value,
    iter_harmonized_values,
    parse_harmonized_mapping,
)


def _package(version: str = "3.0") -> dict:
    return {
        "miniml_schema_version": version,
        "source": {"format": "test"},
        "sample": [
            {
                "iid": "GSM1",
                "channel": [
                    {
                        "source": {
                            "value": "lung",
                            "hz_tissue": "lung",
                            "hz_tissue_id": "UBERON:0002048",
                            "hz_tissue_onto": "uberon",
                            "hz_tissue_hierarchy_depth": 0,
                        },
                        "characteristics": [
                            {"name": "disease", "value": "IPF"},
                            {"name": "hz_disease", "value": "idiopathic pulmonary fibrosis"},
                            {"name": "hz_disease_id", "value": "MONDO:0002771"},
                            {"name": "hz_disease_onto", "value": "mondo"},
                            {"name": "hz_disease_hierarchy_depth", "value": 0},
                            {"name": "hz_disease(1)", "value": "pulmonary fibrosis"},
                            {"name": "hz_disease_id(1)", "value": "MONDO:0003782"},
                            {"name": "hz_disease_onto(1)", "value": "mondo"},
                            {"name": "hz_disease_hierarchy_depth(1)", "value": 1},
                        ],
                    }
                ],
            }
        ],
        "series": {"iid": "GSE1"},
    }


def test_miniml_v3_round_trips_flat_harmonized_values() -> None:
    result = MINiMLCodec().decode(_package(), strict=True)
    encoded = MINiMLCodec.encode(result.package)

    assert MINIML_SCHEMA_VERSION == "3.0"
    assert encoded["sample"] == _package()["sample"]
    assert encoded["series"] == _package()["series"]
    assert "annotations" not in repr(encoded)

    channel = encoded["sample"][0]["channel"][0]
    source_values = list(iter_harmonized_values(channel["source"]))
    characteristic_values = list(
        iter_harmonized_values(channel["characteristics"])
    )
    assert source_values == [
        HarmonizedValue(
            field="tissue",
            value="lung",
            term_accession_number="UBERON:0002048",
            term_source_ref="uberon",
            hierarchy_depth=0,
        )
    ]
    assert [value.index for value in characteristic_values] == [0, 1]
    assert [value.value for value in characteristic_values] == [
        "idiopathic pulmonary fibrosis",
        "pulmonary fibrosis",
    ]


def test_v2_annotations_migrate_to_occurrence_local_hz_values() -> None:
    package = _package("2.0")
    channel = package["sample"][0]["channel"][0]
    channel["source"] = {
        "value": "lung",
        "annotations": [
            {
                "field": "tissue",
                "value": "lung",
                "term_accession_number": "UBERON:0002048",
                "term_source_ref": "uberon",
                "hierarchy_depth": 0,
            }
        ],
    }
    channel["characteristics"] = [
        {
            "name": "disease",
            "value": "IPF",
            "annotations": [
                {
                    "field": "disease",
                    "value": "idiopathic pulmonary fibrosis",
                    "term_accession_number": "MONDO:0002771",
                    "term_source_ref": "mondo",
                    "hierarchy_depth": 0,
                },
                {
                    "field": "disease",
                    "value": "pulmonary fibrosis",
                    "term_accession_number": "MONDO:0003782",
                    "term_source_ref": "mondo",
                    "hierarchy_depth": 1,
                },
            ],
        }
    ]

    decoded = MINiMLCodec().decode(package, strict=True)
    encoded = MINiMLCodec.encode(decoded.package)

    assert encoded["miniml_schema_version"] == "3.0"
    assert encoded["sample"][0]["channel"][0]["source"]["hz_tissue_id"] == (
        "UBERON:0002048"
    )
    rows = encoded["sample"][0]["channel"][0]["characteristics"]
    assert {row["name"] for row in rows} == {
        "disease",
        "hz_disease",
        "hz_disease_id",
        "hz_disease_onto",
        "hz_disease_hierarchy_depth",
        "hz_disease(1)",
        "hz_disease_id(1)",
        "hz_disease_onto(1)",
        "hz_disease_hierarchy_depth(1)",
    }


def test_harmonized_companions_require_a_value_and_aligned_suffix() -> None:
    with pytest.raises(MINiMLModelError, match="without a corresponding value"):
        parse_harmonized_mapping({"hz_disease_id": "MONDO:1"})

    with pytest.raises(MINiMLModelError, match="invalid harmonized field"):
        parse_harmonized_mapping({"hz_disease_1": "fibrosis"})

    package = _package()
    package["sample"][0]["channel"][0]["characteristics"] = [
        {"name": "disease", "value": "IPF"},
        {"name": "hz_disease_id", "value": "MONDO:1"},
    ]
    with pytest.raises(MINiMLModelError, match="without a corresponding value"):
        MINiMLCodec().decode(package)


def test_v3_rejects_annotation_objects() -> None:
    package = _package()
    package["sample"][0]["channel"][0]["source"]["annotations"] = [
        {"field": "tissue", "value": "lung"}
    ]

    with pytest.raises(MINiMLModelError, match="annotations"):
        MINiMLCodec().decode(package)


def test_append_harmonized_value_allocates_aligned_mapping_collision() -> None:
    destination = {
        "value": "IPF",
        "hz_disease": "idiopathic pulmonary fibrosis",
        "hz_disease_id": "MONDO:0002771",
        "hz_disease_onto": "MONDO",
        "hz_disease_hierarchy_depth": 0,
    }

    appended = append_harmonized_value(
        destination,
        HarmonizedValue(
            field="disease",
            value="pulmonary fibrosis",
            term_accession_number="MONDO:0003782",
            term_source_ref="MONDO",
            hierarchy_depth=1,
        ),
    )

    assert appended.index == 1
    assert destination["hz_disease(1)"] == "pulmonary fibrosis"
    assert destination["hz_disease_id(1)"] == "MONDO:0003782"
    assert destination["hz_disease_onto(1)"] == "MONDO"
    assert destination["hz_disease_hierarchy_depth(1)"] == 1


def test_append_harmonized_value_writes_named_rows_and_reuses_exact_identity() -> None:
    destination = [
        {"tag": "disease", "value": "IPF"},
        {"tag": "hz_disease", "value": "idiopathic pulmonary fibrosis"},
        {"tag": "hz_disease_id", "value": "MONDO:0002771"},
        {"tag": "hz_disease_onto", "value": "MONDO"},
        {"tag": "hz_disease_hierarchy_depth", "value": 0},
    ]
    candidate = HarmonizedValue(
        field="disease",
        value="pulmonary fibrosis",
        term_accession_number="MONDO:0003782",
        term_source_ref="MONDO",
        hierarchy_depth=1,
    )

    first = append_harmonized_value(destination, candidate, name_key="tag")
    second = append_harmonized_value(destination, candidate, name_key="tag")

    assert first == second
    assert first.index == 1
    assert (
        destination.count(
            {"tag": "hz_disease(1)", "value": "pulmonary fibrosis"}
        )
        == 1
    )
    assert {row["tag"] for row in destination if row["tag"].startswith("hz_disease")} == {
        "hz_disease",
        "hz_disease_id",
        "hz_disease_onto",
        "hz_disease_hierarchy_depth",
        "hz_disease(1)",
        "hz_disease_id(1)",
        "hz_disease_onto(1)",
        "hz_disease_hierarchy_depth(1)",
    }


def test_append_harmonized_value_rejects_malformed_existing_groups() -> None:
    with pytest.raises(MINiMLModelError, match="without a corresponding value"):
        append_harmonized_value(
            {"value": "IPF", "hz_disease_id": "MONDO:1"},
            HarmonizedValue(field="disease", value="fibrosis"),
        )

    with pytest.raises(MINiMLModelError, match="name_key"):
        append_harmonized_value(
            [{"name": "disease", "value": "IPF"}],
            HarmonizedValue(field="disease", value="fibrosis"),
        )
