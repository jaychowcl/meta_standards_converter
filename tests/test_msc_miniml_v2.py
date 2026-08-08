# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
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
    MINiMLV1Migrator,
)


def package_v2() -> dict:
    return {
        "miniml_schema_version": "2.0",
        "source": {
            "format": "MAGE-TAB",
            "version": "1.1",
            "documents": [
                {"kind": "idf", "name": "study.idf.txt"},
                {"kind": "sdrf", "name": "study.sdrf.txt"},
            ],
        },
        "database": [],
        "organization": [],
        "contributor": [
            {
                "iid": "person-1",
                "person": {"first": "Ada", "last": "Lovelace"},
                "roles": [
                    {
                        "value": "investigator",
                        "term_source_ref": "EFO",
                        "term_accession_number": "EFO:0001739",
                    }
                ],
            }
        ],
        "platform": [],
        "sample": [
            {
                "iid": "sample-1",
                "channel": [
                    {
                        "source": {"value": "lung"},
                        "characteristics": [
                            {
                                "name": "age",
                                "value": "10",
                                "unit": {
                                    "value": "year",
                                    "term_source_ref": "UO",
                                    "term_accession_number": "UO:0000036",
                                    "annotations": [
                                        {
                                            "field": "unit",
                                            "value": "year",
                                            "term_source_ref": "UO",
                                            "term_accession_number": "UO:0000036",
                                            "hierarchy_depth": 0,
                                        }
                                    ],
                                },
                            }
                        ],
                    }
                ],
            }
        ],
        "series": {
            "iid": "E-MTAB-1",
            "sample_ref": [{"ref": "sample-1"}],
            "experiment_date": "2025-01-01",
            "protocols": [
                {
                    "name": "P-extract",
                    "type": {"value": "nucleic acid extraction"},
                    "parameters": ["duration"],
                },
                {"name": "P-library", "type": {"value": "library preparation"}},
            ],
            "assay_paths": [
                {
                    "document": "study.sdrf.txt",
                    "steps": [
                        {"kind": "sample", "name": "sample-1", "sample_ref": "sample-1"},
                        {
                            "kind": "protocol_application",
                            "protocol_ref": "P-extract",
                            "parameter_values": [
                                {
                                    "name": "duration",
                                    "value": "30",
                                    "unit": {"value": "minute"},
                                }
                            ],
                        },
                        {"kind": "protocol_application", "protocol_ref": "P-library"},
                        {"kind": "assay", "name": "assay-1"},
                    ],
                }
            ],
            "quality_controls": [{"value": "biological replicate"}],
            "replicate_types": [{"value": "technical replicate"}],
            "normalization_types": [{"value": "quantile normalization"}],
        },
    }


def test_v2_is_typed_source_neutral_and_preserves_semantic_order() -> None:
    result = MINiMLCodec().decode(package_v2())
    package = result.package

    assert MINIML_SCHEMA_VERSION == "2.0"
    assert package.source.format == "MAGE-TAB"
    assert [protocol.name for protocol in package.series.protocols] == [
        "P-extract",
        "P-library",
    ]
    assert [step.kind for step in package.series.assay_paths[0].steps] == [
        "sample",
        "protocol_application",
        "protocol_application",
        "assay",
    ]
    characteristic = package.samples[0].channels[0].characteristics[0]
    assert characteristic.name == "age"
    assert characteristic.unit.term_accession_number == "UO:0000036"
    assert characteristic.unit.annotations[0].field == "unit"
    assert MINiMLCodec().encode(package) == package_v2()


@pytest.mark.parametrize("version", [None, "1.0", "3.0"])
def test_runtime_codec_rejects_every_non_v2_document(version: str | None) -> None:
    payload = package_v2()
    if version is None:
        payload.pop("miniml_schema_version")
    else:
        payload["miniml_schema_version"] = version

    with pytest.raises(MINiMLModelError, match="requires MSC MINiML schema version '2.0'"):
        MINiMLCodec().decode(payload)


def test_v2_rejects_opaque_mage_tab_and_unknown_sibling_fields() -> None:
    payload = package_v2()
    payload["mage_tab"] = {"model": {}, "roundtrip": {}}
    with pytest.raises(MINiMLModelError, match="unsupported MINiML package field: mage_tab"):
        MINiMLCodec().decode(payload)

    payload = package_v2()
    payload["series"]["idf_layout"] = []
    with pytest.raises(MINiMLModelError, match="unsupported series field: idf_layout"):
        MINiMLCodec().decode(payload)


def test_v2_validates_protocol_references_and_natural_identity() -> None:
    payload = package_v2()
    payload["series"]["protocols"].append({"name": "P-extract"})
    with pytest.raises(MINiMLModelError, match="duplicate protocol name: P-extract"):
        MINiMLCodec().decode(payload)

    payload = package_v2()
    payload["series"]["assay_paths"][0]["steps"][1]["protocol_ref"] = "P-missing"
    with pytest.raises(MINiMLModelError, match="unknown protocol reference: P-missing"):
        MINiMLCodec().decode(payload)


def test_explicit_v1_migrator_folds_mage_tab_and_hz_fields() -> None:
    legacy = {
        "miniml_schema_version": "1.0",
        "version": "magetabv1.1",
        "schema_location": "https://www.ebi.ac.uk/arrayexpress/help/magetab_spec.html",
        "database": [],
        "organization": [],
        "contributor": [],
        "platform": [],
        "sample": [
            {
                "iid": "sample-1",
                "channel": [
                    {
                        "source": "lung",
                        "characteristics": [
                            {"tag": "organism", "value": "Homo sapiens"}
                        ],
                    }
                ],
            }
        ],
        "series": {"iid": "E-MTAB-1", "sample_ref": [{"ref": "sample-1"}]},
        "mage_tab": {
            "model": {
                "schema_version": 1,
                "idf_layout": [],
                "protocols": [
                    {
                        "id": "protocol:1",
                        "position": 0,
                        "name": "P-extract",
                        "type": "extraction",
                    }
                ],
                "declarations": {"quality_control": [], "replicate": [], "normalization": []},
                "assay_paths": [
                    {
                        "id": "study.sdrf.txt:1",
                        "sdrf": "study.sdrf.txt",
                        "row_index": 1,
                        "binding": {"sample_name": "sample-1"},
                        "steps": [
                            {"kind": "node", "name": "Sample Name", "value": "sample-1", "column_index": 0},
                            {"kind": "protocol", "value": "P-extract", "column_index": 1},
                            {
                                "kind": "attribute",
                                "attribute_type": "parameter value",
                                "name": "duration",
                                "value": "30",
                                "unit": "minute",
                                "hz_unit": "minute",
                                "hz_unit_id": "UO:0000031",
                                "hz_unit_onto": "UO",
                                "column_index": 2,
                            },
                            {"kind": "node", "name": "Assay Name", "value": "assay-1", "column_index": 4},
                        ],
                    }
                ],
                "sdrfs": [{"name": "study.sdrf.txt", "columns": []}],
                "investigation_fields": [],
            },
            "roundtrip": {"schema_version": 1, "idf_rows": [], "sdrfs": []},
        },
    }

    result = MINiMLV1Migrator().migrate(legacy)
    rendered = MINiMLCodec().encode(result.package)

    assert rendered["miniml_schema_version"] == "2.0"
    assert "mage_tab" not in rendered
    assert "roundtrip" not in str(rendered)
    assert rendered["series"]["protocols"][0]["name"] == "P-extract"
    application = rendered["series"]["assay_paths"][0]["steps"][1]
    assert application["protocol_ref"] == "P-extract"
    assert application["parameter_values"][0]["unit"]["annotations"][0][
        "term_accession_number"
    ] == "UO:0000031"
    assert any(item.code == "source_layout_dropped" for item in result.diagnostics)
