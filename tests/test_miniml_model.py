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

from meta_standards_converter.miniml import (
    MINIML_SCHEMA_VERSION,
    MINiMLCodec,
    MINiMLModelError,
    MINiMLPackage,
    Series,
)


def complete_package() -> dict:
    return {
        "miniml_schema_version": "2.0",
        "source": {
            "format": "GEO MINiML",
            "version": "0.5.4",
            "schema_location": "https://example.org/MINiML.xsd",
        },
        "database": [
            {
                "iid": "GEO",
                "name": "Gene Expression Omnibus",
                "public_id": "GEO",
                "organization_ref": {"ref": "org-1"},
            }
        ],
        "organization": [
            {
                "iid": "org-1",
                "name": "Example Institute",
                "address": {
                    "line": ["1 Science Road"],
                    "city": "London",
                    "postal_code": "NW1",
                    "country": "UK",
                },
            }
        ],
        "contributor": [
            {
                "iid": "person-1",
                "person": {"first": "Ada", "last": "Lovelace"},
                "email": "ada@example.org",
                "organization_ref": {"ref": "org-1"},
            }
        ],
        "platform": [
            {
                "iid": "GPL1",
                "accession": [{"database": "GEO", "value": "GPL1"}],
                "title": "Sequencing platform",
                "technology": "high-throughput sequencing",
                "distribution": "virtual",
                "organism": [{"taxid": "9606", "value": "Homo sapiens"}],
                "manufacturer": "Example",
                "contact_ref": [{"ref": "person-1", "position": "0"}],
            }
        ],
        "sample": [
            {
                "iid": "GSM1",
                "accession": [{"database": "GEO", "value": "GSM1"}],
                "title": "Sample one",
                "channel_count": "1",
                "channel": [
                    {
                        "source": "lung",
                        "organism": [
                            {"taxid": "9606", "value": "Homo sapiens"}
                        ],
                        "characteristics": [
                            {"name": "disease state", "value": "normal"}
                        ],
                        "molecule": "total RNA",
                    }
                ],
                "platform_ref": {"ref": "GPL1"},
                "raw_data": [
                    {
                        "type": "FASTQ",
                        "checksum": "0123456789abcdef0123456789abcdef",
                        "value": "https://example.org/read.fastq.gz",
                    }
                ],
            }
        ],
        "series": {
            "iid": "GSE1",
            "accession": [{"database": "GEO", "value": "GSE1"}],
            "title": "Example study",
            "summary": "A complete model fixture.",
            "sample_ref": [{"ref": "GSM1", "position": "0"}],
            "variable": [
                {
                    "position": "0",
                    "factor": "disease state",
                    "description": "Case versus control",
                    "sample_ref": [{"ref": "GSM1"}],
                }
            ],
            "relation": [
                {"type": "SubSeries of", "target": "GSE2", "comment": "test"}
            ],
            "data_table": [
                {
                    "title": "Series matrix",
                    "column": [
                        {"position": "0", "name": "ID_REF", "type": "identifier"}
                    ],
                    "external_data": {"rows": "1", "value": "matrix.txt"},
                }
            ],
            "extensions": {"vendor_note": {"value": "preserved"}},
        },
    }


def test_complete_xsd_derived_package_round_trips_through_python_model() -> None:
    model = MINiMLPackage.from_mapping(complete_package())

    assert isinstance(model.series, Series)
    assert model.series.sample_ref[0].ref == "GSM1"
    assert model.samples[0].channels[0].characteristics[0].name == "disease state"
    assert model.series.extras["vendor_note"] == {"value": "preserved"}

    canonical = model.to_mapping()
    assert canonical["miniml_schema_version"] == MINIML_SCHEMA_VERSION
    assert canonical["series"]["extensions"]["vendor_note"] == {"value": "preserved"}
    assert MINiMLCodec().decode(canonical, strict=True).package == model


def test_v2_singletons_are_normalized_without_dropping_extensions() -> None:
    model = MINiMLPackage.from_mapping(
        {
            "miniml_schema_version": "2.0",
            "source": {"format": "test"},
            "contributor": {
                "iid": "contributor-1",
                "address": "University Hospital, London, UK",
            },
            "sample": {"iid": "GSM1", "extensions": {"custom_sample_field": "kept"}},
            "series": {
                "iid": "GSE1",
                "accession": {"value": "GSE1"},
                "sample_ref": {"ref": "GSM1"},
            },
            "extensions": {"custom_root": {"source": "legacy"}},
        }
    )

    canonical = model.to_mapping()
    assert canonical["database"] == []
    assert canonical["contributor"][0]["address"] == "University Hospital, London, UK"
    assert canonical["sample"] == [
        {"iid": "GSM1", "extensions": {"custom_sample_field": "kept"}}
    ]
    assert canonical["series"]["accession"] == [{"value": "GSE1"}]
    assert canonical["series"]["sample_ref"] == [{"ref": "GSM1"}]
    assert canonical["extensions"]["custom_root"] == {"source": "legacy"}


def test_compatibility_validation_reports_xsd_deviations_as_warnings() -> None:
    payload = complete_package()
    payload["platform"][0]["technology"] = "single-cell spatial sequencing"
    payload["sample"][0]["raw_data"][0]["checksum"] = "not-an-md5"
    payload["series"]["sample_ref"].append({"ref": "GSM404"})
    payload["sample"][0]["channel_count"] = "2"

    issues = MINiMLPackage.from_mapping(payload).validate()
    codes = {issue.code for issue in issues}

    assert {
        "xsd_enumeration",
        "xsd_checksum",
        "unresolved_reference",
        "channel_count_mismatch",
    } <= codes
    assert all(issue.severity == "warning" for issue in issues)


@pytest.mark.parametrize(
    ("payload", "message"),
    [
        ({"miniml_schema_version": "2.0", "series": []}, "series must be an object"),
        ({"miniml_schema_version": "2.0", "series": {}}, "series requires iid or accession"),
        (
            {"miniml_schema_version": "2.0", "source": {"format": "test"}, "series": {"iid": "GSE1"}, "sample": ["bad"]},
            r"sample\[0\] must be an object",
        ),
        (
            {
                "miniml_schema_version": "2.0",
                "source": {"format": "test"},
                "series": {"iid": "GSE1"},
                "sample": [{"iid": "GSM1"}, {"iid": "GSM1"}],
            },
            "duplicate sample iid",
        ),
    ],
)
def test_structural_corruption_fails(payload: dict, message: str) -> None:
    with pytest.raises(MINiMLModelError, match=message):
        MINiMLPackage.from_mapping(payload)


def test_dump_and_load_are_deterministic_and_atomic(tmp_path: Path) -> None:
    destination = tmp_path / "package.json"
    model = MINiMLPackage.from_mapping(complete_package())

    model.dump(destination)

    assert MINiMLPackage.load(destination) == model
    assert destination.read_text(encoding="utf-8").endswith("\n")
