# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
import json

import pytest

from meta_standards_converter.converters.json_source import JSONPackageSource


def package(study: str, sample: str) -> dict:
    return {
        "series": {"accession": [{"value": study}]},
        "sample": [{"iid": sample}],
    }


def dataset(
    dataset_id: str,
    status: str,
    metadata: dict,
    *,
    diagnostics: list[dict] | None = None,
) -> dict:
    return {
        "dataset_id": dataset_id,
        "source_repository": "geo",
        "source_ordinal": 0,
        "status": status,
        "metadata": metadata,
        "publication_ids": [],
        "review": None,
        "harmonization": None,
        "diagnostics": diagnostics or [],
    }


def atlas_v2(datasets: list[dict]) -> dict:
    completed = sum(item["status"] == "harmonized" for item in datasets)
    failed = sum(item["status"] == "failed" for item in datasets)
    return {
        "schema_version": "2.0",
        "atlas": {"atlas_id": "atlas-test", "title": "Test", "theme": "test"},
        "run": {
            "run_id": "run-test",
            "created_at": "2026-07-31T00:00:00Z",
            "config": {"queries": [], "metadata_repositories": [], "options": {}},
            "status": "complete",
        },
        "datasets": datasets,
        "publications": [],
        "summary": {
            "dataset_count": len(datasets),
            "publication_count": 0,
            "completed_dataset_count": completed,
            "failed_dataset_count": failed,
        },
    }


def test_miniml_packages_are_grouped_by_study(tmp_path):
    source = tmp_path / "miniml.json"
    source.write_text(
        json.dumps([package("GSE2", "GSM2"), package("GSE1", "GSM1")]),
        encoding="utf-8",
    )

    result = JSONPackageSource().load(source)

    assert [group.dataset_id for group in result.groups] == ["GSE2", "GSE1"]
    assert result.groups[0].packages[0]["sample"][0]["iid"] == "GSM2"
    assert result.warnings == ()


def test_atlas_v2_document_yields_only_harmonized_metadata(tmp_path):
    source = tmp_path / "atlas.json"
    source.write_text(
        json.dumps(
            atlas_v2(
                [
                    dataset("GSE1", "harmonized", package("GSE1", "GSM1")),
                    dataset(
                        "GSE2",
                        "failed",
                        package("GSE2", "GSM2"),
                        diagnostics=[
                            {
                                "code": "harmonization_failed",
                                "message": "lookup failed",
                                "severity": "error",
                                "path": None,
                            }
                        ],
                    ),
                    dataset("GSE3", "collected", {}),
                ]
            )
        ),
        encoding="utf-8",
    )

    result = JSONPackageSource().load(source)

    assert [group.dataset_id for group in result.groups] == ["GSE1"]
    assert result.groups[0].source_accession == "GSE1"
    assert result.warnings == (
        "GSE2: dataset status failed; no convertible metadata; "
        "harmonization_failed: lookup failed",
        "GSE3: dataset status collected; no convertible metadata",
    )


def test_atlas_v2_dataset_expands_multiple_miniml_packages(tmp_path):
    source = tmp_path / "atlas.json"
    source.write_text(
        json.dumps(
            atlas_v2(
                [
                    dataset(
                        "GSE1",
                        "harmonized",
                        {
                            "packages": [
                                package("GSE1", "GSM1"),
                                package("GSE1-related", "GSM2"),
                            ]
                        },
                    )
                ]
            )
        ),
        encoding="utf-8",
    )

    result = JSONPackageSource().load(source)

    assert len(result.groups) == 1
    assert result.groups[0].dataset_id == "GSE1"
    assert [
        item["sample"][0]["iid"] for item in result.groups[0].packages
    ] == ["GSM1", "GSM2"]


def test_conflicting_duplicate_samples_are_rejected(tmp_path):
    first = package("GSE1", "GSM1")
    second = package("GSE1", "GSM1")
    second["sample"][0]["title"] = "different"
    source = tmp_path / "miniml.json"
    source.write_text(json.dumps([first, second]), encoding="utf-8")

    try:
        JSONPackageSource().load(source)
    except ValueError as error:
        assert "GSM1" in str(error)
        assert "conflicting" in str(error)
    else:
        raise AssertionError("expected conflicting duplicate sample failure")


@pytest.mark.parametrize(
    "payload",
    [
        [],
        atlas_v2([]),
        atlas_v2([dataset("GSE1", "collected", {})]),
    ],
)
def test_sources_without_convertible_groups_fail_closed(tmp_path, payload):
    source = tmp_path / "empty.json"
    source.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(
        ValueError, match=r"^JSON source contains no convertible package groups\.$"
    ):
        JSONPackageSource().load(source)


@pytest.mark.parametrize(
    "payload",
    [
        package("GSE1", "GSM1") | {"sample": []},
        atlas_v2(
            [
                dataset(
                    "GSE1",
                    "harmonized",
                    package("GSE1", "GSM1") | {"sample": []},
                )
            ]
        ),
    ],
)
def test_sources_without_convertible_samples_fail_closed(tmp_path, payload):
    source = tmp_path / "samples.json"
    source.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(
        ValueError, match=r"^JSON source contains no convertible samples\.$"
    ):
        JSONPackageSource().load(source)
