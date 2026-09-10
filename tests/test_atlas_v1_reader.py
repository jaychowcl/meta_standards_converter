# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
import ast
import json
from pathlib import Path
import tomllib

import pytest

from meta_standards_converter.atlas_v1 import AtlasV1Error, AtlasV1Reader


ROOT = Path(__file__).parents[1]
FIXTURE = ROOT / "tests" / "fixtures" / "contracts" / "atlas-document-v1.json"


def test_owner_contract_fixture_yields_only_harmonized_dataset_metadata():
    result = AtlasV1Reader().load(FIXTURE)

    assert [dataset.dataset_id for dataset in result.datasets] == ["GSE100"]
    assert result.datasets[0].source_repository == "geo"
    assert result.datasets[0].source_ordinal == 0
    assert result.datasets[0].metadata["sample"][0]["iid"] == "GSM100"
    assert result.warnings == (
        "GSE200: dataset status failed; no convertible metadata; "
        "collection_failed: fixture collection failure",
    )


def test_reader_accepts_and_validates_harmonization_status_contract_v2():
    payload = json.loads(FIXTURE.read_text(encoding="utf-8"))
    payload["datasets"][0]["harmonization"].pop("degraded_stages")
    payload["datasets"][0]["harmonization"]["status"] = {
        "contract_version": "2.0",
        "execution": "succeeded",
        "completeness": "complete",
        "evidence_confidence": "high",
        "validation": "valid",
        "publication": "publishable",
        "terminal_reason": "harmonization_complete",
        "errors": [],
    }

    result = AtlasV1Reader().from_mapping(payload)

    assert [dataset.dataset_id for dataset in result.datasets] == ["GSE100"]


def test_reader_rejects_invalid_harmonization_status_contract_v2():
    payload = json.loads(FIXTURE.read_text(encoding="utf-8"))
    payload["datasets"][0]["harmonization"]["status"] = {
        "contract_version": "2.0",
        "execution": "succeeded",
        "completeness": "partial",
        "evidence_confidence": "high",
        "validation": "valid",
        "publication": "publishable",
        "terminal_reason": "contradictory",
        "errors": [],
    }

    with pytest.raises(AtlasV1Error, match="harmonization.status"):
        AtlasV1Reader().from_mapping(payload)


def test_reader_rejects_legacy_unversioned_envelopes_with_cutover_guidance():
    with pytest.raises(
        AtlasV1Error,
        match=(
            r"^Legacy unversioned accessions envelopes are not supported; "
            r"use the Atlas v1 contract\.$"
        ),
    ):
        AtlasV1Reader().from_mapping({"accessions": []})


@pytest.mark.parametrize(
    ("mutate", "message"),
    [
        (lambda payload: payload.update(schema_version="3.0"), "unsupported Atlas"),
        (
            lambda payload: payload["datasets"].append(
                dict(payload["datasets"][0])
            ),
            "duplicate dataset_id",
        ),
        (
            lambda payload: payload["summary"].update(dataset_count=99),
            "summary dataset_count",
        ),
        (
            lambda payload: payload["datasets"][0].update(metadata=[]),
            "metadata must be an object",
        ),
        (
            lambda payload: payload["run"].update(status="finished"),
            "run.status is not a supported value",
        ),
        (
            lambda payload: payload["summary"].update(extra=1),
            "unknown summary fields",
        ),
        (
            lambda payload: payload["publications"][0].update(extra=1),
            r"unknown publications\[0\] fields",
        ),
        (
            lambda payload: payload["datasets"][1]["diagnostics"][0].update(
                severity="critical"
            ),
            "severity is not a supported value",
        ),
    ],
)
def test_reader_fails_closed_on_invalid_v1_contracts(mutate, message):
    payload = json.loads(FIXTURE.read_text(encoding="utf-8"))
    mutate(payload)

    with pytest.raises(AtlasV1Error, match=message):
        AtlasV1Reader().from_mapping(payload)


@pytest.mark.parametrize(
    ("mutate", "message"),
    [
        (lambda payload: payload["summary"].update(dataset_count=True), "dataset_count"),
        (lambda payload: payload["run"]["config"].update(queries=[1]), "queries\\[0\\]"),
        (
            lambda payload: payload["datasets"][0].update(publication_ids=[{}]),
            "publication_ids\\[0\\]",
        ),
        (
            lambda payload: payload["datasets"][0].update(source_ordinal=True),
            "source_ordinal",
        ),
    ],
)
def test_reader_rejects_python_scalar_aliases_and_malformed_string_lists(
    mutate, message
):
    payload = json.loads(FIXTURE.read_text(encoding="utf-8"))
    mutate(payload)

    with pytest.raises(AtlasV1Error, match=message):
        AtlasV1Reader().from_mapping(payload)


def test_runtime_and_build_metadata_do_not_depend_on_thematicatlases():
    pyproject = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    assert pyproject["project"]["version"] == "7.0.0"
    dependencies = pyproject["project"]["dependencies"]
    assert not any("thematicatlases" in item.lower() for item in dependencies)

    offenders = []
    for path in (ROOT / "src" / "meta_standards_converter").rglob("*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                names = [alias.name for alias in node.names]
            elif isinstance(node, ast.ImportFrom):
                names = [node.module or ""]
            else:
                continue
            if any(name.lower().startswith("thematicatlases") for name in names):
                offenders.append(path.relative_to(ROOT))
    assert offenders == []
