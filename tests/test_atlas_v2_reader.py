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

from meta_standards_converter.atlas_v2 import AtlasV2Error, AtlasV2Reader


ROOT = Path(__file__).parents[1]
FIXTURE = ROOT / "tests" / "fixtures" / "contracts" / "atlas-document-v2.json"


def test_owner_contract_fixture_yields_only_harmonized_dataset_metadata():
    result = AtlasV2Reader().load(FIXTURE)

    assert [dataset.dataset_id for dataset in result.datasets] == ["GSE100"]
    assert result.datasets[0].source_repository == "geo"
    assert result.datasets[0].source_ordinal == 0
    assert result.datasets[0].metadata["sample"][0]["iid"] == "GSM100"
    assert result.warnings == (
        "GSE200: dataset status failed; no convertible metadata; "
        "collection_failed: fixture collection failure",
    )


def test_reader_rejects_legacy_v1_envelopes_with_cutover_guidance():
    with pytest.raises(
        AtlasV2Error,
        match=(
            r"^Legacy Atlas v1 envelopes are not supported; "
            r"use the pinned v1 tools\.$"
        ),
    ):
        AtlasV2Reader().from_mapping({"accessions": []})


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
def test_reader_fails_closed_on_invalid_v2_contracts(mutate, message):
    payload = json.loads(FIXTURE.read_text(encoding="utf-8"))
    mutate(payload)

    with pytest.raises(AtlasV2Error, match=message):
        AtlasV2Reader().from_mapping(payload)


def test_runtime_and_build_metadata_do_not_depend_on_thematicatlases():
    pyproject = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    assert pyproject["project"]["version"] == "4.0.0"
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
