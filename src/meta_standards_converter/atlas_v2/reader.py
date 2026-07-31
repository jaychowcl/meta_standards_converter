# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Validate Atlas v2 JSON and adapt harmonized datasets for conversion."""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any, Mapping


SCHEMA_VERSION = "2.0"
_TOP_LEVEL_FIELDS = {
    "schema_version",
    "atlas",
    "run",
    "datasets",
    "publications",
    "summary",
}
_DATASET_FIELDS = {
    "dataset_id",
    "source_repository",
    "source_ordinal",
    "status",
    "metadata",
    "publication_ids",
    "review",
    "harmonization",
    "diagnostics",
}
_DATASET_STATUSES = {
    "discovered",
    "collected",
    "reviewed",
    "harmonized",
    "failed",
}
_RUN_STATUSES = {"planned", "running", "complete", "partial", "failed"}
_JUDGEMENTS = {"relevant", "not_relevant", "unsure"}
_DIAGNOSTIC_SEVERITIES = {"warning", "error"}


class AtlasV2Error(ValueError):
    """The input is not a valid or supported Atlas v2 document."""


@dataclass(frozen=True)
class AtlasV2Dataset:
    """One harmonized dataset adapted from an Atlas v2 document."""

    dataset_id: str
    source_repository: str
    source_ordinal: int
    metadata: Mapping[str, Any]


@dataclass(frozen=True)
class AtlasV2ReadResult:
    """Convertible datasets plus diagnostics for skipped dataset states."""

    datasets: tuple[AtlasV2Dataset, ...]
    warnings: tuple[str, ...] = ()


def _mapping(value: Any, name: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise AtlasV2Error(f"{name} must be an object")
    return value


def _list(value: Any, name: str) -> list[Any]:
    if not isinstance(value, list):
        raise AtlasV2Error(f"{name} must be an array")
    return value


def _only(value: Mapping[str, Any], allowed: set[str], name: str) -> None:
    unknown = sorted(set(value) - allowed)
    if unknown:
        raise AtlasV2Error(f"unknown {name} fields: {', '.join(unknown)}")


def _nonblank(value: Any, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise AtlasV2Error(f"{name} must be a nonblank string")
    return value


def _optional_string(value: Any, name: str) -> None:
    if value is not None and not isinstance(value, str):
        raise AtlasV2Error(f"{name} must be a string or null")


def _review(value: Any, name: str) -> None:
    if value is None:
        return
    review = _mapping(value, name)
    _only(review, {"judgement", "criteria", "evidence"}, name)
    if review.get("judgement") not in _JUDGEMENTS:
        raise AtlasV2Error(f"{name}.judgement is not a supported value")
    for index, raw in enumerate(_list(review.get("criteria"), f"{name}.criteria")):
        criterion = _mapping(raw, f"{name}.criteria[{index}]")
        _only(
            criterion,
            {"criterion", "passed", "confidence", "evidence"},
            f"{name}.criteria[{index}]",
        )
        _nonblank(criterion.get("criterion"), f"{name}.criteria[{index}].criterion")
        passed = criterion.get("passed")
        if passed is not None and not isinstance(passed, bool):
            raise AtlasV2Error(
                f"{name}.criteria[{index}].passed must be boolean or null"
            )
        _optional_string(
            criterion.get("confidence"), f"{name}.criteria[{index}].confidence"
        )
        _optional_string(
            criterion.get("evidence"), f"{name}.criteria[{index}].evidence"
        )
    for index, raw in enumerate(_list(review.get("evidence"), f"{name}.evidence")):
        _mapping(raw, f"{name}.evidence[{index}]")


def _harmonization(value: Any, name: str) -> None:
    if value is None:
        return
    harmonization = _mapping(value, name)
    _only(harmonization, {"targets", "degraded_stages"}, name)
    for index, raw in enumerate(
        _list(harmonization.get("targets"), f"{name}.targets")
    ):
        target = _mapping(raw, f"{name}.targets[{index}]")
        _only(
            target,
            {"raw_value", "path", "match"},
            f"{name}.targets[{index}]",
        )
        _nonblank(target.get("raw_value"), f"{name}.targets[{index}].raw_value")
        _nonblank(target.get("path"), f"{name}.targets[{index}].path")
        if target.get("match") is not None:
            match = _mapping(target["match"], f"{name}.targets[{index}].match")
            _only(
                match,
                {"label", "accession", "framework", "confidence"},
                f"{name}.targets[{index}].match",
            )
            for field in ("label", "accession", "framework"):
                _nonblank(
                    match.get(field), f"{name}.targets[{index}].match.{field}"
                )
            _optional_string(
                match.get("confidence"),
                f"{name}.targets[{index}].match.confidence",
            )
    degraded = _list(
        harmonization.get("degraded_stages"), f"{name}.degraded_stages"
    )
    for index, stage in enumerate(degraded):
        _nonblank(stage, f"{name}.degraded_stages[{index}]")


class AtlasV2Reader:
    """Read the versioned wire format without importing its producer package."""

    def load(self, path: str | Path) -> AtlasV2ReadResult:
        with Path(path).open(encoding="utf-8") as handle:
            payload = json.load(handle)
        return self.from_mapping(payload)

    def from_mapping(self, value: Mapping[str, Any]) -> AtlasV2ReadResult:
        data = _mapping(value, "AtlasDocument")
        if "accessions" in data and "schema_version" not in data:
            raise AtlasV2Error(
                "Legacy Atlas v1 envelopes are not supported; "
                "use the pinned v1 tools."
            )
        _only(data, _TOP_LEVEL_FIELDS, "AtlasDocument")
        version = str(data.get("schema_version", ""))
        if version != SCHEMA_VERSION:
            raise AtlasV2Error(f"unsupported Atlas schema version: {version!r}")

        atlas = _mapping(data.get("atlas"), "atlas")
        _only(atlas, {"atlas_id", "title", "theme"}, "atlas")
        _nonblank(atlas.get("atlas_id"), "atlas.atlas_id")
        _nonblank(atlas.get("title"), "atlas.title")
        _nonblank(atlas.get("theme"), "atlas.theme")

        run = _mapping(data.get("run"), "run")
        _only(run, {"run_id", "created_at", "config", "status"}, "run")
        _nonblank(run.get("run_id"), "run.run_id")
        _nonblank(run.get("created_at"), "run.created_at")
        if run.get("status") not in _RUN_STATUSES:
            raise AtlasV2Error("run.status is not a supported value")
        config = _mapping(run.get("config"), "run.config")
        _only(
            config,
            {"queries", "metadata_repositories", "options"},
            "run.config",
        )
        _list(config.get("queries"), "run.config.queries")
        _list(
            config.get("metadata_repositories"),
            "run.config.metadata_repositories",
        )
        _mapping(config.get("options"), "run.config.options")

        publications = _list(data.get("publications"), "publications")
        publication_ids: list[str] = []
        for index, raw in enumerate(publications):
            publication = _mapping(raw, f"publications[{index}]")
            _only(
                publication,
                {"publication_id", "title", "abstract", "review"},
                f"publications[{index}]",
            )
            publication_id = _nonblank(
                publication.get("publication_id"),
                f"publications[{index}].publication_id",
            )
            publication_ids.append(publication_id)
            _optional_string(
                publication.get("title"), f"publications[{index}].title"
            )
            _optional_string(
                publication.get("abstract"), f"publications[{index}].abstract"
            )
            _review(publication.get("review"), f"publications[{index}].review")
        if len(publication_ids) != len(set(publication_ids)):
            raise AtlasV2Error("duplicate publication_id in AtlasDocument")

        raw_datasets = _list(data.get("datasets"), "datasets")
        dataset_ids: list[str] = []
        converted: list[AtlasV2Dataset] = []
        warnings: list[str] = []
        failed_count = 0
        harmonized_count = 0
        for index, raw in enumerate(raw_datasets):
            dataset = _mapping(raw, f"datasets[{index}]")
            _only(dataset, _DATASET_FIELDS, f"datasets[{index}]")
            dataset_id = _nonblank(
                dataset.get("dataset_id"), f"datasets[{index}].dataset_id"
            )
            dataset_ids.append(dataset_id)
            repository = _nonblank(
                dataset.get("source_repository"),
                f"datasets[{index}].source_repository",
            )
            ordinal = dataset.get("source_ordinal")
            if isinstance(ordinal, bool) or not isinstance(ordinal, int):
                raise AtlasV2Error(
                    f"datasets[{index}].source_ordinal must be an integer"
                )
            status = dataset.get("status")
            if status not in _DATASET_STATUSES:
                raise AtlasV2Error(
                    f"datasets[{index}].status is not a supported value"
                )
            metadata = _mapping(
                dataset.get("metadata"), f"datasets[{index}].metadata"
            )
            _review(dataset.get("review"), f"datasets[{index}].review")
            _harmonization(
                dataset.get("harmonization"), f"datasets[{index}].harmonization"
            )
            references = _list(
                dataset.get("publication_ids"),
                f"datasets[{index}].publication_ids",
            )
            unknown_publications = sorted(set(references) - set(publication_ids))
            if unknown_publications:
                raise AtlasV2Error(
                    f"dataset {dataset_id!r} references unknown publication: "
                    + ", ".join(unknown_publications)
                )
            diagnostics = _list(
                dataset.get("diagnostics"), f"datasets[{index}].diagnostics"
            )
            diagnostic_messages = []
            for diagnostic_index, raw_diagnostic in enumerate(diagnostics):
                diagnostic = _mapping(
                    raw_diagnostic,
                    f"datasets[{index}].diagnostics[{diagnostic_index}]",
                )
                _only(
                    diagnostic,
                    {"code", "message", "severity", "path"},
                    f"datasets[{index}].diagnostics[{diagnostic_index}]",
                )
                code = _nonblank(
                    diagnostic.get("code"),
                    f"datasets[{index}].diagnostics[{diagnostic_index}].code",
                )
                message = _nonblank(
                    diagnostic.get("message"),
                    f"datasets[{index}].diagnostics[{diagnostic_index}].message",
                )
                severity = diagnostic.get("severity")
                if severity not in _DIAGNOSTIC_SEVERITIES:
                    raise AtlasV2Error(
                        f"datasets[{index}].diagnostics[{diagnostic_index}]."
                        "severity is not a supported value"
                    )
                _optional_string(
                    diagnostic.get("path"),
                    f"datasets[{index}].diagnostics[{diagnostic_index}].path",
                )
                diagnostic_messages.append(f"{code}: {message}")

            if status == "harmonized":
                harmonized_count += 1
                converted.append(
                    AtlasV2Dataset(
                        dataset_id=dataset_id,
                        source_repository=repository,
                        source_ordinal=ordinal,
                        metadata=deepcopy(dict(metadata)),
                    )
                )
            else:
                if status == "failed":
                    failed_count += 1
                warning = (
                    f"{dataset_id}: dataset status {status}; "
                    "no convertible metadata"
                )
                if diagnostic_messages:
                    warning += "; " + "; ".join(diagnostic_messages)
                warnings.append(warning)

        if len(dataset_ids) != len(set(dataset_ids)):
            raise AtlasV2Error("duplicate dataset_id in AtlasDocument")

        summary = _mapping(data.get("summary"), "summary")
        _only(
            summary,
            {
                "dataset_count",
                "publication_count",
                "completed_dataset_count",
                "failed_dataset_count",
            },
            "summary",
        )
        expected_summary = {
            "dataset_count": len(raw_datasets),
            "publication_count": len(publications),
            "completed_dataset_count": harmonized_count,
            "failed_dataset_count": failed_count,
        }
        for field, expected in expected_summary.items():
            if summary.get(field) != expected:
                raise AtlasV2Error(f"summary {field} does not match document")

        return AtlasV2ReadResult(tuple(converted), tuple(warnings))


__all__ = [
    "SCHEMA_VERSION",
    "AtlasV2Dataset",
    "AtlasV2Error",
    "AtlasV2ReadResult",
    "AtlasV2Reader",
]
