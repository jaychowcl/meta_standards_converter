# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Load parsed MINiML or canonical Atlas v2 JSON as study-scoped packages."""

from __future__ import annotations

import json
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

from meta_standards_converter.atlas_v2 import AtlasV2Reader


@dataclass(frozen=True)
class DatasetPackageGroup:
    """One study-sized group of parsed MINiML packages."""

    dataset_id: str
    packages: tuple[dict[str, Any], ...]
    source_accession: str | None = None


@dataclass(frozen=True)
class SourceLoadResult:
    """Groups and non-fatal source diagnostics."""

    groups: tuple[DatasetPackageGroup, ...]
    warnings: tuple[str, ...] = ()


class JSONPackageSource:
    """Recognize native MINiML payloads and canonical Atlas v2 documents."""

    def __init__(self, atlas_reader: AtlasV2Reader | None = None) -> None:
        self._atlas_reader = atlas_reader or AtlasV2Reader()

    def load(self, path: str | Path) -> SourceLoadResult:
        source = Path(path)
        payload = json.loads(source.read_text(encoding="utf-8"))
        if isinstance(payload, Mapping) and (
            "schema_version" in payload or "accessions" in payload
        ):
            return self._validate_result(self._atlas_v2(payload))
        packages = payload if isinstance(payload, list) else [payload]
        if not packages:
            raise ValueError("JSON source contains no convertible package groups.")
        for index, package in enumerate(packages, start=1):
            if not isinstance(package, Mapping):
                raise ValueError(
                    f"Parsed MINiML package {index} must be a JSON object."
                )
        return self._validate_result(
            SourceLoadResult(self._group_packages(packages, source.stem))
        )

    def _validate_result(self, result: SourceLoadResult) -> SourceLoadResult:
        if not result.groups:
            raise ValueError("JSON source contains no convertible package groups.")
        if not any(
            isinstance(sample, Mapping)
            for group in result.groups
            for package in group.packages
            for sample in self._as_list(package.get("sample"))
        ):
            raise ValueError("JSON source contains no convertible samples.")
        return result

    @staticmethod
    def _as_list(value: Any) -> list[Any]:
        if value is None:
            return []
        return value if isinstance(value, list) else [value]

    def _atlas_v2(self, payload: Mapping[str, Any]) -> SourceLoadResult:
        result = self._atlas_reader.from_mapping(payload)
        groups = tuple(
            DatasetPackageGroup(
                dataset.dataset_id,
                self._dedupe_samples(
                    self._atlas_packages(dataset.dataset_id, dataset.metadata)
                ),
                source_accession=dataset.dataset_id,
            )
            for dataset in result.datasets
        )
        return SourceLoadResult(groups, result.warnings)

    @staticmethod
    def _atlas_packages(
        dataset_id: str, metadata: Mapping[str, Any]
    ) -> list[Mapping[str, Any]]:
        if set(metadata) != {"packages"}:
            return [metadata]
        packages = metadata["packages"]
        if not isinstance(packages, list):
            raise ValueError(
                f"Atlas v2 dataset {dataset_id} metadata.packages must be a list."
            )
        if any(not isinstance(package, Mapping) for package in packages):
            raise ValueError(
                f"Atlas v2 dataset {dataset_id} metadata.packages must contain objects."
            )
        return packages

    def _group_packages(
        self, packages: list[Mapping[str, Any]], fallback: str
    ) -> tuple[DatasetPackageGroup, ...]:
        grouped: dict[str, list[Mapping[str, Any]]] = {}
        for package in packages:
            dataset_id = self._dataset_id(package) or fallback
            grouped.setdefault(dataset_id, []).append(package)
        return tuple(
            DatasetPackageGroup(dataset_id, self._dedupe_samples(values))
            for dataset_id, values in grouped.items()
        )

    def _dedupe_samples(
        self, packages: list[Mapping[str, Any]]
    ) -> tuple[dict[str, Any], ...]:
        seen: dict[str, Mapping[str, Any]] = {}
        retained: list[dict[str, Any]] = []
        for package in packages:
            copied = deepcopy(dict(package))
            raw_samples = copied.get("sample", [])
            samples = raw_samples if isinstance(raw_samples, list) else [raw_samples]
            unique = []
            for sample in samples:
                if not isinstance(sample, Mapping):
                    continue
                accession = self._sample_accession(sample)
                if not accession:
                    unique.append(dict(sample))
                    continue
                existing = seen.get(accession)
                if existing is not None:
                    if existing != sample:
                        raise ValueError(
                            f"Sample {accession} has conflicting metadata"
                        )
                    continue
                seen[accession] = sample
                unique.append(dict(sample))
            copied["sample"] = unique
            if unique or not samples:
                retained.append(copied)
        return tuple(retained)

    def _dataset_id(self, package: Mapping[str, Any]) -> str | None:
        series = package.get("series")
        if not isinstance(series, Mapping):
            return None
        accessions = series.get("accession", [])
        values = accessions if isinstance(accessions, list) else [accessions]
        for accession in values:
            value = (
                accession.get("value")
                if isinstance(accession, Mapping)
                else accession
            )
            if value:
                return str(value)
        return str(series.get("iid")) if series.get("iid") else None

    def _sample_accession(self, sample: Mapping[str, Any]) -> str | None:
        accessions = sample.get("accession", [])
        values = accessions if isinstance(accessions, list) else [accessions]
        for accession in values:
            value = (
                accession.get("value")
                if isinstance(accession, Mapping)
                else accession
            )
            if value:
                return str(value)
        return str(sample.get("iid")) if sample.get("iid") else None
