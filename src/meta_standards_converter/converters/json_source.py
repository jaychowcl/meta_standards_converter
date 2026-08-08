# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Load parsed MINiML or canonical Atlas v1 JSON as study-scoped packages."""

from __future__ import annotations

import json
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

from meta_standards_converter.atlas_v1 import AtlasV1Reader
from meta_standards_converter.miniml import MINiMLCodec, MINiMLPackage, MINiMLValidationIssue


@dataclass(frozen=True)
class DatasetPackageGroup:
    """One study-sized group of parsed MINiML packages."""

    dataset_id: str
    packages: tuple[MINiMLPackage, ...]
    source_accession: str | None = None
    harmonization_overrides: Mapping[str, Any] | None = None
    source_packages: tuple[MINiMLPackage, ...] | None = None
    harmonization_resolution: Any | None = None

    def resolved(self, *, enabled: bool) -> "DatasetPackageGroup":
        from .harmonization_overrides import resolve_harmonization_overrides

        codec = MINiMLCodec()
        resolution = resolve_harmonization_overrides(
            tuple(codec.encode(package) for package in self.packages),
            self.harmonization_overrides,
            enabled=enabled,
        )
        return DatasetPackageGroup(
            self.dataset_id,
            codec.decode_many(resolution.packages).packages,
            source_accession=self.source_accession,
            harmonization_overrides=self.harmonization_overrides,
            source_packages=self.packages,
            harmonization_resolution=resolution,
        )


@dataclass(frozen=True)
class SourceLoadResult:
    """Groups and non-fatal source diagnostics."""

    groups: tuple[DatasetPackageGroup, ...]
    warnings: tuple[str, ...] = ()
    diagnostics: tuple[MINiMLValidationIssue, ...] = ()


class JSONPackageSource:
    """Recognize native MINiML payloads and canonical Atlas v1 documents."""

    def __init__(self, atlas_reader: AtlasV1Reader | None = None) -> None:
        self._atlas_reader = atlas_reader or AtlasV1Reader()
        self._codec = MINiMLCodec()

    def load(self, path: str | Path) -> SourceLoadResult:
        source = Path(path)
        payload = json.loads(source.read_text(encoding="utf-8"))
        if isinstance(payload, Mapping) and "miniml_json" in payload:
            return self._validate_result(self._agentic_envelope(payload, source.stem))
        if isinstance(payload, Mapping) and (
            "schema_version" in payload or "accessions" in payload
        ):
            return self._validate_result(self._atlas_v1(payload))
        packages = payload if isinstance(payload, list) else [payload]
        if not packages:
            raise ValueError("JSON source contains no convertible package groups.")
        for index, package in enumerate(packages, start=1):
            if not isinstance(package, Mapping):
                raise ValueError(
                    f"Parsed MINiML package {index} must be a JSON object."
                )
        return self._validate_result(
            self._source_result(self._group_packages(packages, source.stem))
        )

    def _agentic_envelope(
        self, payload: Mapping[str, Any], fallback: str
    ) -> SourceLoadResult:
        document = payload.get("miniml_json")
        packages = document if isinstance(document, list) else [document]
        if not packages or any(not isinstance(package, Mapping) for package in packages):
            raise ValueError("Agentic Curator miniml_json must be an object or non-empty object list.")
        profile = payload.get("harmonization_overrides")
        groups = self._group_packages(packages, fallback)
        return self._source_result(tuple(
            DatasetPackageGroup(
                group.dataset_id,
                group.packages,
                source_accession=group.source_accession,
                harmonization_overrides=deepcopy(profile) if isinstance(profile, Mapping) else profile,
            )
            for group in groups
        ))

    def _validate_result(self, result: SourceLoadResult) -> SourceLoadResult:
        if not result.groups:
            raise ValueError("JSON source contains no convertible package groups.")
        if not any(
            package.samples
            for group in result.groups
            for package in group.packages
        ):
            raise ValueError("JSON source contains no convertible samples.")
        return result

    @staticmethod
    def _as_list(value: Any) -> list[Any]:
        if value is None:
            return []
        return value if isinstance(value, list) else [value]

    def _atlas_v1(self, payload: Mapping[str, Any]) -> SourceLoadResult:
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
        return self._source_result(groups, result.warnings)

    @staticmethod
    def _atlas_packages(
        dataset_id: str, metadata: Mapping[str, Any]
    ) -> list[Mapping[str, Any]]:
        if set(metadata) != {"packages"}:
            return [metadata]
        packages = metadata["packages"]
        if not isinstance(packages, list):
            raise ValueError(
                f"Atlas v1 dataset {dataset_id} metadata.packages must be a list."
            )
        if any(not isinstance(package, Mapping) for package in packages):
            raise ValueError(
                f"Atlas v1 dataset {dataset_id} metadata.packages must contain objects."
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
    ) -> tuple[MINiMLPackage, ...]:
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
        return self._codec.decode_many(retained).packages

    def _source_result(
        self,
        groups: tuple[DatasetPackageGroup, ...],
        warnings: tuple[str, ...] = (),
    ) -> SourceLoadResult:
        diagnostics = tuple(
            issue
            for group in groups
            for package in group.packages
            for issue in package.validate()
        )
        return SourceLoadResult(groups, tuple(warnings), diagnostics)

    def _dataset_id(self, package: Mapping[str, Any]) -> str | None:
        series = package.get("series")
        if not isinstance(series, Mapping):
            return None
        iid = series.get("iid")
        if iid:
            return str(iid)
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
        return None

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
