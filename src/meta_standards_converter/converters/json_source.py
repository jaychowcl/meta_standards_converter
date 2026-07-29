# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Load parsed MINiML or ThematicAtlases JSON as study-scoped packages."""

from __future__ import annotations

import json
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping


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
    """Recognize MINiML payloads and completed ThematicAtlases accessions."""

    def load(self, path: str | Path) -> SourceLoadResult:
        source = Path(path)
        payload = json.loads(source.read_text(encoding="utf-8"))
        if isinstance(payload, Mapping) and isinstance(
            payload.get("accessions"), list
        ):
            return self._atlas(payload)
        packages = payload if isinstance(payload, list) else [payload]
        if not packages:
            raise ValueError(
                "Parsed MINiML JSON must contain a non-empty package object "
                "or list; expected a non-empty list of packages or an Atlas "
                "object"
            )
        for index, package in enumerate(packages, start=1):
            if not isinstance(package, Mapping):
                raise ValueError(
                    f"Parsed MINiML package {index} must be a JSON object."
                )
        return SourceLoadResult(self._group_packages(packages, source.stem))

    def _atlas(self, payload: Mapping[str, Any]) -> SourceLoadResult:
        packages_by_source: list[tuple[Mapping[str, Any], str]] = []
        warnings: list[str] = []
        for record in payload.get("accessions", []):
            if not isinstance(record, Mapping):
                continue
            accession = str(record.get("datalink_id") or "unknown accession")
            status = str(
                record.get("ontology_harmonization_run_status") or "unknown"
            )
            metadata = record.get("accession_metadata")
            if status != "completed" or not isinstance(metadata, (dict, list)):
                error = record.get("ontology_harmonization_error")
                if status == "error" and error:
                    warnings.append(f"{accession}: harmonization error: {error}")
                else:
                    warnings.append(
                        f"{accession}: harmonization status {status}; "
                        "no convertible metadata"
                    )
                continue
            values = metadata if isinstance(metadata, list) else [metadata]
            for package in values:
                if not isinstance(package, Mapping):
                    raise ValueError(
                        f"{accession}: accession_metadata must contain objects"
                    )
                packages_by_source.append((package, accession))

        grouped: dict[str, list[Mapping[str, Any]]] = {}
        sources: dict[str, str] = {}
        for package, accession in packages_by_source:
            dataset_id = self._dataset_id(package) or accession
            grouped.setdefault(dataset_id, []).append(package)
            sources.setdefault(dataset_id, accession)
        groups = tuple(
            DatasetPackageGroup(
                dataset_id,
                self._dedupe_samples(values),
                source_accession=sources[dataset_id],
            )
            for dataset_id, values in grouped.items()
        )
        return SourceLoadResult(groups, tuple(warnings))

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
