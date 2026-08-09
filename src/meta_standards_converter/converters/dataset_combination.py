# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Scientific compatibility evidence for catalogue H5AD datasets."""

from __future__ import annotations

from collections.abc import Callable, Mapping
import re
from typing import Any


class DatasetCompatibilityError(ValueError):
    """A public-safe scientific reason that prevents dataset combination."""


class DatasetCombinationPolicy:
    """Inspect compatibility evidence without claiming matrix integration."""

    def __init__(
        self,
        *,
        scientific_modules: Callable[[], tuple[Any, Any, Any, Any]],
        attach_sample_values: Callable[[Any, dict[str, dict]], None],
        package_version: Callable[[], str],
        metadata_schema_version: str,
    ) -> None:
        self._scientific_modules = scientific_modules
        self._attach_sample_values = attach_sample_values
        self._package_version = package_version
        self._metadata_schema_version = metadata_schema_version

    def combine(
        self,
        adatas: dict[str, object],
        *,
        allow_unverified: bool = False,
    ):
        raise DatasetCompatibilityError(
            "Expression matrix combination is disabled: outer concatenation is "
            "not scientific integration. Publish the per-sample H5AD catalogue "
            "and perform an explicit integration workflow separately."
        )

    def missing_combination_evidence(
        self,
        adatas: Mapping[str, object],
    ) -> dict[str, list[str]]:
        if len(adatas) < 2:
            return {}
        evidence: dict[str, dict[str, str | None]] = {}
        for sample_id, adata in adatas.items():
            organism = None
            if "msc.sample.channel.organism.value" in adata.obs:
                values = {
                    str(value).strip()
                    for value in adata.obs["msc.sample.channel.organism.value"]
                    if str(value).strip()
                }
                if len(values) == 1:
                    organism = next(iter(values))
            provenance = adata.uns.get("meta_standards_converter")
            provenance = provenance if isinstance(provenance, dict) else {}
            reference = str(provenance.get("reference") or "").strip() or None
            raw_modality = str(provenance.get("modality") or "").strip()
            modality = (
                raw_modality
                if raw_modality and raw_modality.casefold() != "unknown"
                else None
            )
            namespace = self.feature_namespace(adata)
            evidence[sample_id] = {
                "organism": organism,
                "reference": reference,
                "modality": modality,
                "feature_namespace": (
                    namespace if namespace != "unknown" else None
                ),
            }
        missing: dict[str, list[str]] = {}
        for dimension in (
            "organism",
            "reference",
            "modality",
            "feature_namespace",
        ):
            absent = [
                sample_id
                for sample_id, values in evidence.items()
                if not values[dimension]
            ]
            if absent:
                missing[dimension] = absent
        return missing

    @staticmethod
    def feature_namespace(adata) -> str:
        values = adata.var.get("gene_ids", adata.var_names)
        values = [
            str(value).split(".")[0].upper()
            for value in list(values)[:100]
            if value
        ]
        if values and all(
            re.fullmatch(r"ENS(?:G|MUSG|RNOG)\d+", value) for value in values
        ):
            return "ensembl"
        if values and all(value.isdecimal() for value in values):
            return "entrez"
        if values and all(
            re.fullmatch(r"[A-Z][A-Z0-9.-]*", value)
            and any(character.isalpha() for character in value)
            for value in values
        ):
            return "symbol"
        return "unknown"

    @staticmethod
    def declared_reference(adata) -> str | None:
        for key in ("genome", "reference_genome", "genome_build"):
            value = adata.uns.get(key)
            if isinstance(value, (str, int, float)) and str(value).strip():
                return str(value).strip()
        existing = adata.uns.get("meta_standards_converter")
        if isinstance(existing, dict):
            for key in ("reference", "genome", "genome_build"):
                value = existing.get(key)
                if value:
                    return str(value).strip()
        if "genome" in adata.var:
            values = {
                str(value).strip()
                for value in adata.var["genome"]
                if str(value).strip()
            }
            if len(values) == 1:
                return next(iter(values))
        return None
