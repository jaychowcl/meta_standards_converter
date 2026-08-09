# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Scientific compatibility and assembly policy for combined H5AD datasets."""

from __future__ import annotations

from collections.abc import Callable, Mapping
from typing import Any


class DatasetCompatibilityError(ValueError):
    """A public-safe scientific reason that prevents dataset combination."""


class DatasetCombinationPolicy:
    """Validate evidence and combine compatible AnnData values deterministically."""

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
        if not adatas:
            raise DatasetCompatibilityError("No sample H5ADs were produced.")
        anndata, _numpy, _pandas, sparse = self._scientific_modules()
        missing_evidence = self.missing_combination_evidence(adatas)
        if missing_evidence and not allow_unverified:
            rendered = ", ".join(
                f"{dimension} ({', '.join(sample_ids)})"
                for dimension, sample_ids in sorted(missing_evidence.items())
            )
            raise DatasetCompatibilityError(
                "Cannot combine samples without positive compatibility evidence "
                f"for every sample: {rendered}. Set allow_unverified_combination "
                "only to publish an explicitly acknowledged partial result."
            )
        organisms = {
            str(adata.obs["msc.sample.channel.organism.value"].iloc[0]).strip()
            for adata in adatas.values()
            if "msc.sample.channel.organism.value" in adata.obs
            and str(adata.obs["msc.sample.channel.organism.value"].iloc[0]).strip()
        }
        if len(organisms) > 1:
            raise DatasetCompatibilityError(
                f"Cannot combine samples with incompatible organisms: {sorted(organisms)}"
            )
        references = {
            str(adata.uns.get("meta_standards_converter", {}).get("reference")).strip()
            for adata in adatas.values()
            if isinstance(adata.uns.get("meta_standards_converter"), dict)
            and adata.uns["meta_standards_converter"].get("reference")
        }
        if len(references) > 1:
            raise DatasetCompatibilityError(
                "Cannot combine samples with incompatible reference builds: "
                f"{sorted(references)}"
            )
        modalities = {
            str(adata.uns.get("meta_standards_converter", {}).get("modality")).strip()
            for adata in adatas.values()
            if isinstance(adata.uns.get("meta_standards_converter"), dict)
            and adata.uns["meta_standards_converter"].get("modality")
            not in (None, "", "unknown")
        }
        if len(modalities) > 1:
            raise DatasetCompatibilityError(
                "Cannot combine incompatible expression modalities: "
                f"{sorted(modalities)}"
            )
        namespaces = {self.feature_namespace(adata) for adata in adatas.values()}
        namespaces.discard("unknown")
        if len(namespaces) > 1:
            raise DatasetCompatibilityError(
                "Cannot combine incompatible feature identifier namespaces: "
                f"{sorted(namespaces)}"
            )
        combined = anndata.concat(
            adatas,
            axis="obs",
            join="outer",
            merge="first",
            label="msc.combination.batch",
            index_unique=None,
            fill_value=0,
        )
        if not sparse.issparse(combined.X):
            combined.X = sparse.csr_matrix(combined.X)
        else:
            combined.X = combined.X.tocsr()
        combined.uns["meta_standards_converter"] = {
            "combined_samples": list(adatas),
            "join": "outer",
            "fill_value": 0,
            "converter_version": self._package_version(),
            "metadata_schema_version": self._metadata_schema_version,
            "path_base": "artifact_parent",
            "combination_compatibility": {
                "verified": not missing_evidence,
                "missing_evidence": missing_evidence,
                "explicitly_acknowledged": bool(
                    missing_evidence and allow_unverified
                ),
            },
            "sample_provenance": {
                sample_id: dict(adata.uns.get("meta_standards_converter", {}))
                for sample_id, adata in adatas.items()
                if isinstance(adata.uns.get("meta_standards_converter"), dict)
            },
        }
        sample_values = {}
        for sample_id, adata in adatas.items():
            values = adata.uns.get("msc_metadata", {}).get("sample_values")
            fields = {}
            if values is not None:
                for row in values.to_dict("records"):
                    fields.setdefault(str(row["field"]), []).append(row["value"])
            fields["msc.combination.batch"] = [sample_id]
            sample_values[sample_id] = fields
        self._attach_sample_values(combined, sample_values)
        return combined

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
            if not any(values[dimension] for values in evidence.values()):
                continue
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
        if values and sum(
            value.startswith(("ENSG", "ENSMUSG", "ENSRNOG")) for value in values
        ) >= len(values) / 2:
            return "ensembl"
        if values and sum(value.isalnum() for value in values) >= len(values) / 2:
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
