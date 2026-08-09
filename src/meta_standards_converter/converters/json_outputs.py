# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Public orchestration for JSON-origin manifest and AnnData outputs."""

from __future__ import annotations

import json
import tempfile
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any, Mapping, Sequence

from meta_standards_converter.artifact_bundle import (
    DurableArtifactBundlePublisher,
    PublishedArtifactBundle,
)

from .json2h5ad import BatchConversionResult, ConversionResult, JSON2H5ADConverter
from .json2tabular import JSON2TSVConverter, TabularMetadataProjector


@dataclass
class AnnDataMetadataExportResult:
    """Aggregated observation metadata and optional catalogue sidecars."""

    study_accession: str
    obs: Any
    obs_path: str
    manifest_path: str
    var: Any | None = None
    var_path: str | None = None
    uns: Mapping[str, Any] | None = None
    uns_path: str | None = None
    warnings: tuple[str, ...] = ()
    errors: tuple[str, ...] = ()
    failures: tuple[str, ...] = ()
    bundle_pointer_path: str | None = None

    @property
    def partial(self) -> bool:
        return bool(self.errors or self.failures)

    def to_dict(self) -> dict[str, Any]:
        artifacts = {"obs": self.obs_path, "manifest": self.manifest_path}
        if self.bundle_pointer_path:
            artifacts["bundle_pointer"] = self.bundle_pointer_path
        if self.var_path:
            artifacts["var"] = self.var_path
        if self.uns_path:
            artifacts["uns"] = self.uns_path
        return {
            "operation": "anndata_metadata",
            "status": "partial" if self.partial else "complete",
            "study_accession": self.study_accession,
            "expression_integration": "none",
            "metadata_aggregation": "observation_rows",
            "artifacts": artifacts,
            "obs": {
                "rows": int(self.obs.shape[0]),
                "columns": ["cell_id", *map(str, self.obs.columns)],
            },
            "var": None
            if self.var is None
            else {
                "rows": int(self.var.shape[0]),
                "columns": ["feature_id", *map(str, self.var.columns)],
            },
            "uns": None if self.uns is None else {"keys": sorted(map(str, self.uns))},
            "warnings": list(self.warnings),
            "errors": list(self.errors),
            "failures": list(self.failures),
        }


@dataclass
class AnnDataMetadataBatchResult:
    conversions: dict[str, AnnDataMetadataExportResult]
    failures: tuple[str, ...] = ()
    warnings: tuple[str, ...] = ()

    @property
    def partial(self) -> bool:
        return bool(self.failures) or any(item.partial for item in self.conversions.values())

    def to_dict(self) -> dict[str, Any]:
        return {
            "operation": "anndata_metadata",
            "status": "partial" if self.partial else "complete",
            "datasets": {
                key: value.to_dict() for key, value in self.conversions.items()
            },
            "warnings": list(self.warnings),
            "failures": list(self.failures),
        }


class JSONDataOutputOrchestrator:
    """Coordinate manifest, H5AD, and AnnData-component JSON workflows."""

    def __init__(
        self,
        *,
        tabular_projectors: Sequence[TabularMetadataProjector] | None = None,
        h5ad_converter: JSON2H5ADConverter | None = None,
    ) -> None:
        self.tabular_projectors = tabular_projectors
        self.h5ad_converter = h5ad_converter or JSON2H5ADConverter()

    def export_manifest(
        self,
        source: str | Path,
        *,
        outdir: str | Path,
        output_format: str = "tsv",
        allow_invalid: bool = False,
        overwrite: bool = False,
        use_harmonization_overrides: bool = False,
    ):
        output_dir = Path(outdir)
        output_dir.parent.mkdir(parents=True, exist_ok=True)
        suffix = ".tsv" if output_format == "tsv" else ".csv"
        table = output_dir / f"{Path(source).stem}{suffix}"
        manifest = output_dir / f"{Path(source).stem}.json2tsv.json"
        existing = [path for path in (table, manifest) if path.exists()]
        if existing and not overwrite:
            raise FileExistsError(f"Output already exists: {existing[0]}")
        converter = JSON2TSVConverter(
            metadata_projectors=self.tabular_projectors,
            output_format=output_format,
        )
        with tempfile.TemporaryDirectory(
            prefix=f".{output_dir.name}.manifest-", dir=output_dir.parent
        ) as temporary:
            staged_table = Path(temporary) / table.name
            staged_manifest = Path(temporary) / manifest.name
            converted = converter.convert_source(
                source,
                staged_table,
                allow_invalid=allow_invalid,
                overwrite=True,
                use_harmonization_overrides=use_harmonization_overrides,
            )
            result = replace(
                converted,
                output_path=str(table),
                manifest_path=str(manifest),
            )
            staged_manifest.write_text(
                json.dumps(result.to_dict(), indent=2, ensure_ascii=False, sort_keys=True)
                + "\n",
                encoding="utf-8",
            )
            published = _publish_bundle(
                {"table": staged_table, "manifest": staged_manifest},
                {"table": table, "manifest": manifest},
                overwrite=overwrite,
            )
        return replace(result, bundle_pointer_path=str(published.pointer_path))

    def export_h5ad(self, source: str | Path, *, outdir: str | Path, **options):
        return self.h5ad_converter.convert(str(source), out=str(outdir), **options)

    def export_anndata_metadata(
        self,
        source: str | Path,
        *,
        outdir: str | Path,
        include_var: bool = False,
        include_uns: bool = False,
        overwrite: bool = False,
        **options,
    ) -> AnnDataMetadataExportResult | AnnDataMetadataBatchResult:
        destination = Path(outdir)
        destination.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(
            prefix=f".{destination.name}.json2obs-",
            dir=destination.parent,
        ) as temporary:
            assembly_root = Path(temporary) / "assembly"
            converted = self.h5ad_converter.convert(
                str(source),
                out=str(assembly_root),
                overwrite=True,
                **options,
            )
            if isinstance(converted, BatchConversionResult):
                results: dict[str, AnnDataMetadataExportResult] = {}
                for dataset_id, conversion in converted.conversions.items():
                    target = destination / dataset_id
                    results[dataset_id] = self._export_components(
                        conversion,
                        target,
                        include_var=include_var,
                        include_uns=include_uns,
                        overwrite=overwrite,
                    )
                return AnnDataMetadataBatchResult(
                    conversions=results,
                    failures=tuple(converted.failures),
                    warnings=tuple(converted.warnings),
                )
            return self._export_components(
                converted,
                destination,
                include_var=include_var,
                include_uns=include_uns,
                overwrite=overwrite,
            )

    def _export_components(
        self,
        conversion: ConversionResult,
        destination: Path,
        *,
        include_var: bool,
        include_uns: bool,
        overwrite: bool,
    ) -> AnnDataMetadataExportResult:
        if not conversion.sample_h5ads:
            raise ValueError(
                "Cannot export observations because no sample H5AD was produced: "
                + "; ".join(conversion.failures)
            )
        anndata, _numpy, pandas, _sparse = self.h5ad_converter._scientific_modules()
        if include_var and len(conversion.sample_h5ads) != 1:
            raise ValueError(
                "A combined var table is not available for a multi-sample "
                "catalogue; export each sample H5AD separately."
            )
        obs_frames = []
        var = None
        sample_uns: dict[str, Any] = {}
        for sample_id, h5ad_path in conversion.sample_h5ads.items():
            adata = anndata.read_h5ad(h5ad_path, backed="r")
            try:
                obs_frames.append(adata.obs.copy())
                if include_var:
                    var = adata.var.copy()
                if include_uns:
                    sample_uns[sample_id] = dict(adata.uns)
            finally:
                adata.file.close()
        obs = pandas.concat(obs_frames, axis="index", join="outer", sort=False)
        uns = None
        if include_uns:
            uns = (
                next(iter(sample_uns.values()))
                if len(sample_uns) == 1
                else {
                    "catalogue_contract_version": "1.0",
                    "expression_integration": "none",
                    "samples": sample_uns,
                }
            )

        destinations = {
            "obs": destination / f"{conversion.study_accession}.obs.csv",
            "manifest": destination / f"{conversion.study_accession}.json2obs.json",
        }
        if include_var:
            destinations["var"] = destination / f"{conversion.study_accession}.var.csv"
        if include_uns:
            destinations["uns"] = destination / f"{conversion.study_accession}.uns.json"
        existing = [path for path in destinations.values() if path.exists()]
        if existing and not overwrite:
            raise FileExistsError(f"Output already exists: {existing[0]}")

        destination.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(
            prefix=f".{destination.name}.publish-", dir=destination.parent
        ) as temporary:
            stage = Path(temporary)
            staged = {name: stage / path.name for name, path in destinations.items()}
            obs.rename_axis("cell_id").to_csv(staged["obs"], lineterminator="\n")
            if var is not None:
                var.rename_axis("feature_id").to_csv(staged["var"], lineterminator="\n")
            if uns is not None:
                staged["uns"].write_text(
                    json.dumps(
                        {"schema_version": "1.0", "values": _json_value(uns)},
                        indent=2,
                        ensure_ascii=False,
                        sort_keys=True,
                    )
                    + "\n",
                    encoding="utf-8",
                )
            result = AnnDataMetadataExportResult(
                study_accession=conversion.study_accession,
                obs=obs,
                obs_path=str(destinations["obs"]),
                manifest_path=str(destinations["manifest"]),
                var=var,
                var_path=str(destinations["var"]) if include_var else None,
                uns=uns,
                uns_path=str(destinations["uns"]) if include_uns else None,
                warnings=tuple(conversion.warnings),
                errors=tuple(conversion.errors),
                failures=tuple(conversion.failures),
            )
            staged["manifest"].write_text(
                json.dumps(result.to_dict(), indent=2, ensure_ascii=False, sort_keys=True)
                + "\n",
                encoding="utf-8",
            )
            published = _publish_bundle(staged, destinations, overwrite=overwrite)
        return replace(result, bundle_pointer_path=str(published.pointer_path))


def _json_value(value: Any) -> Any:
    """Convert AnnData metadata to reconstructable, JSON-safe values."""
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, Mapping):
        return {str(key): _json_value(item) for key, item in value.items()}
    if hasattr(value, "to_dict") and hasattr(value, "columns") and hasattr(value, "index"):
        return {
            "type": "dataframe",
            "columns": [str(item) for item in value.columns],
            "index": [_json_value(item) for item in value.index.tolist()],
            "records": [
                {str(key): _json_value(item) for key, item in row.items()}
                for row in value.to_dict(orient="records")
            ],
        }
    if hasattr(value, "tolist"):
        values = value.tolist()
        if isinstance(values, list):
            return {
                "type": "ndarray",
                "dtype": str(getattr(value, "dtype", "object")),
                "shape": list(getattr(value, "shape", (len(values),))),
                "values": _json_value(values),
            }
        return _json_value(values)
    if isinstance(value, (list, tuple)):
        return [_json_value(item) for item in value]
    raise TypeError(f"Unsupported AnnData uns value: {type(value).__name__}")


def _publish_bundle(
    staged: Mapping[str, Path],
    destinations: Mapping[str, Path],
    *,
    overwrite: bool,
) -> PublishedArtifactBundle:
    """Publish one crash-durable generation plus legacy compatibility views."""

    return DurableArtifactBundlePublisher().publish(
        staged,
        destinations,
        overwrite=overwrite,
    )
