# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Convert parsed MINiML JSON packages into annotated H5AD datasets."""

from __future__ import annotations

import json
import hashlib
import logging
import os
import tempfile
from dataclasses import dataclass, field
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from urllib.parse import urlparse


logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class Asset:
    """One processed or raw data source associated with a sample or study."""

    scope_id: str
    path: str
    kind: str
    role: str = "primary"
    source: str = "json"
    members: tuple[dict, ...] = ()
    features_path: str | None = None
    barcodes_path: str | None = None
    orientation: str = "auto"
    md5: str | None = None


@dataclass
class PipelineRun:
    pipeline: str
    revision: str
    command: list[str]
    work_dir: str
    out_dir: str
    returncode: int | None = None
    log_path: str | None = None


@dataclass
class ConversionResult:
    """Files and diagnostics produced for one parsed GEO study."""

    study_accession: str
    combined_h5ad: str | None = None
    sample_h5ads: dict[str, str] = field(default_factory=dict)
    retained_h5ads: list[str] = field(default_factory=list)
    pipeline_runs: list[PipelineRun] = field(default_factory=list)
    manifest_path: str | None = None
    warnings: list[str] = field(default_factory=list)
    failures: list[str] = field(default_factory=list)

    @property
    def primary_h5ad(self) -> str | None:
        if self.combined_h5ad:
            return self.combined_h5ad
        return next(iter(self.sample_h5ads.values()), None)

    @property
    def partial(self) -> bool:
        return bool(self.failures)

    def __str__(self) -> str:
        return self.primary_h5ad or self.manifest_path or self.study_accession


class SourcePlanner:
    """Discover assets and select the best available source for every sample."""

    SOURCE_RANK = {"json": 0, "cli": 1, "manifest": 2}
    KIND_RANK = {"raw": 0, "matrix": 1, "h5ad": 2}

    def plan(
        self,
        packages: list[dict],
        explicit_assets: list[Asset] | None = None,
        force_reprocess: bool = False,
    ) -> dict[str, Asset]:
        assets = self.discover(packages)
        assets.extend(explicit_assets or [])
        samples = self.samples(packages)
        planned = {}

        for sample_id in samples:
            candidates = [asset for asset in assets if asset.scope_id == sample_id]
            if force_reprocess:
                candidates = [asset for asset in candidates if asset.kind == "raw"]
                if not candidates:
                    raise ValueError(f"{sample_id} has no raw FASTQ files for forced reprocessing.")
            if not candidates:
                raise ValueError(f"{sample_id} has no supported H5AD, matrix, or raw FASTQ source.")
            planned[sample_id] = max(
                candidates,
                key=lambda asset: (
                    self.SOURCE_RANK.get(asset.source, -1),
                    self.KIND_RANK.get(asset.kind, -1),
                ),
            )
        return planned

    def discover(self, packages: list[dict]) -> list[Asset]:
        assets = []
        for package in packages:
            for sample in self._as_list(package.get("sample")):
                if not isinstance(sample, dict):
                    continue
                sample_id = self.sample_accession(sample)
                if not sample_id:
                    continue
                for key in ("supplementary_data", "raw_data"):
                    for entry in self._as_list(sample.get(key)):
                        path = self._value(entry)
                        kind = self.classify(path)
                        if path and kind:
                            assets.append(Asset(sample_id, path, kind))

                fastqs = []
                for run in self._as_list(sample.get("sra_run")):
                    if not isinstance(run, dict):
                        continue
                    for fastq in self._as_list(run.get("fastq_files")):
                        if isinstance(fastq, dict) and (fastq.get("uri") or fastq.get("filename")):
                            fastqs.append(dict(fastq, run=run.get("run")))
                if fastqs:
                    assets.append(
                        Asset(
                            sample_id,
                            fastqs[0].get("uri") or fastqs[0].get("filename"),
                            "raw",
                            members=tuple(fastqs),
                        )
                    )
        return assets

    def samples(self, packages: list[dict]) -> list[str]:
        values = []
        for package in packages:
            for sample in self._as_list(package.get("sample")):
                if isinstance(sample, dict):
                    accession = self.sample_accession(sample)
                    if accession and accession not in values:
                        values.append(accession)
        return values

    def sample_accession(self, sample: dict) -> str | None:
        for accession in self._as_list(sample.get("accession")):
            value = self._value(accession)
            if isinstance(value, str) and value.upper().startswith("GSM"):
                return value.upper()
        value = sample.get("iid")
        return str(value) if value else None

    def classify(self, path: str | None) -> str | None:
        if not path:
            return None
        filename = os.path.basename(urlparse(str(path)).path).lower()
        for suffix in (".gz", ".bz2", ".xz", ".zip"):
            if filename.endswith(suffix):
                filename = filename[: -len(suffix)]
                break
        if filename.endswith(".h5ad"):
            return "h5ad"
        if filename.endswith((".h5", ".mtx", ".csv", ".tsv", ".txt")):
            return "matrix"
        if filename.endswith((".fastq", ".fq")):
            return "raw"
        return None

    def _value(self, value):
        if isinstance(value, dict):
            return value.get("value") or value.get("uri") or value.get("filename")
        return value

    def _as_list(self, value) -> list:
        if value is None:
            return []
        return value if isinstance(value, list) else [value]


class JSON2H5ADConverter:
    """Top-level JSON-to-H5AD conversion orchestrator."""

    def __init__(self, planner: SourcePlanner | None = None):
        self.planner = planner or SourcePlanner()

    def convert(
        self,
        json_path: str,
        out: str | None = None,
        explicit_assets: list[Asset] | None = None,
        force_reprocess: bool = False,
        matrix_orientation: str = "auto",
        overwrite: bool = False,
        **options,
    ) -> ConversionResult:
        if not os.path.exists(json_path):
            raise FileNotFoundError(f"MINiML JSON file not found: {json_path}")
        with open(json_path, encoding="utf-8") as handle:
            packages = json.load(handle)
        if not isinstance(packages, list) or not packages:
            raise ValueError("Parsed MINiML JSON must contain a non-empty list of packages.")

        study_accession = self._study_accession(packages) or Path(json_path).stem
        out_path = Path(out or ".")
        out_path.mkdir(parents=True, exist_ok=True)
        planned = self.planner.plan(
            packages,
            explicit_assets=explicit_assets,
            force_reprocess=force_reprocess,
        )
        sample_lookup = self._sample_lookup(packages)
        result = ConversionResult(study_accession=study_accession)
        adatas = {}

        for sample_id, asset in planned.items():
            if asset.kind == "raw":
                raise NotImplementedError(
                    f"Raw FASTQ processing for {sample_id} requires nf-core orchestration."
                )
            adata = self._read_processed_asset(
                asset,
                orientation=(asset.orientation if asset.orientation != "auto" else matrix_orientation),
            )
            self._normalize(
                adata,
                sample=sample_lookup[sample_id],
                study_accession=study_accession,
                asset=asset,
            )
            sample_path = out_path / f"{sample_id}.h5ad"
            self._write_h5ad(adata, sample_path, overwrite=overwrite)
            result.sample_h5ads[sample_id] = str(sample_path)
            adatas[sample_id] = adata

        try:
            combined = self._combine(adatas)
        except ValueError as exc:
            result.failures.append(str(exc))
        else:
            combined_path = out_path / f"{study_accession}.h5ad"
            self._write_h5ad(combined, combined_path, overwrite=overwrite)
            result.combined_h5ad = str(combined_path)

        result.manifest_path = str(out_path / f"{study_accession}.json2h5ad.json")
        self._write_manifest(result, planned, json_path=json_path, overwrite=overwrite)
        return result

    def _study_accession(self, packages: list[dict]) -> str | None:
        for package in packages:
            series = package.get("series") if isinstance(package, dict) else None
            if not isinstance(series, dict):
                continue
            for accession in SourcePlanner()._as_list(series.get("accession")):
                value = SourcePlanner()._value(accession)
                if isinstance(value, str) and value.upper().startswith("GSE"):
                    return value.upper()
        return None

    def _sample_lookup(self, packages: list[dict]) -> dict[str, dict]:
        lookup = {}
        for package in packages:
            for sample in self.planner._as_list(package.get("sample")):
                if isinstance(sample, dict):
                    accession = self.planner.sample_accession(sample)
                    if accession:
                        lookup[accession] = sample
        return lookup

    def _scientific_modules(self):
        try:
            import anndata
            import numpy
            import pandas
            import scanpy
            from scipy import sparse
        except ImportError as exc:
            raise RuntimeError(
                "json2h5ad requires optional dependencies; install "
                "meta-standards-converter[h5ad]."
            ) from exc
        return anndata, numpy, pandas, scanpy, sparse

    def _read_processed_asset(self, asset: Asset, orientation: str = "auto"):
        anndata, numpy, pandas, scanpy, sparse = self._scientific_modules()
        path = self._local_path(asset.path)
        if asset.kind == "h5ad":
            adata = anndata.read_h5ad(path)
        elif self._underlying_suffix(path) == ".h5":
            adata = scanpy.read_10x_h5(path, gex_only=True)
        elif self._underlying_suffix(path) == ".mtx" or Path(path).is_dir():
            matrix_dir = path if Path(path).is_dir() else str(Path(path).parent)
            adata = scanpy.read_10x_mtx(matrix_dir, var_names="gene_ids", make_unique=True)
        else:
            separator = "," if self._underlying_suffix(path) == ".csv" else "\t"
            frame = pandas.read_csv(path, sep=separator, index_col=0)
            try:
                values = frame.apply(pandas.to_numeric, errors="raise")
            except (TypeError, ValueError) as exc:
                raise ValueError(f"Matrix {asset.path} contains nonnumeric values.") from exc
            raw = values.to_numpy()
            if raw.size == 0:
                raise ValueError(f"Matrix {asset.path} is empty.")
            if not numpy.isfinite(raw).all() or (raw < 0).any():
                raise ValueError(f"Matrix {asset.path} must contain finite nonnegative values.")
            if orientation == "auto":
                raise ValueError(
                    f"Matrix orientation for {asset.path} is ambiguous; specify "
                    "genes-by-observations or observations-by-genes."
                )
            if orientation == "genes-by-observations":
                raw = raw.T
                obs_names = values.columns.astype(str)
                var_names = values.index.astype(str)
            elif orientation == "observations-by-genes":
                obs_names = values.index.astype(str)
                var_names = values.columns.astype(str)
            else:
                raise ValueError(f"Unsupported matrix orientation: {orientation}")
            adata = anndata.AnnData(
                X=sparse.csr_matrix(raw),
                obs=pandas.DataFrame(index=obs_names),
                var=pandas.DataFrame(index=var_names),
            )
        if adata.n_obs == 0 or adata.n_vars == 0:
            raise ValueError(f"Processed asset {asset.path} contains an empty matrix.")
        if not sparse.issparse(adata.X):
            adata.X = sparse.csr_matrix(adata.X)
        else:
            adata.X = adata.X.tocsr()
        adata.var_names_make_unique()
        return adata

    def _normalize(self, adata, sample: dict, study_accession: str, asset: Asset) -> None:
        sample_id = self.planner.sample_accession(sample)
        metadata = self._sample_metadata(sample)
        original_names = [str(value) for value in adata.obs_names]
        adata.obs_names = [
            value if value.endswith(f"-{sample_id}") else f"{value}-{sample_id}"
            for value in original_names
        ]
        adata.obs_names_make_unique()
        annotations = {
            "geo_accession": sample_id,
            "geo_series_accession": study_accession,
            "geo_title": metadata.get("title"),
            "geo_organism": metadata.get("organism"),
            "geo_organism_taxid": metadata.get("organism_taxid"),
            "geo_organism_part": metadata.get("organism_part"),
            "geo_disease": metadata.get("disease"),
            "geo_genotype": metadata.get("genotype"),
            "geo_source_name": metadata.get("source"),
            "geo_source_tier": asset.kind,
            "geo_source_uri": asset.path,
        }
        for key, value in annotations.items():
            adata.obs[key] = "" if value is None else str(value)
        provenance = {
            "study_accession": study_accession,
            "sample_accession": sample_id,
            "source_tier": asset.kind,
            "source_uri": asset.path,
            "source_origin": asset.source,
            "source_sha256": self._sha256(asset.path),
            "converter_version": self._package_version(),
        }
        existing = adata.uns.get("meta_standards_converter")
        if isinstance(existing, dict):
            provenance = {**existing, **provenance}
        adata.uns["meta_standards_converter"] = provenance

    def _sample_metadata(self, sample: dict) -> dict:
        metadata = {"title": self._text(sample.get("title"))}
        channels = [x for x in self.planner._as_list(sample.get("channel")) if isinstance(x, dict)]
        channel = channels[0] if channels else {}
        metadata["source"] = self._text(channel.get("source"))
        organisms = self.planner._as_list(channel.get("organism"))
        organism = organisms[0] if organisms else None
        metadata["organism"] = self._text(organism)
        metadata["organism_taxid"] = organism.get("taxid") if isinstance(organism, dict) else None
        characteristics = {}
        for item in self.planner._as_list(channel.get("characteristics")):
            if isinstance(item, dict) and item.get("tag"):
                characteristics[str(item["tag"]).strip().lower()] = self._text(item)
        metadata["organism_part"] = characteristics.get("organism part") or characteristics.get("tissue")
        metadata["disease"] = characteristics.get("disease")
        metadata["genotype"] = characteristics.get("genotype")
        return metadata

    def _combine(self, adatas: dict[str, object]):
        if not adatas:
            raise ValueError("No sample H5ADs were produced.")
        anndata, _numpy, _pandas, _scanpy, sparse = self._scientific_modules()
        organisms = {
            str(adata.obs["geo_organism"].iloc[0]).strip()
            for adata in adatas.values()
            if "geo_organism" in adata.obs and str(adata.obs["geo_organism"].iloc[0]).strip()
        }
        if len(organisms) > 1:
            raise ValueError(f"Cannot combine samples with incompatible organisms: {sorted(organisms)}")
        namespaces = {self._feature_namespace(adata) for adata in adatas.values()}
        namespaces.discard("unknown")
        if len(namespaces) > 1:
            raise ValueError(f"Cannot combine incompatible feature identifier namespaces: {sorted(namespaces)}")
        combined = anndata.concat(
            adatas,
            axis="obs",
            join="outer",
            merge="first",
            label="geo_batch",
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
        }
        return combined

    def _feature_namespace(self, adata) -> str:
        values = adata.var.get("gene_ids", adata.var_names)
        values = [str(value).split(".")[0].upper() for value in list(values)[:100] if value]
        if values and sum(value.startswith(("ENSG", "ENSMUSG", "ENSRNOG")) for value in values) >= len(values) / 2:
            return "ensembl"
        if values and sum(value.isalnum() for value in values) >= len(values) / 2:
            return "symbol"
        return "unknown"

    def _write_h5ad(self, adata, path: Path, overwrite: bool) -> None:
        if path.exists() and not overwrite:
            raise FileExistsError(f"Output already exists: {path}")
        path.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(dir=path.parent, suffix=".h5ad", delete=False) as handle:
            temporary = Path(handle.name)
        try:
            adata.write_h5ad(temporary, compression="gzip")
            os.replace(temporary, path)
        finally:
            if temporary.exists():
                temporary.unlink()

    def _write_manifest(
        self,
        result: ConversionResult,
        planned: dict[str, Asset],
        json_path: str,
        overwrite: bool,
    ) -> None:
        path = Path(result.manifest_path)
        if path.exists() and not overwrite:
            raise FileExistsError(f"Output already exists: {path}")
        payload = {
            "study_accession": result.study_accession,
            "source_json": os.path.abspath(json_path),
            "combined_h5ad": result.combined_h5ad,
            "sample_h5ads": result.sample_h5ads,
            "retained_h5ads": result.retained_h5ads,
            "warnings": result.warnings,
            "failures": result.failures,
            "assets": {
                sample: {
                    "path": asset.path,
                    "kind": asset.kind,
                    "source": asset.source,
                    "role": asset.role,
                }
                for sample, asset in planned.items()
            },
        }
        temporary = path.with_suffix(path.suffix + ".tmp")
        with open(temporary, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
        os.replace(temporary, path)

    def _local_path(self, value: str) -> str:
        parsed = urlparse(value)
        if parsed.scheme in ("", "file"):
            return parsed.path if parsed.scheme == "file" else value
        raise NotImplementedError(f"Remote processed asset download is not implemented for {value}")

    def _underlying_suffix(self, path: str) -> str:
        name = Path(path).name.lower()
        for suffix in (".gz", ".bz2", ".xz", ".zip"):
            if name.endswith(suffix):
                name = name[: -len(suffix)]
                break
        return Path(name).suffix

    def _text(self, value) -> str | None:
        if isinstance(value, dict):
            value = value.get("value") or value.get("name")
        if value is None:
            return None
        return " ".join(str(value).split()) or None

    def _sha256(self, path: str) -> str | None:
        local = self._local_path(path)
        if not os.path.isfile(local):
            return None
        digest = hashlib.sha256()
        with open(local, "rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        return digest.hexdigest()

    def _package_version(self) -> str:
        try:
            return version("meta-standards-converter")
        except PackageNotFoundError:
            return "development"


class json2h5ad(JSON2H5ADConverter):
    """Backward-compatible converter name used by the existing console script."""
