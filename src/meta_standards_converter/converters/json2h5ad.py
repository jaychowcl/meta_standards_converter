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
import logging
import os
from dataclasses import dataclass, field
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
        **options,
    ) -> ConversionResult:
        if not os.path.exists(json_path):
            raise FileNotFoundError(f"MINiML JSON file not found: {json_path}")
        with open(json_path, encoding="utf-8") as handle:
            packages = json.load(handle)
        if not isinstance(packages, list) or not packages:
            raise ValueError("Parsed MINiML JSON must contain a non-empty list of packages.")

        study_accession = self._study_accession(packages) or Path(json_path).stem
        raise NotImplementedError(
            f"json2h5ad conversion for {study_accession} requires the H5AD backend implementation."
        )

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


class json2h5ad(JSON2H5ADConverter):
    """Backward-compatible converter name used by the existing console script."""

