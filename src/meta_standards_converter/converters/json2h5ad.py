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
import csv
import gzip
import hashlib
import logging
import os
import re
import shutil
import subprocess
import tempfile
import urllib.request
from dataclasses import dataclass, field, replace
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from urllib.parse import urlparse

import requests

from meta_standards_converter.harmonizers.harmonizers import Harmonizer


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
    study_scope: str | None = None
    reference: str | None = None
    annotation_source: str | None = None
    annotation_format: str | None = None
    annotation_sha256: str | None = None
    effective_annotation: str | None = None


class AssetManifest:
    """Load explicit asset mappings from CSV/TSV or compact CLI specifications."""

    def load(self, path: str) -> list[Asset]:
        delimiter = "\t" if Path(path).suffix.lower() in {".tsv", ".tab"} else ","
        with open(path, encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter=delimiter))
        if not rows or not {"scope_id", "path"}.issubset(rows[0]):
            raise ValueError("Asset manifest requires scope_id and path columns.")
        processed = []
        raw_groups = {}
        planner = SourcePlanner()
        for row in rows:
            scope_id = (row.get("scope_id") or "").strip()
            asset_path = (row.get("path") or "").strip()
            if not scope_id or not asset_path:
                raise ValueError("Asset manifest scope_id and path values cannot be blank.")
            kind = (row.get("kind") or planner.classify(asset_path) or "").strip().lower()
            if not kind and os.path.isdir(asset_path):
                kind = "matrix"
            if kind not in {"h5ad", "matrix", "raw"}:
                raise ValueError(f"Unsupported asset kind for {asset_path}: {kind or 'unknown'}")
            role = (row.get("role") or "primary").strip()
            if kind == "raw":
                key = (scope_id, role)
                raw_groups.setdefault(key, []).append({
                    "uri": asset_path,
                    "read": (row.get("read") or "").strip() or None,
                    "lane": (row.get("lane") or "").strip() or None,
                    "run": (row.get("run") or row.get("lane") or "").strip() or None,
                    "md5": (row.get("md5") or "").strip() or None,
                })
                continue
            processed.append(Asset(
                scope_id=scope_id,
                path=asset_path,
                kind=kind,
                role=role,
                source="manifest",
                features_path=(row.get("features_path") or "").strip() or None,
                barcodes_path=(row.get("barcodes_path") or "").strip() or None,
                orientation=(row.get("orientation") or "auto").strip(),
                md5=(row.get("md5") or "").strip() or None,
            ))
        for (scope_id, role), members in raw_groups.items():
            processed.append(Asset(
                scope_id=scope_id,
                path=members[0]["uri"],
                kind="raw",
                role=role,
                source="manifest",
                members=tuple(members),
            ))
        return processed

    def parse_spec(self, spec: str) -> Asset:
        if "=" not in spec:
            raise ValueError("--asset must use ACCESSION=PATH_OR_URL syntax.")
        scope_id, path = (value.strip() for value in spec.split("=", 1))
        if not scope_id or not path:
            raise ValueError("--asset accession and path cannot be blank.")
        kind = SourcePlanner().classify(path)
        if not kind and os.path.isdir(path):
            kind = "matrix"
        if not kind:
            raise ValueError(f"Cannot infer asset kind from {path}; use an asset manifest.")
        if kind == "raw":
            member = {"uri": path, "read": None, "run": None}
            return Asset(scope_id, path, kind, source="cli", members=(member,))
        return Asset(scope_id, path, kind, source="cli")


class AssetDownloader:
    """Stream remote processed assets into a deterministic local cache."""

    def __init__(self, cache_dir: str, session=None, urlopen=None):
        self.cache_dir = Path(cache_dir)
        self.session = session or requests.Session()
        self.urlopen = urlopen or urllib.request.urlopen

    def localize(self, value: str, md5: str | None = None) -> str:
        parsed = urlparse(value)
        if parsed.scheme in ("", "file"):
            return parsed.path if parsed.scheme == "file" else value
        if parsed.scheme not in {"http", "https", "ftp"}:
            raise ValueError(f"Unsupported asset URL scheme: {parsed.scheme}")
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        basename = os.path.basename(parsed.path) or "asset"
        prefix = hashlib.sha256(value.encode("utf-8")).hexdigest()[:12]
        destination = self.cache_dir / f"{prefix}-{basename}"
        if destination.exists():
            self._verify_md5(destination, md5)
            return str(destination)
        with tempfile.NamedTemporaryFile(dir=self.cache_dir, delete=False) as handle:
            temporary = Path(handle.name)
            try:
                if parsed.scheme in {"http", "https"}:
                    response = self.session.get(value, stream=True, timeout=30)
                    response.raise_for_status()
                    for chunk in response.iter_content(chunk_size=1024 * 1024):
                        if chunk:
                            handle.write(chunk)
                else:
                    with self.urlopen(value, timeout=30) as response:
                        shutil.copyfileobj(response, handle, length=1024 * 1024)
            except Exception:
                temporary.unlink(missing_ok=True)
                raise
        try:
            self._verify_md5(temporary, md5)
            os.replace(temporary, destination)
        finally:
            temporary.unlink(missing_ok=True)
        return str(destination)

    def _verify_md5(self, path: Path, expected: str | None) -> None:
        if not expected:
            return
        digest = hashlib.md5(usedforsecurity=False)
        with open(path, "rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        if digest.hexdigest().lower() != expected.lower():
            raise ValueError(f"MD5 checksum mismatch for {path}")


@dataclass
class PipelineRun:
    pipeline: str
    revision: str
    command: list[str]
    work_dir: str
    out_dir: str
    returncode: int | None = None
    log_path: str | None = None
    annotation_source: str | None = None
    annotation_format: str | None = None
    annotation_sha256: str | None = None
    effective_annotation: str | None = None


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


@dataclass
class RawProcessingResult:
    assets: dict[str, Asset]
    retained_h5ads: list[str] = field(default_factory=list)
    runs: list[PipelineRun] = field(default_factory=list)


class ReferenceResolver:
    """Resolve explicit or safely confirmed nf-core reference parameters."""

    GENOME_BY_TAXID = {"9606": "GRCh38", "10090": "GRCm39"}
    GENOME_BY_NAME = {"homo sapiens": "GRCh38", "mus musculus": "GRCm39"}

    def resolve(
        self,
        packages: list[dict],
        genome: str | None = None,
        fasta: str | None = None,
        gtf: str | None = None,
        gff: str | None = None,
        accept_inferred: bool = False,
    ) -> dict:
        if gtf and gff:
            raise ValueError("Reference annotation options gtf and gff are mutually exclusive.")
        annotation = gtf or gff
        if genome:
            return {
                "genome": genome,
                **({"gtf": gtf} if gtf else {}),
                **({"gff": gff} if gff else {}),
            }
        if fasta or annotation:
            if annotation and not fasta:
                raise ValueError("A custom annotation requires either genome or fasta.")
            if fasta and not annotation:
                raise ValueError("A custom fasta reference requires an annotation in gtf or gff format.")
            return {
                "fasta": fasta,
                **({"gtf": gtf} if gtf else {}),
                **({"gff": gff} if gff else {}),
            }

        taxids = set()
        names = set()
        planner = SourcePlanner()
        for package in packages:
            for sample in planner._as_list(package.get("sample")):
                if not isinstance(sample, dict):
                    continue
                for channel in planner._as_list(sample.get("channel")):
                    if not isinstance(channel, dict):
                        continue
                    for organism in planner._as_list(channel.get("organism")):
                        if isinstance(organism, dict):
                            if organism.get("taxid"):
                                taxids.add(str(organism["taxid"]))
                            if organism.get("value") or organism.get("name"):
                                names.add(str(organism.get("value") or organism.get("name")).lower())
                        elif organism:
                            names.add(str(organism).lower())
        candidates = {
            self.GENOME_BY_TAXID[value]
            for value in taxids
            if value in self.GENOME_BY_TAXID
        } | {
            self.GENOME_BY_NAME[value]
            for value in names
            if value in self.GENOME_BY_NAME
        }
        if len(candidates) != 1:
            raise ValueError(
                "Unable to infer one supported genome reference; provide --genome or --fasta and --gtf."
            )
        inferred = next(iter(candidates))
        if not accept_inferred:
            raise ValueError(
                f"Inferred reference {inferred}; rerun with --accept-inferred-reference "
                "or provide an explicit reference."
            )
        return {"genome": inferred, "inferred": True}


class AnnotationConverter:
    """Validate local annotations and normalize GFF3 input to GTF."""

    def __init__(self, command_runner=None, which=None):
        self.command_runner = command_runner or subprocess.run
        self.which = which or shutil.which

    def prepare(self, reference: dict, reference_dir: Path) -> dict:
        prepared = dict(reference)
        if prepared.get("fasta"):
            prepared["fasta"] = str(self._local_file(prepared["fasta"], "FASTA"))

        annotation = prepared.get("gtf") or prepared.get("gff")
        if not annotation:
            return prepared
        source = self._local_file(annotation, "annotation")
        annotation_format = self._annotation_format(source)
        digest = self._sha256(source)
        metadata = {
            "annotation_source": str(source),
            "annotation_format": annotation_format,
            "annotation_sha256": digest,
        }

        if annotation_format == "gtf":
            prepared["gtf"] = str(source)
            prepared["effective_annotation"] = str(source)
            return {**prepared, **metadata}

        reference_dir.mkdir(parents=True, exist_ok=True)
        destination = reference_dir / f"{digest}.gtf"
        temporary = Path(str(destination) + ".tmp")
        if not destination.exists() or destination.stat().st_size == 0:
            if not self.which("gffread"):
                raise RuntimeError("GFF3 annotations require gffread on PATH.")
            command = ["gffread", str(source), "-T", "-o", str(temporary)]
            completed = self.command_runner(
                command,
                text=True,
                capture_output=True,
                check=False,
            )
            if completed.returncode or not temporary.exists() or temporary.stat().st_size == 0:
                temporary.unlink(missing_ok=True)
                detail = (completed.stderr or completed.stdout or "conversion produced no GTF").strip()
                raise RuntimeError(f"gffread annotation conversion failed: {detail}")
            os.replace(temporary, destination)

        prepared.pop("gff", None)
        prepared["gtf"] = str(destination)
        prepared["effective_annotation"] = str(destination)
        return {**prepared, **metadata}

    def _local_file(self, value: str, label: str) -> Path:
        parsed = urlparse(str(value))
        if parsed.scheme not in ("", "file"):
            raise ValueError(f"User-supplied {label} must be a local file: {value}")
        path = Path(parsed.path if parsed.scheme == "file" else value).expanduser().resolve()
        if not path.is_file():
            raise FileNotFoundError(f"User-supplied {label} file not found: {path}")
        return path

    def _annotation_format(self, path: Path) -> str:
        name = path.name.lower()
        if name.endswith((".gtf", ".gtf.gz")):
            return "gtf"
        if name.endswith((".gff", ".gff.gz", ".gff3", ".gff3.gz")):
            return "gff3"
        raise ValueError(f"Unsupported annotation format: {path}")

    def _sha256(self, path: Path) -> str:
        digest = hashlib.sha256()
        with open(path, "rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        return digest.hexdigest()


class NFCoreRunner:
    """Prepare, execute, and inspect pinned nf-core RNA-seq workflows."""

    REVISIONS = {"scrnaseq": "4.2.0", "rnaseq": "3.26.0"}
    HTTPS_ARCHIVE_HOSTS = {"ftp.sra.ebi.ac.uk", "ftp.ncbi.nlm.nih.gov"}

    def __init__(
        self,
        command_runner=None,
        which=None,
        reference_resolver=None,
        runtime_runner=None,
        annotation_converter=None,
    ):
        self.command_runner = command_runner or subprocess.run
        self.which = which or shutil.which
        self.reference_resolver = reference_resolver or ReferenceResolver()
        self.runtime_runner = runtime_runner or subprocess.run
        self.annotation_converter = annotation_converter or AnnotationConverter(which=self.which)

    def process(
        self,
        assets: dict[str, Asset],
        packages: list[dict],
        out: str,
        study_accession: str,
        pipeline: str = "auto",
        genome: str | None = None,
        fasta: str | None = None,
        gtf: str | None = None,
        gff: str | None = None,
        accept_inferred_reference: bool = False,
        profile: str = "docker",
        revision: str | None = None,
        params_file: str | None = None,
        nextflow_config: str | None = None,
        work_dir: str | None = None,
        resume: bool = False,
    ) -> RawProcessingResult:
        if pipeline == "auto":
            groups = self._pipeline_groups(assets, packages)
            if len(groups) > 1:
                combined = RawProcessingResult(assets={})
                for selected, selected_assets in groups.items():
                    result = self.process(
                        selected_assets,
                        packages=packages,
                        out=out,
                        study_accession=study_accession,
                        pipeline=selected,
                        genome=genome,
                        fasta=fasta,
                        gtf=gtf,
                        gff=gff,
                        accept_inferred_reference=accept_inferred_reference,
                        profile=profile,
                        revision=revision,
                        params_file=params_file,
                        nextflow_config=nextflow_config,
                        work_dir=work_dir,
                        resume=resume,
                    )
                    combined.assets.update(result.assets)
                    combined.retained_h5ads.extend(result.retained_h5ads)
                    combined.runs.extend(result.runs)
                return combined
            pipeline = next(iter(groups), "rnaseq")
        if pipeline not in self.REVISIONS:
            raise ValueError(f"Unsupported nf-core pipeline: {pipeline}")
        reference = self.reference_resolver.resolve(
            packages,
            genome=genome,
            fasta=fasta,
            gtf=gtf,
            gff=gff,
            accept_inferred=accept_inferred_reference,
        )
        self._preflight(profile)
        revision = revision or self.REVISIONS[pipeline]
        reference = self.annotation_converter.prepare(
            reference,
            Path(out).resolve() / "nfcore" / study_accession / "reference",
        )
        run_dir = Path(out).resolve() / "nfcore" / study_accession / pipeline
        result_dir = run_dir / "results"
        run_dir.mkdir(parents=True, exist_ok=True)
        samplesheet = run_dir / "samplesheet.csv"
        self._write_samplesheet(samplesheet, assets, pipeline)
        params_path = run_dir / "params.json"
        params = {
            "input": str(samplesheet),
            "outdir": str(result_dir),
            **{
                key: reference[key]
                for key in ("genome", "fasta", "gtf")
                if reference.get(key)
            },
        }
        if params_file:
            with open(params_file, encoding="utf-8") as handle:
                supplied = json.load(handle)
            if not isinstance(supplied, dict):
                raise ValueError("nf-core params file must contain a JSON object.")
            params = {**supplied, **params}
        with open(params_path, "w", encoding="utf-8") as handle:
            json.dump(params, handle, indent=2, sort_keys=True)
            handle.write("\n")

        nextflow_work = Path(work_dir).resolve() if work_dir else run_dir / "work"
        command = [
            "nextflow",
            "run",
            f"nf-core/{pipeline}",
            "-r",
            revision,
            "-profile",
            profile,
            "-params-file",
            str(params_path),
            "-work-dir",
            str(nextflow_work),
        ]
        if nextflow_config:
            command.extend(["-c", nextflow_config])
        if resume:
            command.append("-resume")
        log_path = run_dir / "nextflow.log"
        completed = self.command_runner(
            command,
            cwd=run_dir,
            text=True,
            capture_output=True,
            check=False,
        )
        with open(log_path, "w", encoding="utf-8") as handle:
            handle.write(completed.stdout or "")
            handle.write(completed.stderr or "")
        run = PipelineRun(
            pipeline=pipeline,
            revision=revision,
            command=command,
            work_dir=str(nextflow_work),
            out_dir=str(result_dir),
            returncode=completed.returncode,
            log_path=str(log_path),
            annotation_source=reference.get("annotation_source"),
            annotation_format=reference.get("annotation_format"),
            annotation_sha256=reference.get("annotation_sha256"),
            effective_annotation=reference.get("effective_annotation"),
        )
        if completed.returncode:
            raise RuntimeError(f"Nextflow {pipeline} failed; see {log_path}")
        if pipeline == "scrnaseq":
            processed, retained = self._scrnaseq_assets(result_dir, assets, reference)
        else:
            processed, retained = self._rnaseq_assets(result_dir, assets, reference)
        return RawProcessingResult(assets=processed, retained_h5ads=retained, runs=[run])

    def _preflight(self, profile: str) -> None:
        missing = []
        if not self.which("nextflow"):
            missing.append("nextflow")
        if not self.which("java"):
            missing.append("java")
        runtime = next(
            (name for name in ("docker", "podman", "apptainer", "singularity") if name in profile.split(",")),
            None,
        )
        if runtime and not self.which(runtime):
            missing.append(runtime)
        if missing:
            raise RuntimeError(f"Missing nf-core runtime requirements: {', '.join(missing)}")
        require_rootless = os.environ.get(
            "META_STANDARDS_REQUIRE_ROOTLESS_DOCKER", ""
        ).strip().lower() in {"1", "true", "yes", "on"}
        if runtime == "docker" and require_rootless:
            completed = self.runtime_runner(
                ["docker", "info", "--format", "{{json .SecurityOptions}}"],
                text=True,
                capture_output=True,
                check=False,
            )
            if completed.returncode:
                detail = (completed.stderr or completed.stdout or "unknown error").strip()
                raise RuntimeError(
                    f"Cannot connect to the required rootless Docker daemon: {detail}"
                )
            if "rootless" not in (completed.stdout or "").lower():
                raise RuntimeError(
                    "json2h5ad requires a rootless Docker daemon in this deployment."
                )

    def _write_samplesheet(self, path: Path, assets: dict[str, Asset], pipeline: str) -> None:
        header = ["sample", "fastq_1", "fastq_2"]
        if pipeline == "rnaseq":
            header.append("strandedness")
        rows = []
        for sample_id, asset in assets.items():
            pairs = self._fastq_pairs(sample_id, asset, require_paired=pipeline == "scrnaseq")
            for first, second in pairs:
                row = [
                    sample_id,
                    self._workflow_uri(first),
                    self._workflow_uri(second) if second else "",
                ]
                if pipeline == "rnaseq":
                    row.append("auto")
                rows.append(row)
        with open(path, "w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, lineterminator="\n")
            writer.writerow(header)
            writer.writerows(rows)

    def _workflow_uri(self, value: str) -> str:
        parsed = urlparse(value)
        if parsed.scheme.lower() == "ftp" and (parsed.hostname or "").lower() in self.HTTPS_ARCHIVE_HOSTS:
            return parsed._replace(scheme="https").geturl()
        return value

    def _fastq_pairs(self, sample_id: str, asset: Asset, require_paired: bool) -> list[tuple[str, str | None]]:
        by_run = {}
        for index, member in enumerate(asset.members):
            path = member.get("uri") or member.get("filename")
            if not path:
                continue
            run = member.get("run") or "run1"
            read = member.get("read") or self._read_number(path)
            by_run.setdefault(run, []).append((read, path, index))
        pairs = []
        for run, entries in by_run.items():
            first = next((path for read, path, _ in entries if read == "1"), None)
            second = next((path for read, path, _ in entries if read == "2"), None)
            if not first and len(entries) in (1, 2):
                ordered = [path for _read, path, _index in sorted(entries, key=lambda item: item[2])]
                first = ordered[0]
                second = ordered[1] if len(ordered) == 2 else None
            if not first or (require_paired and not second):
                raise ValueError(f"Ambiguous FASTQ pairing for {sample_id} run {run}.")
            pairs.append((first, second))
        if not pairs:
            raise ValueError(f"No FASTQ files available for {sample_id}.")
        return pairs

    def _read_number(self, path: str) -> str | None:
        name = os.path.basename(urlparse(path).path)
        match = re.search(r"(?:^|[_\.])R?([12])(?:[_\.]|$)", name, re.IGNORECASE)
        return match.group(1) if match else None

    def _pipeline(self, packages: list[dict]) -> str:
        text = json.dumps(packages).lower()
        return "scrnaseq" if any(value in text for value in ("single cell", "single-cell", "10x", "chromium", "visium")) else "rnaseq"

    def _pipeline_groups(self, assets: dict[str, Asset], packages: list[dict]) -> dict[str, dict[str, Asset]]:
        planner = SourcePlanner()
        samples = {}
        for package in packages:
            for sample in planner._as_list(package.get("sample")):
                if isinstance(sample, dict):
                    sample_id = planner.sample_accession(sample)
                    if sample_id:
                        samples[sample_id] = sample
        groups = {}
        for sample_id, asset in assets.items():
            sample_text = json.dumps(samples.get(sample_id, {})).lower()
            selected = (
                "scrnaseq"
                if any(value in sample_text for value in ("single cell", "single-cell", "10x", "chromium", "visium"))
                else "rnaseq"
            )
            groups.setdefault(selected, {})[sample_id] = asset
        return groups

    def _scrnaseq_assets(self, result_dir: Path, inputs: dict[str, Asset], reference: dict):
        paths = sorted(result_dir.glob("**/*.h5ad"), key=lambda path: str(path).lower())
        processed = {}
        for sample_id in inputs:
            candidates = [
                path
                for path in paths
                if sample_id.lower() in str(path.relative_to(result_dir)).lower()
            ]
            if not candidates:
                raise RuntimeError(f"nf-core/scrnaseq produced no H5AD for {sample_id}.")
            chosen = max(candidates, key=self._scrnaseq_rank)
            processed[sample_id] = Asset(
                sample_id,
                str(chosen),
                "h5ad",
                source="nfcore",
                reference=reference.get("genome") or reference.get("fasta"),
                annotation_source=reference.get("annotation_source"),
                annotation_format=reference.get("annotation_format"),
                annotation_sha256=reference.get("annotation_sha256"),
                effective_annotation=reference.get("effective_annotation"),
            )
        return processed, [str(path) for path in paths]

    def _scrnaseq_rank(self, path: Path) -> int:
        location = str(path).lower()
        if "cellbender_filter" in location:
            return 3
        if "filtered" in location:
            return 2
        return 1

    def _rnaseq_assets(self, result_dir: Path, inputs: dict[str, Asset], reference: dict):
        counts = sorted(result_dir.glob("**/*.merged.gene_counts.tsv"))
        if not counts:
            raise RuntimeError("nf-core/rnaseq produced no merged gene-count matrix.")
        counts_path = counts[0]
        tpm_candidate = Path(str(counts_path).replace(".merged.gene_counts.tsv", ".merged.gene_tpm.tsv"))
        tpm_path = str(tpm_candidate) if tpm_candidate.exists() else None
        assets = {
            sample_id: Asset(
                sample_id,
                str(counts_path),
                "matrix",
                role="rnaseq_counts",
                source="nfcore",
                features_path=tpm_path,
                orientation="genes-by-observations",
                reference=reference.get("genome") or reference.get("fasta"),
                annotation_source=reference.get("annotation_source"),
                annotation_format=reference.get("annotation_format"),
                annotation_sha256=reference.get("annotation_sha256"),
                effective_annotation=reference.get("effective_annotation"),
            )
            for sample_id in inputs
        }
        return assets, []


class SourcePlanner:
    """Discover assets and select the best available source for every sample."""

    SOURCE_RANK = {"json": 0, "cli": 1, "manifest": 2, "nfcore": 3}
    KIND_RANK = {"raw": 0, "matrix": 1, "h5ad": 2}

    def plan(
        self,
        packages: list[dict],
        explicit_assets: list[Asset] | None = None,
        force_reprocess: bool = False,
    ) -> dict[str, Asset]:
        assets = self.discover(packages)
        assets.extend(explicit_assets or [])
        assets = self._coalesce_raw_assets(assets)
        samples = self.samples(packages)
        study_by_sample = self._study_by_sample(packages)
        planned = {}

        for sample_id in samples:
            study_id = study_by_sample.get(sample_id)
            candidates = []
            for asset in assets:
                if asset.scope_id == sample_id:
                    candidates.append(asset)
                elif study_id and asset.scope_id == study_id:
                    candidates.append(replace(asset, scope_id=sample_id, study_scope=study_id))
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

    def _coalesce_raw_assets(self, assets: list[Asset]) -> list[Asset]:
        retained = [asset for asset in assets if asset.kind != "raw"]
        groups = {}
        for asset in assets:
            if asset.kind != "raw":
                continue
            key = (asset.scope_id, asset.source, asset.role)
            members = list(asset.members) or [{"uri": asset.path, "md5": asset.md5}]
            groups.setdefault(key, []).extend(members)
        for (scope_id, source, role), members in groups.items():
            deduped = []
            seen = set()
            for member in members:
                path = member.get("uri") or member.get("filename")
                if not path or path in seen:
                    continue
                seen.add(path)
                deduped.append(member)
            retained.append(Asset(
                scope_id=scope_id,
                path=deduped[0].get("uri") or deduped[0].get("filename"),
                kind="raw",
                role=role,
                source=source,
                members=tuple(deduped),
            ))
        return retained

    def _study_by_sample(self, packages: list[dict]) -> dict[str, str]:
        result = {}
        for package in packages:
            series = package.get("series") if isinstance(package, dict) else None
            study_id = None
            if isinstance(series, dict):
                for accession in self._as_list(series.get("accession")):
                    value = self._value(accession)
                    if isinstance(value, str) and value.upper().startswith("GSE"):
                        study_id = value.upper()
                        break
            if not study_id:
                continue
            for sample in self._as_list(package.get("sample")):
                if isinstance(sample, dict):
                    sample_id = self.sample_accession(sample)
                    if sample_id:
                        result[sample_id] = study_id
        return result

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
                            assets.append(Asset(
                                sample_id,
                                path,
                                kind,
                                md5=entry.get("md5") if isinstance(entry, dict) else None,
                            ))

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
    MINIML_SCHEMA_VERSION = "1.0"
    PUBLICATION_POLICY = "citation_metadata_only"
    PUBLICATION_FIELDS = (
        "pubmed_id",
        "doi",
        "title",
        "author_list",
        "status",
        "status_term_source_ref",
        "status_term_accession_number",
    )
    PROTOCOL_PATHS = (
        ("treatment_protocol", "Treatment-Protocol", "channel"),
        ("growth_protocol", "Growth-Protocol", "channel"),
        ("extract_protocol", "Extract-Protocol", "channel"),
        ("label_protocol", "Label-Protocol", "channel"),
        ("hybridization_protocol", "Hybridization-Protocol", "sample"),
        ("scan_protocol", "Scan-Protocol", "sample"),
        ("data_processing", "Data-Processing", "sample"),
    )

    """Top-level JSON-to-H5AD conversion orchestrator."""

    def __init__(
        self,
        planner: SourcePlanner | None = None,
        pipeline_runner: NFCoreRunner | None = None,
        downloader: AssetDownloader | None = None,
    ):
        self.planner = planner or SourcePlanner()
        self.pipeline_runner = pipeline_runner or NFCoreRunner()
        self.downloader = downloader

    def convert(
        self,
        json_path: str,
        out: str | None = None,
        explicit_assets: list[Asset] | None = None,
        asset_manifest: str | None = None,
        asset_specs: list[str] | None = None,
        force_reprocess: bool = False,
        matrix_orientation: str = "auto",
        overwrite: bool = False,
        pipeline: str = "auto",
        genome: str | None = None,
        fasta: str | None = None,
        gtf: str | None = None,
        gff: str | None = None,
        accept_inferred_reference: bool = False,
        profile: str = "docker",
        revision: str | None = None,
        params_file: str | None = None,
        nextflow_config: str | None = None,
        work_dir: str | None = None,
        resume: bool = False,
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
        manifest_handler = AssetManifest()
        supplied_assets = list(explicit_assets or [])
        supplied_assets.extend(manifest_handler.parse_spec(spec) for spec in (asset_specs or []))
        if asset_manifest:
            supplied_assets.extend(manifest_handler.load(asset_manifest))
        if self.downloader is None:
            self.downloader = AssetDownloader(str(out_path / ".cache"))
        planned = self.planner.plan(
            packages,
            explicit_assets=supplied_assets,
            force_reprocess=force_reprocess,
        )
        sample_context = self._sample_context(packages)
        characteristic_columns = self._characteristic_columns(packages)
        source_json = os.path.abspath(json_path)
        source_json_sha256 = self._sha256(json_path)
        result = ConversionResult(study_accession=study_accession)
        adatas = {}

        raw_assets = {sample: asset for sample, asset in planned.items() if asset.kind == "raw"}
        if raw_assets:
            processed = self.pipeline_runner.process(
                raw_assets,
                packages=packages,
                out=str(out_path),
                study_accession=study_accession,
                pipeline=pipeline,
                genome=genome,
                fasta=fasta,
                gtf=gtf,
                gff=gff,
                accept_inferred_reference=accept_inferred_reference,
                profile=profile,
                revision=revision,
                params_file=params_file,
                nextflow_config=nextflow_config,
                work_dir=work_dir,
                resume=resume,
            )
            planned.update(processed.assets)
            result.retained_h5ads.extend(processed.retained_h5ads)
            result.pipeline_runs.extend(processed.runs)

        for sample_id, asset in planned.items():
            if asset.kind == "raw":
                raise RuntimeError(f"nf-core did not replace the raw source for {sample_id}.")
            adata = self._read_processed_asset(
                asset,
                orientation=(asset.orientation if asset.orientation != "auto" else matrix_orientation),
            )
            self._normalize(
                adata,
                sample=sample_context[sample_id][0],
                package=sample_context[sample_id][1],
                study_accession=study_accession,
                asset=asset,
                characteristic_columns=characteristic_columns,
            )
            self._attach_miniml(
                adata,
                packages=packages,
                source_json=source_json,
                source_json_sha256=source_json_sha256,
                sample_id=sample_id,
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
            self._attach_miniml(
                combined,
                packages=packages,
                source_json=source_json,
                source_json_sha256=source_json_sha256,
            )
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
        return {sample_id: context[0] for sample_id, context in self._sample_context(packages).items()}

    def _sample_context(self, packages: list[dict]) -> dict[str, tuple[dict, dict, int]]:
        lookup = {}
        for package_index, package in enumerate(packages):
            for sample in self.planner._as_list(package.get("sample")):
                if isinstance(sample, dict):
                    accession = self.planner.sample_accession(sample)
                    if accession:
                        lookup[accession] = (sample, package, package_index)
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
        path = self._local_path(asset.path, md5=asset.md5)
        if asset.kind == "h5ad":
            adata = self._read_h5ad(anndata, path)
            if asset.study_scope:
                accession_column = next(
                    (
                        column
                        for column in (
                            "msc_accession",
                            "geo_accession",
                            "sample_id",
                            "sample",
                            "gsm_accession",
                        )
                        if column in adata.obs
                    ),
                    None,
                )
                if not accession_column:
                    raise ValueError(
                        f"Study H5AD {asset.path} cannot be mapped to samples; "
                        "obs needs msc_accession, geo_accession, sample_id, sample, or gsm_accession."
                    )
                mask = adata.obs[accession_column].astype(str).str.upper() == asset.scope_id.upper()
                if not mask.any():
                    raise ValueError(f"Study H5AD {asset.path} contains no observations for {asset.scope_id}.")
                adata = adata[mask].copy()
        elif self._underlying_suffix(path) == ".h5":
            adata = scanpy.read_10x_h5(path, gex_only=True)
        elif self._underlying_suffix(path) == ".mtx" or Path(path).is_dir():
            matrix_dir = path if Path(path).is_dir() else str(Path(path).parent)
            adata = scanpy.read_10x_mtx(matrix_dir, var_names="gene_ids", make_unique=True)
        else:
            separator = "," if self._underlying_suffix(path) == ".csv" else "\t"
            frame = pandas.read_csv(path, sep=separator, index_col=0)
            feature_annotations = None
            if asset.role == "rnaseq_counts":
                if asset.scope_id not in frame.columns:
                    raise ValueError(
                        f"RNA-seq count matrix {asset.path} has no column for {asset.scope_id}."
                    )
                if "gene_name" in frame.columns:
                    feature_annotations = frame[["gene_name"]].copy()
                frame = frame[[asset.scope_id]]
            try:
                values = frame.apply(pandas.to_numeric, errors="raise")
            except (TypeError, ValueError) as exc:
                raise ValueError(f"Matrix {asset.path} contains nonnumeric values.") from exc
            raw = values.to_numpy()
            if raw.size == 0:
                raise ValueError(f"Matrix {asset.path} is empty.")
            if not numpy.isfinite(raw).all() or (raw < 0).any():
                raise ValueError(f"Matrix {asset.path} must contain finite nonnegative values.")
            if asset.study_scope and asset.role != "rnaseq_counts":
                if orientation == "genes-by-observations" and asset.scope_id in values.columns:
                    values = values[[asset.scope_id]]
                    raw = values.to_numpy()
                elif orientation == "observations-by-genes" and asset.scope_id in values.index:
                    values = values.loc[[asset.scope_id]]
                    raw = values.to_numpy()
                else:
                    raise ValueError(
                        f"Study matrix {asset.path} cannot be mapped to {asset.scope_id}."
                    )
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
            if feature_annotations is not None:
                adata.var["gene_name"] = (
                    feature_annotations.reindex(var_names)["gene_name"].astype(str).to_numpy()
                )
            if asset.role == "rnaseq_counts" and asset.features_path:
                tpm = pandas.read_csv(asset.features_path, sep="\t", index_col=0)
                if asset.scope_id not in tpm.columns:
                    raise ValueError(
                        f"RNA-seq TPM matrix {asset.features_path} has no column for {asset.scope_id}."
                    )
                tpm_values = tpm[[asset.scope_id]].reindex(values.index)
                tpm_values = tpm_values.apply(pandas.to_numeric, errors="raise")
                if tpm_values.isna().any().any():
                    raise ValueError("RNA-seq TPM features do not align with count features.")
                adata.layers["tpm"] = sparse.csr_matrix(tpm_values.to_numpy().T)
        if adata.n_obs == 0 or adata.n_vars == 0:
            raise ValueError(f"Processed asset {asset.path} contains an empty matrix.")
        if not sparse.issparse(adata.X):
            adata.X = sparse.csr_matrix(adata.X)
        else:
            adata.X = adata.X.tocsr()
        adata.var_names_make_unique()
        return adata

    def _read_h5ad(self, anndata, path: str):
        if not str(path).lower().endswith(".gz"):
            return anndata.read_h5ad(path)
        temporary_path = None
        try:
            with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as temporary:
                temporary_path = temporary.name
                with gzip.open(path, "rb") as compressed:
                    shutil.copyfileobj(compressed, temporary, length=1024 * 1024)
            return anndata.read_h5ad(temporary_path)
        finally:
            if temporary_path:
                Path(temporary_path).unlink(missing_ok=True)

    def _normalize(
        self,
        adata,
        sample: dict,
        package: dict,
        study_accession: str,
        asset: Asset,
        characteristic_columns: list[str],
    ) -> None:
        sample_id = self.planner.sample_accession(sample)
        metadata = self._sample_metadata(sample, package)
        modality = self._sample_modality(sample)
        original_names = [str(value) for value in adata.obs_names]
        adata.obs_names = [
            value if value.endswith(f"-{sample_id}") else f"{value}-{sample_id}"
            for value in original_names
        ]
        adata.obs_names_make_unique()
        annotations = {
            "msc_accession": sample_id,
            "msc_series_accession": study_accession,
            "msc_title": metadata.get("title"),
            "msc_description": metadata.get("description"),
            "msc_organism": metadata.get("organism"),
            "msc_organism_taxid": metadata.get("organism_taxid"),
            "msc_organism_part": metadata.get("organism_part"),
            "msc_developmental_stage": metadata.get("developmental_stage"),
            "msc_disease": metadata.get("disease"),
            "msc_genotype": metadata.get("genotype"),
            "msc_source_name": metadata.get("source"),
            "msc_biomaterial_provider": metadata.get("biomaterial_provider"),
            "msc_material_type": metadata.get("material_type"),
            "msc_molecule": metadata.get("molecule"),
            "msc_platform_accession": metadata.get("platform_accession"),
            "msc_sra_accession": metadata.get("sra_accession"),
            "msc_ena_accession": metadata.get("ena_accession"),
            "msc_biosample_accession": metadata.get("biosample_accession"),
            "msc_sra_run_accessions": metadata.get("sra_run_accessions"),
            "msc_library_strategy": metadata.get("library_strategy"),
            "msc_library_source": metadata.get("library_source"),
            "msc_library_selection": metadata.get("library_selection"),
            "msc_library_layout": metadata.get("library_layout"),
            "msc_instrument_model": metadata.get("instrument_model"),
            "msc_protocol_types": metadata.get("protocol_types"),
            "msc_protocol_term_source_refs": metadata.get("protocol_term_source_refs"),
            "msc_protocol_term_accession_numbers": metadata.get("protocol_term_accession_numbers"),
            "msc_metadata_source": metadata.get("metadata_source"),
            "msc_metadata_source_name": metadata.get("metadata_source_name"),
            "msc_metadata_source_uri": metadata.get("metadata_source_uri"),
            "msc_source_tier": asset.kind,
            "msc_source_uri": asset.path,
            "msc_modality": modality,
        }
        for column in characteristic_columns:
            annotations[f"msc_characteristic_{column}"] = metadata["characteristics"].get(column)
        for key, value in annotations.items():
            adata.obs[key] = "" if value is None else str(value)
        provenance = {
            "study_accession": study_accession,
            "sample_accession": sample_id,
            "source_tier": asset.kind,
            "source_uri": asset.path,
            "source_origin": asset.source,
            "source_sha256": self._sha256(asset.path, md5=asset.md5),
            "converter_version": self._package_version(),
            "modality": modality,
        }
        declared_reference = asset.reference or self._declared_reference(adata)
        if declared_reference:
            provenance["reference"] = declared_reference
        for key in (
            "annotation_source",
            "annotation_format",
            "annotation_sha256",
            "effective_annotation",
        ):
            value = getattr(asset, key)
            if value:
                provenance[key] = value
        existing = adata.uns.get("meta_standards_converter")
        if isinstance(existing, dict):
            provenance = {**existing, **provenance}
        adata.uns["meta_standards_converter"] = provenance

    def _sample_metadata(self, sample: dict, package: dict) -> dict:
        metadata = {
            "title": self._join_values(sample.get("title")),
            "description": self._join_values(sample.get("description")),
        }
        channels = [x for x in self.planner._as_list(sample.get("channel")) if isinstance(x, dict)]
        metadata["source"] = self._join_values(channel.get("source") for channel in channels)
        organisms = [
            organism
            for channel in channels
            for organism in self.planner._as_list(channel.get("organism"))
        ]
        metadata["organism"] = self._join_values(organisms)
        metadata["organism_taxid"] = self._join_values(
            organism.get("taxid") for organism in organisms if isinstance(organism, dict)
        )
        characteristic_values = {}
        for channel in channels:
            for item in self.planner._as_list(channel.get("characteristics")):
                if not isinstance(item, dict) or not item.get("tag"):
                    continue
                slug = self._metadata_slug(item.get("tag"))
                value = self._join_values(item.get("value"))
                if slug and value:
                    characteristic_values.setdefault(slug, []).append(value)
        characteristics = {
            slug: self._join_values(values)
            for slug, values in characteristic_values.items()
        }
        metadata["characteristics"] = characteristics
        metadata["organism_part"] = (
            characteristics.get("organism_part")
            or characteristics.get("tissue")
            or metadata["source"]
        )
        metadata["developmental_stage"] = characteristics.get("developmental_stage")
        metadata["disease"] = characteristics.get("disease")
        metadata["genotype"] = characteristics.get("genotype")

        metadata["biomaterial_provider"] = self._join_values(
            channel.get("biomaterial_provider") for channel in channels
        )
        metadata["molecule"] = self._join_values(channel.get("molecule") for channel in channels)
        material_types = []
        for value in self._values(channel.get("molecule") for channel in channels):
            material_types.append(re.sub(r"^total\s+", "", value, flags=re.IGNORECASE))
        metadata["material_type"] = self._join_values(material_types) or metadata["organism_part"]

        runs = [item for item in self.planner._as_list(sample.get("sra_run")) if isinstance(item, dict)]
        metadata["sra_accession"] = self._join_values(sample.get("sra_accession"))
        metadata["ena_accession"] = self._join_values(sample.get("ena_accession"))
        metadata["biosample_accession"] = self._join_values(run.get("biosample") for run in runs)
        metadata["sra_run_accessions"] = self._join_values(run.get("run") for run in runs)
        metadata["library_strategy"] = self._join_values(
            [sample.get("library_strategy"), *(run.get("library_strategy") for run in runs)]
        )
        metadata["library_source"] = self._join_values(
            [sample.get("library_source"), *(run.get("library_source") for run in runs)]
        )
        metadata["library_selection"] = self._join_values(
            [sample.get("library_selection"), *(run.get("library_selection") for run in runs)]
        )
        metadata["library_layout"] = self._join_values(run.get("library_layout") for run in runs)
        metadata["instrument_model"] = self._join_values(
            [sample.get("instrument_model"), *(run.get("instrument_model") for run in runs)]
        )
        metadata["platform_accession"] = self._platform_accessions(sample, package)

        protocol_types = []
        protocol_sources = []
        protocol_accessions = []
        for field, label, scope in self.PROTOCOL_PATHS:
            containers = channels if scope == "channel" else [sample]
            if not self._join_values(container.get(field) for container in containers):
                continue
            protocol_type, source_ref, accession = Harmonizer().geoprotocols2efo(label)
            protocol_types.append(protocol_type)
            protocol_sources.append(source_ref)
            protocol_accessions.append(accession)
        metadata["protocol_types"] = self._join_values(protocol_types)
        metadata["protocol_term_source_refs"] = self._join_values(protocol_sources)
        metadata["protocol_term_accession_numbers"] = self._join_values(protocol_accessions)

        database = self._metadata_database(package)
        metadata["metadata_source"] = self._join_values(
            database.get("public_id") or database.get("iid") or database.get("name")
        )
        metadata["metadata_source_name"] = self._join_values(database.get("name"))
        metadata["metadata_source_uri"] = self._join_values(database.get("web_link"))
        return metadata

    def _characteristic_columns(self, packages: list[dict]) -> list[str]:
        columns = []
        for package in packages:
            for sample in self.planner._as_list(package.get("sample")):
                if not isinstance(sample, dict):
                    continue
                for channel in self.planner._as_list(sample.get("channel")):
                    if not isinstance(channel, dict):
                        continue
                    for item in self.planner._as_list(channel.get("characteristics")):
                        if not isinstance(item, dict):
                            continue
                        slug = self._metadata_slug(item.get("tag"))
                        if slug and slug not in columns:
                            columns.append(slug)
        return columns

    def _metadata_slug(self, value) -> str:
        value = self._join_values(value).lower()
        return re.sub(r"_+", "_", re.sub(r"[^a-z0-9]+", "_", value)).strip("_")

    def _values(self, values) -> list[str]:
        flattened = []

        def visit(value):
            if value is None:
                return
            if isinstance(value, dict):
                for key in ("value", "name", "predefined", "public_id", "iid"):
                    if value.get(key) is not None:
                        visit(value[key])
                        return
                return
            if isinstance(value, (list, tuple, set)):
                for item in value:
                    visit(item)
                return
            if not isinstance(value, (str, bytes)) and hasattr(value, "__iter__"):
                for item in value:
                    visit(item)
                return
            cleaned = self._text(value)
            if cleaned and cleaned not in flattened:
                flattened.append(cleaned)

        visit(values)
        return flattened

    def _join_values(self, values) -> str:
        return "; ".join(self._values(values))

    def _metadata_database(self, package: dict) -> dict:
        return next(
            (
                item
                for item in self.planner._as_list(package.get("database"))
                if isinstance(item, dict)
            ),
            {},
        )

    def _platform_accessions(self, sample: dict, package: dict) -> str:
        references = set(self._values(sample.get("platform_ref")))
        values = []
        for platform in self.planner._as_list(package.get("platform")):
            if not isinstance(platform, dict):
                continue
            identifiers = {platform.get("iid"), *self._values(platform.get("accession"))}
            if references and not references.intersection(identifier for identifier in identifiers if identifier):
                continue
            values.extend(self._values(platform.get("accession")))
        return self._join_values(values or references)

    def _attach_miniml(
        self,
        adata,
        packages: list[dict],
        source_json: str,
        source_json_sha256: str | None,
        sample_id: str | None = None,
    ) -> None:
        _anndata, _numpy, pandas, _scanpy, _sparse = self._scientific_modules()
        rows = []
        for package_index, package in enumerate(packages):
            if sample_id is not None and not self._package_has_sample(package, sample_id):
                continue
            entities = self._metadata_entities(package, sample_id=sample_id)
            for entity_type, entity_id, entity in entities:
                safe_entity = self._publication_safe(entity)
                self._flatten_metadata(
                    safe_entity,
                    rows=rows,
                    package_index=package_index,
                    entity_type=entity_type,
                    entity_id=entity_id,
                )
        fields = pandas.DataFrame(
            rows,
            columns=(
                "package_index",
                "entity_type",
                "entity_id",
                "path",
                "value",
                "value_type",
            ),
        )
        fields.index = [f"field_{index:06d}" for index in range(len(fields))]
        database = next(
            (
                self._metadata_database(package)
                for package in packages
                if sample_id is None or self._package_has_sample(package, sample_id)
            ),
            {},
        )
        adata.uns["msc_miniml"] = {
            "schema_version": self.MINIML_SCHEMA_VERSION,
            "source_json": source_json,
            "source_sha256": source_json_sha256 or "",
            "publication_policy": self.PUBLICATION_POLICY,
            "metadata_source": self._join_values(
                database.get("public_id") or database.get("iid") or database.get("name")
            ),
            "metadata_source_name": self._join_values(database.get("name")),
            "metadata_source_uri": self._join_values(database.get("web_link")),
            "fields": fields,
        }

    def _package_has_sample(self, package: dict, sample_id: str) -> bool:
        return any(
            self.planner.sample_accession(sample) == sample_id
            for sample in self.planner._as_list(package.get("sample"))
            if isinstance(sample, dict)
        )

    def _metadata_entities(self, package: dict, sample_id: str | None):
        entity_groups = {
            "database": [
                item for item in self.planner._as_list(package.get("database")) if isinstance(item, dict)
            ],
            "contributor": [
                item for item in self.planner._as_list(package.get("contributor")) if isinstance(item, dict)
            ],
            "platform": [
                item for item in self.planner._as_list(package.get("platform")) if isinstance(item, dict)
            ],
            "sample": [
                item for item in self.planner._as_list(package.get("sample")) if isinstance(item, dict)
            ],
        }
        series = package.get("series") if isinstance(package.get("series"), dict) else None
        package_scalars = {
            key: value
            for key, value in package.items()
            if key not in {"database", "contributor", "platform", "sample", "series"}
        }
        entities = [("package", "package", package_scalars)]
        if series is not None:
            entities.append(("series", self._entity_id("series", series, 0), series))

        if sample_id is None:
            for entity_type in ("database", "contributor", "platform", "sample"):
                for index, entity in enumerate(entity_groups[entity_type]):
                    entities.append((entity_type, self._entity_id(entity_type, entity, index), entity))
            return entities

        selected_sample = next(
            (
                sample
                for sample in entity_groups["sample"]
                if self.planner.sample_accession(sample) == sample_id
            ),
            None,
        )
        if selected_sample is None:
            return entities
        selected = {"sample": [selected_sample], "platform": [], "contributor": [], "database": []}
        references = self._metadata_references(selected_sample)
        if series is not None:
            references.update(self._metadata_references(series))
        references.discard("")
        changed = True
        while changed:
            changed = False
            for entity_type in ("platform", "contributor", "database"):
                for index, entity in enumerate(entity_groups[entity_type]):
                    if entity in selected[entity_type]:
                        continue
                    identifiers = self._entity_identifiers(entity_type, entity, index)
                    if not references.intersection(identifiers):
                        continue
                    selected[entity_type].append(entity)
                    references.update(self._metadata_references(entity))
                    changed = True
        database = self._metadata_database(package)
        if database and database not in selected["database"]:
            selected["database"].append(database)
        for entity_type in ("database", "contributor", "platform", "sample"):
            for index, entity in enumerate(selected[entity_type]):
                entities.append((entity_type, self._entity_id(entity_type, entity, index), entity))
        return entities

    def _metadata_references(self, value) -> set[str]:
        references = set()

        def visit(current, key=None):
            if isinstance(current, dict):
                for child_key, child in current.items():
                    if child_key == "sample_ref":
                        continue
                    if child_key.endswith("_ref") or child_key == "database":
                        references.update(self._values(child))
                    else:
                        visit(child, child_key)
            elif isinstance(current, list):
                for child in current:
                    visit(child, key)

        visit(value)
        return references

    def _entity_identifiers(self, entity_type: str, entity: dict, index: int) -> set[str]:
        return {
            value
            for value in (
                self._entity_id(entity_type, entity, index),
                self._join_values(entity.get("iid")),
                *self._values(entity.get("accession")),
                self._join_values(entity.get("public_id")),
                self._join_values(entity.get("name")),
            )
            if value
        }

    def _entity_id(self, entity_type: str, entity: dict, index: int) -> str:
        if entity_type == "sample":
            return self.planner.sample_accession(entity) or self._join_values(entity.get("iid")) or f"sample_{index}"
        prefixes = {"series": "GSE", "platform": "GPL"}
        prefix = prefixes.get(entity_type)
        for value in self._values(entity.get("accession")):
            if not prefix or value.upper().startswith(prefix):
                return value
        return (
            self._join_values(entity.get("public_id"))
            or self._join_values(entity.get("iid"))
            or self._join_values(entity.get("name"))
            or f"{entity_type}_{index}"
        )

    def _publication_safe(self, value):
        if isinstance(value, list):
            return [self._publication_safe(item) for item in value]
        if not isinstance(value, dict):
            return value
        safe = {}
        for key, child in value.items():
            if key == "pubmed_publication":
                publications = []
                for publication in self.planner._as_list(child):
                    if not isinstance(publication, dict):
                        continue
                    publications.append(
                        {
                            field: self._publication_safe(publication[field])
                            for field in self.PUBLICATION_FIELDS
                            if field in publication
                        }
                    )
                safe[key] = publications
            else:
                safe[key] = self._publication_safe(child)
        return safe

    def _flatten_metadata(
        self,
        value,
        rows: list,
        package_index: int,
        entity_type: str,
        entity_id: str,
        path: str = "",
    ) -> None:
        if isinstance(value, dict):
            if not value:
                rows.append((package_index, entity_type, entity_id, path, "", "empty_object"))
                return
            for key, child in value.items():
                child_path = f"{path}.{key}" if path else str(key)
                self._flatten_metadata(
                    child,
                    rows,
                    package_index,
                    entity_type,
                    entity_id,
                    child_path,
                )
            return
        if isinstance(value, list):
            if not value:
                rows.append((package_index, entity_type, entity_id, path, "", "empty_list"))
                return
            for index, child in enumerate(value):
                self._flatten_metadata(
                    child,
                    rows,
                    package_index,
                    entity_type,
                    entity_id,
                    f"{path}[{index}]",
                )
            return
        if value is None:
            value_type, serialized = "null", ""
        elif isinstance(value, bool):
            value_type, serialized = "boolean", "true" if value else "false"
        elif isinstance(value, int):
            value_type, serialized = "integer", str(value)
        elif isinstance(value, float):
            value_type, serialized = "number", repr(value)
        else:
            value_type, serialized = "string", str(value)
        rows.append((package_index, entity_type, entity_id, path, serialized, value_type))

    def _sample_modality(self, sample: dict) -> str:
        text = json.dumps(sample).lower()
        if any(value in text for value in ("single cell", "single-cell", "10x", "chromium", "visium")):
            return "single_cell"
        if any(sample.get(key) for key in ("library_source", "library_strategy", "type", "sra_run")):
            return "bulk"
        return "unknown"

    def _combine(self, adatas: dict[str, object]):
        if not adatas:
            raise ValueError("No sample H5ADs were produced.")
        anndata, _numpy, _pandas, _scanpy, sparse = self._scientific_modules()
        organisms = {
            str(adata.obs["msc_organism"].iloc[0]).strip()
            for adata in adatas.values()
            if "msc_organism" in adata.obs and str(adata.obs["msc_organism"].iloc[0]).strip()
        }
        if len(organisms) > 1:
            raise ValueError(f"Cannot combine samples with incompatible organisms: {sorted(organisms)}")
        references = {
            str(adata.uns.get("meta_standards_converter", {}).get("reference")).strip()
            for adata in adatas.values()
            if isinstance(adata.uns.get("meta_standards_converter"), dict)
            and adata.uns["meta_standards_converter"].get("reference")
        }
        if len(references) > 1:
            raise ValueError(
                f"Cannot combine samples with incompatible reference builds: {sorted(references)}"
            )
        modalities = {
            str(adata.uns.get("meta_standards_converter", {}).get("modality")).strip()
            for adata in adatas.values()
            if isinstance(adata.uns.get("meta_standards_converter"), dict)
            and adata.uns["meta_standards_converter"].get("modality") not in (None, "", "unknown")
        }
        if len(modalities) > 1:
            raise ValueError(f"Cannot combine incompatible expression modalities: {sorted(modalities)}")
        namespaces = {self._feature_namespace(adata) for adata in adatas.values()}
        namespaces.discard("unknown")
        if len(namespaces) > 1:
            raise ValueError(f"Cannot combine incompatible feature identifier namespaces: {sorted(namespaces)}")
        combined = anndata.concat(
            adatas,
            axis="obs",
            join="outer",
            merge="first",
            label="msc_batch",
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
            "sample_provenance": {
                sample_id: dict(adata.uns.get("meta_standards_converter", {}))
                for sample_id, adata in adatas.items()
                if isinstance(adata.uns.get("meta_standards_converter"), dict)
            },
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

    def _declared_reference(self, adata) -> str | None:
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
            values = {str(value).strip() for value in adata.var["genome"] if str(value).strip()}
            if len(values) == 1:
                return next(iter(values))
        return None

    def _write_h5ad(self, adata, path: Path, overwrite: bool) -> None:
        if path.exists() and not overwrite:
            raise FileExistsError(f"Output already exists: {path}")
        path.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(dir=path.parent, suffix=".h5ad", delete=False) as handle:
            temporary = Path(handle.name)
        try:
            adata.write_h5ad(temporary, compression="gzip")
            os.replace(temporary, path)
            path.chmod(0o660)
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
            "pipeline_runs": [
                {
                    "pipeline": run.pipeline,
                    "revision": run.revision,
                    "command": run.command,
                    "work_dir": run.work_dir,
                    "out_dir": run.out_dir,
                    "returncode": run.returncode,
                    "log_path": run.log_path,
                    "annotation_source": run.annotation_source,
                    "annotation_format": run.annotation_format,
                    "annotation_sha256": run.annotation_sha256,
                    "effective_annotation": run.effective_annotation,
                }
                for run in result.pipeline_runs
            ],
            "warnings": result.warnings,
            "failures": result.failures,
            "assets": {
                sample: {
                    "path": asset.path,
                    "kind": asset.kind,
                    "source": asset.source,
                    "role": asset.role,
                    "reference": asset.reference,
                    "annotation_source": asset.annotation_source,
                    "annotation_format": asset.annotation_format,
                    "annotation_sha256": asset.annotation_sha256,
                    "effective_annotation": asset.effective_annotation,
                }
                for sample, asset in planned.items()
            },
        }
        temporary = path.with_suffix(path.suffix + ".tmp")
        with open(temporary, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
        os.replace(temporary, path)

    def _local_path(self, value: str, md5: str | None = None) -> str:
        if self.downloader is None:
            parsed = urlparse(value)
            if parsed.scheme in ("", "file"):
                return parsed.path if parsed.scheme == "file" else value
            raise RuntimeError("Remote asset downloader has not been configured.")
        return self.downloader.localize(value, md5=md5)

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

    def _sha256(self, path: str, md5: str | None = None) -> str | None:
        local = self._local_path(path, md5=md5)
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
