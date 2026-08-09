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
from dataclasses import dataclass, field, replace
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any, Mapping, Protocol, Sequence
from urllib.parse import urlparse
from uuid import uuid4

from meta_standards_converter.harmonizers.harmonizers import Harmonizer
from meta_standards_converter.retrieval import AssetDownloader, RetrievalPolicy
from meta_standards_converter.runtime_contracts import get_resource_profile
from .json_source import JSONPackageSource
from .mage_tab_projection import _parameter_rows, _parameter_summary


logger = logging.getLogger(__name__)

_TENX_MEMBER = re.compile(
    r"^(?P<prefix>.+?)[._](?P<role>"
    r"matrix\.mtx|barcodes\.tsv|genes\.tsv|features\.tsv)$",
    re.IGNORECASE,
)


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


@dataclass(frozen=True)
class MetadataProjectionContext:
    """Read-only conversion context supplied to metadata projectors."""

    sample: Mapping[str, Any]
    package: Mapping[str, Any]
    study_accession: str
    sample_accession: str
    asset: Asset
    base_metadata: Mapping[str, Any]


@dataclass(frozen=True)
class AnnDataMetadataProjection:
    """Metadata additions returned by an :class:`AnnDataMetadataProjector`."""

    obs: Mapping[str, Any] = field(default_factory=dict)
    var: Mapping[str, Any] = field(default_factory=dict)
    uns: Mapping[str, Any] = field(default_factory=dict)
    obs_renames: Mapping[str, str] = field(default_factory=dict)
    obs_drops: tuple[str, ...] = ()
    warnings: tuple[str, ...] = ()
    errors: tuple[str, ...] = ()


class AnnDataProjectionError(ValueError):
    """Raised when an AnnData projector reports invalid projected metadata."""

    def __init__(self, errors: Sequence[str]):
        self.errors = tuple(str(error) for error in errors)
        super().__init__("; ".join(self.errors))


class DatasetBundleRecoveryError(RuntimeError):
    """A dataset bundle failed to publish and could not be fully restored."""

    def __init__(
        self,
        publication_error: BaseException,
        recovery_errors: Sequence[BaseException],
        recovery_paths: Sequence[Path],
    ) -> None:
        self.publication_error = publication_error
        self.recovery_errors = tuple(recovery_errors)
        self.recovery_paths = tuple(recovery_paths)
        locations = ", ".join(str(path) for path in self.recovery_paths) or "none"
        super().__init__(
            "dataset bundle publication and recovery failed; preserved recovery "
            f"paths: {locations}; recovery failed: {self.recovery_errors[0]}"
        )


class AnnDataMetadataProjector(Protocol):
    """Optional extension that adds organization-neutral metadata to AnnData."""

    def project_sample(
        self,
        *,
        adata: Any,
        context: MetadataProjectionContext,
    ) -> AnnDataMetadataProjection:
        """Return metadata additions for one sample AnnData object."""

    def project_combined(
        self,
        *,
        adata: Any,
        contexts: Sequence[MetadataProjectionContext],
    ) -> AnnDataMetadataProjection:
        """Return metadata additions for the combined study AnnData object."""


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
    warnings: list[str] = field(default_factory=list)


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
    errors: list[str] = field(default_factory=list)

    @property
    def primary_h5ad(self) -> str | None:
        if self.combined_h5ad:
            return self.combined_h5ad
        return next(iter(self.sample_h5ads.values()), None)

    @property
    def partial(self) -> bool:
        return bool(self.failures or self.errors)

    def __str__(self) -> str:
        return self.primary_h5ad or self.manifest_path or self.study_accession

    def to_dict(self) -> dict[str, Any]:
        return {
            "operation": "h5ad",
            "status": "partial" if self.partial else "complete",
            "study_accession": self.study_accession,
            "artifacts": {
                "combined_h5ad": self.combined_h5ad,
                "sample_h5ads": dict(self.sample_h5ads),
                "retained_h5ads": list(self.retained_h5ads),
                "manifest": self.manifest_path,
            },
            "warnings": list(self.warnings),
            "errors": list(self.errors),
            "failures": list(self.failures),
        }


@dataclass
class BatchConversionResult:
    """Per-study conversions and diagnostics for one JSON source."""

    conversions: dict[str, ConversionResult] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)
    failures: list[str] = field(default_factory=list)

    @property
    def partial(self) -> bool:
        return bool(self.failures) or any(
            result.partial for result in self.conversions.values()
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "operation": "h5ad",
            "status": "partial" if self.partial else "complete",
            "datasets": {
                key: value.to_dict() for key, value in self.conversions.items()
            },
            "warnings": list(self.warnings),
            "failures": list(self.failures),
        }


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
            warnings=self._extract_warnings(completed.stdout, completed.stderr),
        )
        if completed.returncode:
            raise RuntimeError(f"Nextflow {pipeline} failed; see {log_path}")
        if pipeline == "scrnaseq":
            processed, retained = self._scrnaseq_assets(result_dir, assets, reference)
        else:
            processed, retained = self._rnaseq_assets(result_dir, assets, reference)
        return RawProcessingResult(assets=processed, retained_h5ads=retained, runs=[run])

    def _extract_warnings(self, *streams: str | None) -> list[str]:
        text = "\n".join(stream or "" for stream in streams)
        text = re.sub(r"\x1b\[[0-?]*[ -/]*[@-~]", "", text)
        lines = text.splitlines()
        warnings = []
        index = 0
        while index < len(lines):
            line = lines[index].strip()
            if not line.startswith("WARN:"):
                index += 1
                continue
            message = line.partition(":")[2].strip()
            if message and not set(message) <= {"~"}:
                if message not in warnings:
                    warnings.append(message)
                index += 1
                continue
            block = []
            index += 1
            while index < len(lines):
                continuation = lines[index].strip()
                if continuation and set(continuation) <= {"~"}:
                    index += 1
                    break
                if continuation.startswith("WARN:"):
                    break
                if continuation:
                    block.append(continuation)
                index += 1
            message = " ".join(block)
            if message and message not in warnings:
                warnings.append(message)
        return warnings

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
        assets = self._group_10x_assets(assets)
        assets = self._coalesce_raw_assets(assets)
        samples = self.samples(packages)
        study_by_sample = self._study_by_sample(packages)
        assets_by_scope = self._index_assets_by_scope(assets)
        planned = {}

        for sample_id in samples:
            study_id = study_by_sample.get(sample_id)
            candidates = list(assets_by_scope.get(sample_id, ()))
            if study_id:
                candidates.extend(
                    replace(asset, scope_id=sample_id, study_scope=study_id)
                    for asset in assets_by_scope.get(study_id, ())
                )
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

    @staticmethod
    def _index_assets_by_scope(
        assets: list[Asset],
    ) -> dict[str, list[Asset]]:
        by_scope: dict[str, list[Asset]] = {}
        for asset in assets:
            by_scope.setdefault(asset.scope_id, []).append(asset)
        return by_scope

    def _group_10x_assets(self, assets: list[Asset]) -> list[Asset]:
        groups: dict[
            tuple[str, str, str], dict[str, tuple[int, Asset]]
        ] = {}
        for index, asset in enumerate(assets):
            member = _tenx_member(asset.path)
            if member is None:
                continue
            prefix, role = member
            groups.setdefault(
                (asset.scope_id, prefix, asset.source), {}
            )[role] = (index, asset)
        complete = {
            key: members
            for key, members in groups.items()
            if {"matrix", "barcodes", "features"} <= set(members)
        }
        member_groups = {
            index: key
            for key, members in complete.items()
            for index, _asset in members.values()
        }
        emitted: set[tuple[str, str, str]] = set()
        result: list[Asset] = []
        for index, asset in enumerate(assets):
            key = member_groups.get(index)
            if key is None:
                result.append(asset)
                continue
            if key in emitted:
                continue
            emitted.add(key)
            members = complete[key]
            matrix = members["matrix"][1]
            result.append(
                replace(
                    matrix,
                    role="10x_mtx",
                    barcodes_path=members["barcodes"][1].path,
                    features_path=members["features"][1].path,
                )
            )
        return result

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
            series = package.get("series") if isinstance(package, Mapping) else None
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
    H5AD_METADATA_SCHEMA_VERSION = "1.0"
    PUBLICATION_POLICY = "citation_metadata_only"
    OBS_METADATA_FIELDS = {
        "msc.sample.title": "title",
        "msc.sample.description": "description",
        "msc.sample.channel.organism.value": "organism",
        "msc.sample.channel.organism.taxid": "organism_taxid",
        "msc.sample.channel.organism_part": "organism_part",
        "msc.sample.channel.developmental_stage": "developmental_stage",
        "msc.sample.channel.disease": "disease",
        "msc.sample.channel.genotype": "genotype",
        "msc.sample.channel.source": "source",
        "msc.sample.channel.biomaterial_provider": "biomaterial_provider",
        "msc.sample.channel.material_type": "material_type",
        "msc.sample.channel.molecule": "molecule",
        "msc.platform.accession": "platform_accession",
        "msc.archive.sra_accession": "sra_accession",
        "msc.archive.ena_accession": "ena_accession",
        "msc.archive.biosample_accession": "biosample_accession",
        "msc.archive.sra_run_accessions": "sra_run_accessions",
        "msc.library.strategy": "library_strategy",
        "msc.library.source": "library_source",
        "msc.library.selection": "library_selection",
        "msc.library.layout": "library_layout",
        "msc.instrument.model": "instrument_model",
        "msc.protocol.types": "protocol_types",
        "msc.protocol.term_source_refs": "protocol_term_source_refs",
        "msc.protocol.term_accession_numbers": "protocol_term_accession_numbers",
        "msc.database.identifier": "metadata_source",
        "msc.database.name": "metadata_source_name",
        "msc.database.uri": "metadata_source_uri",
    }
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
        metadata_projectors: Sequence[AnnDataMetadataProjector] | None = None,
        package_source: JSONPackageSource | None = None,
        retrieval_policy: RetrievalPolicy | None = None,
        resource_profile: str = "standard",
        resource_overrides: Mapping[str, int | float] | None = None,
    ):
        if downloader is not None and retrieval_policy is not None:
            raise ValueError(
                "downloader and retrieval_policy are mutually exclusive"
            )
        self.planner = planner or SourcePlanner()
        self.pipeline_runner = pipeline_runner or NFCoreRunner()
        self.downloader = downloader
        self.metadata_projectors = tuple(metadata_projectors or ())
        self.package_source = package_source or JSONPackageSource()
        self.retrieval_policy = retrieval_policy or RetrievalPolicy(
            resource_profile=get_resource_profile(
                resource_profile,
                overrides=resource_overrides,
            )
        )

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
        processed_checkpoint_dir: str | None = None,
        allow_invalid: bool = False,
        allow_unverified_combination: bool = False,
        use_harmonization_overrides: bool = False,
        **options,
    ) -> ConversionResult | BatchConversionResult:
        if not os.path.exists(json_path):
            raise FileNotFoundError(f"MINiML JSON file not found: {json_path}")
        loaded = self.package_source.load(json_path)
        loaded = replace(
            loaded,
            groups=tuple(
                group.resolved(enabled=use_harmonization_overrides)
                for group in loaded.groups
            ),
        )
        if not loaded.groups:
            raise ValueError("JSON source contains no convertible package groups.")
        for group in loaded.groups:
            self._validate_dataset_id(group.dataset_id)
        conversion_options = dict(
            explicit_assets=explicit_assets,
            asset_manifest=asset_manifest,
            asset_specs=asset_specs,
            force_reprocess=force_reprocess,
            matrix_orientation=matrix_orientation,
            overwrite=overwrite,
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
            processed_checkpoint_dir=processed_checkpoint_dir,
            allow_invalid=allow_invalid,
            allow_unverified_combination=allow_unverified_combination,
            **options,
        )
        if len(loaded.groups) > 1:
            return self._convert_groups(
                loaded,
                source_json=json_path,
                out=out,
                **conversion_options,
            )
        group = loaded.groups[0]
        result = self._convert_packages(
            list(group.packages),
            dataset_id=group.dataset_id,
            source_packages=list(group.source_packages or group.packages),
            harmonization_resolution=group.harmonization_resolution,
            source_json=json_path,
            out=out,
            **conversion_options,
        )
        for warning in loaded.warnings:
            if warning not in result.warnings:
                result.warnings.append(warning)
        return result

    def convert_source(
        self,
        json_path: str,
        out: str | None = None,
        allow_invalid: bool = False,
        use_harmonization_overrides: bool = False,
        **options,
    ) -> BatchConversionResult:
        if not os.path.exists(json_path):
            raise FileNotFoundError(f"JSON file not found: {json_path}")
        loaded = self.package_source.load(json_path)
        loaded = replace(
            loaded,
            groups=tuple(
                group.resolved(enabled=use_harmonization_overrides)
                for group in loaded.groups
            ),
        )
        if not loaded.groups:
            raise ValueError("JSON source contains no convertible package groups.")
        for group in loaded.groups:
            self._validate_dataset_id(group.dataset_id)
        return self._convert_groups(
            loaded,
            source_json=json_path,
            out=out,
            allow_invalid=allow_invalid,
            **options,
        )

    def _convert_groups(
        self,
        loaded,
        *,
        source_json: str,
        out: str | None,
        **options,
    ) -> BatchConversionResult:
        result = BatchConversionResult(warnings=list(loaded.warnings))
        root = Path(out or ".")
        multiple = len(loaded.groups) > 1
        for group in loaded.groups:
            self._validate_dataset_id(group.dataset_id)
            group_out = root / group.dataset_id if multiple else root
            try:
                converted = self._convert_packages(
                    list(group.packages),
                    dataset_id=group.dataset_id,
                    source_packages=list(group.source_packages or group.packages),
                    harmonization_resolution=group.harmonization_resolution,
                    source_json=source_json,
                    out=str(group_out),
                    **options,
                )
            except Exception as error:
                result.failures.append(f"{group.dataset_id}: {error}")
                continue
            result.conversions[group.dataset_id] = converted
        return result

    def _convert_packages(
        self,
        packages: list[dict],
        *,
        dataset_id: str | None = None,
        source_packages: list[dict] | None = None,
        harmonization_resolution=None,
        source_json: str,
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
        processed_checkpoint_dir: str | None = None,
        allow_invalid: bool = False,
        allow_unverified_combination: bool = False,
        **options,
    ) -> ConversionResult:

        study_accession = dataset_id or self._study_accession(packages) or Path(source_json).stem
        self._validate_path_component(study_accession, "study_accession")
        out_path = Path(out or ".")
        out_path.mkdir(parents=True, exist_ok=True)
        manifest_handler = AssetManifest()
        supplied_assets = list(explicit_assets or [])
        supplied_assets.extend(manifest_handler.parse_spec(spec) for spec in (asset_specs or []))
        if asset_manifest:
            supplied_assets.extend(manifest_handler.load(asset_manifest))
        if self.downloader is None:
            self.downloader = AssetDownloader(
                str(out_path / ".cache"),
                policy=self.retrieval_policy,
            )
        planned = self.planner.plan(
            packages,
            explicit_assets=supplied_assets,
            force_reprocess=force_reprocess,
        )
        for sample_id in planned:
            self._validate_path_component(sample_id, "sample_id")
        sample_context = self._sample_context(packages)
        characteristic_columns = self._characteristic_columns(packages)
        source_json = os.path.abspath(source_json)
        source_json_sha256 = self._sha256(source_json)
        checkpoint_root = (
            Path(processed_checkpoint_dir) / study_accession
            if processed_checkpoint_dir
            else None
        )
        result = ConversionResult(study_accession=study_accession)
        result.warnings.extend(getattr(harmonization_resolution, "warnings", ()))
        adatas = {}
        combined_adata = None
        projection_contexts: list[MetadataProjectionContext] = []

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
            for run in processed.runs:
                for warning in run.warnings:
                    rendered = f"{run.pipeline}: {warning}"
                    if rendered not in result.warnings:
                        result.warnings.append(rendered)

        for sample_id, asset in planned.items():
            if asset.kind == "raw":
                raise RuntimeError(f"nf-core did not replace the raw source for {sample_id}.")
            orientation = asset.orientation if asset.orientation != "auto" else matrix_orientation
            checkpoint = self._processed_checkpoint(
                checkpoint_root,
                sample_id=sample_id,
                source_json_sha256=source_json_sha256,
                sample=sample_context[sample_id][0],
                asset=asset,
                orientation=orientation,
            )
            restored = self._load_processed_checkpoint(checkpoint) if resume else None
            if restored is not None:
                adata, checkpoint_metadata = restored
                result.warnings.extend(checkpoint_metadata.get("warnings", ()))
                result.errors.extend(checkpoint_metadata.get("errors", ()))
                base_metadata = {}
            else:
                adata = self._read_processed_asset(asset, orientation=orientation)
                base_metadata = self._normalize(
                    adata,
                    sample=sample_context[sample_id][0],
                    package=sample_context[sample_id][1],
                    study_accession=study_accession,
                    asset=asset,
                    characteristic_columns=characteristic_columns,
                    artifact_parent=out_path,
                    harmonization_resolution=harmonization_resolution,
                )
            projection_context = MetadataProjectionContext(
                sample=sample_context[sample_id][0],
                package=sample_context[sample_id][1],
                study_accession=study_accession,
                sample_accession=sample_id,
                asset=asset,
                base_metadata=base_metadata,
            )
            if restored is None:
                warning_start = len(result.warnings)
                error_start = len(result.errors)
                self._project_sample_metadata(
                    adata,
                    projection_context,
                    warnings=result.warnings,
                    errors=result.errors,
                    allow_invalid=allow_invalid,
                )
                self._attach_miniml(
                    adata,
                    packages=source_packages or packages,
                    source_json=source_json,
                    source_json_sha256=source_json_sha256,
                    sample_id=sample_id,
                    artifact_parent=out_path,
                )
                self._attach_harmonization(adata, harmonization_resolution)
                self._write_processed_checkpoint(
                    checkpoint,
                    adata,
                    warnings=result.warnings[warning_start:],
                    errors=result.errors[error_start:],
                )
            projection_contexts.append(projection_context)
            sample_path = out_path / f"{sample_id}.h5ad"
            result.sample_h5ads[sample_id] = str(sample_path)
            adatas[sample_id] = adata

        self._ensure_global_observation_ids(adatas)
        try:
            combined = self._combine(
                adatas,
                allow_unverified=allow_unverified_combination,
            )
        except ValueError as exc:
            result.failures.append(str(exc))
        else:
            combined_adata = combined
            compatibility = combined.uns.get(
                "meta_standards_converter", {}
            ).get("combination_compatibility", {})
            if compatibility.get("explicitly_acknowledged"):
                missing = compatibility.get("missing_evidence", {})
                rendered = ", ".join(
                    f"{dimension} ({', '.join(sample_ids)})"
                    for dimension, sample_ids in sorted(missing.items())
                )
                result.failures.append(
                    "Combined samples with unverified compatibility evidence "
                    f"after it was explicitly acknowledged: {rendered}"
                )
            self._project_combined_metadata(
                combined,
                projection_contexts,
                warnings=result.warnings,
                errors=result.errors,
                allow_invalid=allow_invalid,
            )
            self._attach_miniml(
                combined,
                packages=source_packages or packages,
                source_json=source_json,
                source_json_sha256=source_json_sha256,
                artifact_parent=out_path,
            )
            self._attach_harmonization(combined, harmonization_resolution)
            combined_path = out_path / f"{study_accession}.h5ad"
            result.combined_h5ad = str(combined_path)

        result.manifest_path = str(out_path / f"{study_accession}.json2h5ad.json")
        self._write_dataset_bundle(
            adatas=adatas,
            combined=combined_adata,
            result=result,
            planned=planned,
            json_path=source_json,
            overwrite=overwrite,
        )
        return result

    def _processed_checkpoint(
        self,
        root: Path | None,
        *,
        sample_id: str,
        source_json_sha256: str,
        sample: Mapping[str, Any],
        asset: Asset,
        orientation: str,
    ) -> tuple[Path, Path, str] | None:
        if root is None:
            return None
        payload = {
            "schema_version": "1",
            "converter_version": self._package_version(),
            "source_json_sha256": source_json_sha256,
            "sample_id": sample_id,
            "sample": sample,
            "asset": vars(asset),
            "orientation": orientation,
        }
        fingerprint = hashlib.sha256(
            json.dumps(payload, ensure_ascii=False, sort_keys=True).encode("utf-8")
        ).hexdigest()
        key = hashlib.sha256(sample_id.encode("utf-8")).hexdigest()[:20]
        return root / f"{key}.h5ad", root / f"{key}.json", fingerprint

    def _load_processed_checkpoint(self, checkpoint):
        if checkpoint is None:
            return None
        h5ad_path, manifest_path, fingerprint = checkpoint
        if not h5ad_path.is_file() or not manifest_path.is_file():
            return None
        try:
            metadata = json.loads(manifest_path.read_text(encoding="utf-8"))
            if metadata.get("fingerprint") != fingerprint:
                return None
            anndata = self._scientific_modules()[0]
            return anndata.read_h5ad(h5ad_path), metadata
        except (OSError, ValueError, json.JSONDecodeError):
            return None

    def _write_processed_checkpoint(
        self,
        checkpoint,
        adata,
        *,
        warnings: Sequence[str],
        errors: Sequence[str],
    ) -> None:
        if checkpoint is None:
            return
        h5ad_path, manifest_path, fingerprint = checkpoint
        h5ad_path.parent.mkdir(parents=True, exist_ok=True)
        self._write_h5ad(adata, h5ad_path, overwrite=True)
        with tempfile.NamedTemporaryFile(
            dir=manifest_path.parent,
            suffix=".json",
            mode="w",
            encoding="utf-8",
            delete=False,
        ) as handle:
            temporary = Path(handle.name)
            json.dump(
                {
                    "schema_version": "1",
                    "fingerprint": fingerprint,
                    "warnings": list(warnings),
                    "errors": list(errors),
                },
                handle,
                ensure_ascii=False,
                sort_keys=True,
            )
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        try:
            os.replace(temporary, manifest_path)
        finally:
            temporary.unlink(missing_ok=True)

    @staticmethod
    def _validate_dataset_id(dataset_id: str) -> None:
        JSON2H5ADConverter._validate_path_component(dataset_id, "dataset_id")

    @staticmethod
    def _validate_path_component(value: Any, name: str) -> str:
        rendered = str(value)
        if (
            not rendered
            or rendered in {".", ".."}
            or Path(rendered).name != rendered
            or "/" in rendered
            or "\\" in rendered
        ):
            raise ValueError(f"Unsafe {name} path component: {rendered!r}")
        return rendered

    def _study_accession(self, packages: list[dict]) -> str | None:
        for package in packages:
            series = package.get("series") if isinstance(package, Mapping) else None
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
            from scipy import sparse
        except ImportError as exc:
            raise RuntimeError(
                "json2h5ad requires optional dependencies; install "
                "meta-standards-converter[h5ad]."
            ) from exc
        return anndata, numpy, pandas, sparse

    def _scanpy_module(self):
        try:
            import scanpy
        except ImportError as exc:
            raise RuntimeError(
                "10x matrix input requires Scanpy; install "
                "meta-standards-converter[h5ad]."
            ) from exc
        return scanpy

    def _read_processed_asset(self, asset: Asset, orientation: str = "auto"):
        if (
            self._underlying_suffix(asset.path) == ".mtx"
            and asset.features_path
            and asset.barcodes_path
        ):
            with tempfile.TemporaryDirectory(prefix="msc-10x-") as directory:
                prepared = Path(directory)
                paths = {
                    "matrix": self._local_path(asset.path, md5=asset.md5),
                    "barcodes": self._local_path(asset.barcodes_path),
                    "features": self._local_path(asset.features_path),
                }
                legacy = ".genes.tsv" in os.path.basename(
                    urlparse(str(asset.features_path)).path
                ).casefold()
                for role, source in paths.items():
                    destination = prepared / _tenx_local_name(
                        role,
                        asset.features_path if role == "features" else source,
                        legacy=legacy,
                    )
                    _copy_10x_member(source, destination, decompress=legacy)
                return self._read_processed_asset(
                    replace(
                        asset,
                        path=str(prepared),
                        features_path=None,
                        barcodes_path=None,
                    ),
                    orientation=orientation,
                )
        anndata, numpy, pandas, sparse = self._scientific_modules()
        path = self._local_path(asset.path, md5=asset.md5)
        if asset.kind == "h5ad":
            adata = self._read_h5ad(anndata, path)
            if asset.study_scope:
                accession_column = next(
                    (
                        column
                        for column in (
                            "msc.sample.accession",
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
                        "obs needs msc.sample.accession, geo_accession, "
                        "sample_id, sample, or gsm_accession."
                    )
                mask = adata.obs[accession_column].astype(str).str.upper() == asset.scope_id.upper()
                if not mask.any():
                    raise ValueError(f"Study H5AD {asset.path} contains no observations for {asset.scope_id}.")
                adata = adata[mask].copy()
        elif self._underlying_suffix(path) == ".h5":
            scanpy = self._scanpy_module()
            adata = scanpy.read_10x_h5(path, gex_only=True)
        elif self._underlying_suffix(path) == ".mtx" or Path(path).is_dir():
            scanpy = self._scanpy_module()
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
        artifact_parent: Path,
        harmonization_resolution=None,
    ) -> dict:
        sample_id = self.planner.sample_accession(sample)
        metadata_values = self._sample_metadata_values(sample, package)
        metadata = self._render_sample_metadata(metadata_values)
        modality = self._sample_modality(sample)
        original_names = [str(value) for value in adata.obs_names]
        adata.obs["msc.observation.original_id"] = original_names
        candidates = [
            value
            if self._is_sample_qualified(value, sample_id)
            else f"{value}-{sample_id}"
            for value in original_names
        ]
        adata.obs_names = self._deduplicate_observation_ids(candidates)
        source_uri, source_uri_scope = self._portable_location(asset.path, artifact_parent)
        canonical_values = {
            "msc.sample.accession": (sample_id,),
            "msc.series.accession": (study_accession,),
            **{
                column: metadata_values.get(key, ())
                for column, key in self.OBS_METADATA_FIELDS.items()
            },
            "msc.asset.tier": (asset.kind,),
            "msc.asset.uri": (source_uri,) if source_uri else (),
            "msc.asset.uri_scope": (source_uri_scope,) if source_uri_scope else (),
            "msc.expression.modality": (modality,),
        }
        canonical_values.update(_parameter_summary(package, sample))
        for column in characteristic_columns:
            canonical_values[f"msc.characteristics.{column}"] = metadata_values[
                "characteristics"
            ].get(column, ())
        for item in getattr(harmonization_resolution, "selections", ()):
            if item.sample_accession != sample_id:
                continue
            prefix = f"msc.harmonization.{item.destination}"
            canonical_values.update({
                f"{prefix}.value": (item.value,),
                f"{prefix}.id": (item.identifier,) if item.identifier else (),
                f"{prefix}.ontology": (item.ontology,) if item.ontology else (),
                f"{prefix}.source_field": (item.source_field,),
                f"{prefix}.hierarchy_depth": (
                    (item.hierarchy_depth,) if item.hierarchy_depth is not None else ()
                ),
            })
        for key, values in canonical_values.items():
            adata.obs[key] = self._join_values(values)
        self._attach_sample_values(
            adata,
            {sample_id: canonical_values},
        )
        provenance = {
            "study_accession": study_accession,
            "sample_accession": sample_id,
            "source_tier": asset.kind,
            "source_uri": source_uri,
            "source_uri_scope": source_uri_scope,
            "path_base": "artifact_parent",
            "source_origin": asset.source,
            "source_sha256": self._sha256(asset.path, md5=asset.md5),
            "converter_version": self._package_version(),
            "metadata_schema_version": self.H5AD_METADATA_SCHEMA_VERSION,
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
                if key in {"annotation_source", "effective_annotation"}:
                    portable, scope = self._portable_location(value, artifact_parent)
                    provenance[key] = portable
                    provenance[f"{key}_scope"] = scope
                else:
                    provenance[key] = value
        existing = adata.uns.get("meta_standards_converter")
        if isinstance(existing, dict):
            provenance = {**existing, **provenance}
        adata.uns["meta_standards_converter"] = provenance
        return metadata

    def _project_sample_metadata(
        self,
        adata,
        context: MetadataProjectionContext,
        *,
        warnings: list[str],
        errors: list[str],
        allow_invalid: bool,
    ) -> None:
        for projector in self.metadata_projectors:
            callback = getattr(projector, "project_sample", None)
            if callback is None:
                continue
            projection = callback(adata=adata, context=context)
            self._apply_metadata_projection(adata, projection)
            self._extend_warnings(warnings, projection.warnings)
            self._extend_warnings(errors, projection.errors)
            if projection.errors and not allow_invalid:
                raise AnnDataProjectionError(projection.errors)

    def _project_combined_metadata(
        self,
        adata,
        contexts: Sequence[MetadataProjectionContext],
        *,
        warnings: list[str],
        errors: list[str],
        allow_invalid: bool,
    ) -> None:
        for projector in self.metadata_projectors:
            callback = getattr(projector, "project_combined", None)
            if callback is None:
                continue
            projection = callback(adata=adata, contexts=tuple(contexts))
            self._apply_metadata_projection(adata, projection)
            self._extend_warnings(warnings, projection.warnings)
            self._extend_warnings(errors, projection.errors)
            if projection.errors and not allow_invalid:
                raise AnnDataProjectionError(projection.errors)

    def _apply_metadata_projection(
        self,
        adata,
        projection: AnnDataMetadataProjection,
    ) -> None:
        if not isinstance(projection, AnnDataMetadataProjection):
            raise TypeError(
                "metadata projector must return AnnDataMetadataProjection"
            )
        self._apply_obs_transforms(
            adata.obs,
            renames=projection.obs_renames,
            drops=projection.obs_drops,
        )
        self._apply_axis_projection(
            adata.obs,
            projection.obs,
            axis_name="obs",
            axis_length=adata.n_obs,
        )
        self._apply_axis_projection(
            adata.var,
            projection.var,
            axis_name="var",
            axis_length=adata.n_vars,
        )
        for key, value in projection.uns.items():
            if key in adata.uns:
                raise ValueError(f"uns metadata key {key!r} already exists")
            adata.uns[key] = value

    @staticmethod
    def _apply_obs_transforms(
        frame,
        *,
        renames: Mapping[str, str],
        drops: Sequence[str],
    ) -> None:
        renames = dict(renames)
        drops = tuple(drops)
        drop_set = set(drops)
        if len(drop_set) != len(drops):
            raise ValueError("obs drop columns must be unique")
        overlap = sorted(set(renames) & drop_set)
        if overlap:
            raise ValueError(
                f"obs metadata key {overlap[0]!r} cannot be renamed and dropped"
            )
        for source, target in renames.items():
            if not isinstance(source, str) or not isinstance(target, str):
                raise TypeError("obs rename sources and targets must be strings")
            if source not in frame:
                raise ValueError(f"obs rename source {source!r} does not exist")
            if not target:
                raise ValueError(f"obs rename target for {source!r} must be nonblank")
        missing_drops = [column for column in drops if column not in frame]
        if missing_drops:
            raise ValueError(f"obs drop source {missing_drops[0]!r} does not exist")
        targets = list(renames.values())
        if len(set(targets)) != len(targets):
            raise ValueError("obs rename targets must be unique")
        survivors = set(frame.columns) - set(renames) - drop_set
        collisions = sorted(set(targets) & survivors)
        if collisions:
            raise ValueError(
                f"obs rename target {collisions[0]!r} already exists"
            )
        if drops:
            frame.drop(columns=list(drops), inplace=True)
        if renames:
            frame.rename(columns=renames, inplace=True)

    @staticmethod
    def _apply_axis_projection(
        frame,
        values: Mapping[str, Any],
        *,
        axis_name: str,
        axis_length: int,
    ) -> None:
        for key, value in values.items():
            if key in frame:
                raise ValueError(f"{axis_name} metadata key {key!r} already exists")
            if isinstance(value, Sequence) and not isinstance(
                value, (str, bytes, bytearray)
            ):
                if len(value) != axis_length:
                    raise ValueError(
                        f"{axis_name} metadata key {key!r} must contain "
                        f"{axis_length} values"
                    )
                frame[key] = list(value)
            else:
                frame[key] = "" if value is None else value

    @staticmethod
    def _extend_warnings(target: list[str], values: Sequence[str]) -> None:
        for value in values:
            rendered = str(value)
            if rendered and rendered not in target:
                target.append(rendered)

    def _sample_metadata(self, sample: dict, package: dict) -> dict:
        return self._render_sample_metadata(
            self._sample_metadata_values(sample, package)
        )

    def _sample_metadata_values(self, sample: dict, package: dict) -> dict:
        metadata = {
            "title": tuple(self._values(sample.get("title"))),
            "description": tuple(self._values(sample.get("description"))),
        }
        channels = [x for x in self.planner._as_list(sample.get("channel")) if isinstance(x, dict)]
        metadata["source"] = tuple(
            self._values(channel.get("source") for channel in channels)
        )
        organisms = [
            organism
            for channel in channels
            for organism in self.planner._as_list(channel.get("organism"))
        ]
        organism_values = []
        for channel in channels:
            harmonized = []
            for organism in self.planner._as_list(channel.get("organism")):
                if not isinstance(organism, dict):
                    continue
                for annotation in self.planner._as_list(organism.get("annotations")):
                    if (
                        isinstance(annotation, dict)
                        and annotation.get("field") in {"organism", "species_name"}
                    ):
                        harmonized.extend(self._values(annotation.get("value")))
            organism_values.extend(harmonized or self._values(channel.get("organism")))
        metadata["organism"] = tuple(self._values(organism_values))
        metadata["organism_taxid"] = tuple(
            self._values(
                organism.get("taxid")
                for organism in organisms
                if isinstance(organism, dict)
            )
        )
        characteristic_values = {}
        for channel in channels:
            for annotation in self.planner._as_list(channel.get("annotations")):
                if not isinstance(annotation, dict) or not annotation.get("field"):
                    continue
                annotation_slug = "harmonized_" + self._metadata_slug(annotation["field"])
                characteristic_values.setdefault(annotation_slug, []).extend(
                    self._values(annotation.get("value"))
                )
                if annotation.get("term_accession_number"):
                    characteristic_values.setdefault(f"{annotation_slug}_id", []).extend(
                        self._values(annotation["term_accession_number"])
                    )
                if annotation.get("term_source_ref"):
                    characteristic_values.setdefault(f"{annotation_slug}_onto", []).extend(
                        self._values(annotation["term_source_ref"])
                    )
            for item in self.planner._as_list(channel.get("characteristics")):
                if not isinstance(item, dict) or not item.get("name", item.get("tag")):
                    continue
                slug = self._metadata_slug(item.get("name", item.get("tag")))
                values = self._values(item.get("value"))
                if slug and values:
                    characteristic_values.setdefault(slug, []).extend(values)
                for annotation in self.planner._as_list(item.get("annotations")):
                    if not isinstance(annotation, dict) or not annotation.get("field"):
                        continue
                    annotation_slug = "harmonized_" + self._metadata_slug(
                        annotation["field"]
                    )
                    characteristic_values.setdefault(annotation_slug, []).extend(
                        self._values(annotation.get("value"))
                    )
                    if annotation.get("term_accession_number"):
                        characteristic_values.setdefault(
                            f"{annotation_slug}_id", []
                        ).extend(self._values(annotation["term_accession_number"]))
                    if annotation.get("term_source_ref"):
                        characteristic_values.setdefault(
                            f"{annotation_slug}_onto", []
                        ).extend(self._values(annotation["term_source_ref"]))
        characteristics = {
            slug: tuple(self._values(values))
            for slug, values in characteristic_values.items()
        }
        metadata["characteristics"] = characteristics
        metadata["organism_part"] = (
            characteristics.get("organism_part")
            or characteristics.get("tissue")
            or metadata["source"]
        )
        metadata["developmental_stage"] = characteristics.get("developmental_stage", ())
        metadata["disease"] = characteristics.get("disease", ())
        metadata["genotype"] = characteristics.get("genotype", ())

        metadata["biomaterial_provider"] = tuple(
            self._values(channel.get("biomaterial_provider") for channel in channels)
        )
        metadata["molecule"] = tuple(
            self._values(channel.get("molecule") for channel in channels)
        )
        explicit_material_types = self._values(
            channel.get("material_type") for channel in channels
        )
        material_types = []
        for value in self._values(channel.get("molecule") for channel in channels):
            material_types.append(re.sub(r"^total\s+", "", value, flags=re.IGNORECASE))
        metadata["material_type"] = (
            tuple(explicit_material_types)
            or tuple(self._values(material_types))
            or metadata["organism_part"]
        )

        runs = [item for item in self.planner._as_list(sample.get("sra_run")) if isinstance(item, dict)]
        metadata["sra_accession"] = tuple(self._values(sample.get("sra_accession")))
        metadata["ena_accession"] = tuple(self._values(sample.get("ena_accession")))
        metadata["biosample_accession"] = tuple(
            self._values(run.get("biosample") for run in runs)
        )
        metadata["sra_run_accessions"] = tuple(
            self._values(run.get("run") for run in runs)
        )
        metadata["library_strategy"] = tuple(
            self._values(
                [sample.get("library_strategy"), *(run.get("library_strategy") for run in runs)]
            )
        )
        metadata["library_source"] = tuple(
            self._values(
                [sample.get("library_source"), *(run.get("library_source") for run in runs)]
            )
        )
        metadata["library_selection"] = tuple(
            self._values(
                [sample.get("library_selection"), *(run.get("library_selection") for run in runs)]
            )
        )
        metadata["library_layout"] = tuple(
            self._values(run.get("library_layout") for run in runs)
        )
        metadata["instrument_model"] = tuple(
            self._values(
                [sample.get("instrument_model"), *(run.get("instrument_model") for run in runs)]
            )
        )
        metadata["platform_accession"] = tuple(
            self._platform_accession_values(sample, package)
        )

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
        metadata["protocol_types"] = tuple(self._values(protocol_types))
        metadata["protocol_term_source_refs"] = tuple(self._values(protocol_sources))
        metadata["protocol_term_accession_numbers"] = tuple(self._values(protocol_accessions))

        database = self._metadata_database(package)
        metadata["metadata_source"] = tuple(
            self._values(
                database.get("public_id") or database.get("iid") or database.get("name")
            )
        )
        metadata["metadata_source_name"] = tuple(self._values(database.get("name")))
        metadata["metadata_source_uri"] = tuple(self._values(database.get("web_link")))
        return metadata

    def _render_sample_metadata(self, values: Mapping[str, Any]) -> dict:
        return {
            key: (
                {
                    characteristic: self._join_values(items)
                    for characteristic, items in value.items()
                }
                if key == "characteristics"
                else self._join_values(value)
            )
            for key, value in values.items()
        }

    def _attach_sample_values(
        self,
        adata,
        samples: Mapping[str, Mapping[str, Sequence[Any]]],
    ) -> None:
        _anndata, _numpy, pandas, _sparse = self._scientific_modules()
        rows = []
        for sample_accession, fields in samples.items():
            for field_name, values in fields.items():
                for ordinal, value in enumerate(values):
                    if value is None or str(value) == "":
                        continue
                    rows.append(
                        (
                            str(sample_accession),
                            str(field_name),
                            ordinal,
                            str(value),
                            self._metadata_value_type(value),
                        )
                    )
        frame = pandas.DataFrame(
            rows,
            columns=(
                "sample_accession",
                "field",
                "ordinal",
                "value",
                "value_type",
            ),
        )
        frame.index = [f"value_{index:06d}" for index in range(len(frame))]
        adata.uns["msc_metadata"] = {
            "schema_version": self.H5AD_METADATA_SCHEMA_VERSION,
            "sample_values": frame,
        }

    @staticmethod
    def _metadata_value_type(value: Any) -> str:
        if isinstance(value, bool):
            return "boolean"
        if isinstance(value, int):
            return "integer"
        if isinstance(value, float):
            return "number"
        return "string"

    @staticmethod
    def _is_sample_qualified(observation_id: str, sample_id: str) -> bool:
        token = re.escape(str(sample_id))
        return re.search(
            rf"(?:^|[-_.:]){token}(?:$|[-_.:])",
            str(observation_id),
            flags=re.IGNORECASE,
        ) is not None

    @staticmethod
    def _deduplicate_observation_ids(values: Sequence[str]) -> list[str]:
        totals = {}
        for value in values:
            totals[value] = totals.get(value, 0) + 1
        seen = {}
        rendered = []
        for value in values:
            if totals[value] == 1:
                rendered.append(value)
                continue
            seen[value] = seen.get(value, 0) + 1
            rendered.append(f"{value}-{seen[value]}")
        return rendered

    def _ensure_global_observation_ids(self, adatas: Mapping[str, Any]) -> None:
        locations: dict[str, list[tuple[str, int]]] = {}
        for sample_id, adata in adatas.items():
            for position, value in enumerate(adata.obs_names.astype(str)):
                locations.setdefault(value, []).append((sample_id, position))
        collisions = {
            value: entries
            for value, entries in locations.items()
            if len(entries) > 1
        }
        if not collisions:
            return
        for value, entries in collisions.items():
            for sample_id, position in entries:
                names = list(adatas[sample_id].obs_names.astype(str))
                names[position] = f"{value}-{sample_id}"
                adatas[sample_id].obs_names = names
        all_values = [
            value
            for adata in adatas.values()
            for value in adata.obs_names.astype(str)
        ]
        if len(all_values) != len(set(all_values)):
            raise ValueError("Observation identifiers remain non-unique after sample qualification.")

    def _characteristic_columns(self, packages: list[dict]) -> list[str]:
        columns = []
        for package in packages:
            for sample in self.planner._as_list(package.get("sample")):
                if not isinstance(sample, dict):
                    continue
                for channel in self.planner._as_list(sample.get("channel")):
                    if not isinstance(channel, dict):
                        continue
                    for annotation in self.planner._as_list(channel.get("annotations")):
                        if not isinstance(annotation, dict) or not annotation.get("field"):
                            continue
                        annotation_slug = "harmonized_" + self._metadata_slug(annotation["field"])
                        for candidate in (annotation_slug, f"{annotation_slug}_id", f"{annotation_slug}_onto"):
                            if candidate not in columns:
                                columns.append(candidate)
                    for item in self.planner._as_list(channel.get("characteristics")):
                        if not isinstance(item, dict):
                            continue
                        slug = self._metadata_slug(item.get("name", item.get("tag")))
                        if slug and slug not in columns:
                            columns.append(slug)
                        for annotation in self.planner._as_list(item.get("annotations")):
                            if not isinstance(annotation, dict) or not annotation.get("field"):
                                continue
                            annotation_slug = "harmonized_" + self._metadata_slug(
                                annotation["field"]
                            )
                            for candidate in (
                                annotation_slug,
                                f"{annotation_slug}_id",
                                f"{annotation_slug}_onto",
                            ):
                                if candidate not in columns:
                                    columns.append(candidate)
        return columns

    def _metadata_slug(self, value) -> str:
        value = self._join_values(value).lower()
        return re.sub(r"_+", "_", re.sub(r"[^a-z0-9]+", "_", value)).strip("_")

    def _values(self, values) -> list[str]:
        flattened = []
        seen = set()

        def visit(value):
            if value is None:
                return
            if isinstance(value, dict):
                for key in ("value", "name", "predefined", "public_id", "iid", "ref"):
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
            normalized = cleaned.casefold()
            if cleaned and normalized not in seen:
                flattened.append(cleaned)
                seen.add(normalized)

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
        return self._join_values(self._platform_accession_values(sample, package))

    def _platform_accession_values(self, sample: dict, package: dict) -> list[str]:
        references = set(self._values(sample.get("platform_ref")))
        values = []
        for platform in self.planner._as_list(package.get("platform")):
            if not isinstance(platform, dict):
                continue
            identifiers = {platform.get("iid"), *self._values(platform.get("accession"))}
            if references and not references.intersection(identifier for identifier in identifiers if identifier):
                continue
            values.extend(self._values(platform.get("accession")))
        return self._values(values or references)

    def _attach_miniml(
        self,
        adata,
        packages: list[dict],
        source_json: str,
        source_json_sha256: str | None,
        sample_id: str | None = None,
        artifact_parent: Path | None = None,
    ) -> None:
        _anndata, _numpy, pandas, _sparse = self._scientific_modules()
        rows = []
        transported_packages = []
        for package_index, package in enumerate(packages):
            if sample_id is not None and not self._package_has_sample(package, sample_id):
                continue
            to_mapping = getattr(package, "to_mapping", None)
            transported_packages.append(
                to_mapping() if callable(to_mapping) else package
            )
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
        portable_source, source_scope = self._portable_location(
            source_json, artifact_parent or Path.cwd()
        )
        adata.uns["msc_miniml"] = {
            "schema_version": self.MINIML_SCHEMA_VERSION,
            "packages_json": json.dumps(
                transported_packages, sort_keys=True, ensure_ascii=False
            ),
            "source_json": portable_source,
            "source_json_scope": source_scope,
            "path_base": "artifact_parent",
            "source_sha256": source_json_sha256 or "",
            "publication_policy": self.PUBLICATION_POLICY,
            "metadata_source": self._join_values(
                database.get("public_id") or database.get("iid") or database.get("name")
            ),
            "metadata_source_name": self._join_values(database.get("name")),
            "metadata_source_uri": self._join_values(database.get("web_link")),
            "fields": fields,
        }
        parameter_occurrences = []
        for package in packages:
            if sample_id is None:
                parameter_occurrences.extend(_parameter_rows(package))
                continue
            sample = next(
                (
                    item for item in self.planner._as_list(package.get("sample"))
                    if isinstance(item, dict)
                    and self.planner.sample_accession(item) == sample_id
                ),
                None,
            )
            if sample is not None:
                parameter_occurrences.extend(_parameter_rows(package, sample=sample))
        if parameter_occurrences:
            parameters = pandas.DataFrame(parameter_occurrences)
            parameters.index = [
                f"parameter_{index:06d}" for index in range(len(parameters))
            ]
            adata.uns["msc_assay"] = {
                "schema_version": "2.0",
                "parameters": parameters,
            }

    def _attach_harmonization(self, adata, resolution) -> None:
        if resolution is None or not resolution.enabled:
            return
        _anndata, _numpy, pandas, _sparse = self._scientific_modules()
        rows = [
            {
                "sample_accession": item.sample_accession,
                "destination": item.destination,
                "value": item.value,
                "id": item.identifier or "",
                "ontology": item.ontology or "",
                "source_field": item.source_field,
                "hierarchy_depth": (
                    item.hierarchy_depth if item.hierarchy_depth is not None else -1
                ),
                "status": item.status,
            }
            for item in resolution.selections
        ]
        selections = pandas.DataFrame(rows, columns=(
            "sample_accession", "destination", "value", "id", "ontology",
            "source_field", "hierarchy_depth", "status",
        ))
        selections.index = [f"selection_{index:06d}" for index in range(len(selections))]
        adata.uns["msc_harmonization"] = {
            "schema_version": "1.0",
            "enabled": True,
            "applied": bool(resolution.applied),
            "profile": dict(resolution.profile or {}),
            "selections": selections,
            "warnings": list(resolution.warnings),
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

    def _combine(
        self,
        adatas: dict[str, object],
        *,
        allow_unverified: bool = False,
    ):
        if not adatas:
            raise ValueError("No sample H5ADs were produced.")
        anndata, _numpy, _pandas, sparse = self._scientific_modules()
        missing_evidence = self._missing_combination_evidence(adatas)
        if missing_evidence and not allow_unverified:
            rendered = ", ".join(
                f"{dimension} ({', '.join(sample_ids)})"
                for dimension, sample_ids in sorted(missing_evidence.items())
            )
            raise ValueError(
                "Cannot combine samples without positive compatibility evidence "
                f"for every sample: {rendered}. Set allow_unverified_combination "
                "only to publish an explicitly acknowledged partial result."
            )
        organisms = {
            str(
                adata.obs["msc.sample.channel.organism.value"].iloc[0]
            ).strip()
            for adata in adatas.values()
            if "msc.sample.channel.organism.value" in adata.obs
            and str(
                adata.obs["msc.sample.channel.organism.value"].iloc[0]
            ).strip()
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
            "metadata_schema_version": self.H5AD_METADATA_SCHEMA_VERSION,
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

    def _missing_combination_evidence(
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
                    for value in adata.obs[
                        "msc.sample.channel.organism.value"
                    ]
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
            namespace = self._feature_namespace(adata)
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

    def _write_dataset_bundle(
        self,
        *,
        adatas: Mapping[str, Any],
        combined: Any,
        result: ConversionResult,
        planned: dict[str, Asset],
        json_path: str,
        overwrite: bool,
    ) -> None:
        destinations = [Path(path) for path in result.sample_h5ads.values()]
        if result.combined_h5ad:
            destinations.append(Path(result.combined_h5ad))
        destinations.append(Path(result.manifest_path))
        for destination in destinations:
            if destination.exists() and not overwrite:
                raise FileExistsError(f"Output already exists: {destination}")

        output_dir = Path(result.manifest_path).parent
        output_dir.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(
            prefix=".json2h5ad-staging-", dir=output_dir
        ) as temporary_dir:
            staging = Path(temporary_dir)
            staged: list[tuple[Path, Path]] = []
            for sample_id, adata in adatas.items():
                destination = Path(result.sample_h5ads[sample_id])
                source = staging / destination.name
                self._write_h5ad(adata, source, overwrite=True)
                staged.append((source, destination))
            if combined is not None and result.combined_h5ad:
                destination = Path(result.combined_h5ad)
                source = staging / destination.name
                self._write_h5ad(combined, source, overwrite=True)
                staged.append((source, destination))
            manifest_destination = Path(result.manifest_path)
            manifest_source = staging / manifest_destination.name
            self._write_manifest(
                result,
                planned,
                json_path=json_path,
                overwrite=True,
                output_path=manifest_source,
            )
            staged.append((manifest_source, manifest_destination))
            self._commit_dataset_bundle(staged, overwrite=overwrite, staging=staging)

    def _commit_dataset_bundle(
        self,
        staged: Sequence[tuple[Path, Path]],
        *,
        overwrite: bool,
        staging: Path,
    ) -> None:
        output_dir = staged[0][1].parent
        backup_dir = output_dir / f".json2h5ad-recovery-{uuid4().hex}"
        backups: list[tuple[Path, Path]] = []
        installed: list[Path] = []
        recovery_failed = False
        try:
            if overwrite:
                backup_dir.mkdir()
                for index, (_source, destination) in enumerate(staged):
                    if destination.exists():
                        backup = backup_dir / f"{index}-{destination.name}"
                        os.replace(destination, backup)
                        backups.append((backup, destination))
                self._fsync_directory(output_dir)
            for source, destination in staged:
                os.replace(source, destination)
                installed.append(destination)
                if destination.suffix == ".h5ad":
                    destination.chmod(0o660)
            self._fsync_directory(output_dir)
        except BaseException as publication_error:
            recovery_errors: list[BaseException] = []
            for destination in reversed(installed):
                try:
                    destination.unlink(missing_ok=True)
                except BaseException as error:
                    recovery_errors.append(error)
            for backup, destination in reversed(backups):
                try:
                    if backup.exists():
                        os.replace(backup, destination)
                except BaseException as error:
                    recovery_errors.append(error)
            try:
                self._fsync_directory(output_dir)
            except BaseException as error:
                recovery_errors.append(error)
            if recovery_errors:
                recovery_failed = True
                recovery_paths = [
                    path for path, _destination in backups if path.exists()
                ]
                if backup_dir.exists():
                    recovery_paths.append(backup_dir)
                raise DatasetBundleRecoveryError(
                    publication_error, recovery_errors, recovery_paths
                ) from publication_error
            raise
        finally:
            if backup_dir.exists() and not recovery_failed:
                shutil.rmtree(backup_dir)

    @staticmethod
    def _fsync_directory(path: Path) -> None:
        descriptor = os.open(
            path, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0)
        )
        try:
            os.fsync(descriptor)
        finally:
            os.close(descriptor)

    def _write_manifest(
        self,
        result: ConversionResult,
        planned: dict[str, Asset],
        json_path: str,
        overwrite: bool,
        output_path: Path | None = None,
    ) -> None:
        logical_path = Path(result.manifest_path)
        path = output_path or logical_path
        if logical_path.exists() and not overwrite:
            raise FileExistsError(f"Output already exists: {logical_path}")
        base = logical_path.parent.resolve()
        source_json, source_json_scope = self._portable_location(json_path, base)
        combined_h5ad, combined_h5ad_scope = self._portable_location(result.combined_h5ad, base)
        sample_h5ads = {}
        sample_h5ad_scopes = {}
        for sample, value in result.sample_h5ads.items():
            sample_h5ads[sample], sample_h5ad_scopes[sample] = self._portable_location(value, base)
        retained_h5ads = []
        retained_h5ad_scopes = []
        for value in result.retained_h5ads:
            portable, scope = self._portable_location(value, base)
            retained_h5ads.append(portable)
            retained_h5ad_scopes.append(scope)
        payload = {
            "path_base": "artifact_parent",
            "h5ad_metadata_schema_version": self.H5AD_METADATA_SCHEMA_VERSION,
            "study_accession": result.study_accession,
            "source_json": source_json,
            "source_json_scope": source_json_scope,
            "combined_h5ad": combined_h5ad,
            "combined_h5ad_scope": combined_h5ad_scope,
            "sample_h5ads": sample_h5ads,
            "sample_h5ad_scopes": sample_h5ad_scopes,
            "retained_h5ads": retained_h5ads,
            "retained_h5ad_scopes": retained_h5ad_scopes,
            "pipeline_runs": [
                self._portable_pipeline_run(run, base)
                for run in result.pipeline_runs
            ],
            "warnings": result.warnings,
            "failures": result.failures,
            "errors": result.errors,
            "partial": result.partial,
            "assets": {
                sample: self._portable_asset(asset, base)
                for sample, asset in planned.items()
            },
        }
        temporary = path.with_suffix(path.suffix + ".tmp")
        with open(temporary, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
        os.replace(temporary, path)

    def _portable_pipeline_run(self, run: PipelineRun, base: Path) -> dict:
        work_dir, work_dir_scope = self._portable_location(run.work_dir, base)
        out_dir, out_dir_scope = self._portable_location(run.out_dir, base)
        log_path, log_path_scope = self._portable_location(run.log_path, base)
        annotation_source, annotation_source_scope = self._portable_location(
            run.annotation_source, base
        )
        effective_annotation, effective_annotation_scope = self._portable_location(
            run.effective_annotation, base
        )
        return {
            "pipeline": run.pipeline,
            "revision": run.revision,
            "command": [self._portable_command_argument(value, base) for value in run.command],
            "work_dir": work_dir,
            "work_dir_scope": work_dir_scope,
            "out_dir": out_dir,
            "out_dir_scope": out_dir_scope,
            "returncode": run.returncode,
            "log_path": log_path,
            "log_path_scope": log_path_scope,
            "annotation_source": annotation_source,
            "annotation_source_scope": annotation_source_scope,
            "annotation_format": run.annotation_format,
            "annotation_sha256": run.annotation_sha256,
            "effective_annotation": effective_annotation,
            "effective_annotation_scope": effective_annotation_scope,
            "warnings": run.warnings,
        }

    def _portable_asset(self, asset: Asset, base: Path) -> dict:
        path, path_scope = self._portable_location(asset.path, base)
        annotation_source, annotation_source_scope = self._portable_location(
            asset.annotation_source, base
        )
        effective_annotation, effective_annotation_scope = self._portable_location(
            asset.effective_annotation, base
        )
        return {
            "path": path,
            "path_scope": path_scope,
            "kind": asset.kind,
            "source": asset.source,
            "role": asset.role,
            "reference": asset.reference,
            "annotation_source": annotation_source,
            "annotation_source_scope": annotation_source_scope,
            "annotation_format": asset.annotation_format,
            "annotation_sha256": asset.annotation_sha256,
            "effective_annotation": effective_annotation,
            "effective_annotation_scope": effective_annotation_scope,
        }

    def _portable_command_argument(self, value: str, base: Path) -> str:
        if isinstance(value, str) and os.path.isabs(value):
            return self._portable_location(value, base)[0]
        return value

    def _portable_location(self, value: str | None, base: Path) -> tuple[str | None, str | None]:
        if value in (None, ""):
            return value, None
        rendered = str(value)
        parsed = urlparse(rendered)
        if parsed.scheme and parsed.scheme.lower() != "file":
            return rendered, "remote"
        local = parsed.path if parsed.scheme.lower() == "file" else rendered
        absolute = Path(local)
        if not absolute.is_absolute():
            absolute = Path.cwd() / absolute
        relative = os.path.relpath(absolute.resolve(strict=False), Path(base).resolve(strict=False))
        scope = "external" if relative == ".." or relative.startswith(f"..{os.sep}") else "internal"
        return relative, scope

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


def _tenx_member(value: str) -> tuple[str, str] | None:
    filename = os.path.basename(urlparse(str(value)).path)
    for suffix in (".gz", ".bz2", ".xz", ".zip"):
        if filename.casefold().endswith(suffix):
            filename = filename[: -len(suffix)]
            break
    match = _TENX_MEMBER.fullmatch(filename)
    if match is None:
        return None
    raw_role = match.group("role").casefold()
    role = (
        "matrix"
        if raw_role == "matrix.mtx"
        else "barcodes"
        if raw_role == "barcodes.tsv"
        else "features"
    )
    return match.group("prefix").casefold(), role


def _tenx_local_name(role: str, source: str, *, legacy: bool) -> str:
    if legacy:
        return {
            "matrix": "matrix.mtx",
            "barcodes": "barcodes.tsv",
            "features": "genes.tsv",
        }[role]
    filename = os.path.basename(urlparse(str(source)).path).casefold()
    compressed = ".gz" if filename.endswith(".gz") else ""
    if role == "matrix":
        return f"matrix.mtx{compressed}"
    if role == "barcodes":
        return f"barcodes.tsv{compressed}"
    stem = "features" if ".features.tsv" in filename else "genes"
    return f"{stem}.tsv{compressed}"


def _copy_10x_member(source: str, destination: Path, *, decompress: bool) -> None:
    if decompress and str(source).casefold().endswith(".gz"):
        with gzip.open(source, "rb") as input_stream:
            with destination.open("wb") as output_stream:
                shutil.copyfileobj(input_stream, output_stream)
        return
    shutil.copyfile(source, destination)
