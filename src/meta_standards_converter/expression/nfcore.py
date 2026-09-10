# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from __future__ import annotations
import json
import csv
import logging
import os
import re
import shutil
import subprocess
from pathlib import Path
from urllib.parse import urlparse
from meta_standards_converter.expression.assets import Asset
from meta_standards_converter.expression.planning import SourcePlanner
from meta_standards_converter.expression.references import ReferenceResolver
from meta_standards_converter.expression.references import AnnotationConverter
from meta_standards_converter.expression.catalogue import PipelineRun
from meta_standards_converter.expression.catalogue import RawProcessingResult

logger = logging.getLogger(__name__)

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
