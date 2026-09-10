# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from __future__ import annotations
from meta_standards_converter.runtime_contracts import SafeErrorEnvelope
import json
import os
import shutil
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping, Sequence
from urllib.parse import urlparse
from uuid import uuid4
from meta_standards_converter.expression.assets import Asset


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
    memory_report: list[dict[str, Any]] = field(default_factory=list)

    @property
    def primary_h5ad(self) -> str | None:
        return next(iter(self.sample_h5ads.values()), self.combined_h5ad)

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
            "memory_report": list(self.memory_report),
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


class CataloguePublisher:
    MINIML_SCHEMA_VERSION = "1.0"

    H5AD_METADATA_SCHEMA_VERSION = "2.0"

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

    def publish(
        self,
        *,
        sample_sources: Mapping[str, Path],
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
            for sample_id, checkpoint_path in sample_sources.items():
                destination = Path(result.sample_h5ads[sample_id])
                source = staging / destination.name
                shutil.copy2(checkpoint_path, source)
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
            "artifact_kind": "per_sample_h5ad_catalogue",
            "expression_integration": "none",
            "sample_count": len(sample_h5ads),
            "combination_compatibility": {
                "assessed": False,
                "verified": False,
                "reason": "catalogue_only_no_expression_matrix_combination",
            },
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
            "memory_report": result.memory_report,
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


def _safe_group_failure(error: BaseException, dataset_id: str) -> str:
    envelope = SafeErrorEnvelope.from_exception(
        error,
        provider="meta_standards_converter",
        stage="json2h5ad_group_conversion",
        item_id=dataset_id,
    )
    return (
        f"{dataset_id}: {envelope.error_type} "
        f"[correlation_id={envelope.correlation_id}]"
    )
