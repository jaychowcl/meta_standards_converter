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
import heapq
import logging
import os
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field, replace
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any, Callable, Mapping, Protocol, Sequence
from urllib.parse import urlparse
from uuid import uuid4

from meta_standards_converter.retrieval import AssetDownloader, RetrievalPolicy
from meta_standards_converter.runtime_contracts import (
    SafeErrorEnvelope,
    get_resource_profile,
)
from meta_standards_converter.miniml import (
    harmonized_value_mappings,
    iter_harmonization_operations,
    iter_harmonization_patches,
)
from meta_standards_converter.sources.json import JSONPackageSource
from .dataset_combination import (
    DatasetCombinationPolicy,
    DatasetCompatibilityError as _DatasetCompatibilityError,
)
from meta_standards_converter.metadata.projection.assay import _parameter_rows, _parameter_summary
from meta_standards_converter.metadata.provenance import patch_provenance_columns
from meta_standards_converter.metadata.interpretation import MINiMLMetadataProvider, MINiMLMetadataService


from meta_standards_converter.expression.assets import Asset
from meta_standards_converter.expression.assets import AssetManifest
from meta_standards_converter.expression.planning import SourcePlanner
from meta_standards_converter.expression.references import ReferenceResolver
from meta_standards_converter.expression.references import AnnotationConverter
from meta_standards_converter.expression.nfcore import NFCoreRunner
from meta_standards_converter.expression.catalogue import PipelineRun
from meta_standards_converter.expression.catalogue import ConversionResult
from meta_standards_converter.expression.catalogue import BatchConversionResult
from meta_standards_converter.expression.catalogue import RawProcessingResult
from meta_standards_converter.expression.catalogue import DatasetBundleRecoveryError
from meta_standards_converter.expression.catalogue import _safe_group_failure
from meta_standards_converter.metadata.projection.anndata import MetadataProjectionContext
from meta_standards_converter.metadata.projection.anndata import AnnDataMetadataProjection
from meta_standards_converter.metadata.projection.anndata import AnnDataProjectionError
from meta_standards_converter.metadata.projection.anndata import AnnDataMetadataProjector
from meta_standards_converter.expression.memory import _available_memory_bytes
from meta_standards_converter.expression.memory import _path_size_bytes
from meta_standards_converter.expression.memory import _estimate_asset_memory_bytes
from meta_standards_converter.expression.readers import ProcessedAssetReader, underlying_suffix

from meta_standards_converter.expression.normalization import AnnDataNormalizer
from meta_standards_converter.expression.checkpoints import ProcessedCheckpointStore
from meta_standards_converter.expression.catalogue import CataloguePublisher
from meta_standards_converter.metadata.projection.anndata import AnnDataProjectorRunner

logger = logging.getLogger(__name__)







_TENX_MEMBER = re.compile(
    r"^(?P<prefix>.+?)[._](?P<role>"
    r"matrix\.mtx|barcodes\.tsv|genes\.tsv|features\.tsv)$",
    re.IGNORECASE,
)


































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
    """Top-level JSON-to-H5AD conversion orchestrator."""

    def __init__(
        self,
        planner: SourcePlanner | None = None,
        reader=None,
        asset_cache_dir=None,
        pipeline_runner: NFCoreRunner | None = None,
        downloader: AssetDownloader | None = None,
        metadata_projectors: Sequence[AnnDataMetadataProjector] | None = None,
        package_source: JSONPackageSource | None = None,
        retrieval_policy: RetrievalPolicy | None = None,
        resource_profile: str = "standard",
        resource_overrides: Mapping[str, int | float] | None = None,
        metadata_service: MINiMLMetadataProvider | None = None,
        combination_policy: DatasetCombinationPolicy | None = None,
        available_memory: Callable[[], int] | None = None,
        memory_estimator: Callable[[str, Asset], int] | None = None,
    ):
        if downloader is not None and retrieval_policy is not None:
            raise ValueError(
                "downloader and retrieval_policy are mutually exclusive"
            )
        self.planner = planner or SourcePlanner()
        self.reader = reader or ProcessedAssetReader()
        self.asset_cache_dir = asset_cache_dir
        self.pipeline_runner = pipeline_runner or NFCoreRunner()
        self.downloader = downloader
        self.metadata_projectors = tuple(metadata_projectors or ())
        self.package_source = package_source or JSONPackageSource()
        self.metadata_service = metadata_service or MINiMLMetadataService()
        self.normalizer = AnnDataNormalizer(metadata_service=self.metadata_service,
            planner=self.planner, localize=self._local_path, package_version=self._package_version,
            combination_policy=combination_policy)
        self.projectors = AnnDataProjectorRunner(self.metadata_projectors)
        self.checkpoints = ProcessedCheckpointStore(self._package_version)
        self.publisher = CataloguePublisher()
        resolved_resource_profile = (
            retrieval_policy.resource_profile
            if retrieval_policy is not None
            else get_resource_profile(
                resource_profile,
                overrides=resource_overrides,
            )
        )
        self.retrieval_policy = retrieval_policy or RetrievalPolicy(
            resource_profile=resolved_resource_profile
        )
        self.resource_profile = resolved_resource_profile
        self.available_memory = available_memory or _available_memory_bytes
        self.memory_estimator = memory_estimator or _estimate_asset_memory_bytes

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
        force_memory: bool = False,
        processed_checkpoint_dir: str | None = None,
        allow_invalid: bool = False,
        allow_unverified_combination: bool = False,
        use_harmonization_overrides: bool = False,
        **options,
    ) -> ConversionResult | BatchConversionResult:
        if not os.path.exists(json_path):
            raise FileNotFoundError(f"MINiML JSON file not found: {json_path}")
        if force_memory and not resume:
            raise ValueError("force_memory requires resume=True")
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
            force_memory=force_memory,
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
                result.failures.append(_safe_group_failure(error, group.dataset_id))
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
        force_memory: bool = False,
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
                str(self.asset_cache_dir or out_path / ".cache"),
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
        characteristic_columns = self.normalizer.characteristic_columns(packages)
        source_json = os.path.abspath(source_json)
        source_json_sha256 = self._sha256(source_json)
        checkpoint_root = (
            Path(processed_checkpoint_dir)
            if processed_checkpoint_dir
            else out_path / ".processed"
        ) / study_accession
        result = ConversionResult(study_accession=study_accession)
        result.warnings.extend(getattr(harmonization_resolution, "warnings", ()))
        checkpoint_sources: dict[str, Path] = {}
        used_observation_ids: set[str] = set()

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
            checkpoint = self.checkpoints.key(
                checkpoint_root,
                sample_id=sample_id,
                source_json_sha256=source_json_sha256,
                sample=sample_context[sample_id][0],
                asset=asset,
                orientation=orientation,
            )
            checkpoint_metadata = (
                self.checkpoints.metadata(checkpoint) if resume else None
            )
            if checkpoint_metadata is not None:
                checkpoint_path = checkpoint[0]
                names = self.checkpoints.observation_ids(checkpoint_path)
                if names is not None and not (used_observation_ids & names):
                    result.warnings.extend(checkpoint_metadata.get("warnings", ()))
                    result.errors.extend(checkpoint_metadata.get("errors", ()))
                    used_observation_ids.update(names)
                    sample_path = out_path / f"{sample_id}.h5ad"
                    result.sample_h5ads[sample_id] = str(sample_path)
                    checkpoint_sources[sample_id] = checkpoint_path
                    continue

            local_path = self._local_path(asset.path, md5=asset.md5)
            available_bytes = int(self.available_memory())
            estimated_peak_bytes = int(self.memory_estimator(local_path, asset))
            if available_bytes <= 0 or estimated_peak_bytes <= 0:
                raise ValueError("Memory admission inputs must be positive byte counts.")
            fixed_limit = int(self.resource_profile.max_in_memory_matrix_bytes)
            limit_bytes = (
                int(available_bytes * self.resource_profile.force_memory_fraction)
                if force_memory
                else min(
                    fixed_limit,
                    int(available_bytes * self.resource_profile.available_memory_fraction),
                )
            )
            admitted = estimated_peak_bytes <= limit_bytes
            result.memory_report.append(
                {
                    "sample_id": sample_id,
                    "estimated_peak_bytes": estimated_peak_bytes,
                    "available_memory_bytes": available_bytes,
                    "fixed_profile_bytes": fixed_limit,
                    "limit_bytes": limit_bytes,
                    "force_memory": bool(force_memory),
                    "decision": "admitted" if admitted else "skipped",
                    "reason": (
                        "within_memory_limit"
                        if admitted
                        else "estimated_peak_exceeds_memory_limit"
                    ),
                }
            )
            if not admitted:
                result.failures.append(
                    f"{sample_id}: skipped by memory admission; estimated peak "
                    f"{estimated_peak_bytes} bytes exceeds {limit_bytes} bytes"
                )
                continue

            adata = self.reader.read(asset, orientation=orientation, localize=self._local_path)
            base_metadata = self.normalizer.normalize(
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
            warning_start = len(result.warnings)
            error_start = len(result.errors)
            self.projectors.project_sample(
                adata,
                projection_context,
                warnings=result.warnings,
                errors=result.errors,
                allow_invalid=allow_invalid,
            )
            self.normalizer.attach_miniml(
                adata,
                packages=source_packages or packages,
                source_json=source_json,
                source_json_sha256=source_json_sha256,
                sample_id=sample_id,
                artifact_parent=out_path,
            )
            self.normalizer.attach_harmonization(
                adata,
                harmonization_resolution,
                packages=source_packages or packages,
                sample_id=sample_id,
            )
            self.normalizer.ensure_observation_ids(
                adata, sample_id, used_observation_ids
            )
            self.checkpoints.write(
                checkpoint,
                adata,
                warnings=result.warnings[warning_start:],
                errors=result.errors[error_start:],
            )
            sample_path = out_path / f"{sample_id}.h5ad"
            result.sample_h5ads[sample_id] = str(sample_path)
            checkpoint_sources[sample_id] = checkpoint[0]
            del adata

        if allow_unverified_combination:
            result.warnings.append(
                "allow_unverified_combination is ignored because expression "
                "matrices are published only as a per-sample catalogue."
            )

        result.manifest_path = str(out_path / f"{study_accession}.json2h5ad.json")
        self.publisher.publish(
            sample_sources=checkpoint_sources,
            result=result,
            planned=planned,
            json_path=source_json,
            overwrite=overwrite,
        )
        return result






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
        return self.metadata_service.study_accession(packages)

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

















































    def _local_path(self, value: str, md5: str | None = None) -> str:
        if self.downloader is None:
            parsed = urlparse(value)
            if parsed.scheme in ("", "file"):
                return parsed.path if parsed.scheme == "file" else value
            raise RuntimeError("Remote asset downloader has not been configured.")
        return self.downloader.localize(value, md5=md5)


    def _text(self, value) -> str | None:
        return self.metadata_service.text(value)

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
