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

import tempfile

from dataclasses import dataclass, replace

from pathlib import Path

from typing import Any, Mapping, Sequence

from meta_standards_converter.artifact_bundle import (
    DurableArtifactBundlePublisher,
    PublishedArtifactBundle,
)

from .json2h5ad import BatchConversionResult, ConversionResult, JSON2H5ADConverter




from meta_standards_converter.expression.components import (AnnDataComponentExporter, AnnDataMetadataExportResult, AnnDataMetadataBatchResult)

class JSON2OBSConverter:
    def __init__(self, h5ad_converter=None, components=None):
        self.h5ad_converter = h5ad_converter or JSON2H5ADConverter()
        self.components = components or AnnDataComponentExporter()

    def convert(
        self,
        source: str | Path,
        *,
        outdir: str | Path,
        include_var: bool = False,
        include_uns: bool = False,
        overwrite: bool = False,
        replacement_profile: Mapping[str, Any] | None = None,
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
                replacement_profile=replacement_profile,
                **options,
            )
            if isinstance(converted, BatchConversionResult):
                results: dict[str, AnnDataMetadataExportResult] = {}
                for dataset_id, conversion in converted.conversions.items():
                    target = destination / dataset_id
                    results[dataset_id] = self.components.export(
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
            return self.components.export(
                converted,
                destination,
                include_var=include_var,
                include_uns=include_uns,
                overwrite=overwrite,
            )
