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
from dataclasses import replace
from meta_standards_converter.expression.components import _publish_bundle

import csv

from dataclasses import dataclass

from pathlib import Path

from typing import Any, Mapping, Protocol, Sequence

from meta_standards_converter.sources.json import JSONPackageSource

from meta_standards_converter.metadata.provenance import patch_provenance_columns

from meta_standards_converter.metadata.projection.assay import _parameter_summary

from meta_standards_converter.metadata.interpretation import MINiMLMetadataProvider, MINiMLMetadataService


from meta_standards_converter.metadata.projection.tabular import (TabularMetadataContext, TabularMetadataProjection, TabularMetadataProjector, TabularConversionResult, TabularProjectionError, MSCMetadataProjector)
class JSON2DelimitedConverter:
    """Host one or more injected sample projectors and write one table."""

    delimiter = "\t"

    def __init__(
        self,
        metadata_projectors: Sequence[TabularMetadataProjector] | None = None,
        package_source: JSONPackageSource | None = None,
        metadata_service: MINiMLMetadataProvider | None = None,
    ) -> None:
        self.metadata_projectors = tuple(
            [MSCMetadataProjector()]
            if metadata_projectors is None
            else metadata_projectors
        )
        self.package_source = package_source or JSONPackageSource()
        self.metadata_service = metadata_service or MINiMLMetadataService()

    def convert_source(
        self,
        source: str | Path,
        destination: str | Path,
        *,
        allow_invalid: bool = False,
        overwrite: bool = False,
        use_harmonization_overrides: bool = False,
    ) -> TabularConversionResult:
        loaded = self.package_source.load(source)
        records: list[dict[str, Any]] = []
        preferred: list[str] = []
        warnings = list(loaded.warnings)
        errors: list[str] = []
        for original_group in loaded.groups:
            group = original_group.resolved(enabled=use_harmonization_overrides)
            resolution = group.harmonization_resolution
            warnings.extend(getattr(resolution, "warnings", ()))
            for package in group.packages:
                study_accession = (
                    self.metadata_service.study_accession([package])
                    or group.dataset_id
                )
                for sample in self.metadata_service.samples(package):
                    sample_accession = (
                        self.metadata_service.sample_accession(sample) or ""
                    )
                    base = self.metadata_service.sample_metadata(
                        sample, package
                    )
                    base = {
                        **base,
                        "modality": self.metadata_service.sample_modality(sample),
                        "harmonization": [
                            vars(item)
                            for item in getattr(resolution, "selections", ())
                            if item.sample_accession == sample_accession
                        ],
                    }
                    context = TabularMetadataContext(
                        package=package,
                        sample=sample,
                        dataset_id=group.dataset_id,
                        study_accession=study_accession,
                        sample_accession=sample_accession,
                        base_metadata=base,
                    )
                    row: dict[str, Any] = {}
                    for projector in self.metadata_projectors:
                        projection = projector.project_sample(context=context)
                        if not isinstance(projection, TabularMetadataProjection):
                            raise TypeError(
                                "tabular projector must return "
                                "TabularMetadataProjection"
                            )
                        collisions = set(row).intersection(projection.values)
                        if collisions:
                            raise ValueError(
                                "Tabular projector column collision: "
                                + ", ".join(sorted(collisions))
                            )
                        row.update(projection.values)
                        for column in projection.columns:
                            if column not in preferred:
                                preferred.append(column)
                        warnings.extend(
                            f"{sample_accession}: {item}"
                            for item in projection.warnings
                        )
                        errors.extend(
                            f"{sample_accession}: {item}"
                            for item in projection.errors
                        )
                    records.append(row)
        if errors and not allow_invalid:
            raise TabularProjectionError("; ".join(errors))
        extra = sorted(
            {
                key
                for record in records
                for key in record
                if key not in preferred
            }
        )
        columns = (*preferred, *extra)
        output = Path(destination)
        if output.exists() and not overwrite:
            raise FileExistsError(f"Output already exists: {output}")
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(
                stream,
                fieldnames=columns,
                delimiter=self.delimiter,
                lineterminator="\n",
                extrasaction="ignore",
            )
            writer.writeheader()
            writer.writerows(records)
        return TabularConversionResult(
            row_count=len(records),
            columns=columns,
            dataset_ids=tuple(group.dataset_id for group in loaded.groups),
            warnings=tuple(dict.fromkeys(warnings)),
            errors=tuple(dict.fromkeys(errors)),
            output_path=str(output),
        )

class JSON2TSVConverter(JSON2DelimitedConverter):
    def __init__(
        self,
        metadata_projectors: Sequence[TabularMetadataProjector] | None = None,
        package_source: JSONPackageSource | None = None,
        metadata_service: MINiMLMetadataProvider | None = None,
        *,
        output_format: str = "tsv",
    ) -> None:
        if output_format not in {"tsv", "csv"}:
            raise ValueError("output_format must be 'tsv' or 'csv'")
        self.output_format = output_format
        self.delimiter = "\t" if output_format == "tsv" else ","
        super().__init__(
            metadata_projectors=metadata_projectors,
            package_source=package_source,
            metadata_service=metadata_service,
        )


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
            metadata_projectors=self.metadata_projectors,
            package_source=self.package_source,
            metadata_service=self.metadata_service,
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
