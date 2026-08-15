# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Project parsed JSON samples into organization-neutral TSV or CSV tables."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Protocol, Sequence

from .json_source import JSONPackageSource
from .harmonization_provenance import patch_provenance_columns
from .mage_tab_projection import _parameter_summary
from .miniml_metadata import MINiMLMetadataProvider, MINiMLMetadataService


@dataclass(frozen=True)
class TabularMetadataContext:
    """Read-only sample context supplied to a tabular projector."""

    package: Mapping[str, Any]
    sample: Mapping[str, Any]
    dataset_id: str
    study_accession: str
    sample_accession: str
    base_metadata: Mapping[str, Any]


@dataclass(frozen=True)
class TabularMetadataProjection:
    """One projector's columns, values, and validation diagnostics."""

    values: Mapping[str, Any]
    columns: tuple[str, ...] = ()
    warnings: tuple[str, ...] = ()
    errors: tuple[str, ...] = ()


class TabularMetadataProjector(Protocol):
    def project_sample(
        self, *, context: TabularMetadataContext
    ) -> TabularMetadataProjection:
        """Project one parsed sample into tabular values."""


@dataclass(frozen=True)
class TabularConversionResult:
    """Result and compatibility paths for one tabular artifact generation."""

    row_count: int
    columns: tuple[str, ...]
    dataset_ids: tuple[str, ...]
    warnings: tuple[str, ...] = ()
    errors: tuple[str, ...] = ()
    output_path: str | None = None
    manifest_path: str | None = None
    bundle_pointer_path: str | None = None

    @property
    def partial(self) -> bool:
        return bool(self.errors)

    def to_dict(self) -> dict[str, Any]:
        return {
            "operation": "manifest",
            "status": "partial" if self.partial else "complete",
            "row_count": self.row_count,
            "columns": list(self.columns),
            "dataset_ids": list(self.dataset_ids),
            "warnings": list(self.warnings),
            "errors": list(self.errors),
            "artifacts": {
                "table": self.output_path,
                "manifest": self.manifest_path,
                "bundle_pointer": self.bundle_pointer_path,
            },
        }


class TabularProjectionError(ValueError):
    pass


class MSCMetadataProjector:
    """Default canonical tabular view of MINiML sample metadata."""

    COLUMNS = (
        "msc.sample.accession",
        "msc.series.accession",
        "msc.sample.title",
        "msc.sample.description",
        "msc.sample.channel.organism.value",
        "msc.sample.channel.organism.taxid",
        "msc.sample.channel.organism_part",
        "msc.sample.channel.developmental_stage",
        "msc.sample.channel.disease",
        "msc.sample.channel.genotype",
        "msc.sample.channel.source",
        "msc.sample.channel.biomaterial_provider",
        "msc.sample.channel.material_type",
        "msc.sample.channel.molecule",
        "msc.platform.accession",
        "msc.archive.sra_accession",
        "msc.archive.ena_accession",
        "msc.archive.biosample_accession",
        "msc.archive.sra_run_accessions",
        "msc.library.strategy",
        "msc.library.source",
        "msc.library.selection",
        "msc.library.layout",
        "msc.instrument.model",
        "msc.protocol.types",
        "msc.protocol.term_source_refs",
        "msc.protocol.term_accession_numbers",
        "msc.database.identifier",
        "msc.database.name",
        "msc.database.uri",
        "msc.expression.modality",
    )
    VALUE_MAP = {
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

    def project_sample(
        self, *, context: TabularMetadataContext
    ) -> TabularMetadataProjection:
        values = {
            "msc.sample.accession": context.sample_accession,
            "msc.series.accession": context.study_accession,
            **{
                column: context.base_metadata.get(key)
                for column, key in self.VALUE_MAP.items()
            },
            "msc.expression.modality": context.base_metadata.get("modality"),
        }
        characteristics = context.base_metadata.get("characteristics", {})
        if isinstance(characteristics, Mapping):
            values.update(
                {
                    f"msc.characteristics.{key}": value
                    for key, value in characteristics.items()
                }
            )
        harmonization = context.base_metadata.get("harmonization", ())
        for item in harmonization if isinstance(harmonization, Sequence) else ():
            if not isinstance(item, Mapping):
                continue
            prefix = f"msc.harmonization.{item['destination']}"
            values.update({
                f"{prefix}.value": item.get("value"),
                f"{prefix}.id": item.get("identifier"),
                f"{prefix}.ontology": item.get("ontology"),
                f"{prefix}.source_field": item.get("source_field"),
                f"{prefix}.hierarchy_depth": item.get("hierarchy_depth"),
            })
        values.update(
            patch_provenance_columns(
                context.package,
                context.sample,
                occupied=set(values),
            )
        )
        values.update({
            key: "; ".join(str(item) for item in items)
            for key, items in _parameter_summary(context.package, context.sample).items()
        })
        return TabularMetadataProjection(values=values, columns=self.COLUMNS)


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


class json2tsv(JSON2TSVConverter):
    pass
