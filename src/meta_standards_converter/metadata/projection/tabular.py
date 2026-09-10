# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from __future__ import annotations
from dataclasses import dataclass
from typing import Any, Mapping, Protocol, Sequence
from meta_standards_converter.metadata.provenance import patch_provenance_columns
from meta_standards_converter.metadata.projection.assay import _parameter_summary


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
