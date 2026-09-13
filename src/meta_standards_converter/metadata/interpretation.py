# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Format-neutral MINiML sample metadata projection services."""

from __future__ import annotations

import json
import re
from typing import Any, Mapping, Protocol, Sequence

from meta_standards_converter.metadata.projection.assay import (
    _material_types_for_sample,
    _protocols_for_sample,
)
from meta_standards_converter.metadata.ontology_mappings import Harmonizer
from meta_standards_converter.miniml import harmonized_value_mappings


class MINiMLMetadataProvider(Protocol):
    """Public dependency contract shared by tabular and AnnData exporters."""

    def study_accession(self, packages: Sequence[Mapping[str, Any]]) -> str | None:
        """Return the first GEO series accession represented by ``packages``."""

    def samples(self, package: Mapping[str, Any]) -> tuple[Mapping[str, Any], ...]:
        """Return mapping-valued samples from one MINiML package."""

    def sample_accession(self, sample: Mapping[str, Any]) -> str | None:
        """Return the canonical sample accession when present."""

    def sample_metadata(
        self,
        sample: Mapping[str, Any],
        package: Mapping[str, Any],
    ) -> dict[str, Any]:
        """Render canonical sample metadata for an output adapter."""

    def sample_modality(self, sample: Mapping[str, Any]) -> str:
        """Classify the expression modality using the established policy."""

    def sample_metadata_values(
        self,
        sample: Mapping[str, Any],
        package: Mapping[str, Any],
    ) -> dict[str, Any]:
        """Return lossless tuple-valued metadata for AnnData transport."""

    def render_sample_metadata(
        self, values: Mapping[str, Any]
    ) -> dict[str, Any]:
        """Render lossless values into the canonical delimited representation."""

    def metadata_slug(self, value: Any) -> str:
        """Return the established canonical metadata column slug."""

    def values(self, values: Any) -> list[str]:
        """Flatten and case-insensitively de-duplicate metadata values."""

    def join_values(self, values: Any) -> str:
        """Render flattened values using the canonical delimiter."""

    def metadata_database(
        self, package: Mapping[str, Any]
    ) -> Mapping[str, Any]:
        """Return the package's primary metadata database declaration."""

    def platform_accession_values(
        self,
        sample: Mapping[str, Any],
        package: Mapping[str, Any],
    ) -> list[str]:
        """Resolve platform references to platform accessions."""

    def text(self, value: Any) -> str | None:
        """Normalize scalar metadata whitespace."""


class MINiMLMetadataService:
    """Own canonical, output-format-neutral MINiML metadata interpretation."""

    PROTOCOL_PATHS = (
        ("treatment_protocol", "Treatment-Protocol", "channel"),
        ("growth_protocol", "Growth-Protocol", "channel"),
        ("extract_protocol", "Extract-Protocol", "channel"),
        ("label_protocol", "Label-Protocol", "channel"),
        ("hybridization_protocol", "Hybridization-Protocol", "sample"),
        ("scan_protocol", "Scan-Protocol", "sample"),
        ("data_processing", "Data-Processing", "sample"),
    )

    def study_accession(
        self, packages: Sequence[Mapping[str, Any]]
    ) -> str | None:
        for package in packages:
            series = package.get("series")
            if not isinstance(series, Mapping):
                continue
            if package.get("source", {}).get("format") in {"SRA", "ENA"} and series.get("iid"):
                return str(series["iid"])
            for accession in self._as_list(series.get("accession")):
                value = self._value(accession)
                if isinstance(value, str) and value.upper().startswith("GSE"):
                    return value.upper()
        return None

    def samples(self, package: Mapping[str, Any]) -> tuple[Mapping[str, Any], ...]:
        return tuple(
            sample
            for sample in self._as_list(package.get("sample"))
            if isinstance(sample, Mapping)
        )

    def sample_accession(self, sample: Mapping[str, Any]) -> str | None:
        iid = sample.get("iid")
        if isinstance(iid, str) and iid.startswith(("SRS", "ERS", "DRS", "SAMN", "SAMEA", "SAMD")):
            return iid
        for accession in self._as_list(sample.get("accession")):
            value = self._value(accession)
            if isinstance(value, str) and value.upper().startswith("GSM"):
                return value.upper()
        value = sample.get("iid")
        return str(value) if value else None

    def sample_metadata(
        self,
        sample: Mapping[str, Any],
        package: Mapping[str, Any],
    ) -> dict[str, Any]:
        return self.render_sample_metadata(
            self.sample_metadata_values(sample, package)
        )

    def sample_metadata_values(
        self,
        sample: Mapping[str, Any],
        package: Mapping[str, Any],
    ) -> dict[str, Any]:
        metadata: dict[str, Any] = {
            "title": tuple(self.values(sample.get("title"))),
            "description": tuple(self.values(sample.get("description"))),
        }
        channels = [
            item
            for item in self._as_list(sample.get("channel"))
            if isinstance(item, Mapping)
        ]
        metadata["source"] = tuple(
            self.values(channel.get("source") for channel in channels)
        )
        organisms = [
            organism
            for channel in channels
            for organism in self._as_list(channel.get("organism"))
        ]
        organism_values: list[str] = []
        for channel in channels:
            harmonized: list[str] = []
            for organism in self._as_list(channel.get("organism")):
                if not isinstance(organism, Mapping):
                    continue
                for annotation in harmonized_value_mappings(organism):
                    if (
                        isinstance(annotation, Mapping)
                        and annotation.get("field") in {"organism", "species_name"}
                    ):
                        harmonized.extend(self.values(annotation.get("value")))
            organism_values.extend(
                harmonized or self.values(channel.get("organism"))
            )
        metadata["organism"] = tuple(self.values(organism_values))
        metadata["organism_taxid"] = tuple(
            self.values(
                organism.get("taxid")
                for organism in organisms
                if isinstance(organism, Mapping)
            )
        )

        characteristic_values: dict[str, list[str]] = {}
        for channel in channels:
            for annotation in harmonized_value_mappings(channel):
                if not isinstance(annotation, Mapping) or not annotation.get("field"):
                    continue
                annotation_slug = "hz_" + self.metadata_slug(
                    annotation["field"]
                )
                characteristic_values.setdefault(annotation_slug, []).extend(
                    self.values(annotation.get("value"))
                )
                if annotation.get("term_accession_number"):
                    characteristic_values.setdefault(
                        f"{annotation_slug}_id", []
                    ).extend(self.values(annotation["term_accession_number"]))
                if annotation.get("term_source_ref"):
                    characteristic_values.setdefault(
                        f"{annotation_slug}_onto", []
                    ).extend(self.values(annotation["term_source_ref"]))
            characteristic_rows = self._as_list(channel.get("characteristics"))
            for item in characteristic_rows:
                if not isinstance(item, Mapping) or not item.get(
                    "name", item.get("tag")
                ):
                    continue
                if str(item.get("name", item.get("tag"))).startswith("hz_"):
                    continue
                slug = self.metadata_slug(item.get("name", item.get("tag")))
                item_values = self.values(item.get("value"))
                if slug and item_values:
                    characteristic_values.setdefault(slug, []).extend(item_values)
            for annotation in harmonized_value_mappings(characteristic_rows):
                annotation_slug = "hz_" + self.metadata_slug(
                    annotation["field"]
                )
                characteristic_values.setdefault(annotation_slug, []).extend(
                    self.values(annotation.get("value"))
                )
                if annotation.get("term_accession_number"):
                    characteristic_values.setdefault(
                        f"{annotation_slug}_id", []
                    ).extend(self.values(annotation["term_accession_number"]))
                if annotation.get("term_source_ref"):
                    characteristic_values.setdefault(
                        f"{annotation_slug}_onto", []
                    ).extend(self.values(annotation["term_source_ref"]))
        characteristics = {
            slug: tuple(self.values(items))
            for slug, items in characteristic_values.items()
        }
        metadata["characteristics"] = characteristics
        metadata["organism_part"] = (
            characteristics.get("organism_part")
            or characteristics.get("tissue")
            or metadata["source"]
        )
        metadata["developmental_stage"] = characteristics.get(
            "developmental_stage", ()
        )
        metadata["disease"] = characteristics.get("disease", ())
        metadata["genotype"] = characteristics.get("genotype", ())

        metadata["biomaterial_provider"] = tuple(
            self.values(channel.get("biomaterial_provider") for channel in channels)
        )
        metadata["molecule"] = tuple(
            self.values(channel.get("molecule") for channel in channels)
        )
        assay_material_types = self.values(
            _material_types_for_sample(package, sample)
        )
        characteristic_material_types = self.values(
            characteristics.get("material_type")
        )
        legacy_material_types = self.values(
            channel.get("material_type") for channel in channels
        )
        material_types = [
            re.sub(r"^total\s+", "", value, flags=re.IGNORECASE)
            for value in self.values(
                channel.get("molecule") for channel in channels
            )
        ]
        metadata["material_type"] = (
            tuple(assay_material_types)
            or tuple(characteristic_material_types)
            or tuple(legacy_material_types)
            or tuple(self.values(material_types))
            or metadata["organism_part"]
        )

        runs = [
            item
            for item in self._as_list(sample.get("sra_run"))
            if isinstance(item, Mapping)
        ]
        metadata["sra_accession"] = tuple(
            self.values(sample.get("sra_accession"))
        )
        metadata["ena_accession"] = tuple(
            self.values(sample.get("ena_accession"))
        )
        metadata["biosample_accession"] = tuple(
            self.values(run.get("biosample") for run in runs)
        )
        metadata["sra_run_accessions"] = tuple(
            self.values(run.get("run") for run in runs)
        )
        metadata["library_strategy"] = tuple(
            self.values(
                [
                    sample.get("library_strategy"),
                    *(run.get("library_strategy") for run in runs),
                ]
            )
        )
        metadata["library_source"] = tuple(
            self.values(
                [
                    sample.get("library_source"),
                    *(run.get("library_source") for run in runs),
                ]
            )
        )
        metadata["library_selection"] = tuple(
            self.values(
                [
                    sample.get("library_selection"),
                    *(run.get("library_selection") for run in runs),
                ]
            )
        )
        metadata["library_layout"] = tuple(
            self.values(run.get("library_layout") for run in runs)
        )
        metadata["instrument_model"] = tuple(
            self.values(
                [
                    sample.get("instrument_model"),
                    *(run.get("instrument_model") for run in runs),
                ]
            )
        )
        metadata["platform_accession"] = tuple(
            self.platform_accession_values(sample, package)
        )

        protocol_types: list[str] = []
        protocol_sources: list[str] = []
        protocol_accessions: list[str] = []
        typed_protocols = _protocols_for_sample(package, sample)
        if typed_protocols:
            for protocol in typed_protocols:
                protocol_type = protocol.get("type")
                if isinstance(protocol_type, Mapping):
                    protocol_types.extend(self.values(protocol_type.get("value")))
                    protocol_sources.extend(
                        self.values(protocol_type.get("term_source_ref"))
                    )
                    protocol_accessions.extend(
                        self.values(protocol_type.get("term_accession_number"))
                    )
                else:
                    protocol_types.extend(self.values(protocol_type))
        else:
            for field, label, scope in self.PROTOCOL_PATHS:
                containers = channels if scope == "channel" else [sample]
                if not self.join_values(
                    container.get(field) for container in containers
                ):
                    continue
                protocol_type, source_ref, accession = Harmonizer().geoprotocols2efo(
                    label
                )
                protocol_types.append(protocol_type)
                protocol_sources.append(source_ref)
                protocol_accessions.append(accession)
        metadata["protocol_types"] = tuple(self.values(protocol_types))
        metadata["protocol_term_source_refs"] = tuple(
            self.values(protocol_sources)
        )
        metadata["protocol_term_accession_numbers"] = tuple(
            self.values(protocol_accessions)
        )

        database = self.metadata_database(package)
        metadata["metadata_source"] = tuple(
            self.values(
                database.get("public_id")
                or database.get("iid")
                or database.get("name")
            )
        )
        metadata["metadata_source_name"] = tuple(
            self.values(database.get("name"))
        )
        metadata["metadata_source_uri"] = tuple(
            self.values(database.get("web_link"))
        )
        return metadata

    def render_sample_metadata(self, values: Mapping[str, Any]) -> dict[str, Any]:
        return {
            key: (
                {
                    characteristic: self.join_values(items)
                    for characteristic, items in value.items()
                }
                if key == "characteristics"
                else self.join_values(value)
            )
            for key, value in values.items()
        }

    def sample_modality(self, sample: Mapping[str, Any]) -> str:
        text = json.dumps(sample).lower()
        if any(
            value in text
            for value in (
                "single cell",
                "single-cell",
                "10x",
                "chromium",
                "visium",
            )
        ):
            return "single_cell"
        if any(
            sample.get(key)
            for key in ("library_source", "library_strategy", "type", "sra_run")
        ):
            return "bulk"
        return "unknown"

    def metadata_slug(self, value: Any) -> str:
        rendered = self.join_values(value).lower()
        return re.sub(
            r"_+", "_", re.sub(r"[^a-z0-9]+", "_", rendered)
        ).strip("_")

    def values(self, values: Any) -> list[str]:
        flattened: list[str] = []
        seen: set[str] = set()

        def visit(value: Any) -> None:
            if value is None:
                return
            if isinstance(value, Mapping):
                for key in (
                    "value",
                    "name",
                    "predefined",
                    "public_id",
                    "iid",
                    "ref",
                ):
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
            cleaned = self.text(value)
            if not cleaned:
                return
            normalized = cleaned.casefold()
            if normalized not in seen:
                flattened.append(cleaned)
                seen.add(normalized)

        visit(values)
        return flattened

    def join_values(self, values: Any) -> str:
        return "; ".join(self.values(values))

    def metadata_database(self, package: Mapping[str, Any]) -> Mapping[str, Any]:
        return next(
            (
                item
                for item in self._as_list(package.get("database"))
                if isinstance(item, Mapping)
            ),
            {},
        )

    def platform_accession_values(
        self,
        sample: Mapping[str, Any],
        package: Mapping[str, Any],
    ) -> list[str]:
        references = set(self.values(sample.get("platform_ref")))
        values: list[str] = []
        for platform in self._as_list(package.get("platform")):
            if not isinstance(platform, Mapping):
                continue
            identifiers = {
                platform.get("iid"),
                *self.values(platform.get("accession")),
            }
            if references and not references.intersection(
                identifier for identifier in identifiers if identifier
            ):
                continue
            values.extend(self.values(platform.get("accession")))
        return self.values(values or references)

    @staticmethod
    def text(value: Any) -> str | None:
        if isinstance(value, Mapping):
            value = value.get("value") or value.get("name")
        if value is None:
            return None
        return " ".join(str(value).split()) or None

    @staticmethod
    def _value(value: Any) -> Any:
        if isinstance(value, Mapping):
            return value.get("value") or value.get("uri") or value.get("filename")
        return value

    @staticmethod
    def _as_list(value: Any) -> list[Any]:
        if value is None:
            return []
        return value if isinstance(value, list) else [value]
