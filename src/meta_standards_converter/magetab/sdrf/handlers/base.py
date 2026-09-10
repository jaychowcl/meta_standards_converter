# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from __future__ import annotations
from meta_standards_converter.magetab.sdrf.renderer import SDRFRenderer
from dataclasses import dataclass, field

from collections import OrderedDict

import requests

import xml.etree.ElementTree as ET

from meta_standards_converter.magetab.protocols import ProtocolRegistry
from meta_standards_converter.magetab.technology import (
    detect_ae_technology,
    has_array_files,
    normalized_extension,
)

from meta_standards_converter.sources.insdc import INSDCWebfetcher

from meta_standards_converter.helpers.json_helper import JSONHandler

from meta_standards_converter.magetab.sdrf.model import SDRFAttr
from meta_standards_converter.magetab.sdrf.model import SDRFNode
from meta_standards_converter.magetab.sdrf.model import SDRFEdge
from meta_standards_converter.magetab.sdrf.model import SDRFPath
from meta_standards_converter.magetab.sdrf.model import ColumnGroup
from meta_standards_converter.magetab.sdrf.model import SDRFAudit

def classify_file(path: str) -> str:
    ext = normalized_extension(path)
    if ext in {".fastq", ".fq", ".bam", ".sam", ".cram"}:
        return "sequencing_raw"
    if ext in {".cel", ".gpr", ".idat", ".exp", ".rpt", ".cab", ".tif", ".tiff"}:
        return "array_raw"
    if ext in {".txt", ".tsv", ".csv", ".mtx", ".h5", ".h5ad"}:
        return "matrix_or_derived"
    return "supplementary"


class _BaseSDRFHandler(SDRFRenderer):
    def __init__(self, parent, data: dict, protocol_registry=None):
        self.parent = parent
        self.data = data
        self.insdc_handler = getattr(parent, "insdc_fetcher", None) or INSDCWebfetcher()
        self.sra_cache = {}
        self.samples = [x for x in self._as_list(data.get("sample")) if isinstance(x, dict)]
        self.platforms = {
            platform.get("iid"): platform
            for platform in self._as_list(data.get("platform"))
            if isinstance(platform, dict) and platform.get("iid")
        }
        self.sample_by_id = {
            sample.get("iid"): sample
            for sample in self.samples
            if sample.get("iid")
        }
        self.series_accession = self._series_accession()
        if protocol_registry is None:
            protocol_registry = ProtocolRegistry(series_accession=self.series_accession)
        self.protocol_registry = protocol_registry
        self.audit = SDRFAudit()
        self.factor_tags = self._factor_tags()

    def _as_list(self, value):
        if value is None:
            return []
        if isinstance(value, list):
            return value
        return [value]

    def build(self) -> list:
        self.preregister_protocols()
        paths = self.build_paths()
        columns = self.plan_columns(paths=paths)
        return self.render_paths(columns=columns, paths=paths)

    def preregister_protocols(self) -> None:
        for sample in self.ordered_samples():
            self.preregister_data_processing(sample=sample)

    def build_paths(self) -> list[SDRFPath]:
        paths = []
        for sample in self.ordered_samples():
            for channel in self.channels(sample=sample):
                path = SDRFPath(parts=[self.source_node(sample=sample, channel=channel)])
                path.parts.append(self.sample_collection_edge())
                processing_edge = self.data_processing_edge(sample=sample)
                if processing_edge:
                    path.parts.append(processing_edge)
                path.parts.extend(self.factor_nodes(sample=sample, channel=channel))
                paths.append(path)
        return paths or [SDRFPath(parts=[SDRFNode(kind="Source Name", key="source", value=None)])]


    def ordered_samples(self) -> list:
        refs = [
            ref.get("ref")
            for series in self._as_list(self.data.get("series"))
            if isinstance(series, dict)
            for ref in self._as_list(series.get("sample_ref"))
            if isinstance(ref, dict) and ref.get("ref")
        ]
        ordered = [self.sample_by_id[ref] for ref in refs if ref in self.sample_by_id]
        seen = {sample.get("iid") for sample in ordered}
        ordered.extend(sample for sample in self.samples if sample.get("iid") not in seen)
        return ordered

    def channels(self, sample: dict) -> list[dict]:
        channels = [channel for channel in self._as_list(sample.get("channel")) if isinstance(channel, dict)]
        if len(channels) > 1:
            self.audit.warnings.append(f"Sample {self.sample_accession(sample=sample)} has {len(channels)} channels; emitted {len(channels)} channel paths.")
        return channels or [{}]

    def source_node(self, sample: dict, channel: dict) -> SDRFNode:
        accession = self.sample_accession(sample=sample)
        value = accession
        attrs = []
        attrs.extend(SDRFAttr(label="Comment[BioSD_SAMPLE]", value=x) for x in self.biosample_accessions(sample=sample))
        attrs.extend(self.sample_comment_attrs(sample=sample, channel=channel))
        attrs.extend(self.characteristic_attrs(channel=channel))
        provider = self.provider(channel=channel)
        if provider:
            attrs.append(SDRFAttr(label="Provider", value=provider))
        attrs.append(SDRFAttr(label="Material Type", value=self.material_type(channel=channel)))
        # Legacy greedy GEO fallback:
        # attrs.extend(_GEOFallbackComments(self).geo_fallback_attrs(sample=sample, channel=channel))
        attrs.extend(self.extra_source_attrs(sample=sample, channel=channel))
        return SDRFNode(kind="Source Name", key=f"source:{accession}", value=value, attrs=attrs)

    def sample_comment_attrs(self, sample: dict, channel: dict) -> list[SDRFAttr]:
        attrs = []
        for label, value in (
            ("Comment[Sample_description]", sample.get("description")),
            ("Comment[Sample_source_name]", channel.get("source")),
            ("Comment[Sample_title]", sample.get("title")),
        ):
            value = self.clean(value)
            if value:
                attrs.append(SDRFAttr(label=label, value=value))
        return attrs

    def characteristic_attrs(self, channel: dict) -> list[SDRFAttr]:
        required_attrs = OrderedDict(
            (
                (tag, SDRFAttr(label=f"Characteristics[{tag}]", value=None, required=True))
                for tag in ("organism", "organism part", "developmental stage", "disease", "genotype")
            )
        )
        extra_attrs = []

        organism_values = []
        for organism in self._as_list(channel.get("organism")):
            organism_value = self.organism_value(organism=organism)
            if organism_value:
                organism_values.append(organism_value)
        if organism_values:
            required_attrs["organism"].value = organism_values[0]
            first_organism = next(
                (item for item in self._as_list(channel.get("organism")) if isinstance(item, dict)),
                {},
            )
            required_attrs["organism"].attrs = self.ontology_companions(first_organism)
            for organism_value in organism_values[1:]:
                extra_attrs.append(SDRFAttr(label="Characteristics[organism]", value=organism_value))

        organism_part_value = self.organism_part_value(channel=channel)
        required_attrs["organism part"].value = organism_part_value

        seen_tags = {}
        first_organism_part_preserved = False
        for characteristic in channel.get("characteristics", []) or []:
            if not isinstance(characteristic, dict):
                continue
            tag = characteristic.get("name") or characteristic.get("tag")
            if not tag:
                continue
            if str(tag).startswith("hz_"):
                continue
            if tag.lower() == "organism part":
                value = self.clean(characteristic.get("value"))
                if not first_organism_part_preserved and value == organism_part_value:
                    required_attrs["organism part"].attrs = self.ontology_companions(characteristic)
                    first_organism_part_preserved = True
                    continue
            seen_tags[tag] = seen_tags.get(tag, 0) + 1
            if seen_tags[tag] > 1:
                self.audit.warnings.append(f"Repeated Characteristics[{tag}] preserved as {seen_tags[tag]} columns.")
            lower_tag = tag.lower()
            value = self.clean(characteristic.get("value"))
            if lower_tag in required_attrs and required_attrs[lower_tag].value is None:
                required_attrs[lower_tag].value = value
                required_attrs[lower_tag].attrs = self.ontology_companions(characteristic)
                continue
            extra_attrs.append(SDRFAttr(
                label=f"Characteristics[{tag}]",
                value=value,
                attrs=self.ontology_companions(characteristic),
            ))
        return list(required_attrs.values()) + extra_attrs

    def ontology_companions(self, item: dict) -> list[SDRFAttr]:
        companions = []
        if item.get("term_source_ref"):
            companions.append(SDRFAttr(label="Term Source REF", value=self.clean(item["term_source_ref"])))
        if item.get("term_accession_number"):
            companions.append(SDRFAttr(label="Term Accession Number", value=self.clean(item["term_accession_number"])))
        return companions

    def organism_part_value(self, channel: dict):
        organism_part = self.characteristic_values(channel=channel, tag="organism part")
        if organism_part:
            return organism_part[0]

        tissue = self.characteristic_values(channel=channel, tag="tissue")
        if tissue:
            return tissue[0]

        return self.clean(channel.get("source"))

    def provider(self, channel: dict):
        providers = [self.clean(x) for x in self._as_list(channel.get("biomaterial_provider")) if self.clean(x)]
        return "; ".join(providers) if providers else None

    def organism_value(self, organism):
        if not isinstance(organism, dict):
            return None
        return self.clean(organism.get("name") or organism.get("value")) or None

    def material_type(self, channel: dict):
        if channel.get("material_type"):
            return self.clean(channel.get("material_type"))
        molecule = channel.get("molecule")
        if molecule:
            return self.clean(molecule).replace("total ", "")
        return "organism part"

    def factor_nodes(self, sample: dict, channel: dict) -> list[SDRFNode]:
        nodes = []
        for tag in self.factor_tags:
            value = self.factor_value(sample=sample, channel=channel, tag=tag)
            if value:
                nodes.append(SDRFNode(kind=f"Factor Value[{tag}]", key=f"factor:{tag}", value=value))
        nodes.extend(self.extra_factor_nodes(sample=sample, channel=channel))
        return nodes

    def extra_source_attrs(self, sample: dict, channel: dict) -> list[SDRFAttr]:
        return []

    def extra_factor_nodes(self, sample: dict, channel: dict) -> list[SDRFNode]:
        return []

    def extraction_edges(self, sample: dict, channel: dict) -> list[SDRFEdge]:
        edges = []
        for kind, text in (
            ("treatment", channel.get("treatment_protocol")),
            ("growth", channel.get("growth_protocol")),
        ):
            ref = self.protocol_registry.get_ref(kind=kind, text=text)
            if ref:
                edges.append(SDRFEdge(protocol_ref=ref))

        extract_ref = self.protocol_registry.get_ref(kind="extraction", text=channel.get("extract_protocol"))
        if not extract_ref:
            self.audit.warnings.append(f"Protocol text missing for extraction in sample {self.sample_accession(sample=sample)}; Protocol REF left blank.")
        edges.append(SDRFEdge(protocol_ref=extract_ref))
        return edges

    def preregister_extraction_protocols(self, sample: dict, channel: dict) -> None:
        for kind, text in (
            ("treatment", channel.get("treatment_protocol")),
            ("growth", channel.get("growth_protocol")),
            ("extraction", channel.get("extract_protocol")),
        ):
            self.protocol_registry.get_ref(kind=kind, text=text)

    def sample_collection_edge(self) -> SDRFEdge:
        ref = self.protocol_registry.ensure_required(
            kind="sample collection",
            label="Sample-Collection-Protocol",
        )
        return SDRFEdge(protocol_ref=ref)

    def nucleic_acid_sequencing_edge(self) -> SDRFEdge:
        ref = self.protocol_registry.ensure_required(
            kind="nucleic acid sequencing",
            label="Nucleic-Acid-Sequencing-Protocol",
        )
        return SDRFEdge(protocol_ref=ref)

    def protocol_edge(self, kind: str, text: str | None, sample: dict, required: bool = False) -> SDRFEdge | None:
        ref = self.protocol_registry.get_ref(kind=kind, text=text)
        if not ref and required:
            self.audit.warnings.append(f"Protocol text missing for {kind} in sample {self.sample_accession(sample=sample)}; Protocol REF left blank.")
            return SDRFEdge(protocol_ref=None)
        if not ref:
            return None
        return SDRFEdge(protocol_ref=ref)

    def preregister_data_processing(self, sample: dict) -> None:
        self.protocol_registry.get_ref(
            kind="data processing",
            text=sample.get("data_processing"),
        )

    def data_processing_edge(self, sample: dict) -> SDRFEdge | None:
        return self.protocol_edge(
            kind="data processing",
            text=sample.get("data_processing"),
            sample=sample,
        )

    def sample_accession(self, sample: dict):
        for accession in self._as_list(sample.get("accession")):
            if isinstance(accession, dict) and accession.get("value"):
                return self.clean(accession.get("value"))
        return self.clean(sample.get("iid"))

    def biosample_accessions(self, sample: dict) -> list:
        values = []
        for relation in self._as_list(sample.get("relation")):
            if not isinstance(relation, dict):
                continue
            if (relation.get("type") or "").lower() != "biosample":
                continue
            target = relation.get("target") or ""
            values.append(self.clean(target.rstrip("/").split("/")[-1] if "/" in target else target))
        return [x for x in values if x]

    def biosample_accession(self, sample: dict):
        accessions = self.biosample_accessions(sample=sample)
        return accessions[0] if accessions else None

    def platform(self, sample: dict) -> dict:
        platform_ref = sample.get("platform_ref") or {}
        return self.platforms.get(platform_ref.get("ref"), {})

    def platform_accession(self, sample: dict):
        platform = self.platform(sample=sample)
        for accession in self._as_list(platform.get("accession")):
            if isinstance(accession, dict) and accession.get("value"):
                return self.clean(accession.get("value"))
        return self.clean(platform.get("iid"))

    def instrument_model(self, sample: dict):
        instrument_model = sample.get("instrument_model") or {}
        if isinstance(instrument_model, dict):
            return self.clean(instrument_model.get("predefined") or instrument_model.get("other"))
        return self.clean(instrument_model)

    def supplementary_files(self, sample: dict) -> list:
        files = []
        for key in ("supplementary_data", "raw_data"):
            files.extend(
                data_file.get("value")
                for data_file in sample.get(key, []) or []
                if isinstance(data_file, dict) and data_file.get("value")
            )

        platform = self.platform(sample=sample)
        files.extend(
            data_file.get("value")
            for data_file in platform.get("supplementary_data", []) or []
            if isinstance(data_file, dict) and data_file.get("value")
        )
        for series in self._as_list(self.data.get("series")):
            if not isinstance(series, dict):
                continue
            files.extend(
                data_file.get("value")
                for data_file in series.get("supplementary_data", []) or []
                if isinstance(data_file, dict) and data_file.get("value")
            )
        return [self.clean(x) for x in files if self.clean(x)]

    def raw_files(self, sample: dict) -> list:
        return [
            data_file
            for data_file in self.supplementary_files(sample=sample)
            if classify_file(data_file) in {"sequencing_raw", "array_raw"}
        ]

    def derived_files(self, sample: dict) -> list:
        return [
            data_file
            for data_file in self.supplementary_files(sample=sample)
            if classify_file(data_file) not in {"sequencing_raw", "array_raw"}
        ]

    def arrayexpress_ftp(self, value):
        value = self.clean(value)
        if value and value.startswith(("ftp://", "http://", "https://")):
            return value
        return None

    def file_node(self, kind: str, value: str | None, attrs: list[SDRFAttr] | None = None) -> SDRFNode:
        attrs = attrs or []
        ftp = self.arrayexpress_ftp(value)
        if ftp and not any(attr.label in {"Comment[ArrayExpress FTP file]", "Comment[Derived ArrayExpress FTP file]"} for attr in attrs):
            label = "Comment[Derived ArrayExpress FTP file]" if kind.startswith("Derived") else "Comment[ArrayExpress FTP file]"
            attrs.append(SDRFAttr(label=label, value=ftp))
        return SDRFNode(kind=kind, key=f"file:{kind}:{value}", value=value, attrs=attrs)

    def sra_runs(self, sample: dict) -> list:
        if "sra_run" in sample:
            return [
                run
                for run in self._as_list(sample.get("sra_run"))
                if isinstance(run, dict)
            ]

        accessions = []
        for relation in self._as_list(sample.get("relation")):
            if not isinstance(relation, dict):
                continue
            if (relation.get("type") or "").lower() != "sra":
                continue
            accessions.extend(self.insdc_handler.extract_sra_accessions(relation.get("target") or ""))

        runs = []
        for accession in dict.fromkeys(accessions):
            if accession not in self.sra_cache:
                self.sra_cache[accession] = self.parent._lookup_sra(sra=accession)
            runs.extend(self.sra_cache[accession])
        return runs

    def characteristic_values(self, channel: dict, tag: str) -> list:
        values = []
        lower_tag = tag.lower()
        if lower_tag == "organism":
            for organism in self._as_list(channel.get("organism")):
                organism_value = self.organism_value(organism=organism)
                if organism_value:
                    values.append(organism_value)
        for characteristic in channel.get("characteristics", []) or []:
            if not isinstance(characteristic, dict):
                continue
            name = characteristic.get("name") or characteristic.get("tag") or ""
            if str(name).startswith("hz_"):
                continue
            if str(name).lower() == lower_tag:
                values.append(self.clean(characteristic.get("value")))
        return [x for x in values if x]

    def characteristic_value(self, sample: dict, tag: str):
        for channel in self.channels(sample=sample):
            values = self.characteristic_values(channel=channel, tag=tag)
            if values:
                return values[0]
        return None

    def clean(self, value):
        if value is None:
            return None
        if isinstance(value, dict) and "value" in value:
            value = value["value"]
        return " ".join(str(value).replace("\t", " ").replace("\n", " ").split())

    def _series_accession(self):
        for series in self._as_list(self.data.get("series")):
            if not isinstance(series, dict):
                continue
            for accession in self._as_list(series.get("accession")):
                if isinstance(accession, dict) and accession.get("value"):
                    return accession.get("value")
        return "GEO"

    def _factor_tags(self) -> list:
        variable_tags = []
        for series in self._as_list(self.data.get("series")):
            if not isinstance(series, dict):
                continue
            for variable in self._as_list(series.get("variable")):
                if not isinstance(variable, dict):
                    continue
                tag = variable.get("factor") or variable.get("name") or variable.get("tag")
                if tag and tag not in variable_tags:
                    variable_tags.append(tag)
        if variable_tags:
            return variable_tags

        values_by_tag = {}
        for sample in self.samples:
            for channel in self.channels(sample=sample):
                for attr in self.characteristic_attrs(channel=channel):
                    if not attr.label.startswith("Characteristics[") or not attr.value:
                        continue
                    tag = attr.label.removeprefix("Characteristics[").removesuffix("]")
                    values_by_tag.setdefault(tag, set()).add(attr.value)

        return [
            tag
            for tag, values in values_by_tag.items()
            if tag != "organism" and len(values) > 1
        ]

    def factor_value(self, sample: dict, channel: dict, tag: str):
        values = self.characteristic_values(channel=channel, tag=tag)
        return values[0] if values else None
