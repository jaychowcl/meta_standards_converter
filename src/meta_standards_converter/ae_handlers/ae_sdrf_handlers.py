# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
'''
Constructor class for ae MAGETAB sdrf
'''
from dataclasses import dataclass, field
from collections import OrderedDict
from urllib.parse import urlparse
import os
import requests
import xml.etree.ElementTree as ET

from meta_standards_converter.insdc_handlers.insdc_webfetcher import INSDCWebfetcher
from meta_standards_converter.helpers.json_helper import JSONHandler


@dataclass
class SDRFAttr:
    label: str
    value: str | None
    attrs: list["SDRFAttr"] = field(default_factory=list)
    required: bool = False


@dataclass
class SDRFNode:
    kind: str
    key: str
    value: str | None
    attrs: list[SDRFAttr] = field(default_factory=list)


@dataclass
class SDRFEdge:
    protocol_ref: str | None
    attrs: list[SDRFAttr] = field(default_factory=list)


@dataclass
class SDRFPath:
    parts: list[SDRFNode | SDRFEdge] = field(default_factory=list)


@dataclass
class ColumnGroup:
    main_key: str
    main_label: str
    companions: list["ColumnGroup"] = field(default_factory=list)


@dataclass
class SDRFAudit:
    warnings: list[str] = field(default_factory=list)
    dropped_values: list[str] = field(default_factory=list)
    validation_errors: list[str] = field(default_factory=list)


class SDRFConstructor():
    def __init__(self, insdc_fetcher=None):
        self.insdc_fetcher = insdc_fetcher or INSDCWebfetcher()

    def _add_sdrf_to_idf(self, idf: list, data: dict) -> list:
        """
        Appends the generated SDRF payload to an IDF row list.
        """
        idf.append(["SDRF File", self._miniml2sdrf(data=data)])
        return idf

    def _miniml2sdrf(self, data: dict, protocol_registry=None, technology_type=None):
        """
        converts miniml json to magetab sdrf.
        """
        tech_type = technology_type or self._detect_sdrf_technology(data=data)
        handler_class = {
            "plate_single_cell_sequencing": _PlateSingleCellSequencingSDRFHandler,
            "droplet_single_cell_sequencing": _DropletSingleCellSequencingSDRFHandler,
            "tenx_v2_droplet_single_cell_sequencing": _TenXV2DropletSingleCellSequencingSDRFHandler,
            "tenx_v3_droplet_single_cell_sequencing": _TenXV3DropletSingleCellSequencingSDRFHandler,
            "single_cell_sequencing": _SingleCellSequencingSDRFHandler,
            "spatial_sequencing": _SpatialSequencingSDRFHandler,
            "bulk_sequencing": _BulkSequencingSDRFHandler,
            "sequencing": _SequencingSDRFHandler,
            "array": _ArraySDRFHandler,
        }.get(tech_type, _GenericSDRFHandler)

        handler = handler_class(parent=self, data=data, protocol_registry=protocol_registry)
        sdrf = handler.build()
        self.last_sdrf_audit = handler.audit
        return sdrf

    def _detect_sdrf_technology(self, data: dict) -> str:
        from meta_standards_converter.ae_handlers.ae_constructor import AEConstructor

        return AEConstructor()._detect_ae_technology(data=data)

    def _has_array_files(self, data: dict) -> bool:
        from meta_standards_converter.ae_handlers.ae_constructor import AEConstructor

        return AEConstructor()._has_array_files(data=data)

    def _lookup_sra(self, sra: str) -> list:
        '''
        take sra accession to return srr info
        '''
        try:
            return self.insdc_fetcher.fetch_sra_runs(accession=sra)
        except (requests.RequestException, ET.ParseError):
            return []


def normalized_extension(path: str) -> str:
    parsed = urlparse(str(path))
    basename = os.path.basename(parsed.path or str(path)).lower()
    for suffix in (".gz", ".zip", ".bz2", ".xz"):
        if basename.endswith(suffix):
            basename = basename[:-len(suffix)]
            break
    return os.path.splitext(basename)[1]


def classify_file(path: str) -> str:
    ext = normalized_extension(path)
    if ext in {".fastq", ".fq", ".bam", ".sam", ".cram"}:
        return "sequencing_raw"
    if ext in {".cel", ".gpr", ".idat", ".exp", ".rpt", ".cab", ".tif", ".tiff"}:
        return "array_raw"
    if ext in {".txt", ".tsv", ".csv", ".mtx", ".h5", ".h5ad"}:
        return "matrix_or_derived"
    return "supplementary"


class _BaseSDRFHandler():
    def __init__(self, parent: SDRFConstructor, data: dict, protocol_registry=None):
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
            from meta_standards_converter.ae_handlers.ae_constructor import ProtocolRegistry

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

    def plan_columns(self, paths: list[SDRFPath]) -> list[ColumnGroup]:
        columns = OrderedDict()
        for path in paths:
            for group, values in self.path_groups(path=path):
                self.merge_column_group(columns=columns, group=group)
        return list(columns.values())

    def merge_column_group(self, columns: OrderedDict, group: ColumnGroup) -> None:
        if group.main_key not in columns:
            columns[group.main_key] = ColumnGroup(
                main_key=group.main_key,
                main_label=group.main_label,
                companions=[],
            )
        target = columns[group.main_key]
        companion_columns = OrderedDict((x.main_key, x) for x in target.companions)
        for companion in group.companions:
            self.merge_column_group(columns=companion_columns, group=companion)
        target.companions = list(companion_columns.values())

    def render_paths(self, columns: list[ColumnGroup], paths: list[SDRFPath]) -> list:
        header = []
        for column in columns:
            header.extend(self.column_labels(column=column))

        rows = [header]
        for path in paths:
            values = {}
            for group, group_values in self.path_groups(path=path):
                values.update(group_values)
            row = []
            for column in columns:
                row.extend(self.column_values(column=column, values=values))
            rows.append(row)
        return rows

    def column_labels(self, column: ColumnGroup) -> list:
        labels = [column.main_label]
        for companion in column.companions:
            labels.extend(self.column_labels(column=companion))
        return labels

    def column_values(self, column: ColumnGroup, values: dict) -> list:
        rendered = [self.render_value(values.get(column.main_key))]
        for companion in column.companions:
            rendered.extend(self.column_values(column=companion, values=values))
        return rendered

    def path_groups(self, path: SDRFPath) -> list[tuple[ColumnGroup, dict]]:
        groups = []
        counts = {}
        for part in path.parts:
            if isinstance(part, SDRFEdge):
                key = self.occurrence_key(counts=counts, label="Protocol REF")
                groups.append(self.group_with_values(
                    key=key,
                    label="Protocol REF",
                    value=part.protocol_ref,
                    attrs=part.attrs,
                ))
                continue

            key = self.occurrence_key(counts=counts, label=part.kind)
            groups.append(self.group_with_values(
                key=key,
                label=part.kind,
                value=part.value,
                attrs=part.attrs,
            ))
        return groups

    def group_with_values(self, key: str, label: str, value, attrs: list[SDRFAttr]) -> tuple[ColumnGroup, dict]:
        companions, values = self.attr_columns(parent_key=key, attrs=attrs)
        values[key] = value
        return ColumnGroup(main_key=key, main_label=label, companions=companions), values

    def attr_columns(self, parent_key: str, attrs: list[SDRFAttr]) -> tuple[list[ColumnGroup], dict]:
        groups = []
        values = {}
        counts = {}
        for attr in attrs:
            if attr.value is None and not attr.attrs and not attr.required:
                continue
            key = f"{parent_key}|{self.occurrence_key(counts=counts, label=attr.label)}"
            companions, companion_values = self.attr_columns(parent_key=key, attrs=attr.attrs)
            groups.append(ColumnGroup(
                main_key=key,
                main_label=attr.label,
                companions=companions,
            ))
            values[key] = attr.value
            values.update(companion_values)
        return groups, values

    def occurrence_key(self, counts: dict, label: str) -> str:
        counts[label] = counts.get(label, 0) + 1
        return f"{label}#{counts[label]}"

    def render_value(self, value):
        if value is None:
            return ""
        return value

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
            for organism_value in organism_values[1:]:
                extra_attrs.append(SDRFAttr(label="Characteristics[organism]", value=organism_value))

        organism_part_value = self.organism_part_value(channel=channel)
        required_attrs["organism part"].value = organism_part_value

        seen_tags = {}
        first_organism_part_preserved = False
        for characteristic in channel.get("characteristics", []) or []:
            if not isinstance(characteristic, dict):
                continue
            tag = characteristic.get("tag")
            if not tag:
                continue
            if tag.lower() == "organism part":
                value = self.clean(characteristic.get("value"))
                if not first_organism_part_preserved and value == organism_part_value:
                    first_organism_part_preserved = True
                    continue
            seen_tags[tag] = seen_tags.get(tag, 0) + 1
            if seen_tags[tag] > 1:
                self.audit.warnings.append(f"Repeated Characteristics[{tag}] preserved as {seen_tags[tag]} columns.")
            lower_tag = tag.lower()
            value = self.clean(characteristic.get("value"))
            if lower_tag in required_attrs and required_attrs[lower_tag].value is None:
                required_attrs[lower_tag].value = value
                continue
            extra_attrs.append(SDRFAttr(label=f"Characteristics[{tag}]", value=value))
        return list(required_attrs.values()) + extra_attrs

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
            accessions.extend(self.insdc_handler._extract_sra(relation.get("target") or ""))

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
            if (characteristic.get("tag") or "").lower() == lower_tag:
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


class _SequencingSDRFHandler(_BaseSDRFHandler):
    def preregister_protocols(self) -> None:
        for sample in self.ordered_samples():
            self.preregister_data_processing(sample=sample)
            self.protocol_registry.get_ref(kind="scan", text=sample.get("scan_protocol"))
            runs = self.sra_runs(sample=sample) or [None]
            for channel in self.channels(sample=sample):
                self.preregister_extraction_protocols(sample=sample, channel=channel)
                for run in runs:
                    self.protocol_registry.get_ref(
                        kind="library construction",
                        text=self.library_protocol_text(sample=sample, channel=channel, run=run),
                    )

    def build_paths(self) -> list[SDRFPath]:
        paths = []
        for sample in self.ordered_samples():
            runs = self.sra_runs(sample=sample) or [None]
            for channel in self.channels(sample=sample):
                for run in runs:
                    path = SDRFPath(parts=[self.source_node(sample=sample, channel=channel)])
                    path.parts.append(self.sample_collection_edge())
                    path.parts.extend(self.extraction_edges(sample=sample, channel=channel))
                    path.parts.append(self.extract_node(sample=sample, channel=channel, run=run))
                    assay_edge = self.protocol_edge(
                        kind="library construction",
                        text=self.library_protocol_text(sample=sample, channel=channel, run=run),
                        sample=sample,
                        required=True,
                    )
                    if assay_edge:
                        path.parts.append(assay_edge)
                    path.parts.append(self.assay_node(sample=sample, run=run))
                    path.parts.append(self.nucleic_acid_sequencing_edge())
                    scan_edge = self.protocol_edge(
                        kind="scan",
                        text=sample.get("scan_protocol"),
                        sample=sample,
                    )
                    if scan_edge:
                        path.parts.append(scan_edge)
                    path.parts.append(self.scan_node(sample=sample, run=run))
                    processing_edge = self.data_processing_edge(sample=sample)
                    if processing_edge:
                        path.parts.append(processing_edge)
                    path.parts.extend(self.factor_nodes(sample=sample, channel=channel))
                    paths.append(path)
        return paths

    def extract_node(self, sample: dict, channel: dict, run: dict | None) -> SDRFNode:
        accession = self.sample_accession(sample=sample)
        attrs = [SDRFAttr(label="Material Type", value=self.material_type(channel=channel))]
        attrs.extend(self.library_attrs(sample=sample, run=run))
        attrs.extend(self.extra_library_attrs(sample=sample, channel=channel, run=run))
        # Legacy greedy SRA library fallback:
        # attrs.extend(_SRAFallbackComments(self).sra_library_fallback_attrs(run=run))
        return SDRFNode(kind="Extract Name", key=f"extract:{accession}", value=accession, attrs=attrs)

    def library_attrs(self, sample: dict, run: dict | None) -> list[SDRFAttr]:
        attrs = []
        values = {
            "Comment[LIBRARY_LAYOUT]": self.geo_first_value(sample=sample, run=run, field="library_layout"),
            "Comment[LIBRARY_SELECTION]": self.geo_first_value(sample=sample, run=run, field="library_selection"),
            "Comment[LIBRARY_SOURCE]": self.geo_first_value(sample=sample, run=run, field="library_source"),
            "Comment[LIBRARY_STRATEGY]": self.geo_first_value(sample=sample, run=run, field="library_strategy"),
        }
        for label, value in values.items():
            value = self.clean(value)
            if label == "Comment[LIBRARY_SOURCE]" and value:
                value = value.upper()
            if value:
                attrs.append(SDRFAttr(label=label, value=value))
        return attrs

    def geo_first_value(self, sample: dict, run: dict | None, field: str):
        geo_value = self.clean(sample.get(field))
        sra_value = self.clean(run.get(field) if run else None)
        if geo_value and sra_value and geo_value != sra_value:
            self.audit.warnings.append(
                f"Sample {self.sample_accession(sample=sample)} GEO {field} differs from SRA {field}; using GEO value."
            )
        return geo_value or sra_value

    def extra_library_attrs(self, sample: dict, channel: dict, run: dict | None) -> list[SDRFAttr]:
        return []

    def library_protocol_text(self, sample: dict, channel: dict, run: dict | None) -> str | None:
        values = [
            sample.get("library_layout"),
            sample.get("library_strategy"),
            sample.get("library_source"),
            sample.get("library_selection"),
        ]
        if run:
            values.extend([
                run.get("library_layout"),
                run.get("library_strategy"),
                run.get("library_source"),
                run.get("library_selection"),
            ])
        return " | ".join(self.clean(x) for x in values if self.clean(x)) or None

    def assay_node(self, sample: dict, run: dict | None) -> SDRFNode:
        accession = self.sample_accession(sample=sample)
        attrs = [SDRFAttr(label="Technology Type", value="sequencing assay")]
        instrument_model = self.geo_first_instrument_model(sample=sample, run=run)
        for label, value in (
            ("Comment[ENA_SAMPLE]", run.get("sample") if run else None),
            ("Comment[ENA_EXPERIMENT]", run.get("experiment") if run else None),
            ("Comment[ENA_RUN]", run.get("run") if run else None),
            ("Comment[SUBMITTED_FILE_NAME]", run.get("submitted_file_name") if run else None),
            ("Comment[MD5]", run.get("md5") if run else None),
            ("Comment[INSTRUMENT_MODEL]", instrument_model),
        ):
            value = self.clean(value)
            if value:
                attrs.append(SDRFAttr(label=label, value=value))
        attrs.extend(self.extra_assay_attrs(sample=sample, run=run))
        # Legacy greedy SRA run fallback:
        # attrs.extend(_SRAFallbackComments(self).sra_assay_fallback_attrs(sample=sample, run=run))
        value = accession
        if run and self.clean(run.get("geo_sample")) and self.clean(run.get("geo_sample")) != accession:
            self.audit.warnings.append(
                f"Sample {accession} differs from SRA geo_sample {self.clean(run.get('geo_sample'))}; using GEO accession."
            )
        return SDRFNode(kind="Assay Name", key=f"assay:{value}", value=value, attrs=attrs)

    def geo_first_instrument_model(self, sample: dict, run: dict | None):
        geo_value = self.instrument_model(sample=sample)
        sra_value = self.clean(run.get("instrument_model") if run else None)
        if geo_value and sra_value and geo_value != sra_value:
            self.audit.warnings.append(
                f"Sample {self.sample_accession(sample=sample)} GEO instrument_model differs from SRA instrument_model; using GEO value."
            )
        return geo_value or sra_value

    def extra_assay_attrs(self, sample: dict, run: dict | None) -> list[SDRFAttr]:
        return []

    def scan_node(self, sample: dict, run: dict | None) -> SDRFNode:
        accession = self.sample_accession(sample=sample)
        return SDRFNode(
            kind="Scan Name",
            key=f"scan:{accession}",
            value=(run.get("scan_name") if run else None) or accession,
            attrs=self.sequencing_file_attrs(sample=sample, run=run),
        )

    def sequencing_file_attrs(self, sample: dict, run: dict | None) -> list[SDRFAttr]:
        attrs = []
        fastqs = run.get("fastq_files", []) if run else []
        for index, fastq in enumerate(fastqs, start=1):
            filename = self.clean(fastq.get("filename") or fastq.get("uri"))
            for label, value in (
                (f"Comment[read{index} file]", fastq.get("filename")),
                ("Comment[FASTQ_URI]", fastq.get("uri")),
                ("Comment[MD5]", fastq.get("md5")),
            ):
                value = self.clean(value)
                if value:
                    attrs.append(SDRFAttr(label=label, value=value))
            # Legacy greedy SRA FASTQ fallback:
            # attrs.extend(_SRAFallbackComments(self).sra_fastq_fallback_attrs(fastq=fastq))
            if filename and not any(attr.label == f"Comment[read{index} file]" for attr in attrs):
                attrs.append(SDRFAttr(label=f"Comment[read{index} file]", value=filename))

        if len(fastqs) > 2:
            self.audit.warnings.append(f"Sample {self.sample_accession(sample=sample)} has {len(fastqs)} FASTQ files; emitted {len(fastqs)} read file comments.")

        if not fastqs:
            for index, data_file in enumerate(self.raw_files(sample=sample), start=1):
                if classify_file(data_file) == "sequencing_raw":
                    attrs.append(SDRFAttr(label=f"Comment[read{index} file]", value=data_file))
                    if self.arrayexpress_ftp(data_file):
                        attrs.append(SDRFAttr(label="Comment[FASTQ_URI]", value=data_file))

        for data_file in self.derived_files(sample=sample):
            attrs.append(SDRFAttr(label="Derived Array Data File", value=data_file))
        return attrs


class _BulkSequencingSDRFHandler(_SequencingSDRFHandler):
    def build_paths(self) -> list[SDRFPath]:
        paths = []
        for sample in self.ordered_samples():
            runs = self.sra_runs(sample=sample) or [None]
            for channel in self.channels(sample=sample):
                for run in runs:
                    fastqs = self.bulk_fastq_files(sample=sample, run=run) or [None]
                    for fastq in fastqs:
                        path = SDRFPath(parts=[self.source_node(sample=sample, channel=channel)])
                        path.parts.append(self.sample_collection_edge())
                        path.parts.extend(self.extraction_edges(sample=sample, channel=channel))
                        path.parts.append(self.extract_node(sample=sample, channel=channel, run=run))
                        assay_edge = self.protocol_edge(
                            kind="library construction",
                            text=self.library_protocol_text(sample=sample, channel=channel, run=run),
                            sample=sample,
                            required=True,
                        )
                        if assay_edge:
                            path.parts.append(assay_edge)
                        path.parts.append(self.assay_node(sample=sample, run=run))
                        path.parts.append(self.nucleic_acid_sequencing_edge())
                        scan_edge = self.protocol_edge(
                            kind="scan",
                            text=sample.get("scan_protocol"),
                            sample=sample,
                        )
                        if scan_edge:
                            path.parts.append(scan_edge)
                        path.parts.append(self.bulk_scan_node(sample=sample, run=run, fastq=fastq))
                        processing_edge = self.data_processing_edge(sample=sample)
                        if processing_edge:
                            path.parts.append(processing_edge)
                        path.parts.extend(self.factor_nodes(sample=sample, channel=channel))
                        paths.append(path)
        return paths

    def bulk_fastq_files(self, sample: dict, run: dict | None) -> list:
        fastqs = run.get("fastq_files", []) if run else []
        if fastqs:
            return fastqs

        return [
            {
                "filename": data_file,
                "uri": data_file,
                "md5": None,
            }
            for data_file in self.raw_files(sample=sample)
            if classify_file(data_file) == "sequencing_raw"
        ]

    def bulk_scan_node(self, sample: dict, run: dict | None, fastq: dict | None) -> SDRFNode:
        accession = self.sample_accession(sample=sample)
        return SDRFNode(
            kind="Scan Name",
            key=f"scan:{accession}",
            value=(run.get("scan_name") if run else None) or accession,
            attrs=self.bulk_file_attrs(sample=sample, run=run, fastq=fastq),
        )

    def bulk_file_attrs(self, sample: dict, run: dict | None, fastq: dict | None) -> list[SDRFAttr]:
        attrs = []
        if fastq:
            fastq_uri = self.clean(fastq.get("uri") or fastq.get("filename"))
            md5 = self.clean(fastq.get("md5") or (run.get("md5") if run else None))
            if fastq_uri:
                attrs.append(SDRFAttr(label="Comment[FASTQ_URI]", value=fastq_uri))
            if md5:
                attrs.append(SDRFAttr(label="Comment[MD5]", value=md5))

        for data_file in self.derived_files(sample=sample):
            attrs.append(SDRFAttr(label="Derived Array Data File", value=data_file))
        return attrs


class _SingleCellSequencingSDRFHandler(_SequencingSDRFHandler):
    def extra_library_attrs(self, sample: dict, channel: dict, run: dict | None) -> list[SDRFAttr]:
        attrs = []
        values = {
            "Comment[library construction]": self.library_construction(sample=sample),
        }
        read_lengths = run.get("read_lengths") if run else None
        if read_lengths:
            values["Comment[cdna read size]"] = read_lengths[-1]
        for label, value in values.items():
            value = self.clean(value)
            if value:
                attrs.append(SDRFAttr(label=label, value=value))
        return attrs

    def extra_assay_attrs(self, sample: dict, run: dict | None) -> list[SDRFAttr]:
        biosample = self.biosample_accession(sample=sample)
        return [SDRFAttr(label="Comment[technical replicate group]", value=biosample)] if biosample else []

    def library_construction(self, sample: dict):
        return None

    def study_text(self, sample: dict) -> str:
        values = [
            sample.get("title"),
            sample.get("description"),
            sample.get("data_processing"),
        ]
        for channel in self._as_list(sample.get("channel")):
            if isinstance(channel, dict):
                values.extend([
                    channel.get("extract_protocol"),
                    channel.get("growth_protocol"),
                    channel.get("treatment_protocol"),
                ])
        for series in self._as_list(self.data.get("series")):
            if not isinstance(series, dict):
                continue
            values.extend([
                series.get("title"),
                series.get("summary"),
                series.get("overall_design"),
            ])
        return " ".join(str(value).lower() for value in values if value)


class _DropletSingleCellSequencingSDRFHandler(_SingleCellSequencingSDRFHandler):
    def extra_library_attrs(self, sample: dict, channel: dict, run: dict | None) -> list[SDRFAttr]:
        attrs = super().extra_library_attrs(sample=sample, channel=channel, run=run)
        text = self.study_text(sample=sample)
        values = {
            "Comment[cdna read]": "read2",
            "Comment[cell barcode read]": "read1",
            "Comment[cell barcode offset]": "0",
            "Comment[cell barcode size]": "16",
            "Comment[umi barcode read]": "read1",
            "Comment[umi barcode offset]": "16",
            "Comment[umi barcode size]": "12",
            "Comment[sample barcode read]": "index1",
            "Comment[single cell isolation]": self.single_cell_isolation(sample=sample),
        }
        for label, value in values.items():
            value = self.clean(value)
            if value:
                attrs.append(SDRFAttr(label=label, value=value))
        return attrs

    def library_construction(self, sample: dict):
        text = self.study_text(sample=sample)
        if "v3" in text or "3' v3" in text or "3 v3" in text:
            return "10xV3"
        if "10x" in text or "chromium" in text or "droplet" in text:
            return "10x technology"
        return None

    def single_cell_isolation(self, sample: dict):
        text = self.study_text(sample=sample)
        if "10x" in text or "chromium" in text:
            return "10x technology"
        if "droplet" in text:
            return "droplet"
        return None


class _TenXDropletSingleCellSequencingSDRFHandler(_DropletSingleCellSequencingSDRFHandler):
    def tenx_library_attrs(
        self,
        library_construction: str,
        cdna_read_size: str,
        umi_barcode_size: str,
    ) -> list[SDRFAttr]:
        values = {
            "Comment[cdna read]": "read2",
            "Comment[cdna read offset]": "0",
            "Comment[cdna read size]": cdna_read_size,
            "Comment[cell barcode offset]": "0",
            "Comment[cell barcode read]": "read1",
            "Comment[cell barcode size]": "16",
            "Comment[end bias]": "3 prime tag",
            "Comment[input molecule]": "polyA RNA",
            "Comment[library construction]": library_construction,
            "Comment[primer]": "oligo-dT",
            "Comment[LIBRARY_STRAND]": "not applicable",
            "Comment[sample barcode offset]": "0",
            "Comment[sample barcode read]": "index1",
            "Comment[sample barcode size]": "8",
            "Comment[single cell isolation]": "10x technology",
            "Comment[spike in]": "",
            "Comment[umi barcode offset]": "16",
            "Comment[umi barcode read]": "read1",
            "Comment[umi barcode size]": umi_barcode_size,
        }
        return [
            SDRFAttr(label=label, value=self.clean(value), required=value == "")
            for label, value in values.items()
        ]


class _TenXV2DropletSingleCellSequencingSDRFHandler(_TenXDropletSingleCellSequencingSDRFHandler):
    def extra_library_attrs(self, sample: dict, channel: dict, run: dict | None) -> list[SDRFAttr]:
        return self.tenx_library_attrs(
            library_construction="10xV2",
            cdna_read_size="98",
            umi_barcode_size="10",
        )


class _TenXV3DropletSingleCellSequencingSDRFHandler(_TenXDropletSingleCellSequencingSDRFHandler):
    def extra_library_attrs(self, sample: dict, channel: dict, run: dict | None) -> list[SDRFAttr]:
        return self.tenx_library_attrs(
            library_construction="10xV3",
            cdna_read_size="91",
            umi_barcode_size="12",
        )


class _PlateSingleCellSequencingSDRFHandler(_BulkSequencingSDRFHandler):
    def extra_source_attrs(self, sample: dict, channel: dict) -> list[SDRFAttr]:
        barcode = self.clean(sample.get("barcode"))
        description = self.clean(sample.get("description"))
        attrs = []
        if barcode:
            attrs.append(SDRFAttr(label="Comment[index]", value=barcode))
        if description:
            attrs.append(SDRFAttr(label="Description", value=description))
        return attrs


class _SpatialSequencingSDRFHandler(_SingleCellSequencingSDRFHandler):
    def extra_library_attrs(self, sample: dict, channel: dict, run: dict | None) -> list[SDRFAttr]:
        attrs = super().extra_library_attrs(sample=sample, channel=channel, run=run)
        values = {
            "Comment[cdna read]": "read2",
            "Comment[cell barcode read]": "read1",
            "Comment[cell barcode offset]": "0",
            "Comment[cell barcode size]": "16",
            "Comment[umi barcode read]": "read1",
            "Comment[umi barcode offset]": "16",
            "Comment[umi barcode size]": "12",
            "Comment[sample barcode read]": "index1",
        }
        for label, value in values.items():
            attrs.append(SDRFAttr(label=label, value=value))
        return attrs

    def library_construction(self, sample: dict):
        text = self.study_text(sample=sample)
        if "visium" in text:
            return "10x Visium"
        return None

    def extra_assay_attrs(self, sample: dict, run: dict | None) -> list[SDRFAttr]:
        attrs = super().extra_assay_attrs(sample=sample, run=run)
        filename = self.clean(run.get("submitted_file_name") if run else None) or ""
        lower_filename = filename.lower()
        if "_i1_" in lower_filename or "index" in lower_filename:
            attrs.extend([
                SDRFAttr(label="Comment[read_type]", value="sample_barcode"),
                SDRFAttr(label="Comment[read_index]", value="index1"),
            ])
        elif "_r1_" in lower_filename:
            attrs.extend([
                SDRFAttr(label="Comment[read_type]", value="cell_barcode"),
                SDRFAttr(label="Comment[read_index]", value="read1"),
            ])
        elif "_r2_" in lower_filename:
            attrs.extend([
                SDRFAttr(label="Comment[read_type]", value="cdna"),
                SDRFAttr(label="Comment[read_index]", value="read2"),
            ])
        return attrs


class _ArraySDRFHandler(_BaseSDRFHandler):
    def preregister_protocols(self) -> None:
        for sample in self.ordered_samples():
            self.preregister_data_processing(sample=sample)
            for channel in self.channels(sample=sample):
                self.preregister_extraction_protocols(sample=sample, channel=channel)
                self.protocol_registry.get_ref(kind="labeling", text=channel.get("label_protocol"))
                self.protocol_registry.get_ref(kind="hybridization", text=sample.get("hybridization_protocol"))
                self.protocol_registry.get_ref(kind="scan", text=sample.get("scan_protocol"))

    def build_paths(self) -> list[SDRFPath]:
        paths = []
        for sample in self.ordered_samples():
            for channel in self.channels(sample=sample):
                path = SDRFPath(parts=[self.source_node(sample=sample, channel=channel)])
                path.parts.append(self.sample_collection_edge())
                path.parts.extend(self.extraction_edges(sample=sample, channel=channel))
                path.parts.append(self.array_extract_node(sample=sample, channel=channel))
                label_edge = self.protocol_edge(kind="labeling", text=channel.get("label_protocol"), sample=sample, required=bool(channel.get("label")))
                if label_edge:
                    path.parts.append(label_edge)
                if channel.get("label"):
                    path.parts.append(self.labeled_extract_node(sample=sample, channel=channel))
                hybridization_edge = self.protocol_edge(kind="hybridization", text=sample.get("hybridization_protocol"), sample=sample, required=True)
                if hybridization_edge:
                    path.parts.append(hybridization_edge)
                path.parts.append(self.array_assay_node(sample=sample))
                scan_edge = self.protocol_edge(kind="scan", text=sample.get("scan_protocol"), sample=sample, required=False)
                if scan_edge:
                    path.parts.append(scan_edge)
                path.parts.append(SDRFNode(kind="Scan Name", key=f"scan:{self.sample_accession(sample=sample)}", value=self.sample_accession(sample=sample)))
                path.parts.extend(self.array_file_nodes(sample=sample))
                processing_edge = self.data_processing_edge(sample=sample)
                if processing_edge:
                    path.parts.append(processing_edge)
                path.parts.extend(self.factor_nodes(sample=sample, channel=channel))
                paths.append(path)
        return paths

    def array_extract_node(self, sample: dict, channel: dict) -> SDRFNode:
        accession = self.sample_accession(sample=sample)
        return SDRFNode(
            kind="Extract Name",
            key=f"extract:{accession}",
            value=accession,
            attrs=[SDRFAttr(label="Material Type", value=self.material_type(channel=channel))],
        )

    def labeled_extract_node(self, sample: dict, channel: dict) -> SDRFNode:
        accession = self.sample_accession(sample=sample)
        label = self.clean(channel.get("label"))
        return SDRFNode(
            kind="Labeled Extract Name",
            key=f"labeled_extract:{accession}:{label}",
            value=f"{accession}:{label}" if label else accession,
            attrs=[SDRFAttr(label="Label", value=label)] if label else [],
        )

    def array_assay_node(self, sample: dict) -> SDRFNode:
        accession = self.sample_accession(sample=sample)
        array_design = self.platform_accession(sample=sample)
        attrs = [SDRFAttr(label="Technology Type", value="array assay")]
        if array_design:
            attrs.append(SDRFAttr(
                label="Array Design REF",
                value=array_design,
                attrs=[SDRFAttr(label="Term Source REF", value="ArrayExpress")],
            ))
        attrs.extend(self.extra_assay_attrs(sample=sample))
        return SDRFNode(kind="Assay Name", key=f"assay:{accession}", value=accession, attrs=attrs)

    def extra_assay_attrs(self, sample: dict) -> list[SDRFAttr]:
        return []

    def array_file_nodes(self, sample: dict) -> list[SDRFNode]:
        nodes = []
        files = self.supplementary_files(sample=sample)
        raw_count = 0
        derived_count = 0
        for data_file in files:
            file_class = classify_file(data_file)
            if file_class == "array_raw":
                raw_count += 1
                kind = "Image File" if normalized_extension(data_file) in {".tif", ".tiff"} else "Array Data File"
                nodes.append(self.file_node(kind=kind, value=data_file))
            elif file_class == "matrix_or_derived":
                derived_count += 1
                nodes.append(self.file_node(kind="Derived Array Data Matrix File", value=data_file))
            elif file_class == "sequencing_raw":
                nodes.append(self.file_node(kind="Array Data File", value=data_file))
            else:
                derived_count += 1
                nodes.append(self.file_node(kind="Derived Array Data File", value=data_file))

        if raw_count > 1:
            self.audit.warnings.append(f"Sample {self.sample_accession(sample=sample)} has {raw_count} raw files; emitted {raw_count} array raw file columns.")
        if derived_count > 1:
            self.audit.warnings.append(f"Sample {self.sample_accession(sample=sample)} has {derived_count} derived files; emitted {derived_count} derived file columns.")
        return nodes


class _GenericSDRFHandler(_BaseSDRFHandler):
    pass


# Legacy greedy fallback reference code.
# These classes are intentionally commented out and have no runtime effect.
# To restore GEO fallback comments at the old source-node location:
# attrs.extend(_GEOFallbackComments(self).geo_fallback_attrs(sample=sample, channel=channel))
# To restore SRA fallback comments at the old sequencing locations:
# attrs.extend(_SRAFallbackComments(self).sra_library_fallback_attrs(run=run))
# attrs.extend(_SRAFallbackComments(self).sra_assay_fallback_attrs(sample=sample, run=run))
# attrs.extend(_SRAFallbackComments(self).sra_fastq_fallback_attrs(fastq=fastq))
#
# class _FallbackCommentMixin:
#     def __init__(self, handler):
#         self.handler = handler
#
#     def fallback_comment_attrs(self, prefix: str, payload, mapped_paths: set[str]) -> list[SDRFAttr]:
#         attrs = []
#         for path, value in self.flatten_payload(payload=payload):
#             if self.is_mapped_path(path=path, mapped_paths=mapped_paths):
#                 continue
#             rendered = self.flatten_comment_value(value)
#             if rendered is None:
#                 continue
#             attrs.append(SDRFAttr(label=f"Comment[{prefix}_{self.comment_path(path)}]", value=rendered))
#         return attrs
#
#     def flatten_payload(self, payload, path: str = "") -> list[tuple[str, object]]:
#         if isinstance(payload, dict):
#             flattened = []
#             for key in sorted(payload):
#                 child_path = f"{path}.{key}" if path else str(key)
#                 flattened.extend(self.flatten_payload(payload=payload[key], path=child_path))
#             return flattened
#
#         if isinstance(payload, list):
#             flattened = []
#             for index, value in enumerate(payload, start=1):
#                 child_path = f"{path}.{index}" if path else str(index)
#                 flattened.extend(self.flatten_payload(payload=value, path=child_path))
#             return flattened
#
#         return [(path, payload)]
#
#     def is_mapped_path(self, path: str, mapped_paths: set[str]) -> bool:
#         normalized = ".".join(part for part in path.split(".") if not part.isdigit())
#         return normalized in mapped_paths or any(
#             normalized == mapped_path or normalized.startswith(f"{mapped_path}.")
#             for mapped_path in mapped_paths
#         )
#
#     def flatten_comment_value(self, value):
#         if value in (None, ""):
#             return None
#         return self.handler.clean(value)
#
#     def comment_path(self, path: str) -> str:
#         return "_".join(part for part in path.split(".") if not part.isdigit())
#
#
# class _GEOFallbackComments(_FallbackCommentMixin):
#     def geo_fallback_attrs(self, sample: dict, channel: dict) -> list[SDRFAttr]:
#         attrs = []
#         attrs.extend(self.fallback_comment_attrs(
#             prefix="GEO_sample",
#             payload=sample,
#             mapped_paths=self.mapped_geo_sample_paths(),
#         ))
#         attrs.extend(self.fallback_comment_attrs(
#             prefix="GEO_channel",
#             payload=channel,
#             mapped_paths=self.mapped_geo_channel_paths(),
#         ))
#         platform = self.handler.platform(sample=sample)
#         if platform:
#             attrs.extend(self.fallback_comment_attrs(
#                 prefix="GEO_platform",
#                 payload=platform,
#                 mapped_paths=self.mapped_geo_platform_paths(),
#             ))
#         return attrs
#
#     def mapped_geo_sample_paths(self) -> set[str]:
#         return {
#             "iid",
#             "title",
#             "description",
#             "accession",
#             "accession.value",
#             "platform_ref",
#             "platform_ref.ref",
#             "library_strategy",
#             "library_source",
#             "library_selection",
#             "instrument_model",
#             "instrument_model.predefined",
#             "instrument_model.other",
#             "hybridization_protocol",
#             "scan_protocol",
#             "supplementary_data.value",
#             "raw_data.value",
#             "data_processing",
#             "relation.target",
#             "channel",
#         }
#
#     def mapped_geo_channel_paths(self) -> set[str]:
#         return {
#             "position",
#             "source",
#             "organism.name",
#             "characteristics.tag",
#             "characteristics.value",
#             "biomaterial_provider",
#             "treatment_protocol",
#             "growth_protocol",
#             "molecule",
#             "extract_protocol",
#             "label",
#             "label_protocol",
#         }
#
#     def mapped_geo_platform_paths(self) -> set[str]:
#         return {
#             "iid",
#             "accession",
#             "accession.value",
#             "supplementary_data.value",
#         }
#
#
# class _SRAFallbackComments(_FallbackCommentMixin):
#     def sra_library_fallback_attrs(self, run: dict | None) -> list[SDRFAttr]:
#         if not run:
#             return []
#         payload = {
#             key: value
#             for key, value in run.items()
#             if key.startswith("library_") and key not in self.mapped_sra_library_paths()
#         }
#         return self.fallback_comment_attrs(
#             prefix="SRA_library",
#             payload=payload,
#             mapped_paths=set(),
#         )
#
#     def sra_assay_fallback_attrs(self, sample: dict, run: dict | None) -> list[SDRFAttr]:
#         if not run:
#             return []
#         mapped = self.mapped_sra_assay_paths()
#         if self.handler.biosample_accession(sample=sample):
#             mapped = {*mapped, "biosample"}
#         payload = {
#             key: value
#             for key, value in run.items()
#             if key not in mapped and not key.startswith("library_") and key != "fastq_files"
#         }
#         return self.fallback_comment_attrs(
#             prefix="SRA_run",
#             payload=payload,
#             mapped_paths=set(),
#         )
#
#     def sra_fastq_fallback_attrs(self, fastq: dict) -> list[SDRFAttr]:
#         return self.fallback_comment_attrs(
#             prefix="SRA_fastq",
#             payload=fastq,
#             mapped_paths={"filename", "uri", "md5"},
#         )
#
#     def mapped_sra_library_paths(self) -> set[str]:
#         return {
#             "library_layout",
#             "library_selection",
#             "library_source",
#             "library_strategy",
#             "sample",
#             "biosample",
#             "geo_sample",
#             "experiment",
#             "run",
#             "scan_name",
#             "instrument_model",
#             "fastq_files",
#             "submitted_file_name",
#             "md5",
#             "read_lengths",
#         }
#
#     def mapped_sra_assay_paths(self) -> set[str]:
#         return {
#             "library_layout",
#             "library_selection",
#             "library_source",
#             "library_strategy",
#             "sample",
#             "experiment",
#             "run",
#             "scan_name",
#             "instrument_model",
#             "fastq_files",
#             "submitted_file_name",
#             "md5",
#             "read_lengths",
#         }
