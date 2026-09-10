# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from __future__ import annotations
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
from meta_standards_converter.magetab.sdrf.model import SDRFPath
from meta_standards_converter.magetab.sdrf.handlers.base import _BaseSDRFHandler
from meta_standards_converter.magetab.sdrf.handlers.base import classify_file

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
