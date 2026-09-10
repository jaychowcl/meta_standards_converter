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
