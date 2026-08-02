# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Neutral protocol state and technology detection shared by AE constructors."""

from __future__ import annotations

import os
import re
from urllib.parse import urlparse

from meta_standards_converter.harmonizers.harmonizers import Harmonizer
from meta_standards_converter.helpers.json_helper import JSONHandler


class ProtocolRegistry:
    """Allocate stable protocol references shared by IDF and SDRF construction."""

    LABEL_BY_KIND = {
        "manufacture": "Manufacture-Protocol",
        "treatment": "Treatment-Protocol",
        "growth": "Growth-Protocol",
        "extraction": "Extract-Protocol",
        "extract": "Extract-Protocol",
        "library construction": "Library-Construction-Protocol",
        "labeling": "Label-Protocol",
        "label": "Label-Protocol",
        "hybridization": "Hybridization-Protocol",
        "scan": "Scan-Protocol",
        "data processing": "Data-Processing",
        "sample collection": "Sample-Collection-Protocol",
        "nucleic acid sequencing": "Nucleic-Acid-Sequencing-Protocol",
    }

    def __init__(self, series_accession: str):
        self.series_accession = series_accession
        self.by_key: dict[tuple[str, str], dict] = {}

    def get_ref(self, kind: str, text: str | None, label: str | None = None) -> str | None:
        text = self.clean(text)
        if not text:
            return None
        label = label or self.LABEL_BY_KIND.get(kind, kind)
        key = (kind, text)
        if key not in self.by_key:
            self.by_key[key] = {
                "ref": f"P-{self.series_accession}-{len(self.by_key) + 1}",
                "kind": kind,
                "label": label,
                "text": text,
            }
        return self.by_key[key]["ref"]

    def ensure_required(self, kind: str, label: str | None = None) -> str:
        label = label or self.LABEL_BY_KIND.get(kind, kind)
        required_type = Harmonizer().geoprotocols2efo(protocol_type=label)[0]
        for record in self.records():
            record_type = Harmonizer().geoprotocols2efo(protocol_type=record["label"])[0]
            if record_type == required_type:
                return record["ref"]
        key = (kind, "")
        if key not in self.by_key:
            self.by_key[key] = {
                "ref": f"P-{self.series_accession}-{len(self.by_key) + 1}",
                "kind": kind,
                "label": label,
                "text": "",
                "required": True,
            }
        return self.by_key[key]["ref"]

    def records(self) -> list[dict]:
        return list(self.by_key.values())

    @staticmethod
    def clean(value):
        if value is None:
            return None
        return " ".join(str(value).replace("\t", " ").replace("\n", " ").split())


def normalized_extension(path: str) -> str:
    parsed = urlparse(str(path))
    basename = os.path.basename(parsed.path or str(path)).lower()
    for suffix in (".gz", ".zip", ".bz2", ".xz"):
        if basename.endswith(suffix):
            basename = basename[: -len(suffix)]
            break
    return os.path.splitext(basename)[1]


def has_array_files(data: dict) -> bool:
    handler = JSONHandler()
    values = []
    for path in (
        "platform.*.supplementary_data.*.value",
        "sample.*.supplementary_data.*.value",
        "sample.*.raw_data.*.value",
        "series.supplementary_data.*.value",
    ):
        values.extend(x for x in handler._from_path(data, path) if x)
    extensions = (".cel", ".gpr", ".idat", ".chp", ".txt", ".tif", ".tiff", ".exp", ".rpt", ".cab")
    return any(normalized_extension(value) in extensions for value in values)


def _has_tenx_version(text: str, version: str) -> bool:
    if "10x" not in text and "chromium" not in text:
        return False
    return re.search(rf"(?<![a-z0-9])v{version}(?![a-z0-9])", text) is not None


def detect_ae_technology(data: dict) -> str:
    """Select the shared platform-handler key without importing either constructor."""

    handler = JSONHandler()
    values = lambda path: (str(x).lower() for x in handler._from_path(data, path) if x)
    platform_tech = " ".join(values("platform.*.technology"))
    library_source = " ".join(values("sample.*.library_source"))
    library_strategy = " ".join(values("sample.*.library_strategy"))
    sample_type = " ".join(values("sample.*.type"))
    text_paths = (
        "series.title", "series.summary", "series.overall_design", "series.type.*",
        "sample.*.description", "sample.*.data_processing",
        "sample.*.channel.*.extract_protocol", "sample.*.channel.*.growth_protocol",
        "sample.*.channel.*.treatment_protocol", "sample.*.channel.*.molecule",
        "sample.*.channel.*.characteristics.*.tag",
        "sample.*.channel.*.characteristics.*.value",
        "sample.*.supplementary_data.*.value", "sample.*.raw_data.*.value",
        "series.supplementary_data.*.value",
    )
    text = " ".join(value for path in text_paths for value in values(path))
    relations = [
        x for x in handler._from_path(data, "sample.*.relation.*") if isinstance(x, dict)
    ]
    has_sra = any((relation.get("type") or "").lower() == "sra" for relation in relations)
    if "high-throughput sequencing" in platform_tech or has_sra or library_strategy or sample_type == "sra":
        if "single cell" in library_source or "single-cell" in text or "single cell" in text or "10x" in text:
            if "visium" in text or "spatial" in text:
                return "spatial_sequencing"
            if "10x" not in text and "droplet" not in text and "chromium" not in text:
                return "plate_single_cell_sequencing"
            if _has_tenx_version(text, "3"):
                return "tenx_v3_droplet_single_cell_sequencing"
            if _has_tenx_version(text, "2"):
                return "tenx_v2_droplet_single_cell_sequencing"
            return "droplet_single_cell_sequencing"
        return "bulk_sequencing"
    if "array" in platform_tech or has_array_files(data):
        return "array"
    return "generic"


__all__ = ["ProtocolRegistry", "detect_ae_technology", "has_array_files", "normalized_extension"]
