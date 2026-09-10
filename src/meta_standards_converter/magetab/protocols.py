# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from __future__ import annotations

import os

import re

from urllib.parse import urlparse

from meta_standards_converter.metadata.ontology_mappings import Harmonizer

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
