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
