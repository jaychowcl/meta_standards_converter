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

from meta_standards_converter.magetab.sdrf.handlers.base import _BaseSDRFHandler

class _GenericSDRFHandler(_BaseSDRFHandler):
    pass
