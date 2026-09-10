# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
'''
Constructor class for ae MAGETAB idf and sdrf
'''
from meta_standards_converter.magetab.idf import IDFConstructor
from meta_standards_converter.magetab.semantics import overlay_miniml_semantics
from meta_standards_converter.magetab.protocols import ProtocolRegistry
from meta_standards_converter.magetab.technology import (
    _has_tenx_version,
    detect_ae_technology,
    has_array_files,
    normalized_extension,
    series_identity,
)

import copy
import csv
import os


PLATFORM_HANDLER_KEYS = (
    "plate_single_cell_sequencing",
    "droplet_single_cell_sequencing",
    "tenx_v2_droplet_single_cell_sequencing",
    "tenx_v3_droplet_single_cell_sequencing",
    "single_cell_sequencing",
    "spatial_sequencing",
    "bulk_sequencing",
    "sequencing",
    "array",
    "generic",
)


def validate_platform_handler(value: str) -> str:
    if value not in PLATFORM_HANDLER_KEYS:
        raise ValueError(
            f"Unsupported platform handler: {value}. "
            f"Choose one of: {', '.join(PLATFORM_HANDLER_KEYS)}"
        )
    return value


from meta_standards_converter.magetab.sdrf.constructor import SDRFConstructor
from meta_standards_converter.miniml import MINiMLCodec, MINiMLPackage


class AEConstructor:
    def __init__(self, idf_constructor=None, sdrf_constructor=None):
        self.idf_constructor = idf_constructor or IDFConstructor()
        self.sdrf_constructor = sdrf_constructor or SDRFConstructor()

    def miniml2magetab(self, data: MINiMLPackage, platform_handler: str | None = None) -> list:
        """
        converts miniml json to magetab idf. Walks through sections of idf to extract from miniml
        """
        data = MINiMLCodec().encode(MINiMLCodec().decode(data).package)
        forced = platform_handler is not None
        if forced:
            technology_type = validate_platform_handler(platform_handler)
        else:
            technology_type = self._detect_ae_technology(data=data)
        protocol_registry = ProtocolRegistry(series_accession=self._series_accession(data=data))
        sdrf = self.sdrf_constructor._miniml2sdrf(
            data=data,
            protocol_registry=protocol_registry,
            technology_type=technology_type,
        )
        idf = self.idf_constructor.miniml2idf(
            data=data,
            protocol_registry=protocol_registry,
            technology_type=technology_type,
        )
        sdrf_index = self._sdrf_row_index(rows=idf)
        if sdrf_index is None:
            raise ValueError("IDF does not contain an SDRF File row.")
        idf[sdrf_index] = ["SDRF File", sdrf, *idf[sdrf_index][2:]]
        return overlay_miniml_semantics(data, idf)

    def _detect_ae_technology(self, data: dict) -> str:
        return detect_ae_technology(data)

    def _has_tenx_version(self, text: str, version: str) -> bool:
        return _has_tenx_version(text, version)

    def _has_array_files(self, data: dict) -> bool:
        return has_array_files(data)

    def _series_accession(self, data: dict):
        return series_identity(data) or "GEO"



    def _sdrf_row_index(self, rows: list):
        for index, row in enumerate(rows):
            if row and str(row[0]).strip().lower() == "sdrf file":
                return index
        return None
