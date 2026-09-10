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

from meta_standards_converter.magetab.sdrf.handlers.sequencing import _SequencingSDRFHandler
from meta_standards_converter.magetab.sdrf.handlers.sequencing import _BulkSequencingSDRFHandler
from meta_standards_converter.magetab.sdrf.handlers.single_cell import _SingleCellSequencingSDRFHandler
from meta_standards_converter.magetab.sdrf.handlers.single_cell import _DropletSingleCellSequencingSDRFHandler
from meta_standards_converter.magetab.sdrf.handlers.single_cell import _TenXV2DropletSingleCellSequencingSDRFHandler
from meta_standards_converter.magetab.sdrf.handlers.single_cell import _TenXV3DropletSingleCellSequencingSDRFHandler
from meta_standards_converter.magetab.sdrf.handlers.single_cell import _PlateSingleCellSequencingSDRFHandler
from meta_standards_converter.magetab.sdrf.handlers.spatial import _SpatialSequencingSDRFHandler
from meta_standards_converter.magetab.sdrf.handlers.array import _ArraySDRFHandler
from meta_standards_converter.magetab.sdrf.handlers.generic import _GenericSDRFHandler

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
        return detect_ae_technology(data)

    def _has_array_files(self, data: dict) -> bool:
        return has_array_files(data)

    def _lookup_sra(self, sra: str) -> list:
        '''
        take sra accession to return srr info
        '''
        try:
            return self.insdc_fetcher.fetch_sra_runs(accession=sra)
        except (requests.RequestException, ET.ParseError):
            return []
