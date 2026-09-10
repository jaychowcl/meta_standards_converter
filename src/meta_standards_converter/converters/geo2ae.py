# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from meta_standards_converter.magetab.writer import MAGETabWriter
"""
Converter for GEO XML to ArrayExpress MAGETAB format.
"""

import logging

from meta_standards_converter.sources.geo import GEOWebFetcher, GEOSource
from meta_standards_converter.miniml.geo_parser import GEOParser

from meta_standards_converter.magetab.constructor import AEConstructor
from meta_standards_converter.metadata.enrichment import MINiMLEnricher
from meta_standards_converter.helpers.json_helper import JSONHandler
from meta_standards_converter.runtime_contracts import (
   ResourceProfile,
   get_resource_profile,
)

logger = logging.getLogger(__name__)

class GEO2AEConverter(JSONHandler):
   def metrics(self):
      from meta_standards_converter.sources.contracts import request_metrics
      return request_metrics(self.geo_fetcher, self.enricher)

   def __init__(
      self,
      enricher=None,
      geo_fetcher=None,
      parser=None,
      ae_constructor=None,
      resource_profile: str | ResourceProfile = "standard",
      resource_overrides=None,
   ):
      self.resource_profile = get_resource_profile(
         resource_profile, overrides=resource_overrides
      )
      self.enricher = enricher or MINiMLEnricher(
         resource_profile=self.resource_profile, resource_overrides=None
      )
      self.geo_fetcher = geo_fetcher or GEOWebFetcher(
         resource_profile=self.resource_profile, resource_overrides=None
      )
      self.parser = parser or GEOSource(fetcher=self.geo_fetcher, resource_profile=self.resource_profile)
      self.ae_constructor = ae_constructor or AEConstructor()

   def convert(
      self,
      gse: str,
      related_series: bool = False,
      remove_empty: bool = True,
      out: str = None,
      platform_handler: str | None = None,
   ):
      """
      fetches MINIML from GEO using gse accession, parses into meta_json, then writes via AEConstructor.
      """
      # get gse miniml from GEO
      logger.info("%s: fetching GEO MINiML", gse)
      miniml = self.geo_fetcher.fetch_gse_miniml(gse=gse)
      logger.debug("%s: fetched GEO MINiML with %d characters", gse, len(miniml))

      # parse gse miniml to get json
      logger.info("%s: parsing GEO MINiML", gse)
      meta_jsons = self.parser.parse(
         miniml=miniml,
         remove_empty=remove_empty,
         related_series=related_series,
      )
      logger.debug("%s: parsed %d MINiML package(s)", gse, len(meta_jsons))

      # convert to MAGETAB
      constructor = self.ae_constructor
      magetab_dfs = []
      for index, meta_json in enumerate(meta_jsons, start=1):
         logger.info("%s: enriching parsed package %d", gse, index)
         enriched_json = self.enricher.enrich(data=meta_json)
         logger.info("%s: building MAGE-TAB package %d", gse, index)
         if platform_handler is None:
            magetab = constructor.miniml2magetab(data=enriched_json)
         else:
            magetab = constructor.miniml2magetab(
               data=enriched_json,
               platform_handler=platform_handler,
            )
         magetab_dfs.append(magetab)

      #write to outfile if given
      if out:
         for index, magetab in enumerate(magetab_dfs, start=1):
            logger.info("%s: writing MAGE-TAB package %d to %s", gse, index, out)
            MAGETabWriter().write(magetab=magetab, out = out)

      logger.info("%s: conversion produced %d MAGE-TAB package(s)", gse, len(magetab_dfs))
      return magetab_dfs
