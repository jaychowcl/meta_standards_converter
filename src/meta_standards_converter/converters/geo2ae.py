# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""
Converter for GEO XML to ArrayExpress MAGETAB format.
"""

import logging

from meta_standards_converter.geo_handlers.geo_webfetcher import GEOWebFetcher
from meta_standards_converter.geo_handlers.geo_parser import GEOParser

from meta_standards_converter.ae_handlers.ae_constructor import AEConstructor
from meta_standards_converter.enrichers.miniml_enricher import MINiMLEnricher
from meta_standards_converter.helpers.json_helper import JSONHandler

logger = logging.getLogger(__name__)

class geo2ae(JSONHandler):
   def __init__(self, enricher=None, geo_fetcher=None, parser=None):
      self.enricher = enricher or MINiMLEnricher()
      self.geo_fetcher = geo_fetcher or GEOWebFetcher()
      self.parser = parser or GEOParser(geo_fetcher=self.geo_fetcher)

   def convert(
      self,
      gse: str,
      related_series: bool = False,
      remove_empty: bool = True,
      out: str = None,
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
      constructor = AEConstructor()
      magetab_dfs = []
      for index, meta_json in enumerate(meta_jsons, start=1):
         logger.info("%s: enriching parsed package %d", gse, index)
         enriched_json = self.enricher.enrich(data=meta_json)
         logger.info("%s: building MAGE-TAB package %d", gse, index)
         magetab_dfs.append(constructor.miniml2magetab(data=enriched_json))

      #write to outfile if given 
      if out:
         for index, magetab in enumerate(magetab_dfs, start=1):
            logger.info("%s: writing MAGE-TAB package %d to %s", gse, index, out)
            constructor.magetab2file(magetab=magetab, out = out)
         
      logger.info("%s: conversion produced %d MAGE-TAB package(s)", gse, len(magetab_dfs))
      return magetab_dfs
