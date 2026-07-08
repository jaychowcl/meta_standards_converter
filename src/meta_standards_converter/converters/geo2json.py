# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""
Converter for GEO XML to enriched MINiML JSON packages.
"""

import json
import logging
import os

from meta_standards_converter.enrichers.miniml_enricher import MINiMLEnricher
from meta_standards_converter.geo_handlers.geo_parser import GEOParser
from meta_standards_converter.geo_handlers.geo_webfetcher import GEOWebFetcher
from meta_standards_converter.helpers.json_helper import JSONHandler

logger = logging.getLogger(__name__)


class geo2json(JSONHandler):
   def __init__(self, enricher=None, geo_fetcher=None, parser=None):
      self.enricher = enricher or MINiMLEnricher()
      self.geo_fetcher = geo_fetcher or GEOWebFetcher()
      self.parser = parser or GEOParser(geo_fetcher=self.geo_fetcher)

   def convert(
      self,
      gse: str,
      related_series: bool = False,
      remove_empty: bool = True,
      enrich: bool = True,
      out: str = None,
   ) -> list[dict]:
      """
      Fetches GEO MINiML, parses it to JSON packages, optionally enriches, and optionally writes JSON.
      """
      logger.info("%s: fetching GEO MINiML", gse)
      miniml = self.geo_fetcher.fetch_gse_miniml(gse=gse)
      logger.debug("%s: fetched GEO MINiML with %d characters", gse, len(miniml))

      logger.info("%s: parsing GEO MINiML", gse)
      meta_jsons = self.parser.parse(
         miniml=miniml,
         remove_empty=remove_empty,
         related_series=related_series,
      )
      logger.debug("%s: parsed %d MINiML package(s)", gse, len(meta_jsons))

      if enrich:
         packages = []
         for index, meta_json in enumerate(meta_jsons, start=1):
            logger.info("%s: enriching parsed package %d", gse, index)
            packages.append(self.enricher.enrich(data=meta_json))
      else:
         logger.info("%s: skipping JSON enrichment", gse)
         packages = meta_jsons

      if out:
         self.json2file(gse=gse, packages=packages, out=out)

      logger.info("%s: conversion produced %d JSON package(s)", gse, len(packages))
      return packages

   def json2file(self, gse: str, packages: list[dict], out: str) -> str:
      os.makedirs(out, exist_ok=True)
      path = os.path.join(out, f"{gse}.json")
      logger.info("%s: writing JSON package list to %s", gse, path)
      with open(path, "w", encoding="utf-8") as handle:
         json.dump(packages, handle, indent=2, ensure_ascii=False)
         handle.write("\n")
      return path
