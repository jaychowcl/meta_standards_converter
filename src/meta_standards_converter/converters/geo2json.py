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
import re
from collections.abc import Mapping, Sequence

from meta_standards_converter.enrichers.miniml_enricher import MINiMLEnricher
from meta_standards_converter.geo_handlers.geo_parser import GEOParser
from meta_standards_converter.geo_handlers.geo_webfetcher import GEOWebFetcher
from meta_standards_converter.helpers.json_helper import JSONHandler
from meta_standards_converter.miniml import MINiMLCodec, MINiMLPackage
from meta_standards_converter.runtime_contracts import (
   ResourceProfile,
   get_resource_profile,
)

logger = logging.getLogger(__name__)

_GSE_ACCESSION = re.compile(r"\bGSE\d+\b", re.IGNORECASE)


def _as_list(value):
   if value is None:
      return []
   if isinstance(value, (list, tuple)):
      return list(value)
   return [value]


def _gse_values(value) -> list[str]:
   values = []
   for item in _as_list(value):
      if isinstance(item, Mapping):
         candidates = (item.get("value"), item.get("target"), item.get("comment"))
      else:
         candidates = (item,)
      for candidate in candidates:
         if not isinstance(candidate, str):
            continue
         for match in _GSE_ACCESSION.findall(candidate):
            normalized = match.upper()
            if normalized not in values:
               values.append(normalized)
   return values


def _series_accessions(series: Mapping) -> list[str]:
   return _gse_values([series.get("iid"), *_as_list(series.get("accession"))])


def _relation_gses(series: Mapping, relation_type: str) -> list[tuple[str, str]]:
   matches = []
   expected = relation_type.strip().casefold()
   for relation in _as_list(series.get("relation")):
      if not isinstance(relation, Mapping):
         continue
      rendered_type = str(relation.get("type") or "").strip()
      if rendered_type.casefold() != expected:
         continue
      for accession in _gse_values(
         [relation.get("target"), relation.get("comment")]
      ):
         item = (accession, rendered_type)
         if item not in matches:
            matches.append(item)
   return matches


class geo2json(JSONHandler):
   def __init__(
      self,
      enricher=None,
      geo_fetcher=None,
      parser=None,
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
      self.parser = parser or GEOParser(geo_fetcher=self.geo_fetcher)

   def convert(
      self,
      gse: str,
      related_series: bool = False,
      remove_empty: bool = True,
      enrich: bool = True,
      out: str = None,
   ) -> list[MINiMLPackage]:
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
         codec = MINiMLCodec()
         meta_jsons = self._inherit_parent_publications(
            packages=meta_jsons,
            remove_empty=remove_empty,
         )
         packages = []
         for index, meta_json in enumerate(meta_jsons, start=1):
            logger.info("%s: enriching parsed package %d", gse, index)
            packages.append(codec.decode(self.enricher.enrich(data=meta_json)).package)
      else:
         logger.info("%s: skipping JSON enrichment", gse)
         packages = meta_jsons

      if out:
         self.json2file(gse=gse, packages=packages, out=out)

      logger.info("%s: conversion produced %d JSON package(s)", gse, len(packages))
      return packages

   def _inherit_parent_publications(
      self,
      *,
      packages: Sequence[MINiMLPackage],
      remove_empty: bool,
   ) -> list[MINiMLPackage]:
      """Hydrate one unambiguous direct parent publication without adding packages."""

      codec = MINiMLCodec()
      canonical = [codec.encode(codec.decode(package).package) for package in packages]
      available = {}
      for package in canonical:
         series = package.get("series")
         if not isinstance(series, Mapping):
            continue
         for accession in _series_accessions(series):
            available.setdefault(accession, package)

      hydrated = []
      for package in canonical:
         series = package.get("series")
         if not isinstance(series, dict):
            hydrated.append(codec.decode(package).package)
            continue
         if _as_list(series.get("pubmed_id")) or _as_list(
            series.get("pubmed_publication")
         ):
            hydrated.append(codec.decode(package).package)
            continue
         child_accessions = _series_accessions(series)
         parent_relations = _relation_gses(series, "SubSeries of")
         parent_accessions = list(dict.fromkeys(item[0] for item in parent_relations))
         if len(child_accessions) != 1 or len(parent_accessions) != 1:
            hydrated.append(codec.decode(package).package)
            continue

         child_accession = child_accessions[0]
         parent_accession = parent_accessions[0]
         parent_package = available.get(parent_accession)
         if parent_package is None:
            try:
               parent_miniml = self.geo_fetcher.fetch_gse_miniml(
                  gse=parent_accession
               )
               parsed_parents = self.parser.parse(
                  miniml=parent_miniml,
                  remove_empty=remove_empty,
                  related_series=False,
               )
               parent_candidates = []
               for parsed_parent in parsed_parents:
                  candidate = codec.encode(codec.decode(parsed_parent).package)
                  candidate_series = candidate.get("series")
                  if (
                     isinstance(candidate_series, Mapping)
                     and parent_accession in _series_accessions(candidate_series)
                  ):
                     parent_candidates.append(candidate)
               if len(parent_candidates) == 1:
                  parent_package = parent_candidates[0]
            except Exception as error:
               logger.warning(
                  "%s: parent-series publication lookup degraded "
                  "parent=%s error_type=%s",
                  child_accession,
                  parent_accession,
                  type(error).__name__,
               )

         parent_series = (
            parent_package.get("series")
            if isinstance(parent_package, Mapping)
            else None
         )
         if not isinstance(parent_series, Mapping):
            hydrated.append(codec.decode(package).package)
            continue
         reciprocal_children = {
            accession
            for accession, _relation_type in _relation_gses(
               parent_series, "SuperSeries of"
            )
         }
         if child_accession not in reciprocal_children:
            hydrated.append(codec.decode(package).package)
            continue
         pubmed_ids = list(
            dict.fromkeys(
               str(value).strip()
               for value in _as_list(parent_series.get("pubmed_id"))
               if str(value).strip()
            )
         )
         if len(pubmed_ids) != 1:
            hydrated.append(codec.decode(package).package)
            continue

         relation_type = next(
            relation_type
            for accession, relation_type in parent_relations
            if accession == parent_accession
         )
         extensions = series.get("extensions")
         if extensions is None:
            extensions = {}
         if not isinstance(extensions, Mapping):
            hydrated.append(codec.decode(package).package)
            continue
         extensions = dict(extensions)
         extensions["publication_inheritance"] = {
            "relation_type": relation_type,
            "source_series": parent_accession,
            "pubmed_ids": pubmed_ids,
         }
         series["pubmed_id"] = pubmed_ids
         series["extensions"] = extensions
         hydrated.append(codec.decode(package).package)
      return hydrated

   def json2file(self, gse: str, packages: list[MINiMLPackage], out: str) -> str:
      os.makedirs(out, exist_ok=True)
      path = os.path.join(out, f"{gse}.json")
      logger.info("%s: writing JSON package list to %s", gse, path)
      with open(path, "w", encoding="utf-8") as handle:
         json.dump(MINiMLCodec().encode_many(packages), handle, indent=2, ensure_ascii=False)
         handle.write("\n")
      return path
