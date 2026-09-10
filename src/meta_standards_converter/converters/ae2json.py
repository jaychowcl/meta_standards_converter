# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Convert MAGE-TAB IDF/SDRF metadata to parsed MINiML-compatible JSON."""

import json
import logging
import os

from meta_standards_converter.ae_handlers.ae_parser import AEParser
from meta_standards_converter.sources.magetab import AEWebFetcher
from meta_standards_converter.miniml import MINiMLCodec, MINiMLPackage
from meta_standards_converter.retrieval import RetrievalPolicy
from meta_standards_converter.runtime_contracts import (
    ResourceProfile,
    get_resource_profile,
)


logger = logging.getLogger(__name__)


class AE2JSONConverter:
    def metrics(self):
        from meta_standards_converter.sources.contracts import request_metrics
        return request_metrics(self.fetcher)

    def __init__(
        self,
        fetcher=None,
        parser=None,
        resource_profile: str | ResourceProfile = "standard",
        resource_overrides=None,
        source_hosts=(),
    ):
        profile = get_resource_profile(
            resource_profile,
            overrides=resource_overrides,
        )
        self.fetcher = fetcher or AEWebFetcher(
            resource_profile=profile,
            retrieval_policy=RetrievalPolicy(
                resource_profile=profile,
                allowed_hosts=frozenset(source_hosts),
                allowed_schemes=frozenset({"https"}),
            ),
        )
        self.parser = parser or AEParser()

    def convert(
        self,
        source: str,
        out: str | None = None,
        sdrf_sources: list[str] | None = None,
    ) -> list[MINiMLPackage]:
        logger.info("%s: resolving MAGE-TAB metadata", source)
        resolved = self.fetcher.resolve(source, sdrf_sources=sdrf_sources)
        logger.info("%s: parsing IDF and %d SDRF file(s)", source, len(resolved.sdrfs))
        package = self.parser.parse(resolved)
        packages = [package]
        if out:
            os.makedirs(out, exist_ok=True)
            accession = self._study_accession(package)
            path = os.path.join(out, f"{self._safe_filename(accession)}.json")
            with open(path, "w", encoding="utf-8") as handle:
                json.dump(MINiMLCodec().encode_many(packages), handle, indent=2, ensure_ascii=False)
                handle.write("\n")
            logger.info("%s: wrote parsed JSON to %s", source, path)
        return packages

    def _study_accession(self, package):
        accessions = package.get("series", {}).get("accession", [])
        for item in accessions:
            if isinstance(item, dict) and item.get("database") == "ArrayExpress" and item.get("value"):
                return str(item["value"])
        for item in accessions:
            if isinstance(item, dict) and item.get("value"):
                return str(item["value"])
        return "MAGE-TAB"

    def _safe_filename(self, value):
        token = str(value).replace(os.sep, "_")
        if os.altsep:
            token = token.replace(os.altsep, "_")
        return token or "MAGE-TAB"
