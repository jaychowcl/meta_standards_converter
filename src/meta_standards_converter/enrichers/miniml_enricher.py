# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""
Adds remote metadata lookups to parsed MINiML JSON packages.
"""

import xml.etree.ElementTree as ET
import logging
import time

import requests

from meta_standards_converter.insdc_handlers.insdc_webfetcher import INSDCWebfetcher
from meta_standards_converter.pubmed_handlers.pubmed_webfetcher import PubmedWebFetcher


logger = logging.getLogger(__name__)


class MINiMLEnricher:
    def __init__(self, pubmed_fetcher=None, insdc_fetcher=None):
        self.pubmed_fetcher = pubmed_fetcher or PubmedWebFetcher()
        self.insdc_fetcher = insdc_fetcher or INSDCWebfetcher()

    def enrich(self, data: dict) -> dict:
        started = time.monotonic()
        self._pubmed_failures = 0
        self._sra_failures = 0
        self.enrich_pubmed(data=data)
        self.enrich_sra(data=data)
        series = data.get("series") if isinstance(data.get("series"), dict) else {}
        pubmed_ids = self._dedupe(self._as_list(series.get("pubmed_id")))
        samples = [item for item in self._as_list(data.get("sample")) if isinstance(item, dict)]
        sra_accessions = sum(len(self._as_list(item.get("sra_accession"))) for item in samples)
        sra_runs = sum(len(self._as_list(item.get("sra_run"))) for item in samples)
        logger.info(
            "MINiML enrichment stats samples=%s pubmed_ids=%s pubmed_failures=%s sra_accessions=%s sra_runs=%s sra_failures=%s elapsed_seconds=%.3f",
            len(samples),
            len(pubmed_ids),
            self._pubmed_failures,
            sra_accessions,
            sra_runs,
            self._sra_failures,
            time.monotonic() - started,
        )
        return data

    def enrich_pubmed(self, data: dict) -> dict:
        series = data.get("series")
        if not isinstance(series, dict):
            return data

        pubmed_ids = self._dedupe(self._as_list(series.get("pubmed_id")))
        if not pubmed_ids:
            return data

        series["pubmed_publication"] = [
            self._pubmed_publication(pubmed_id=pubmed_id)
            for pubmed_id in pubmed_ids
        ]
        return data

    def enrich_sra(self, data: dict) -> dict:
        for sample in self._as_list(data.get("sample")):
            if not isinstance(sample, dict):
                continue

            accessions = []
            for relation in self._as_list(sample.get("relation")):
                if not isinstance(relation, dict):
                    continue
                if (relation.get("type") or "").lower() != "sra":
                    continue
                accessions.extend(self.insdc_fetcher._extract_sra(relation.get("target") or ""))

            accessions = self._dedupe(accessions)
            if not accessions:
                continue

            sample["sra_accession"] = accessions
            runs = []
            for accession in accessions:
                try:
                    runs.extend(self.insdc_fetcher.fetch_sra_runs(accession=accession))
                except (requests.RequestException, ET.ParseError):
                    self._sra_failures = getattr(self, "_sra_failures", 0) + 1
                    continue
            sample["sra_run"] = runs
            ena_accessions = self._dedupe(
                run.get("study")
                for run in runs
                if isinstance(run, dict)
            )
            if ena_accessions:
                sample["ena_accession"] = ena_accessions

        return data

    def _pubmed_publication(self, pubmed_id: str) -> dict:
        try:
            doi, authors, title, status, source_ref, accession = self.pubmed_fetcher.pubmed_summary(
                pubmed_id=pubmed_id
            )
        except (requests.RequestException, ET.ParseError):
            self._pubmed_failures = getattr(self, "_pubmed_failures", 0) + 1
            doi, authors, title, status, source_ref, accession = (None, None, None, None, None, None)

        return {
            "pubmed_id": pubmed_id,
            "doi": doi,
            "author_list": authors,
            "title": title,
            "status": status,
            "status_term_source_ref": source_ref,
            "status_term_accession_number": accession,
        }

    def _as_list(self, value):
        if value is None:
            return []
        if isinstance(value, list):
            return value
        return [value]

    def _dedupe(self, values) -> list:
        deduped = []
        for value in values:
            if value and value not in deduped:
                deduped.append(value)
        return deduped
