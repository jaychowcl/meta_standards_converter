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

from meta_standards_converter.sources.insdc import INSDCWebfetcher
from meta_standards_converter.miniml import MINiMLCodec, MINiMLPackage
from meta_standards_converter.sources.pubmed import PubmedWebFetcher
from meta_standards_converter.runtime_contracts import (
    ResourceProfile,
    get_resource_profile,
)


logger = logging.getLogger(__name__)


from typing import Protocol


class MetadataEnrichment(Protocol):
    def enrich(self, data: MINiMLPackage) -> MINiMLPackage: ...


class MINiMLEnricher:

    def metrics(self):
        from meta_standards_converter.sources.contracts import request_metrics
        return request_metrics(self.pubmed_fetcher, self.insdc_fetcher)

    def __init__(
        self,
        pubmed_fetcher=None,
        insdc_fetcher=None,
        resource_profile: str | ResourceProfile = "standard",
        resource_overrides=None,
    ):
        profile = get_resource_profile(
            resource_profile,
            overrides=resource_overrides,
        )
        self.pubmed_fetcher = pubmed_fetcher or PubmedWebFetcher(
            resource_profile=profile
        )
        self.insdc_fetcher = insdc_fetcher or INSDCWebfetcher(
            resource_profile=profile
        )
        self.publication_issues = []

    def enrich(self, data: MINiMLPackage) -> MINiMLPackage:
        started = time.monotonic()
        codec = MINiMLCodec()
        package = codec.decode(data).package
        mutable = codec.encode(package)
        self._pubmed_failures = 0
        self._sra_failures = 0
        self.publication_issues = []
        if mutable.get("source", {}).get("format") in {"SRA", "ENA"}:
            from meta_standards_converter.miniml.archive_entities import declare_ontologies
            from meta_standards_converter.miniml.archive_residuals import finalize, source_records
            self.enrich_pubmed(data=mutable, fill_missing=True)
            declare_ontologies(mutable)
            return finalize(mutable, source_records(package))
        self.enrich_pubmed(data=mutable)
        self.enrich_sra(data=mutable)
        series = mutable.get("series") if isinstance(mutable.get("series"), dict) else {}
        pubmed_ids = self._dedupe(self._as_list(series.get("pubmed_id")))
        samples = [item for item in self._as_list(mutable.get("sample")) if isinstance(item, dict)]
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
        return codec.decode(mutable).package

    def enrich_pubmed(self, data: dict, *, fill_missing: bool = False) -> dict:
        series = data.get("series")
        if not isinstance(series, dict):
            return data

        pubmed_ids = self._dedupe(self._as_list(series.get("pubmed_id")))
        if fill_missing:
            publications = series.setdefault('pubmed_publication', [])
            pubmed_ids = self._dedupe([*pubmed_ids, *[p.get('pubmed_id') for p in publications],
                *[r.get('target') for r in series.get('relation', []) if str(r.get('type', '')).lower() == 'pubmed']])
            if not pubmed_ids:
                return data
            series['pubmed_id'] = pubmed_ids
            fields = ('doi', 'author_list', 'title')
            status_fields = ('status', 'status_term_source_ref', 'status_term_accession_number')
            for pmid in pubmed_ids:
                existing = [p for p in publications if str(p.get('pubmed_id')) == str(pmid)]
                if not existing:
                    existing = [{'pubmed_id': pmid}]
                    publications.extend(existing)
                if all(all(p.get(k) for k in (*fields, *status_fields)) for p in existing):
                    continue
                fetched = self._pubmed_publication(pubmed_id=pmid)
                for publication in existing:
                    for field in fields:
                        if not publication.get(field) and fetched.get(field):
                            publication[field] = fetched[field]
                    compatible = all(not publication.get(k) or publication[k] == fetched.get(k) for k in status_fields)
                    if compatible:
                        for field in status_fields:
                            if not publication.get(field) and fetched.get(field):
                                publication[field] = fetched[field]
            return data
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
                accessions.extend(self.insdc_fetcher.extract_sra_accessions(relation.get("target") or ""))

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
            if not any((doi, authors, title, status)):
                raise ValueError('PubMed response contains no citation metadata')
        except (requests.RequestException, ET.ParseError, ValueError) as error:
            self._pubmed_failures = getattr(self, "_pubmed_failures", 0) + 1
            issue = f'PubMed {pubmed_id}: {type(error).__name__}'
            self.publication_issues.append(issue)
            logger.warning('%s', issue)
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


class MAGETabEvidenceResolver:
    """Resolve only evidence requested by MAGE-TAB construction, per operation."""

    def __init__(self, pubmed_client=None, insdc_client=None):
        self.pubmed = pubmed_client if pubmed_client is not None else PubmedWebFetcher()
        self.insdc = insdc_client if insdc_client is not None else INSDCWebfetcher()

    def publications(self, data):
        if data.get("source", {}).get("format") in {"SRA", "ENA"}:
            return []
        from meta_standards_converter.helpers.json_helper import JSONHandler
        handler = JSONHandler()
        if any(isinstance(p, dict) for p in handler._from_path(data, "series.pubmed_publication.*")):
            return []
        return [self.pubmed.pubmed_summary(pubmed_id=value)
                for value in handler._from_path(data, "series.pubmed_id.*") if value]

    def sample_runs(self, handler, technology_type):
        import requests
        import xml.etree.ElementTree as ET
        sequencing = technology_type not in {"array", "generic"}
        if not sequencing:
            return {}
        cache, result = {}, {}
        for sample in handler.ordered_samples():
            if "sra_run" in sample:
                continue
            accessions = []
            for relation in handler._as_list(sample.get("relation")):
                if isinstance(relation, dict) and (relation.get("type") or "").lower() == "sra":
                    accessions.extend(self.insdc.extract_sra_accessions(relation.get("target") or ""))
            runs = []
            for accession in dict.fromkeys(accessions):
                if accession not in cache:
                    cache[accession] = self.fetch_runs(accession)
                runs.extend(cache[accession])
            result[id(sample)] = runs
        return result

    def fetch_runs(self, accession):
        import requests
        import xml.etree.ElementTree as ET
        try:
            return self.insdc.fetch_sra_runs(accession=accession)
        except (requests.RequestException, ET.ParseError):
            return []
