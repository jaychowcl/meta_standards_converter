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
from copy import deepcopy

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
        *, publication_identifier_resolver=None,
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
        self.profile = profile
        self.publication_identifier_resolver = publication_identifier_resolver

    def enrich(self, data: MINiMLPackage) -> MINiMLPackage:
        started = time.monotonic()
        codec = MINiMLCodec()
        package = codec.decode(data).package
        mutable = codec.encode(package)
        self._pubmed_failures = 0
        self._sra_failures = 0
        self.publication_issues = []
        self._identifier_cache = {}
        self._identifier_http = None
        from ..miniml.publication_identifiers import clean_publication_identifiers
        invalid_records = clean_publication_identifiers(mutable)
        if mutable.get("source", {}).get("format") in {"SRA", "ENA"}:
            self._publication_cache = {}
            from meta_standards_converter.miniml.archive_entities import declare_ontologies
            from meta_standards_converter.miniml.archive_residuals import finalize, source_records
            self.enrich_pubmed(data=mutable, fill_missing=True)
            for sample in mutable.get('sample', []):
                for entity in [sample, *sample.get('sra_run', [])]:
                    if entity.get('pubmed_publication') or entity.get('pubmed_id'):
                        self.enrich_pubmed({'series':entity}, fill_missing=True)
            owners = [mutable.get('series', {})]
            owners.extend(entity for sample in mutable.get('sample', []) for entity in [sample, *sample.get('sra_run', [])])
            for owner in owners:
                for relation in owner.get('relation', []):
                    if relation.get('publication'):
                        self.enrich_pubmed({'series':{'pubmed_publication':[relation['publication']]}}, fill_missing=True)
            del self._publication_cache
            declare_ontologies(mutable)
            return finalize(mutable, source_records(package) + invalid_records)
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

        from ..miniml.publication_identifiers import valid_pubmed_ids
        from ..sources.archive_publications import citation_identifier
        for publication in series.get('pubmed_publication', []):
            if valid_pubmed_ids([publication.get('pubmed_id')]):
                continue
            supplied = [(kind, identifier) for kind in ('pmcid','doi')
                        if (identifier := citation_identifier(kind, publication.get(kind)).get(kind))]
            resolved = {pmid for kind, value in supplied
                        if (pmid := self._resolve_publication_identifier(kind, value))}
            if len(resolved) > 1:
                message = f"Conflicting publication identifiers for {series.get('iid')}: {supplied} resolve to different PMIDs; retaining the source citation without hydration."
                if message not in self.publication_issues:
                    self.publication_issues.append(message)
                    logger.warning(message)
            elif resolved:
                publication['pubmed_id'] = resolved.pop()
        pubmed_ids = valid_pubmed_ids(self._as_list(series.get("pubmed_id")))
        if not fill_missing and series.get('pubmed_publication'):
            # Refresh a saved citation through the same conservative hydration
            # contract used by native imports, preserving legacy ID occurrences.
            declared = deepcopy(series.get('pubmed_id'))
            self.enrich_pubmed(data, fill_missing=True)
            if declared is not None:
                series['pubmed_id'] = declared
            return data
        if fill_missing:
            publications = series.setdefault('pubmed_publication', [])
            pubmed_ids = valid_pubmed_ids([*pubmed_ids, *[p.get('pubmed_id') for p in publications],
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

    def _resolve_publication_identifier(self, kind, value):
        from ..sources.archive_publications import resolve_identifier, citation_identifier
        from ..sources.archive_support import ArchiveHTTP
        cache = getattr(self, '_identifier_cache', None)
        if cache is None:
            self._identifier_cache = cache = {}
        key = (kind, value.casefold() if kind == 'doi' else value)
        if key not in cache:
            try:
                if self.publication_identifier_resolver is not None:
                    resolved = self.publication_identifier_resolver(kind, value)
                else:
                    if getattr(self, '_identifier_http', None) is None:
                        self._identifier_http = ArchiveHTTP('ncbi_eutils', resource_profile=self.profile)
                    resolved = resolve_identifier(kind, value, self._identifier_http)
                cache[key] = citation_identifier('pubmed', resolved).get('pubmed_id')
                if resolved and not cache[key]:
                    raise ValueError('invalid resolved PMID')
            except (requests.RequestException, ET.ParseError, ValueError, KeyError) as error:
                cache[key] = None
                issue = f'{kind} {value}: {type(error).__name__}'
                self.publication_issues.append(issue)
                logger.warning('%s', issue)
        return cache[key]

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

            sample["sra_accession"] = self._dedupe([*self._as_list(sample.get('sra_accession')), *accessions])
            runs = deepcopy(self._as_list(sample.get('sra_run')))
            fetched = []
            for accession in accessions:
                try:
                    fetched.extend(self.insdc_fetcher.fetch_sra_runs(accession=accession))
                except (requests.RequestException, ET.ParseError) as error:
                    self._sra_failures = getattr(self, "_sra_failures", 0) + 1
                    logger.warning('SRA refresh %s failed (%s); retaining saved metadata', accession, type(error).__name__)
                    continue
            self._merge_saved_runs(runs, fetched)
            sample["sra_run"] = runs
            ena_accessions = self._dedupe(
                run.get("study")
                for run in runs
                if isinstance(run, dict)
            )
            if ena_accessions:
                sample["ena_accession"] = self._dedupe([*self._as_list(sample.get('ena_accession')), *ena_accessions])

        return data

    @staticmethod
    def _merge_saved_runs(saved, fetched):
        """Complete sample-local records without selecting conflicting versions."""
        empty = (None, '', [], {})
        groups = {}
        for position, record in enumerate(fetched):
            groups.setdefault(record.get('run') or ('unidentified', position), []).append(record)
        for identity, incoming in groups.items():
            matches = [r for r in saved if r.get('run') and r['run'] == identity]
            owners = [*matches, *incoming]
            if len(matches) > 1 or any(len({r[k] for r in owners if r.get(k)}) > 1
                                      for k in ('experiment', 'sample', 'biosample', 'study')):
                logger.warning('SRA refresh has ambiguous or conflicting identity for %s; retaining saved record', identity)
                continue
            target = matches[0] if matches else deepcopy(incoming[0])
            for key in dict.fromkeys(k for record in incoming for k in record):
                values = [r[key] for r in incoming if r.get(key) not in empty]
                if key == 'fastq_files':
                    if not values:
                        target.setdefault(key, [])
                        continue
                    if target.get(key) in empty:
                        target[key] = []
                    files = target[key]
                    MINiMLEnricher._merge_saved_files(files, [f for group in values for f in group])
                elif target.get(key) in empty and values:
                    if all(v == values[0] for v in values):
                        target[key] = deepcopy(values[0])
                    else:
                        logger.warning('SRA refresh has conflicting %s values for %s; retaining saved field', key, identity)
            if not matches and target not in saved:
                saved.append(target)

    @staticmethod
    def _merge_saved_files(saved, fetched):
        """Use whole-response compatibility before completing a saved URI."""
        empty = (None, '', [], {})
        original = deepcopy(saved)

        def compatible(a, b):
            return all(a.get(k) in empty or v in empty or a[k] == v for k, v in b.items())

        for file in fetched:
            if file in saved:
                continue
            candidates = [i for i, old in enumerate(original) if file.get('uri') and old.get('uri') == file['uri'] and compatible(old, file)]
            if len(candidates) == 1:
                index = candidates[0]
                neighbors = [f for f in fetched if f.get('uri') == file['uri'] and compatible(original[index], f)]
                if all(compatible(a, b) for a in neighbors for b in neighbors):
                    for field, fact in file.items():
                        if saved[index].get(field) in empty:
                            saved[index][field] = deepcopy(fact)
                    continue
            # Ambiguous/conflicting versions keep complete independent facts;
            # filename similarity and response order never select an identity.
            saved.append(deepcopy(file))

    def _pubmed_publication(self, pubmed_id: str) -> dict:
        from ..miniml.publication_identifiers import valid_pubmed_ids
        if not valid_pubmed_ids([pubmed_id]):
            return {'pubmed_id': ''}
        cache = getattr(self, '_publication_cache', None)
        if cache is not None and pubmed_id in cache: return dict(cache[pubmed_id])
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

        result = {
            "pubmed_id": pubmed_id,
            "doi": doi,
            "author_list": authors,
            "title": title,
            "status": status,
            "status_term_source_ref": source_ref,
            "status_term_accession_number": accession,
        }
        if cache is not None: cache[pubmed_id] = result
        return result

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
        from ..miniml.publication_identifiers import valid_pubmed_ids
        return [self.pubmed.pubmed_summary(pubmed_id=value)
                for value in valid_pubmed_ids(handler._from_path(data, "series.pubmed_id.*"))]

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
