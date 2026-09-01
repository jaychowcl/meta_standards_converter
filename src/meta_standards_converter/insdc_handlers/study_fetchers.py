# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Study-scoped NCBI SRA and ENA retrieval with bounded pagination."""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping
import json
import logging
import re
from typing import Any
from urllib.parse import urlencode

from meta_standards_converter.helpers.request_helper import (
    NCBIApplicationIdentity,
    RateLimitedRequester,
    RequestSettings,
)
from meta_standards_converter.runtime_contracts import ResourceProfile, get_resource_profile
from meta_standards_converter.xml_safety import (
    XMLSizeLimitError,
    parse_xml,
    read_limited_response,
)

from .study_models import ProviderDocument, StudyFetchResult


logger = logging.getLogger(__name__)

NCBI_ENUMERATION_PAGE_SIZE = 10_000
NCBI_EFETCH_BATCH_SIZE = 100
ENA_LINK_PAGE_SIZE = 1_000
ENA_BROWSER_BATCH_SIZE = 500
LINKED_BATCH_SIZE = 200
MAX_UNIQUE_RUNS_HARD = 1_000_000
MAX_PAGES_HARD = 100_000
MAX_PROVIDER_REQUESTS_HARD = 100_000

_SUPPORTED_ACCESSION = re.compile(
    r"^(?:PRJ[A-Z]+\d+|SAM[A-Z]+\d+|[SED]RP\d+|[SED]RS\d+|[SED]RX\d+|[SED]RR\d+)$",
    re.IGNORECASE,
)
_STUDY_ACCESSION = re.compile(r"^[SED]RP\d+$", re.IGNORECASE)
_PROJECT_ACCESSION = re.compile(r"^PRJ[A-Z]+\d+$", re.IGNORECASE)


class StudyRetrievalError(RuntimeError):
    """A provider could not construct a complete core study graph."""


class StudyNotFoundError(StudyRetrievalError):
    def __init__(self, provider: str, accession: str) -> None:
        self.provider = provider
        self.accession = accession
        super().__init__(f"{provider} provider could not resolve accession {accession}")


class PaginationSafetyError(StudyRetrievalError):
    """Pagination violated a termination or forward-progress invariant."""


class ProviderEmergencyLimitError(StudyRetrievalError):
    """A provider crossed a high runaway-protection circuit breaker."""


def _chunks(values: Iterable[str], size: int) -> Iterable[list[str]]:
    batch: list[str] = []
    for value in values:
        batch.append(value)
        if len(batch) == size:
            yield batch
            batch = []
    if batch:
        yield batch


def _public_uri(url: str, params: Mapping[str, Any]) -> str:
    """Render a reproducible response URI without serializing credentials."""

    safe = {key: value for key, value in params.items() if key != "api_key"}
    return f"{url}?{urlencode(safe)}" if safe else url


def collect_accession_pages(
    fetch_page: Callable[[int, int], Iterable[str]],
    *,
    page_size: int,
    max_pages: int = MAX_PAGES_HARD,
    max_unique: int = MAX_UNIQUE_RUNS_HARD,
) -> list[str]:
    """Collect stable unique accessions while enforcing pagination invariants."""

    if page_size < 1 or max_pages < 1 or max_unique < 1:
        raise ValueError("pagination limits must be positive")
    ordered: list[str] = []
    seen: set[str] = set()
    signatures: set[tuple[str, ...]] = set()
    offset = 0
    for _page_index in range(max_pages):
        page = [str(value) for value in fetch_page(offset, page_size) if str(value)]
        if not page:
            return ordered
        signature = tuple(page)
        if signature in signatures:
            raise PaginationSafetyError("provider returned a repeated page signature")
        signatures.add(signature)
        additions = [value for value in page if value not in seen]
        if not additions:
            raise PaginationSafetyError("provider page made no new accessions progress")
        for value in additions:
            seen.add(value)
            ordered.append(value)
        if len(ordered) > max_unique:
            raise ProviderEmergencyLimitError(
                f"provider returned more than {max_unique} unique accessions"
            )
        if len(page) < page_size:
            return ordered
        offset += page_size
    raise ProviderEmergencyLimitError(
        f"provider exceeded the internal {max_pages}-page safety guard"
    )


class _BaseStudyFetcher:
    provider = "provider"

    def __init__(
        self,
        *,
        resource_profile: str | ResourceProfile = "standard",
        resource_overrides=None,
    ) -> None:
        self.resource_profile = get_resource_profile(
            resource_profile, overrides=resource_overrides
        )
        self._request_count = 0
        self._downloaded_bytes = 0

    @staticmethod
    def validate_accession(accession: str) -> str:
        normalized = str(accession).strip().upper()
        if not _SUPPORTED_ACCESSION.fullmatch(normalized):
            raise ValueError(
                "expected a supported INSDC accession: PRJ*, SAM*, "
                "[SED]RP, [SED]RS, [SED]RX, or [SED]RR"
            )
        return normalized

    def _register_document(self, document: ProviderDocument) -> ProviderDocument:
        self._downloaded_bytes += len(document.content)
        if self._downloaded_bytes > self.resource_profile.max_aggregate_download_bytes:
            raise ProviderEmergencyLimitError(
                "provider metadata exceeded the configured aggregate download limit"
            )
        return document

    def _before_request(self) -> None:
        self._request_count += 1
        if self._request_count > MAX_PROVIDER_REQUESTS_HARD:
            raise ProviderEmergencyLimitError(
                f"provider exceeded the internal {MAX_PROVIDER_REQUESTS_HARD}-request safety guard"
            )

    def _read(self, response) -> bytes:
        response.raise_for_status()
        return read_limited_response(
            response,
            max_bytes=self.resource_profile.max_xml_bytes,
        )


class SRAStudyFetcher(_BaseStudyFetcher):
    """Resolve any supported accession and fetch complete NCBI SRA packages."""

    provider = "sra"
    eutils = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

    def __init__(
        self,
        *,
        requester=None,
        ncbi_identity: NCBIApplicationIdentity | None = None,
        resource_profile: str | ResourceProfile = "standard",
        resource_overrides=None,
    ) -> None:
        super().__init__(
            resource_profile=resource_profile, resource_overrides=resource_overrides
        )
        self.identity = ncbi_identity or NCBIApplicationIdentity()
        self.requester = requester or RateLimitedRequester(
            service="ncbi_eutils",
            settings=RequestSettings.from_resource_profile(
                self.resource_profile, request_delay=0.5
            ),
        )

    @staticmethod
    def search_field(accession: str) -> str:
        accession = SRAStudyFetcher.validate_accession(accession)
        if accession.startswith("PRJ"):
            return "GPRJ"
        if accession.startswith("SAM"):
            return "BSPL"
        return "ACCN"

    def _get_bytes(self, endpoint: str, params: dict[str, Any]) -> tuple[bytes, str]:
        url = f"{self.eutils}/{endpoint}"
        merged = {**params, **self.identity.params()}
        self._before_request()
        response = self.requester.get(url, params=merged, stream=True)
        return self._read(response), _public_uri(url, merged)

    def _search(
        self,
        accession: str,
        *,
        field: str | None = None,
        label: str = "resolve",
    ) -> dict[str, Any]:
        content, uri = self._get_bytes(
            "esearch.fcgi",
            {
                "db": "sra",
                "term": f"{accession}[{field or self.search_field(accession)}]",
                "retmode": "json",
                "retmax": NCBI_ENUMERATION_PAGE_SIZE,
                "usehistory": "y",
            },
        )
        payload = json.loads(content)
        result = payload.get("esearchresult") or {}
        if int(result.get("count") or 0) == 0:
            raise StudyNotFoundError(self.provider, accession)
        self._last_search_document = self._register_document(
            ProviderDocument(
                "sra_esearch",
                f"{label}-{accession}-sra-esearch.json",
                uri,
                "application/json",
                content,
            )
        )
        return result

    def _history_documents(
        self,
        search: dict[str, Any],
        *,
        prefix: str,
    ) -> tuple[ProviderDocument, ...]:
        count = int(search.get("count") or 0)
        documents: list[ProviderDocument] = []
        signatures: set[tuple[str, ...]] = set()
        seen_accessions: set[str] = set()
        for offset in range(0, count, NCBI_EFETCH_BATCH_SIZE):
            content, uri = self._get_bytes(
                "efetch.fcgi",
                {
                    "db": "sra",
                    "query_key": search["querykey"],
                    "WebEnv": search["webenv"],
                    "retmode": "xml",
                    "retstart": offset,
                    "retmax": NCBI_EFETCH_BATCH_SIZE,
                },
            )
            root = parse_xml(content, max_bytes=self.resource_profile.max_xml_bytes)
            signature = tuple(
                sorted(
                    {
                        node.get("accession")
                        for path in (".//EXPERIMENT", ".//RUN")
                        for node in root.findall(path)
                        if node.get("accession")
                    }
                )
            )
            if signature in signatures:
                raise PaginationSafetyError(
                    "provider returned a repeated NCBI EFetch batch"
                )
            signatures.add(signature)
            additions = set(signature) - seen_accessions
            if offset and not additions:
                raise PaginationSafetyError(
                    "provider NCBI EFetch batch made no accession progress"
                )
            seen_accessions.update(additions)
            documents.append(
                self._register_document(
                    ProviderDocument(
                        "sra_efetch",
                        f"{prefix}-sra-efetch-{offset:09d}.xml",
                        uri,
                        "application/xml",
                        content,
                    )
                )
            )
        return tuple(documents)

    @staticmethod
    def _study_accessions(documents: Iterable[ProviderDocument]) -> list[str]:
        studies: set[str] = set()
        for document in documents:
            root = parse_xml(document.content, max_bytes=max(len(document.content), 1))
            for node in root.findall(".//STUDY") + root.findall(".//STUDY_REF"):
                value = node.get("accession") or node.get("refcenter")
                if value and _STUDY_ACCESSION.fullmatch(value):
                    studies.add(value.upper())
            for node in root.findall(".//STUDY/IDENTIFIERS/SECONDARY_ID"):
                value = (node.text or "").strip().upper()
                if _STUDY_ACCESSION.fullmatch(value):
                    studies.add(value)
        return sorted(studies)

    @staticmethod
    def _linked_accessions(documents: Iterable[ProviderDocument]) -> dict[str, list[str]]:
        values = {"biosample": [], "bioproject": [], "pubmed": []}
        patterns = {
            "biosample": re.compile(r"^SAM[A-Z]+\d+$", re.I),
            "bioproject": _PROJECT_ACCESSION,
            "pubmed": re.compile(r"^\d+$"),
        }
        for document in documents:
            if document.media_type != "application/xml":
                continue
            root = parse_xml(document.content, max_bytes=max(len(document.content), 1))
            for node in root.iter():
                text = (node.text or "").strip()
                namespace = (node.get("namespace") or node.get("db") or "").casefold()
                for key, pattern in patterns.items():
                    if not pattern.fullmatch(text):
                        continue
                    if key == "pubmed" and namespace not in {"pubmed", "pmid"}:
                        continue
                    if text not in values[key]:
                        values[key].append(text)
            for xref in root.findall(".//XREF_LINK"):
                database = (xref.findtext("./DB") or "").strip().casefold()
                identifier = (xref.findtext("./ID") or "").strip()
                if database in {"pubmed", "pmid"} and identifier.isdigit():
                    if identifier not in values["pubmed"]:
                        values["pubmed"].append(identifier)
        return values

    def _linked_documents(
        self, db: str, accessions: list[str], *, kind: str, media_type: str
    ) -> tuple[ProviderDocument, ...]:
        documents: list[ProviderDocument] = []
        for index, batch in enumerate(_chunks(accessions, LINKED_BATCH_SIZE)):
            endpoint = "esummary.fcgi" if db == "pubmed" else "efetch.fcgi"
            content, uri = self._get_bytes(
                endpoint,
                {"db": db, "id": ",".join(batch), "retmode": "xml"},
            )
            if media_type == "application/xml":
                parse_xml(content, max_bytes=self.resource_profile.max_xml_bytes)
            documents.append(
                self._register_document(
                    ProviderDocument(
                        kind,
                        f"{kind}-{index:06d}.xml",
                        uri,
                        media_type,
                        content,
                    )
                )
            )
        return tuple(documents)

    def fetch(self, accession: str) -> list[StudyFetchResult]:
        accession = self.validate_accession(accession)
        self._request_count = 0
        self._downloaded_bytes = 0
        resolution_search = self._search(accession, label="resolve")
        resolution_search_document = self._last_search_document
        resolution = self._history_documents(
            resolution_search, prefix=f"resolve-{accession}"
        )
        studies = self._study_accessions(resolution)
        if not studies:
            raise StudyNotFoundError(self.provider, accession)

        results: list[StudyFetchResult] = []
        for study in studies:
            study_search = self._search(study, field="ACCN", label="study")
            study_search_document = self._last_search_document
            documents = [
                resolution_search_document,
                *resolution,
                study_search_document,
                *self._history_documents(
                    study_search, prefix=study
                ),
            ]
            run_accessions = {
                node.get("accession")
                for document in documents
                if document.kind == "sra_efetch"
                for node in parse_xml(
                    document.content, max_bytes=max(len(document.content), 1)
                ).findall(".//RUN")
                if node.get("accession")
            }
            if len(run_accessions) > MAX_UNIQUE_RUNS_HARD:
                raise ProviderEmergencyLimitError(
                    f"provider returned more than {MAX_UNIQUE_RUNS_HARD} unique runs"
                )
            linked = self._linked_accessions(documents)
            optional_warnings: list[dict] = []
            for db, kind in (
                ("biosample", "biosample_efetch"),
                ("bioproject", "bioproject_efetch"),
                ("pubmed", "pubmed_esummary"),
            ):
                if not linked[db]:
                    continue
                try:
                    documents.extend(
                        self._linked_documents(
                            db, linked[db], kind=kind, media_type="application/xml"
                        )
                    )
                except Exception as error:  # optional linked metadata degrades only
                    optional_warnings.append(
                        {
                            "code": f"missing_{db}_metadata",
                            "message": f"Optional {db} metadata could not be retrieved.",
                            "error_type": type(error).__name__,
                        }
                    )
            results.append(
                StudyFetchResult(
                    self.provider,
                    accession,
                    study,
                    tuple(documents),
                    tuple(optional_warnings),
                )
            )
        return results


class ENAStudyFetcher(_BaseStudyFetcher):
    """Resolve and retrieve study graphs from ENA Portal and Browser APIs."""

    provider = "ena"
    portal = "https://www.ebi.ac.uk/ena/portal/api"
    browser = "https://www.ebi.ac.uk/ena/browser/api"

    def __init__(
        self,
        *,
        requester=None,
        ncbi_requester=None,
        ncbi_identity: NCBIApplicationIdentity | None = None,
        resource_profile: str | ResourceProfile = "standard",
        resource_overrides=None,
    ) -> None:
        super().__init__(
            resource_profile=resource_profile, resource_overrides=resource_overrides
        )
        self.requester = requester or RateLimitedRequester(
            service="ena_portal",
            settings=RequestSettings.from_resource_profile(
                self.resource_profile, request_delay=1.0
            ),
        )
        self.ncbi_requester = ncbi_requester or RateLimitedRequester(
            service="ncbi_eutils",
            settings=RequestSettings.from_resource_profile(
                self.resource_profile, request_delay=0.5
            ),
        )
        self.ncbi_identity = ncbi_identity or NCBIApplicationIdentity()
        self._aux_documents: list[ProviderDocument] = []

    def _request(self, method: str, url: str, **kwargs) -> tuple[bytes, str]:
        self._before_request()
        caller = getattr(self.requester, method)
        response = caller(url, stream=True, **kwargs)
        content = self._read(response)
        params = kwargs.get("params") or {}
        uri = _public_uri(url, params)
        return content, uri

    def _portal_rows(self, result: str, query: str, fields: str) -> list[dict[str, Any]]:
        content, uri = self._request(
            "get",
            f"{self.portal}/search",
            params={
                "result": result,
                "query": query,
                "fields": fields,
                "format": "json",
                # Portal search uses ``limit=0`` for an untruncated result set.
                # Study hierarchy enumeration is separately paginated through
                # /links/study below.
                "limit": 0,
            },
        )
        payload = json.loads(content)
        self._aux_documents.append(
            self._register_document(
                ProviderDocument(
                    "ena_portal_search",
                    f"ena-portal-search-{len(self._aux_documents):06d}.json",
                    uri,
                    "application/json",
                    content,
                )
            )
        )
        return payload if isinstance(payload, list) else []

    @staticmethod
    def _field_values(rows: Iterable[dict[str, Any]], *names: str) -> list[str]:
        values: list[str] = []
        for row in rows:
            for name in names:
                for value in str(row.get(name) or "").split(";"):
                    normalized = value.strip().upper()
                    if normalized and normalized not in values:
                        values.append(normalized)
        return values

    def _resolve(self, accession: str) -> list[tuple[str, str]]:
        if _PROJECT_ACCESSION.fullmatch(accession):
            rows = self._portal_rows(
                "study", f'study_accession="{accession}"',
                "study_accession,secondary_study_accession",
            )
        elif _STUDY_ACCESSION.fullmatch(accession):
            rows = self._portal_rows(
                "study", f'secondary_study_accession="{accession}"',
                "study_accession,secondary_study_accession",
            )
        elif re.fullmatch(r"SAM[A-Z]+\d+|[SED]RS\d+", accession):
            field = "sample_accession" if accession.startswith("SAM") else "secondary_sample_accession"
            rows = self._portal_rows(
                "read_experiment", f'{field}="{accession}"',
                "study_accession,secondary_study_accession",
            )
        elif re.fullmatch(r"[SED]RX\d+", accession):
            rows = self._portal_rows(
                "read_experiment", f'experiment_accession="{accession}"',
                "study_accession,secondary_study_accession",
            )
        else:
            rows = self._portal_rows(
                "read_run", f'run_accession="{accession}"',
                "study_accession,secondary_study_accession",
            )
        pairs: list[tuple[str, str]] = []
        for row in rows:
            projects = self._field_values([row], "study_accession")
            studies = self._field_values([row], "secondary_study_accession")
            if _STUDY_ACCESSION.fullmatch(accession):
                studies = [study for study in studies if study == accession]
            for study in studies:
                pair = (study, projects[0] if projects else "")
                if pair not in pairs:
                    pairs.append(pair)
        return pairs

    def _links(self, project: str, result: str) -> list[str]:
        def page(offset: int, limit: int) -> Iterable[str]:
            content, uri = self._request(
                "get",
                f"{self.portal}/links/study",
                params={
                    "accession": project,
                    "result": result,
                    "format": "json",
                    "offset": offset,
                    "limit": limit,
                },
            )
            payload = json.loads(content)
            self._aux_documents.append(
                self._register_document(
                    ProviderDocument(
                        "ena_portal_links",
                        f"ena-{result}-links-{offset:09d}.json",
                        uri,
                        "application/json",
                        content,
                    )
                )
            )
            if not isinstance(payload, list):
                return []
            candidates: list[str] = []
            for row in payload:
                if isinstance(row, str):
                    candidates.append(row)
                elif isinstance(row, dict):
                    candidates.extend(
                        self._field_values(
                            [row],
                            "accession",
                            "experiment_accession",
                            "run_accession",
                            "analysis_accession",
                            "assembly_accession",
                        )
                    )
            return candidates

        return collect_accession_pages(
            page,
            page_size=ENA_LINK_PAGE_SIZE,
            max_unique=MAX_UNIQUE_RUNS_HARD,
        )

    def _browser_batch(
        self, accessions: list[str], *, kind: str
    ) -> list[ProviderDocument]:
        if not accessions:
            return []
        try:
            content, uri = self._request(
                "post",
                f"{self.browser}/xml",
                json={"accessions": accessions, "includeLinks": True, "set": True},
            )
        except XMLSizeLimitError:
            if len(accessions) == 1:
                raise
            midpoint = len(accessions) // 2
            return [
                *self._browser_batch(accessions[:midpoint], kind=kind),
                *self._browser_batch(accessions[midpoint:], kind=kind),
            ]
        parse_xml(content, max_bytes=self.resource_profile.max_xml_bytes)
        name = f"{kind}-{'-'.join(accessions[:1])}-{len(accessions):06d}.xml"
        return [
            self._register_document(
                ProviderDocument(kind, name, uri, "application/xml", content)
            )
        ]

    def _browser_documents(
        self, accessions: list[str], *, kind: str
    ) -> list[ProviderDocument]:
        documents: list[ProviderDocument] = []
        for batch in _chunks(accessions, ENA_BROWSER_BATCH_SIZE):
            documents.extend(self._browser_batch(batch, kind=kind))
        return documents

    def _simple_document(
        self,
        url: str,
        *,
        kind: str,
        name: str,
        params: dict[str, Any] | None = None,
    ) -> ProviderDocument:
        content, uri = self._request("get", url, params=params or {})
        return self._register_document(
            ProviderDocument(
                kind,
                name,
                uri,
                "application/json" if name.endswith(".json") else "application/xml",
                content,
            )
        )

    def _pubmed_documents(self, accessions: list[str]) -> list[ProviderDocument]:
        documents: list[ProviderDocument] = []
        for index, batch in enumerate(_chunks(accessions, LINKED_BATCH_SIZE)):
            url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
            params = {
                "db": "pubmed",
                "id": ",".join(batch),
                "retmode": "xml",
                **self.ncbi_identity.params(),
            }
            self._before_request()
            response = self.ncbi_requester.get(url, params=params, stream=True)
            content = self._read(response)
            parse_xml(content, max_bytes=self.resource_profile.max_xml_bytes)
            documents.append(
                self._register_document(
                    ProviderDocument(
                        "pubmed_esummary",
                        f"pubmed-esummary-{index:06d}.xml",
                        _public_uri(url, params),
                        "application/xml",
                        content,
                    )
                )
            )
        return documents

    @staticmethod
    def _taxids(documents: Iterable[ProviderDocument]) -> list[str]:
        values: list[str] = []
        for document in documents:
            if document.kind != "ena_sample":
                continue
            root = parse_xml(document.content, max_bytes=max(len(document.content), 1))
            for node in root.findall(".//TAXON_ID"):
                value = (node.text or "").strip()
                if value.isdigit() and value not in values:
                    values.append(value)
        return values

    @staticmethod
    def _pmids(documents: Iterable[ProviderDocument]) -> list[str]:
        values: list[str] = []
        for document in documents:
            if document.kind == "ena_study":
                root = parse_xml(document.content, max_bytes=max(len(document.content), 1))
                for xref in root.findall(".//XREF_LINK"):
                    database = (xref.findtext("./DB") or "").strip().casefold()
                    identifier = (xref.findtext("./ID") or "").strip()
                    if database in {"pubmed", "pmid"} and identifier.isdigit():
                        if identifier not in values:
                            values.append(identifier)
            elif document.kind == "ena_xref":
                payload = json.loads(document.content)
                for row in payload if isinstance(payload, list) else []:
                    if not isinstance(row, dict):
                        continue
                    rendered = " ".join(str(item) for item in row.values())
                    if "pubmed" not in rendered.casefold():
                        continue
                    for identifier in re.findall(r"\b\d{5,9}\b", rendered):
                        if identifier not in values:
                            values.append(identifier)
        return values

    @staticmethod
    def _document_study_accessions(document: ProviderDocument) -> set[str]:
        root = parse_xml(document.content, max_bytes=max(len(document.content), 1))
        values = set()
        for node in root.findall(".//STUDY_REF"):
            value = node.get("accession")
            if value:
                values.add(value.upper())
        return values

    @staticmethod
    def _experiment_sample_accessions(document: ProviderDocument) -> list[str]:
        root = parse_xml(document.content, max_bytes=max(len(document.content), 1))
        values: list[str] = []
        for node in root.findall(".//SAMPLE_DESCRIPTOR"):
            value = node.get("accession")
            if value and value not in values:
                values.append(value)
        return values

    def fetch(self, accession: str) -> list[StudyFetchResult]:
        accession = self.validate_accession(accession)
        self._request_count = 0
        self._downloaded_bytes = 0
        self._aux_documents = []
        resolved = self._resolve(accession)
        if not resolved:
            raise StudyNotFoundError(self.provider, accession)
        resolution_documents = tuple(self._aux_documents)
        results: list[StudyFetchResult] = []
        for study, project in sorted(set(resolved)):
            documents: list[ProviderDocument] = list(resolution_documents)
            warnings: list[dict] = []
            documents.extend(self._browser_documents([study], kind="ena_study"))
            link_document_start = len(self._aux_documents)
            experiment_accessions = self._links(project, "read_experiment") if project else []
            run_accessions = self._links(project, "read_run") if project else []
            provider_object_accessions: dict[str, list[str]] = {}
            for result_kind in ("analysis", "assembly"):
                try:
                    provider_object_accessions[result_kind] = (
                        self._links(project, result_kind) if project else []
                    )
                except Exception as error:
                    provider_object_accessions[result_kind] = []
                    warnings.append(
                        {
                            "code": f"missing_ena_{result_kind}_links",
                            "message": f"Optional ENA {result_kind} links could not be retrieved.",
                            "error_type": type(error).__name__,
                        }
                    )
            documents.extend(self._aux_documents[link_document_start:])
            if len(run_accessions) > MAX_UNIQUE_RUNS_HARD:
                raise ProviderEmergencyLimitError(
                    f"provider returned more than {MAX_UNIQUE_RUNS_HARD} unique runs"
                )
            project_experiments = self._browser_documents(
                experiment_accessions, kind="ena_experiment"
            )
            study_experiments = [
                item
                for item in project_experiments
                if not self._document_study_accessions(item)
                or study in self._document_study_accessions(item)
            ]
            documents.extend(study_experiments)
            # Run XML is retained for exact aliases/statistics; parser filters joins.
            documents.extend(self._browser_documents(run_accessions, kind="ena_run"))
            for result_kind, accessions in provider_object_accessions.items():
                try:
                    documents.extend(
                        self._browser_documents(
                            accessions, kind=f"ena_{result_kind}"
                        )
                    )
                except Exception as error:
                    warnings.append(
                        {
                            "code": f"missing_ena_{result_kind}_metadata",
                            "message": f"Optional ENA {result_kind} metadata could not be retrieved.",
                            "error_type": type(error).__name__,
                        }
                    )
            samples: list[str] = []
            for document in study_experiments:
                for sample in self._experiment_sample_accessions(document):
                    if sample not in samples:
                        samples.append(sample)
            documents.extend(self._browser_documents(samples, kind="ena_sample"))

            optional = (
                (
                    f"{self.portal}/filereport",
                    "ena_file_report",
                    f"{study}-file-report.json",
                    {
                        "accession": study,
                        "result": "read_run",
                        "format": "json",
                        "fields": "run_accession,experiment_accession,submitted_ftp,submitted_md5,submitted_bytes,fastq_ftp,fastq_md5,fastq_bytes,sra_ftp,sra_md5,sra_bytes",
                    },
                ),
                (
                    "https://www.ebi.ac.uk/ena/xref/rest/json/search",
                    "ena_xref",
                    f"{study}-xref.json",
                    {"accession": study},
                ),
            )
            for url, kind, name, params in optional:
                try:
                    documents.append(
                        self._simple_document(url, kind=kind, name=name, params=params)
                    )
                except Exception as error:
                    warnings.append(
                        {
                            "code": f"missing_{kind}",
                            "message": f"Optional ENA {kind} response could not be retrieved.",
                            "error_type": type(error).__name__,
                        }
                    )
            for taxid in self._taxids(documents):
                try:
                    documents.append(
                        self._simple_document(
                            f"https://www.ebi.ac.uk/ena/taxonomy/rest/tax-id/{taxid}",
                            kind="ena_taxonomy",
                            name=f"taxonomy-{taxid}.json",
                        )
                    )
                except Exception as error:
                    warnings.append(
                        {
                            "code": "missing_ena_taxonomy",
                            "message": f"Optional ENA taxonomy response for {taxid} could not be retrieved.",
                            "error_type": type(error).__name__,
                        }
                    )
            pmids = self._pmids(documents)
            if pmids:
                try:
                    documents.extend(self._pubmed_documents(pmids))
                except Exception as error:
                    warnings.append(
                        {
                            "code": "missing_pubmed_metadata",
                            "message": "Optional PubMed summaries could not be retrieved.",
                            "error_type": type(error).__name__,
                        }
                    )
            results.append(
                StudyFetchResult(
                    self.provider,
                    accession,
                    study,
                    tuple(documents),
                    tuple(warnings),
                )
            )
        return results


__all__ = [
    "ENA_BROWSER_BATCH_SIZE",
    "ENA_LINK_PAGE_SIZE",
    "LINKED_BATCH_SIZE",
    "MAX_PAGES_HARD",
    "MAX_PROVIDER_REQUESTS_HARD",
    "MAX_UNIQUE_RUNS_HARD",
    "NCBI_EFETCH_BATCH_SIZE",
    "NCBI_ENUMERATION_PAGE_SIZE",
    "ENAStudyFetcher",
    "PaginationSafetyError",
    "ProviderEmergencyLimitError",
    "SRAStudyFetcher",
    "StudyNotFoundError",
    "StudyRetrievalError",
    "collect_accession_pages",
]
