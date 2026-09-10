# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Resolve local, remote, and BioStudies MAGE-TAB metadata sources."""

from __future__ import annotations

import csv
from dataclasses import dataclass
import io
import json
import os
from pathlib import Path
from urllib.parse import quote, urljoin, urlparse

from meta_standards_converter.helpers.request_helper import (
    RateLimitedRequester,
    RequestSettings,
)
from meta_standards_converter.retrieval import RetrievalPolicy, RetrievalSecurityError
from meta_standards_converter.runtime_contracts import (
    ResourceProfile,
    get_resource_profile,
)
from meta_standards_converter.xml_safety import (
    XMLSizeLimitError,
    read_limited_response,
)


@dataclass(frozen=True)
class TextResource:
    name: str
    text: str
    origin: str


@dataclass(frozen=True)
class MAGETabInput:
    idf: TextResource
    sdrfs: tuple[TextResource, ...]
    source: str
    source_kind: str


class AEWebFetcher:
    """Load an IDF and its SDRFs without persisting remote metadata files."""

    def metrics(self):
        from meta_standards_converter.sources.contracts import request_metrics
        return request_metrics(self.requester)


    API_ROOT = "https://www.ebi.ac.uk/biostudies/api/v1"
    FILE_PAGE_SIZE = 100
    REDIRECT_STATUSES = frozenset({301, 302, 303, 307, 308})

    def __init__(
        self,
        requester=None,
        request_settings=None,
        resource_profile: str | ResourceProfile = "standard",
        resource_overrides=None,
        retrieval_policy: RetrievalPolicy | None = None,
    ):
        self.resource_profile = get_resource_profile(
            resource_profile,
            overrides=resource_overrides,
        )
        self.retrieval_policy = retrieval_policy or RetrievalPolicy(
            resource_profile=self.resource_profile,
            allowed_schemes=frozenset({"https"}),
        )
        self.requester = requester or RateLimitedRequester(
            service="biostudies",
            settings=request_settings
            or RequestSettings.from_resource_profile(
                self.resource_profile,
                request_delay=1.0,
            ),
        )
        self._downloaded_bytes = 0

    def resolve(self, source: str, sdrf_sources: list[str] | None = None) -> MAGETabInput:
        self._downloaded_bytes = 0
        if os.path.exists(source):
            return self._resolve_local(source, sdrf_sources=sdrf_sources)
        if self._is_http(source):
            return self._resolve_http(source, sdrf_sources=sdrf_sources)
        if self._looks_like_missing_file(source):
            raise FileNotFoundError(f"MAGE-TAB IDF file not found: {source}")
        if sdrf_sources:
            raise ValueError("Explicit SDRF overrides cannot be used with an accession source.")
        return self._resolve_accession(source)

    def _resolve_local(self, source: str, sdrf_sources: list[str] | None) -> MAGETabInput:
        idf_path = Path(source).resolve()
        idf = TextResource(idf_path.name, self._read_local(idf_path), str(idf_path))
        references = sdrf_sources or self._sdrf_references(idf.text)
        if not references:
            raise ValueError(f"MAGE-TAB IDF {source} does not reference an SDRF.")
        resources = []
        for reference in references:
            if self._is_http(reference):
                resources.append(self._fetch_http(reference))
                continue
            path = Path(reference)
            if not path.is_absolute():
                path = idf_path.parent / path
            path = path.resolve()
            resources.append(TextResource(path.name, self._read_local(path), str(path)))
        return MAGETabInput(idf, tuple(resources), str(idf_path), "path")

    def _resolve_http(self, source: str, sdrf_sources: list[str] | None) -> MAGETabInput:
        idf = self._fetch_http(source)
        references = sdrf_sources or self._sdrf_references(idf.text)
        if not references:
            raise ValueError(f"MAGE-TAB IDF {source} does not reference an SDRF.")
        resources = [
            self._fetch_http(reference if self._is_http(reference) else urljoin(source, reference))
            for reference in references
        ]
        return MAGETabInput(idf, tuple(resources), source, "url")

    def _resolve_accession(self, accession: str) -> MAGETabInput:
        accession = accession.strip().upper()
        files_url = f"{self.API_ROOT}/files/{quote(accession, safe='')}"
        rows = self._biostudies_file_rows(files_url)
        idf_rows = [row for row in rows if self._file_kind(row) == "idf"]
        sdrf_rows = [row for row in rows if self._file_kind(row) == "sdrf"]
        if len(idf_rows) != 1:
            raise ValueError(
                f"BioStudies accession {accession} must expose exactly one IDF; found {len(idf_rows)}."
            )
        if not sdrf_rows:
            raise ValueError(f"BioStudies accession {accession} exposes no SDRF files.")

        info_url = f"{self.API_ROOT}/studies/{quote(accession, safe='')}/info"
        info = self._fetch_json(info_url)
        base_url = info.get("httpLink") if isinstance(info, dict) else None
        if not base_url or not self._is_http(base_url):
            raise ValueError(f"BioStudies accession {accession} has no HTTP download link.")

        idf = self._fetch_api_file(base_url, idf_rows[0])
        sdrfs = tuple(self._fetch_api_file(base_url, row) for row in sdrf_rows)
        return MAGETabInput(idf, sdrfs, accession, "accession")

    def _biostudies_file_rows(self, files_url: str) -> list[dict]:
        rows: list[dict] = []
        while True:
            payload = self._fetch_json(
                files_url,
                params={"start": len(rows), "length": self.FILE_PAGE_SIZE},
            )
            page = payload.get("data", []) if isinstance(payload, dict) else []
            page = [row for row in page if isinstance(row, dict)]
            rows.extend(page)

            total = payload.get("recordsFiltered") if isinstance(payload, dict) else None
            try:
                total = int(total) if total is not None else None
            except (TypeError, ValueError):
                total = None
            if not page or total is None or len(rows) >= total:
                return rows

    def _fetch_api_file(self, base_url: str, row: dict) -> TextResource:
        path = str(row.get("path") or row.get("Name") or "").lstrip("/")
        if not path:
            raise ValueError("BioStudies file metadata has no path.")
        url = f"{base_url.rstrip('/')}/Files/{quote(path, safe='/')}"
        resource = self._fetch_http(url)
        name = str(row.get("Name") or os.path.basename(path))
        return TextResource(name, resource.text, resource.origin)

    def _fetch_http(self, url: str) -> TextResource:
        response, resolved_url = self._request(url)
        try:
            response.raise_for_status()
            content = self._read_response(response)
        finally:
            response.close()
        text = content.decode("utf-8-sig")
        name = os.path.basename(urlparse(resolved_url).path) or "metadata.txt"
        return TextResource(name, text, resolved_url)

    def _fetch_json(self, url: str, **kwargs):
        response, _ = self._request(url, **kwargs)
        try:
            response.raise_for_status()
            content = self._read_response(response)
        finally:
            response.close()
        try:
            return json.loads(content.decode("utf-8-sig"))
        except (UnicodeDecodeError, json.JSONDecodeError) as error:
            raise ValueError("BioStudies returned invalid UTF-8 JSON metadata.") from error

    def _request(self, url: str, **kwargs):
        current = url
        request_kwargs = dict(kwargs)
        for redirect_count in range(self.resource_profile.max_redirects + 1):
            self.retrieval_policy.validate_url(current)
            response = self.requester.get(
                current,
                stream=True,
                allow_redirects=False,
                **request_kwargs,
            )
            if response.status_code not in self.REDIRECT_STATUSES:
                return response, current
            location = response.headers.get("Location")
            response.close()
            if not location:
                raise RetrievalSecurityError(
                    "MAGE-TAB metadata redirect omitted Location."
                )
            if redirect_count >= self.resource_profile.max_redirects:
                raise RetrievalSecurityError(
                    "MAGE-TAB metadata exceeded the redirect limit."
                )
            target = urljoin(current, location)
            if urlparse(target).scheme.casefold() != urlparse(current).scheme.casefold():
                raise RetrievalSecurityError(
                    "Cross-scheme MAGE-TAB metadata redirects are not allowed."
                )
            self.retrieval_policy.validate_url(target)
            current = target
            request_kwargs = {}
        raise RetrievalSecurityError("MAGE-TAB metadata exceeded the redirect limit.")

    def _read_response(self, response) -> bytes:
        content = read_limited_response(
            response,
            max_bytes=self.resource_profile.max_xml_bytes,
        )
        self._consume_bytes(len(content))
        return content

    def _consume_bytes(self, byte_count: int) -> None:
        aggregate = self._downloaded_bytes + byte_count
        if aggregate > self.resource_profile.max_aggregate_download_bytes:
            raise XMLSizeLimitError(
                "MAGE-TAB metadata exceeded the aggregate run download limit."
            )
        self._downloaded_bytes = aggregate

    def _read_local(self, path: Path) -> str:
        if not path.is_file():
            raise FileNotFoundError(f"MAGE-TAB metadata file not found: {path}")
        declared = path.stat().st_size
        if declared > self.resource_profile.max_xml_bytes:
            raise XMLSizeLimitError(
                f"MAGE-TAB metadata declares {declared} bytes and exceeds the "
                f"{self.resource_profile.max_xml_bytes} byte file limit."
            )
        with path.open("rb") as handle:
            content = handle.read(self.resource_profile.max_xml_bytes + 1)
        if len(content) > self.resource_profile.max_xml_bytes:
            raise XMLSizeLimitError(
                "MAGE-TAB metadata exceeded the configured file limit."
            )
        return content.decode("utf-8-sig")

    def _sdrf_references(self, idf_text: str) -> list[str]:
        references = []
        for row in csv.reader(io.StringIO(idf_text), delimiter="\t"):
            if row and self._normalized_label(row[0]) == "sdrffile":
                references.extend(value.strip() for value in row[1:] if value.strip())
        return references

    def _file_kind(self, row: dict) -> str | None:
        file_type = str(row.get("Type") or "").casefold()
        name = str(row.get("Name") or row.get("path") or "").casefold()
        if "idf file" in file_type or name.endswith((".idf.txt", "_idf.txt")):
            return "idf"
        if "sdrf file" in file_type or name.endswith((".sdrf.txt", "_sdrf.txt")):
            return "sdrf"
        return None

    def _normalized_label(self, value: str) -> str:
        return "".join(str(value).split()).casefold()

    def _is_http(self, value: str) -> bool:
        return urlparse(str(value)).scheme.casefold() in {"http", "https"}

    def _looks_like_missing_file(self, value: str) -> bool:
        lowered = str(value).casefold()
        return os.path.sep in value or lowered.endswith((".idf.txt", "_idf.txt"))
