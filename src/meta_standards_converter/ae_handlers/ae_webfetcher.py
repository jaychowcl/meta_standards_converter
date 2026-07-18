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
import os
from pathlib import Path
from urllib.parse import quote, urljoin, urlparse

from meta_standards_converter.helpers.request_helper import RateLimitedRequester


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

    API_ROOT = "https://www.ebi.ac.uk/biostudies/api/v1"

    def __init__(self, requester=None, request_settings=None):
        self.requester = requester or RateLimitedRequester(
            service="biostudies",
            settings=request_settings,
        )

    def resolve(self, source: str, sdrf_sources: list[str] | None = None) -> MAGETabInput:
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
        files_response = self.requester.get(files_url)
        files_response.raise_for_status()
        payload = files_response.json()
        rows = payload.get("data", []) if isinstance(payload, dict) else []
        idf_rows = [row for row in rows if self._file_kind(row) == "idf"]
        sdrf_rows = [row for row in rows if self._file_kind(row) == "sdrf"]
        if len(idf_rows) != 1:
            raise ValueError(
                f"BioStudies accession {accession} must expose exactly one IDF; found {len(idf_rows)}."
            )
        if not sdrf_rows:
            raise ValueError(f"BioStudies accession {accession} exposes no SDRF files.")

        info_url = f"{self.API_ROOT}/studies/{quote(accession, safe='')}/info"
        info_response = self.requester.get(info_url)
        info_response.raise_for_status()
        info = info_response.json()
        base_url = info.get("httpLink") if isinstance(info, dict) else None
        if not base_url or not self._is_http(base_url):
            raise ValueError(f"BioStudies accession {accession} has no HTTP download link.")

        idf = self._fetch_api_file(base_url, idf_rows[0])
        sdrfs = tuple(self._fetch_api_file(base_url, row) for row in sdrf_rows)
        return MAGETabInput(idf, sdrfs, accession, "accession")

    def _fetch_api_file(self, base_url: str, row: dict) -> TextResource:
        path = str(row.get("path") or row.get("Name") or "").lstrip("/")
        if not path:
            raise ValueError("BioStudies file metadata has no path.")
        url = f"{base_url.rstrip('/')}/Files/{quote(path, safe='/')}"
        resource = self._fetch_http(url)
        name = str(row.get("Name") or os.path.basename(path))
        return TextResource(name, resource.text, url)

    def _fetch_http(self, url: str) -> TextResource:
        response = self.requester.get(url)
        response.raise_for_status()
        content = getattr(response, "content", None)
        if content:
            text = content.decode("utf-8-sig")
        else:
            text = response.text
            if text.startswith("\ufeff"):
                text = text.lstrip("\ufeff")
        name = os.path.basename(urlparse(url).path) or "metadata.txt"
        return TextResource(name, text, url)

    def _read_local(self, path: Path) -> str:
        if not path.is_file():
            raise FileNotFoundError(f"MAGE-TAB metadata file not found: {path}")
        return path.read_text(encoding="utf-8-sig")

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
