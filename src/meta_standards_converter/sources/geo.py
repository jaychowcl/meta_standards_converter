# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""
Fetches data from GEO
"""

import logging
import tarfile
import tempfile
import time
from pathlib import PurePosixPath

from meta_standards_converter.helpers.request_helper import (
    RateLimitedRequester,
    RequestSettings,
)
from meta_standards_converter.runtime_contracts import (
    get_resource_profile,
    require_disk_headroom,
)
from meta_standards_converter.xml_safety import parse_xml, stream_limited_response


logger = logging.getLogger(__name__)
MAX_GEO_ARCHIVE_MEMBERS = 10_000


def _normalise_archive_member_name(name: str) -> str:
    """Return a safe canonical POSIX member name or reject it."""

    if not name or name.startswith("/") or "\\" in name or "\x00" in name:
        raise ValueError(f"GEO archive contains an unsafe path: {name!r}.")
    parts = name.split("/")
    while parts and parts[0] == ".":
        parts.pop(0)
    if not parts or any(part in {"", ".", ".."} for part in parts):
        raise ValueError(f"GEO archive contains an unsafe path: {name!r}.")
    return PurePosixPath(*parts).as_posix()


class GEOWebFetcher:

    def metrics(self):
        from meta_standards_converter.sources.contracts import request_metrics
        return request_metrics(self.requester)


    def __init__(
        self,
        requester=None,
        request_settings=None,
        resource_profile: str = "standard",
        resource_overrides=None,
    ):
        self.resource_profile = get_resource_profile(
            resource_profile,
            overrides=resource_overrides,
        )
        self.requester = requester or RateLimitedRequester(
            service="geo_ftp",
            settings=request_settings
            or RequestSettings.from_resource_profile(
                self.resource_profile,
                request_delay=1.0,
            ),
        )

    def url_gse_miniml(self, gse: str) -> str:
        """
        creates url from gse accession for fetching gse mininml and returns url as string.
        """
        # checks gse valid by checking gse prefix
        if gse[:3].lower() != "gse":
            raise ValueError(f"GSE accession {gse} is not valid. Must start with GSE.")

        # gets gse_nnn for url
        digits = gse[3:]
        if len(digits) <= 3:
            gse_nnn = gse[:3] + "nnn"
        else:
            gse_nnn = gse[:-3] + "nnn"

        # build url
        url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{gse_nnn}/{gse}/miniml/{gse}_family.xml.tgz"
        return url

    def fetch_gse_miniml(self, gse) -> str:
        """
        creates url from gse accession, fetches miniml file, returns miniml as string.
        """
        # create url for fetching
        url = self.url_gse_miniml(gse=gse)
        started = time.monotonic()
        logger.info("GEO MINiML fetch started accession=%s", gse)

        # use url to fetch miniml file
        response = self.requester.get(url, stream=True)
        response.raise_for_status()
        raw_length = response.headers.get("Content-Length")
        if raw_length in (None, ""):
            raise ValueError("GEO archive response requires Content-Length.")
        try:
            declared_bytes = int(raw_length)
        except (TypeError, ValueError) as error:
            raise ValueError("GEO archive Content-Length is invalid.") from error
        require_disk_headroom(
            tempfile.gettempdir(),
            required_bytes=declared_bytes,
            headroom_fraction=self.resource_profile.disk_headroom_fraction,
        )

        # Stream the archive to disk, then inspect every member without
        # extracting paths onto the filesystem.
        with tempfile.NamedTemporaryFile(suffix=".tgz") as archive:
            archive_bytes = stream_limited_response(
                response,
                archive,
                max_bytes=self.resource_profile.max_compressed_archive_bytes,
            )
            archive.flush()
            with tarfile.open(name=archive.name, mode="r|gz") as tar:
                expected_name = f"{gse}_family.xml"
                xml_limit = min(
                    self.resource_profile.max_xml_bytes,
                    self.resource_profile.max_expanded_archive_bytes,
                )
                seen_names: set[str] = set()
                expanded_bytes = 0
                member_count = 0
                encoded: bytes | None = None
                for member in tar:
                    member_count += 1
                    if member_count > MAX_GEO_ARCHIVE_MEMBERS:
                        raise ValueError(
                            "GEO archive exceeds the member-count limit."
                        )
                    member_name = _normalise_archive_member_name(member.name)
                    if member_name in seen_names:
                        raise ValueError(
                            f"GEO archive contains duplicate member {member_name!r}."
                        )
                    seen_names.add(member_name)

                    if not (member.isdir() or member.isfile()):
                        raise ValueError(
                            "GEO archive contains a link or unsafe member type."
                        )
                    if member.isdir():
                        continue
                    if member.size < 0:
                        raise ValueError("GEO archive member has an invalid size.")
                    expanded_bytes += member.size
                    if (
                        expanded_bytes
                        > self.resource_profile.max_expanded_archive_bytes
                    ):
                        raise ValueError(
                            "GEO archive exceeds the expanded-byte limit."
                        )
                    if member_name.casefold().endswith(".xml"):
                        if member_name != expected_name:
                            raise ValueError(
                                "GEO archive contains an unexpected XML member."
                            )
                        if member.size > xml_limit:
                            raise ValueError(
                                "GEO XML member exceeds the "
                                f"{xml_limit} byte expanded limit."
                            )
                        miniml_file = tar.extractfile(member)
                        if miniml_file is None:
                            raise ValueError(
                                "GEO archive XML member could not be read."
                            )
                        encoded = miniml_file.read(xml_limit + 1)
                        if len(encoded) > xml_limit:
                            raise ValueError(
                                "GEO XML member exceeds the "
                                f"{xml_limit} byte expanded limit."
                            )
                if encoded is None:
                    raise ValueError(
                        "GEO archive is missing the expected XML member."
                    )
        parse_xml(encoded, max_bytes=self.resource_profile.max_xml_bytes)
        miniml = encoded.decode("utf-8")

        logger.info(
            "GEO MINiML fetch completed accession=%s archive_bytes=%s xml_characters=%s elapsed_seconds=%.3f",
            gse,
            archive_bytes,
            len(miniml),
            time.monotonic() - started,
        )

        return miniml


from meta_standards_converter.miniml.geo_parser import GEOParser
from meta_standards_converter.miniml import MINiMLV1Migrator
from collections import deque
from meta_standards_converter.runtime_contracts import (OperationStatusV2, ExecutionStatus, CompletenessStatus, EvidenceConfidence, ValidationStatus, PublicationDisposition, SafeErrorEnvelope)

class RelatedSeriesParseResult(list[dict]):
    """List-compatible related-series result with explicit completeness."""

    def __init__(
        self,
        packages,
        *,
        status: OperationStatusV2,
        attempted_accessions,
        failed_accessions,
    ) -> None:
        super().__init__(packages)
        self.status = status
        self.attempted_accessions = tuple(attempted_accessions)
        self.failed_accessions = tuple(failed_accessions)

    def summary_dict(self) -> dict:
        return {
            "status": self.status.to_dict(),
            "attempted_accessions": list(self.attempted_accessions),
            "failed_accessions": list(self.failed_accessions),
            "package_count": len(self),
        }


class GEOSource:
    """Coordinate GEO retrieval and related-series collection around a pure parser."""

    def metrics(self):
        from meta_standards_converter.sources.contracts import request_metrics
        return request_metrics(self.geo_fetcher)


    def __init__(self, fetcher=None, parser=None, resource_profile="standard"):
        self.geo_fetcher = fetcher or GEOWebFetcher(resource_profile=resource_profile)
        self.parser = parser or GEOParser(resource_profile=resource_profile)
        self.resource_profile = self.parser.resource_profile

    def parse(self, miniml, remove_empty=False, related_series=False):
        if not related_series:
            return self.parser.parse(miniml=miniml, remove_empty=remove_empty)
        parsed = self._parse_with_related_series(self.parser.parse_mapping(miniml))
        if remove_empty:
            parsed = [self.parser.remove_empty_fields(p) for p in parsed]
        return [MINiMLV1Migrator().migrate(p).package for p in parsed]

    def parse_related_series(
        self,
        miniml: str,
        remove_empty: bool = False,
        strict: bool = True,
    ) -> RelatedSeriesParseResult:
        root_parsed = self.parser.parse_mapping(miniml=miniml)
        related_parsed = []
        attempted_accessions = []
        failed_accessions = []
        errors = []
        seen_gses = set(self.parser.series_accessions(root_parsed))
        pending_gses = deque()

        for gse in self.parser.related_accessions(root_parsed):
            if gse not in seen_gses:
                seen_gses.add(gse)
                pending_gses.append(gse)

        while pending_gses:
            gse = pending_gses.popleft()
            attempted_accessions.append(gse)
            logger.info(
                "Related-series progress accession=%s pending=%s seen=%s",
                gse,
                len(pending_gses),
                len(seen_gses),
            )
            try:
                related_miniml = self.geo_fetcher.fetch_gse_miniml(gse=gse)
                parsed = self.parser.parse_mapping(miniml=related_miniml)
            except Exception as error:
                if strict:
                    raise
                safe_error = SafeErrorEnvelope.from_exception(
                    error,
                    provider="ncbi_geo",
                    stage="related_series",
                    item_id=gse,
                )
                failed_accessions.append(gse)
                errors.append(safe_error)
                logger.warning(
                    "Related-series collection degraded accession=%s "
                    "error_type=%s correlation_id=%s",
                    gse,
                    safe_error.error_type,
                    safe_error.correlation_id,
                )
                continue

            related_parsed.extend(parsed)
            for related_gse in self.parser.related_accessions(parsed):
                if related_gse not in seen_gses:
                    seen_gses.add(related_gse)
                    pending_gses.append(related_gse)

        if remove_empty:
            related_parsed = [
                self.parser.remove_empty_fields(series_package)
                for series_package in related_parsed
            ]

        if errors:
            status = OperationStatusV2(
                execution=ExecutionStatus.DEGRADED,
                completeness=CompletenessStatus.PARTIAL,
                evidence_confidence=EvidenceConfidence.NOT_ASSESSED,
                validation=ValidationStatus.VALID,
                publication=PublicationDisposition.REVIEW_REQUIRED,
                terminal_reason="related_series_partial",
                errors=tuple(errors),
            )
        elif related_parsed:
            status = OperationStatusV2(
                execution=ExecutionStatus.SUCCEEDED,
                completeness=CompletenessStatus.COMPLETE,
                evidence_confidence=EvidenceConfidence.NOT_ASSESSED,
                validation=ValidationStatus.VALID,
                publication=PublicationDisposition.PUBLISHABLE,
                terminal_reason="related_series_complete",
            )
        else:
            status = OperationStatusV2(
                execution=ExecutionStatus.SUCCEEDED,
                completeness=CompletenessStatus.EMPTY,
                evidence_confidence=EvidenceConfidence.NOT_ASSESSED,
                validation=ValidationStatus.VALID,
                publication=PublicationDisposition.NOT_REQUESTED,
                terminal_reason="no_related_series",
            )

        return RelatedSeriesParseResult(
            related_parsed,
            status=status,
            attempted_accessions=attempted_accessions,
            failed_accessions=failed_accessions,
        )


    def _parse_with_related_series(self, parsed: list[dict]) -> list[dict]:
        all_series = list(parsed)
        seen_gses = set(self.parser.series_accessions(parsed))
        pending_gses = deque()

        for gse in self.parser.related_accessions(parsed):
            if gse not in seen_gses:
                seen_gses.add(gse)
                pending_gses.append(gse)

        while pending_gses:
            gse = pending_gses.popleft()
            related_miniml = self.geo_fetcher.fetch_gse_miniml(gse=gse)
            related_parsed = self.parser.parse_mapping(miniml=related_miniml)
            all_series.extend(related_parsed)

            for related_gse in self.parser.related_accessions(related_parsed):
                if related_gse not in seen_gses:
                    seen_gses.add(related_gse)
                    pending_gses.append(related_gse)

        return all_series


