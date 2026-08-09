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


class GEOWebFetcher:

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

        # Stream the archive to disk, then validate the exact single member.
        with tempfile.NamedTemporaryFile(suffix=".tgz") as archive:
            archive_bytes = stream_limited_response(
                response,
                archive,
                max_bytes=self.resource_profile.max_compressed_archive_bytes,
            )
            archive.flush()
            with tarfile.open(name=archive.name, mode="r:gz") as tar:
                members = tar.getmembers()
                expected_name = f"{gse}_family.xml"
                if (
                    len(members) != 1
                    or not members[0].isfile()
                    or members[0].name != expected_name
                ):
                    raise ValueError(
                        "GEO archive must contain exactly one expected XML member."
                    )
                member = members[0]
                xml_limit = min(
                    self.resource_profile.max_xml_bytes,
                    self.resource_profile.max_expanded_archive_bytes,
                )
                if member.size > xml_limit:
                    raise ValueError(
                        f"GEO XML member exceeds the {xml_limit} byte expanded limit."
                    )
                miniml_file = tar.extractfile(member)
                if miniml_file is None:
                    raise ValueError("GEO archive XML member could not be read.")
                encoded = miniml_file.read(xml_limit + 1)
                if len(encoded) > xml_limit:
                    raise ValueError(
                        f"GEO XML member exceeds the {xml_limit} byte expanded limit."
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
