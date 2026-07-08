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

import io
import tarfile

from meta_standards_converter.helpers.request_helper import RateLimitedRequester


class GEOWebFetcher:

    def __init__(self, requester=None, request_settings=None):
        self.requester = requester or RateLimitedRequester(
            service="geo_ftp",
            settings=request_settings,
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

        # use url to fetch miniml file
        response = self.requester.get(url)
        response.raise_for_status()

        # read .tgz content and extract miniml as string
        fileobj = io.BytesIO(response.content)
        with tarfile.open(fileobj=fileobj, mode="r:gz") as tar:
            miniml_file = tar.extractfile(f"{gse}_family.xml")
            miniml = miniml_file.read().decode("utf-8")

        return miniml
