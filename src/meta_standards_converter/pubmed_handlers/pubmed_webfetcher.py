# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""
Fetches and parses PubMed summary metadata.
"""

import xml.etree.ElementTree as ET

from meta_standards_converter.harmonizers.harmonizers import Harmonizer
from meta_standards_converter.helpers.request_helper import (
    NCBIApplicationIdentity,
    RateLimitedRequester,
    RequestSettings,
)
from meta_standards_converter.runtime_contracts import get_resource_profile
from meta_standards_converter.xml_safety import parse_xml, read_limited_response


class PubmedWebFetcher:
    def __init__(
        self,
        requester=None,
        request_settings=None,
        ncbi_identity: NCBIApplicationIdentity | None = None,
        resource_profile: str = "standard",
        resource_overrides=None,
    ):
        self.resource_profile = get_resource_profile(
            resource_profile,
            overrides=resource_overrides,
        )
        self.ncbi_identity = ncbi_identity or NCBIApplicationIdentity()
        self.requester = requester or RateLimitedRequester(
            service="ncbi_eutils",
            settings=request_settings
            or RequestSettings.from_resource_profile(
                self.resource_profile,
                request_delay=0.5,
            ),
        )

    def fetch_pubmed_summary(self, pubmed_id: str) -> ET.Element:
        url = "https://www.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
        response = self.requester.get(
            url,
            params={
                "db": "pubmed",
                "id": pubmed_id,
                **self.ncbi_identity.params(),
            },
            stream=True,
        )
        response.raise_for_status()
        content = read_limited_response(
            response,
            max_bytes=self.resource_profile.max_xml_bytes,
        )
        return parse_xml(content, max_bytes=self.resource_profile.max_xml_bytes)

    def pubmed_summary(self, pubmed_id: str) -> tuple:
        root = self.fetch_pubmed_summary(pubmed_id=pubmed_id)

        doi = root.findtext(".//Item[@Name='DOI']")
        author_string = ", ".join(
            author.text
            for author in root.findall(".//Item[@Name='AuthorList']/Item")
            if author.text
        ) or None
        title = root.findtext(".//Item[@Name='Title']")
        status, status_term_source_ref, status_term_accession_number = Harmonizer().pubstatus2efo(
            root.findtext(".//Item[@Name='PubStatus']")
        )

        return (
            doi,
            author_string,
            title,
            status,
            status_term_source_ref,
            status_term_accession_number,
        )
