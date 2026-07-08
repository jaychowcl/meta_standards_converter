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
from meta_standards_converter.helpers.request_helper import RateLimitedRequester


class PubmedWebFetcher:
    def __init__(self, requester=None, request_settings=None):
        self.requester = requester or RateLimitedRequester(
            service="ncbi_eutils",
            settings=request_settings,
        )

    def fetch_pubmed_summary(self, pubmed_id: str) -> ET.Element:
        url = f"http://www.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id={pubmed_id}"
        response = self.requester.get(url)
        response.raise_for_status()
        return ET.fromstring(response.content)

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
