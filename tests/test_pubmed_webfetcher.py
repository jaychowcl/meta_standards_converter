# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
import os
import sys
import unittest
from unittest.mock import Mock


ROOT = os.path.dirname(os.path.dirname(__file__))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from meta_standards_converter.ae_handlers.ae_idf_handlers import IDFConstructor  # noqa: E402
from meta_standards_converter.harmonizers.harmonizers import Harmonizer  # noqa: E402
from meta_standards_converter.pubmed_handlers.pubmed_webfetcher import PubmedWebFetcher  # noqa: E402


class TestPubmedStatusHarmonizer(unittest.TestCase):
    def test_pubstatus2efo_maps_known_status(self):
        self.assertEqual(
            ["published", "EFO", "EFO_0001796"],
            Harmonizer().pubstatus2efo("ppublish"),
        )

    def test_pubstatus2efo_returns_empty_fields_for_blank_status(self):
        self.assertEqual([None, None, None], Harmonizer().pubstatus2efo(None))
        self.assertEqual([None, None, None], Harmonizer().pubstatus2efo(""))

    def test_pubstatus2efo_returns_original_label_for_unknown_status(self):
        self.assertEqual(
            ["collection", None, None],
            Harmonizer().pubstatus2efo("collection"),
        )

    def test_pubstatus2efo_returns_first_unknown_composite_status(self):
        self.assertEqual(
            ["collection", None, None],
            Harmonizer().pubstatus2efo("collection+pubmed"),
        )

    def test_pubstatus2efo_maps_first_known_composite_status(self):
        self.assertEqual(
            ["published", "EFO", "EFO_0001796"],
            Harmonizer().pubstatus2efo("ppublish+pubmed"),
        )


class TestPubmedWebFetcher(unittest.TestCase):
    def test_pubmed_summary_parses_esummary_xml(self):
        xml = b"""<?xml version="1.0"?>
<eSummaryResult>
  <DocSum>
    <Item Name="DOI" Type="String">10.1000/example</Item>
    <Item Name="AuthorList" Type="List">
      <Item Name="Author" Type="String">Ada Lovelace</Item>
      <Item Name="Author" Type="String">Grace Hopper</Item>
    </Item>
    <Item Name="Title" Type="String">Example paper</Item>
    <Item Name="PubStatus" Type="String">ppublish</Item>
  </DocSum>
</eSummaryResult>
"""
        response = Mock(content=xml)
        response.raise_for_status = Mock()

        requester = Mock()
        requester.get.return_value = response

        summary = PubmedWebFetcher(requester=requester).pubmed_summary(pubmed_id="12345")

        requester.get.assert_called_once_with(
            "http://www.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id=12345"
        )
        response.raise_for_status.assert_called_once()
        self.assertEqual(
            (
                "10.1000/example",
                "Ada Lovelace, Grace Hopper",
                "Example paper",
                "published",
                "EFO",
                "EFO_0001796",
            ),
            summary,
        )

    def test_pubmed_summary_uses_original_unknown_pubstatus_label(self):
        xml = b"""<?xml version="1.0"?>
<eSummaryResult>
  <DocSum>
    <Item Name="Title" Type="String">Example paper</Item>
    <Item Name="PubStatus" Type="String">collection</Item>
  </DocSum>
</eSummaryResult>
"""
        response = Mock(content=xml)
        response.raise_for_status = Mock()

        requester = Mock()
        requester.get.return_value = response

        summary = PubmedWebFetcher(requester=requester).pubmed_summary(pubmed_id="12345")

        self.assertEqual(
            (
                None,
                None,
                "Example paper",
                "collection",
                None,
                None,
            ),
            summary,
        )

    def test_idf_constructor_delegates_pubmed_lookup_to_fetcher(self):
        fetcher = Mock()
        fetcher.pubmed_summary.return_value = (
            "doi",
            "authors",
            "title",
            "published",
            "EFO",
            "EFO_0001796",
        )

        result = IDFConstructor(pubmed_fetcher=fetcher)._lookup_pubmed_id("123")

        fetcher.pubmed_summary.assert_called_once_with(pubmed_id="123")
        self.assertEqual(fetcher.pubmed_summary.return_value, result)


if __name__ == "__main__":
    unittest.main()
