# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
import json
import os
import sys
import tempfile
import unittest
from unittest.mock import call, patch


ROOT = os.path.dirname(os.path.dirname(__file__))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from meta_standards_converter.converters.geo2json import geo2json  # noqa: E402


class TestGeo2JSONConverter(unittest.TestCase):
    @patch("meta_standards_converter.converters.geo2json.MINiMLEnricher")
    @patch("meta_standards_converter.converters.geo2json.GEOParser")
    @patch("meta_standards_converter.converters.geo2json.GEOWebFetcher")
    def test_convert_fetches_parses_and_enriches_by_default(
        self,
        fetcher_mock,
        parser_mock,
        enricher_mock,
    ):
        fetcher_mock.return_value.fetch_gse_miniml.return_value = "<MINiML />"
        primary_json = {"series": {"accession": "GSE1"}}
        related_json = {"series": {"accession": "GSE2"}}
        enriched_primary = {"series": {"accession": "GSE1", "pubmed_publication": []}}
        enriched_related = {"series": {"accession": "GSE2", "pubmed_publication": []}}
        parser_mock.return_value.parse.return_value = [primary_json, related_json]
        enricher_mock.return_value.enrich.side_effect = [enriched_primary, enriched_related]

        result = geo2json().convert(
            gse="GSE1",
            related_series=True,
            remove_empty=False,
            out=None,
        )

        fetcher_mock.return_value.fetch_gse_miniml.assert_called_once_with(gse="GSE1")
        parser_mock.return_value.parse.assert_called_once_with(
            miniml="<MINiML />",
            remove_empty=False,
            related_series=True,
        )
        self.assertEqual(
            [call(data=primary_json), call(data=related_json)],
            enricher_mock.return_value.enrich.call_args_list,
        )
        self.assertEqual([enriched_primary, enriched_related], result)

    @patch("meta_standards_converter.converters.geo2json.MINiMLEnricher")
    @patch("meta_standards_converter.converters.geo2json.GEOParser")
    @patch("meta_standards_converter.converters.geo2json.GEOWebFetcher")
    def test_convert_can_skip_enrichment(
        self,
        fetcher_mock,
        parser_mock,
        enricher_mock,
    ):
        fetcher_mock.return_value.fetch_gse_miniml.return_value = "<MINiML />"
        parsed_json = {"series": {"accession": "GSE1"}}
        parser_mock.return_value.parse.return_value = [parsed_json]

        result = geo2json().convert(gse="GSE1", enrich=False, out=None)

        enricher_mock.return_value.enrich.assert_not_called()
        self.assertEqual([parsed_json], result)

    @patch("meta_standards_converter.converters.geo2json.MINiMLEnricher")
    @patch("meta_standards_converter.converters.geo2json.GEOParser")
    @patch("meta_standards_converter.converters.geo2json.GEOWebFetcher")
    def test_convert_writes_json_list_file(
        self,
        fetcher_mock,
        parser_mock,
        enricher_mock,
    ):
        fetcher_mock.return_value.fetch_gse_miniml.return_value = "<MINiML />"
        parsed_json = {"series": {"accession": "GSE1", "title": "Börsch"}}
        parser_mock.return_value.parse.return_value = [parsed_json]
        enricher_mock.return_value.enrich.return_value = parsed_json

        with tempfile.TemporaryDirectory() as tmpdir:
            result = geo2json().convert(gse="GSE1", out=tmpdir)
            with open(os.path.join(tmpdir, "GSE1.json"), encoding="utf-8") as handle:
                written = json.load(handle)

        self.assertEqual([parsed_json], result)
        self.assertEqual([parsed_json], written)

    @patch("meta_standards_converter.converters.geo2json.MINiMLEnricher")
    @patch("meta_standards_converter.converters.geo2json.GEOParser")
    @patch("meta_standards_converter.converters.geo2json.GEOWebFetcher")
    def test_convert_emits_stage_logs_without_payload_dump(
        self,
        fetcher_mock,
        parser_mock,
        enricher_mock,
    ):
        fetcher_mock.return_value.fetch_gse_miniml.return_value = "<MINiML><Series /></MINiML>"
        parsed_json = {"series": {"accession": "GSE1"}}
        enriched_json = {"series": {"accession": "GSE1", "secret": "do-not-log"}}
        parser_mock.return_value.parse.return_value = [parsed_json]
        enricher_mock.return_value.enrich.return_value = enriched_json

        with self.assertLogs("meta_standards_converter.converters.geo2json", level="DEBUG") as logs:
            result = geo2json().convert(gse="GSE1", related_series=True, out=None)

        log_output = "\n".join(logs.output)
        self.assertEqual([enriched_json], result)
        self.assertIn("INFO:meta_standards_converter.converters.geo2json:GSE1: fetching GEO MINiML", log_output)
        self.assertIn("DEBUG:meta_standards_converter.converters.geo2json:GSE1: fetched GEO MINiML", log_output)
        self.assertIn("INFO:meta_standards_converter.converters.geo2json:GSE1: parsing GEO MINiML", log_output)
        self.assertIn("INFO:meta_standards_converter.converters.geo2json:GSE1: enriching parsed package 1", log_output)
        self.assertNotIn("do-not-log", log_output)


if __name__ == "__main__":
    unittest.main()
