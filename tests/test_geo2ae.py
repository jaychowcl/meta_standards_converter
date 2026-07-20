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
from unittest.mock import call, patch


ROOT = os.path.dirname(os.path.dirname(__file__))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from meta_standards_converter.converters.geo2ae import geo2ae  # noqa: E402


class TestGeo2AEConverter(unittest.TestCase):
    @patch("builtins.open")
    @patch("meta_standards_converter.converters.geo2ae.AEConstructor")
    @patch("meta_standards_converter.converters.geo2ae.MINiMLEnricher")
    @patch("meta_standards_converter.converters.geo2ae.GEOParser")
    @patch("meta_standards_converter.converters.geo2ae.GEOWebFetcher")
    def test_convert_uses_parser_related_series_option(
        self,
        fetcher_mock,
        parser_mock,
        enricher_mock,
        constructor_mock,
        open_mock,
    ):
        fetcher_mock.return_value.fetch_gse_miniml.return_value = "<MINiML />"
        primary_json = {"series": {"accession": "GSE1"}}
        related_json = {"series": {"accession": "GSE2"}}
        parser_mock.return_value.parse.return_value = [primary_json, related_json]
        constructor_mock.return_value.miniml2magetab.side_effect = [
            "primary-magetab",
            "related-magetab",
        ]
        enricher_mock.return_value.enrich.side_effect = lambda data: data

        converter = geo2ae()
        result = converter.convert(
            gse="GSE1",
            related_series=True,
            out=None,
        )

        parser_mock.return_value.parse.assert_called_once_with(
            miniml="<MINiML />",
            remove_empty=True,
            related_series=True,
        )
        parser_mock.return_value.parse_related_series.assert_not_called()
        self.assertEqual(["primary-magetab", "related-magetab"], result)
        self.assertEqual(
            [call(data=primary_json), call(data=related_json)],
            enricher_mock.return_value.enrich.call_args_list,
        )
        self.assertEqual(
            [call(data=primary_json), call(data=related_json)],
            constructor_mock.return_value.miniml2magetab.call_args_list,
        )
        constructor_mock.return_value.magetab2file.assert_not_called()
        open_mock.assert_not_called()

    @patch("meta_standards_converter.converters.geo2ae.AEConstructor")
    @patch("meta_standards_converter.converters.geo2ae.MINiMLEnricher")
    @patch("meta_standards_converter.converters.geo2ae.GEOParser")
    @patch("meta_standards_converter.converters.geo2ae.GEOWebFetcher")
    def test_convert_forwards_remove_empty_false(
        self,
        fetcher_mock,
        parser_mock,
        enricher_mock,
        constructor_mock,
    ):
        fetcher_mock.return_value.fetch_gse_miniml.return_value = "<MINiML />"
        meta_json = {"series": {"accession": "GSE1"}}
        enriched_json = {"series": {"accession": "GSE1", "pubmed_publication": []}}
        parser_mock.return_value.parse.return_value = [meta_json]
        enricher_mock.return_value.enrich.return_value = enriched_json
        constructor_mock.return_value.miniml2magetab.return_value = "magetab"

        converter = geo2ae()
        result = converter.convert(
            gse="GSE1",
            related_series=False,
            remove_empty=False,
            out=None,
        )

        parser_mock.return_value.parse.assert_called_once_with(
            miniml="<MINiML />",
            remove_empty=False,
            related_series=False,
        )
        parser_mock.return_value.parse_related_series.assert_not_called()
        enricher_mock.return_value.enrich.assert_called_once_with(data=meta_json)
        constructor_mock.return_value.miniml2magetab.assert_called_once_with(data=enriched_json)
        self.assertEqual(["magetab"], result)

    @patch("meta_standards_converter.converters.geo2ae.AEConstructor")
    @patch("meta_standards_converter.converters.geo2ae.MINiMLEnricher")
    @patch("meta_standards_converter.converters.geo2ae.GEOParser")
    @patch("meta_standards_converter.converters.geo2ae.GEOWebFetcher")
    def test_convert_forwards_forced_platform_handler(
        self,
        fetcher_mock,
        parser_mock,
        enricher_mock,
        constructor_mock,
    ):
        fetcher_mock.return_value.fetch_gse_miniml.return_value = "<MINiML />"
        package = {"series": {"accession": "GSE1"}}
        parser_mock.return_value.parse.return_value = [package]
        enricher_mock.return_value.enrich.return_value = package
        constructor_mock.return_value.miniml2magetab.return_value = "magetab"

        result = geo2ae().convert(gse="GSE1", platform_handler="array")

        self.assertEqual(["magetab"], result)
        constructor_mock.return_value.miniml2magetab.assert_called_once_with(
            data=package,
            platform_handler="array",
        )

    @patch("meta_standards_converter.converters.geo2ae.AEConstructor")
    @patch("meta_standards_converter.converters.geo2ae.MINiMLEnricher")
    @patch("meta_standards_converter.converters.geo2ae.GEOParser")
    @patch("meta_standards_converter.converters.geo2ae.GEOWebFetcher")
    def test_convert_emits_stage_logs_without_payload_dump(
        self,
        fetcher_mock,
        parser_mock,
        enricher_mock,
        constructor_mock,
    ):
        fetcher_mock.return_value.fetch_gse_miniml.return_value = "<MINiML><Series /></MINiML>"
        meta_json = {"series": {"accession": "GSE1"}}
        enriched_json = {"series": {"accession": "GSE1", "secret": "do-not-log"}}
        parser_mock.return_value.parse.return_value = [meta_json]
        enricher_mock.return_value.enrich.return_value = enriched_json
        constructor_mock.return_value.miniml2magetab.return_value = "magetab"

        with self.assertLogs("meta_standards_converter.converters.geo2ae", level="DEBUG") as logs:
            result = geo2ae().convert(gse="GSE1", related_series=True, out=".dev")

        log_output = "\n".join(logs.output)
        self.assertEqual(["magetab"], result)
        self.assertIn("INFO:meta_standards_converter.converters.geo2ae:GSE1: fetching GEO MINiML", log_output)
        self.assertIn("DEBUG:meta_standards_converter.converters.geo2ae:GSE1: fetched GEO MINiML", log_output)
        self.assertIn("INFO:meta_standards_converter.converters.geo2ae:GSE1: parsing GEO MINiML", log_output)
        self.assertIn("INFO:meta_standards_converter.converters.geo2ae:GSE1: enriching parsed package 1", log_output)
        self.assertIn("INFO:meta_standards_converter.converters.geo2ae:GSE1: building MAGE-TAB package 1", log_output)
        self.assertIn("INFO:meta_standards_converter.converters.geo2ae:GSE1: writing MAGE-TAB package 1 to .dev", log_output)
        self.assertNotIn("do-not-log", log_output)


if __name__ == "__main__":
    unittest.main()
