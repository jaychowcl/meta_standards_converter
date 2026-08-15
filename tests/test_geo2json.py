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

import requests


ROOT = os.path.dirname(os.path.dirname(__file__))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from meta_standards_converter.converters.geo2json import geo2json  # noqa: E402
from meta_standards_converter.miniml import MINiMLPackage  # noqa: E402
from meta_standards_converter.runtime_contracts import get_resource_profile  # noqa: E402


def _package(series: dict) -> MINiMLPackage:
    return MINiMLPackage.from_mapping({
        "miniml_schema_version": "2.0",
        "source": {"format": "test"},
        "series": series,
    })


class TestGeo2JSONConverter(unittest.TestCase):
    @patch("meta_standards_converter.converters.geo2json.MINiMLEnricher")
    @patch("meta_standards_converter.converters.geo2json.GEOParser")
    @patch("meta_standards_converter.converters.geo2json.GEOWebFetcher")
    def test_typed_resource_profile_is_shared_by_default_network_collaborators(
        self, fetcher_mock, parser_mock, enricher_mock
    ):
        profile = get_resource_profile(
            "standard", overrides={"max_xml_bytes": 4096}
        )

        converter = geo2json(resource_profile=profile)

        fetcher_mock.assert_called_once_with(
            resource_profile=profile, resource_overrides=None
        )
        enricher_mock.assert_called_once_with(
            resource_profile=profile, resource_overrides=None
        )
        parser_mock.assert_called_once_with(geo_fetcher=converter.geo_fetcher)

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
        primary_json = _package({"iid": "GSE1"})
        related_json = _package({"iid": "GSE2"})
        enriched_primary = primary_json
        enriched_related = related_json
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
        parsed_json = _package({"iid": "GSE1"})
        parser_mock.return_value.parse.return_value = [parsed_json]

        result = geo2json().convert(gse="GSE1", enrich=False, out=None)

        enricher_mock.return_value.enrich.assert_not_called()
        self.assertEqual([parsed_json], result)

    @patch("meta_standards_converter.converters.geo2json.MINiMLEnricher")
    @patch("meta_standards_converter.converters.geo2json.GEOParser")
    @patch("meta_standards_converter.converters.geo2json.GEOWebFetcher")
    def test_convert_inherits_one_reciprocal_parent_pubmed_id_before_enrichment(
        self,
        fetcher_mock,
        parser_mock,
        enricher_mock,
    ):
        child = _package(
            {
                "iid": "GSE1",
                "relation": [{"type": "SubSeries of", "target": "GSE2"}],
            }
        )
        parent = _package(
            {
                "iid": "GSE2",
                "pubmed_id": ["12345"],
                "relation": [{"type": "SuperSeries of", "target": "GSE1"}],
            }
        )
        fetcher_mock.return_value.fetch_gse_miniml.side_effect = [
            "<MINiML>child</MINiML>",
            "<MINiML>parent</MINiML>",
        ]
        parser_mock.return_value.parse.side_effect = [[child], [parent]]
        enricher_mock.return_value.enrich.side_effect = lambda data: data

        result = geo2json().convert(gse="GSE1", out=None)

        self.assertEqual(1, len(result))
        self.assertEqual(["12345"], result[0]["series"]["pubmed_id"])
        self.assertEqual(
            {
                "relation_type": "SubSeries of",
                "source_series": "GSE2",
                "pubmed_ids": ["12345"],
            },
            result[0]["extensions"]["publication_inheritance"],
        )
        self.assertEqual(
            [call(gse="GSE1"), call(gse="GSE2")],
            fetcher_mock.return_value.fetch_gse_miniml.call_args_list,
        )
        self.assertEqual(
            [
                call(
                    miniml="<MINiML>child</MINiML>",
                    remove_empty=True,
                    related_series=False,
                ),
                call(
                    miniml="<MINiML>parent</MINiML>",
                    remove_empty=True,
                    related_series=False,
                ),
            ],
            parser_mock.return_value.parse.call_args_list,
        )
        enriched_input = enricher_mock.return_value.enrich.call_args.kwargs["data"]
        self.assertEqual(["12345"], enriched_input["series"]["pubmed_id"])

    @patch("meta_standards_converter.converters.geo2json.MINiMLEnricher")
    @patch("meta_standards_converter.converters.geo2json.GEOParser")
    @patch("meta_standards_converter.converters.geo2json.GEOWebFetcher")
    def test_convert_preserves_direct_publication_without_parent_fetch(
        self,
        fetcher_mock,
        parser_mock,
        enricher_mock,
    ):
        child = _package(
            {
                "iid": "GSE1",
                "pubmed_id": ["999"],
                "relation": [{"type": "SubSeries of", "target": "GSE2"}],
            }
        )
        fetcher_mock.return_value.fetch_gse_miniml.return_value = "<MINiML>child</MINiML>"
        parser_mock.return_value.parse.return_value = [child]
        enricher_mock.return_value.enrich.side_effect = lambda data: data

        result = geo2json().convert(gse="GSE1", out=None)

        self.assertEqual(["999"], result[0]["series"]["pubmed_id"])
        self.assertNotIn("extensions", result[0]["series"])
        fetcher_mock.return_value.fetch_gse_miniml.assert_called_once_with(gse="GSE1")

    @patch("meta_standards_converter.converters.geo2json.MINiMLEnricher")
    @patch("meta_standards_converter.converters.geo2json.GEOParser")
    @patch("meta_standards_converter.converters.geo2json.GEOWebFetcher")
    def test_convert_skips_parent_publication_lookup_when_enrichment_is_disabled(
        self,
        fetcher_mock,
        parser_mock,
        enricher_mock,
    ):
        child = _package(
            {
                "iid": "GSE1",
                "relation": [{"type": "SubSeries of", "target": "GSE2"}],
            }
        )
        fetcher_mock.return_value.fetch_gse_miniml.return_value = "<MINiML>child</MINiML>"
        parser_mock.return_value.parse.return_value = [child]

        result = geo2json().convert(gse="GSE1", enrich=False, out=None)

        self.assertEqual([child], result)
        fetcher_mock.return_value.fetch_gse_miniml.assert_called_once_with(gse="GSE1")
        enricher_mock.return_value.enrich.assert_not_called()

    @patch("meta_standards_converter.converters.geo2json.MINiMLEnricher")
    @patch("meta_standards_converter.converters.geo2json.GEOParser")
    @patch("meta_standards_converter.converters.geo2json.GEOWebFetcher")
    def test_convert_leaves_publication_blank_when_parent_evidence_is_ambiguous(
        self,
        fetcher_mock,
        parser_mock,
        enricher_mock,
    ):
        child = _package(
            {
                "iid": "GSE1",
                "relation": [{"type": "SubSeries of", "target": "GSE2"}],
            }
        )
        parent = _package(
            {
                "iid": "GSE2",
                "pubmed_id": ["123", "456"],
                "relation": [{"type": "SuperSeries of", "target": "GSE1"}],
            }
        )
        fetcher_mock.return_value.fetch_gse_miniml.side_effect = ["child", "parent"]
        parser_mock.return_value.parse.side_effect = [[child], [parent]]
        enricher_mock.return_value.enrich.side_effect = lambda data: data

        result = geo2json().convert(gse="GSE1", out=None)

        self.assertNotIn("pubmed_id", result[0]["series"])
        self.assertNotIn("extensions", result[0]["series"])

    @patch("meta_standards_converter.converters.geo2json.MINiMLEnricher")
    @patch("meta_standards_converter.converters.geo2json.GEOParser")
    @patch("meta_standards_converter.converters.geo2json.GEOWebFetcher")
    def test_convert_leaves_publication_blank_for_multiple_direct_parents(
        self,
        fetcher_mock,
        parser_mock,
        enricher_mock,
    ):
        child = _package(
            {
                "iid": "GSE1",
                "relation": [
                    {"type": "SubSeries of", "target": "GSE2"},
                    {"type": "SubSeries of", "target": "GSE3"},
                ],
            }
        )
        fetcher_mock.return_value.fetch_gse_miniml.return_value = "child"
        parser_mock.return_value.parse.return_value = [child]
        enricher_mock.return_value.enrich.side_effect = lambda data: data

        result = geo2json().convert(gse="GSE1", out=None)

        self.assertNotIn("pubmed_id", result[0]["series"])
        fetcher_mock.return_value.fetch_gse_miniml.assert_called_once_with(gse="GSE1")

    @patch("meta_standards_converter.converters.geo2json.MINiMLEnricher")
    @patch("meta_standards_converter.converters.geo2json.GEOParser")
    @patch("meta_standards_converter.converters.geo2json.GEOWebFetcher")
    def test_convert_preserves_direct_publication_details_without_parent_fetch(
        self,
        fetcher_mock,
        parser_mock,
        enricher_mock,
    ):
        child = _package(
            {
                "iid": "GSE1",
                "relation": [{"type": "SubSeries of", "target": "GSE2"}],
                "pubmed_publication": [
                    {"pubmed_id": "999", "title": "Direct publication"}
                ],
            }
        )
        fetcher_mock.return_value.fetch_gse_miniml.return_value = "child"
        parser_mock.return_value.parse.return_value = [child]
        enricher_mock.return_value.enrich.side_effect = lambda data: data

        result = geo2json().convert(gse="GSE1", out=None)

        self.assertEqual(
            "Direct publication",
            result[0]["series"]["pubmed_publication"][0]["title"],
        )
        fetcher_mock.return_value.fetch_gse_miniml.assert_called_once_with(gse="GSE1")

    @patch("meta_standards_converter.converters.geo2json.MINiMLEnricher")
    @patch("meta_standards_converter.converters.geo2json.GEOParser")
    @patch("meta_standards_converter.converters.geo2json.GEOWebFetcher")
    def test_convert_leaves_publication_blank_for_nonreciprocal_or_failed_parent(
        self,
        fetcher_mock,
        parser_mock,
        enricher_mock,
    ):
        child = _package(
            {
                "iid": "GSE1",
                "relation": [{"type": "SubSeries of", "target": "GSE2"}],
            }
        )
        nonreciprocal = _package({"iid": "GSE2", "pubmed_id": ["123"]})
        enricher_mock.return_value.enrich.side_effect = lambda data: data

        with self.subTest("nonreciprocal"):
            fetcher_mock.return_value.fetch_gse_miniml.side_effect = ["child", "parent"]
            parser_mock.return_value.parse.side_effect = [[child], [nonreciprocal]]
            result = geo2json().convert(gse="GSE1", out=None)
            self.assertNotIn("pubmed_id", result[0]["series"])

        fetcher_mock.return_value.fetch_gse_miniml.reset_mock()
        parser_mock.return_value.parse.reset_mock()
        with self.subTest("fetch failure"):
            fetcher_mock.return_value.fetch_gse_miniml.side_effect = [
                "child",
                requests.RequestException("parent unavailable"),
            ]
            parser_mock.return_value.parse.side_effect = [[child]]
            with self.assertLogs(
                "meta_standards_converter.converters.geo2json", level="WARNING"
            ) as logs:
                result = geo2json().convert(gse="GSE1", out=None)
            self.assertNotIn("pubmed_id", result[0]["series"])
            self.assertNotIn("parent unavailable", "\n".join(logs.output))

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
        parsed_json = _package({"iid": "GSE1", "title": "Börsch"})
        parser_mock.return_value.parse.return_value = [parsed_json]
        enricher_mock.return_value.enrich.return_value = parsed_json

        with tempfile.TemporaryDirectory() as tmpdir:
            result = geo2json().convert(gse="GSE1", out=tmpdir)
            with open(os.path.join(tmpdir, "GSE1.json"), encoding="utf-8") as handle:
                written = json.load(handle)

        self.assertEqual([parsed_json], result)
        self.assertEqual([parsed_json.to_mapping()], written)

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
        parsed_json = _package({"iid": "GSE1"})
        enriched_json = _package({"iid": "GSE1", "extensions": {"secret": "do-not-log"}})
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
