# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
import copy
import json
import os
import sys
import tempfile
import unittest
from unittest.mock import MagicMock, call


ROOT = os.path.dirname(os.path.dirname(__file__))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from meta_standards_converter.ae_handlers.ae_constructor import AEConstructor  # noqa: E402
from meta_standards_converter.ae_handlers.ae_idf_handlers import IDFConstructor  # noqa: E402
from meta_standards_converter.ae_handlers.ae_sdrf_handlers import SDRFConstructor  # noqa: E402
from meta_standards_converter.converters.json2ae import json2ae  # noqa: E402
from meta_standards_converter.geo_handlers.geo_parser import GEOParser  # noqa: E402


def package(accession="GSE1"):
    return {
        "series": {
            "accession": [{"value": accession, "database": "GEO"}],
            "title": f"Study {accession}",
        },
        "sample": [],
        "platform": [],
    }


class TestJSON2AEConverter(unittest.TestCase):
    def write_json(self, directory, payload, name="input.json"):
        path = os.path.join(directory, name)
        with open(path, "w", encoding="utf-8") as handle:
            json.dump(payload, handle)
        return path

    def test_convert_loads_list_enriches_and_builds_each_package_in_order(self):
        enricher = MagicMock()
        constructor = MagicMock()
        first = package("GSE1")
        second = package("GSE2")
        enriched_first = {**first, "enriched": True}
        enriched_second = {**second, "enriched": True}
        enricher.enrich.side_effect = [enriched_first, enriched_second]
        constructor.miniml2magetab.side_effect = ["first-magetab", "second-magetab"]

        with tempfile.TemporaryDirectory() as tmpdir:
            path = self.write_json(tmpdir, [first, second])
            result = json2ae(enricher=enricher, ae_constructor=constructor).convert(path)

        self.assertEqual(["first-magetab", "second-magetab"], result)
        self.assertEqual([call(data=first), call(data=second)], enricher.enrich.call_args_list)
        self.assertEqual(
            [call(data=enriched_first), call(data=enriched_second)],
            constructor.miniml2magetab.call_args_list,
        )
        constructor.magetab2file.assert_not_called()

    def test_convert_accepts_one_package_object_and_can_skip_enrichment(self):
        enricher = MagicMock()
        constructor = MagicMock()
        constructor.miniml2magetab.return_value = "magetab"
        payload = package()

        with tempfile.TemporaryDirectory() as tmpdir:
            path = self.write_json(tmpdir, payload)
            result = json2ae(enricher=enricher, ae_constructor=constructor).convert(
                path,
                enrich=False,
            )

        self.assertEqual(["magetab"], result)
        enricher.enrich.assert_not_called()
        constructor.miniml2magetab.assert_called_once_with(data=payload)

    def test_convert_writes_each_magetab_when_out_is_supplied(self):
        constructor = MagicMock()
        constructor.miniml2magetab.side_effect = ["first", "second"]

        with tempfile.TemporaryDirectory() as tmpdir:
            path = self.write_json(tmpdir, [package("GSE1"), package("GSE2")])
            result = json2ae(
                enricher=MagicMock(enrich=lambda data: data),
                ae_constructor=constructor,
            ).convert(path, out=tmpdir)

        self.assertEqual(["first", "second"], result)
        self.assertEqual(
            [call(magetab="first", out=tmpdir), call(magetab="second", out=tmpdir)],
            constructor.magetab2file.call_args_list,
        )

    def test_convert_rejects_missing_file(self):
        with self.assertRaisesRegex(FileNotFoundError, "MINiML JSON file not found"):
            json2ae().convert("missing.json")

    def test_convert_rejects_empty_package_list(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = self.write_json(tmpdir, [])
            with self.assertRaisesRegex(ValueError, "non-empty package object or list"):
                json2ae().convert(path)

    def test_convert_rejects_non_object_package_before_enrichment(self):
        enricher = MagicMock()
        with tempfile.TemporaryDirectory() as tmpdir:
            path = self.write_json(tmpdir, [package(), "invalid"])
            with self.assertRaisesRegex(ValueError, "package 2 must be a JSON object"):
                json2ae(enricher=enricher).convert(path)
        enricher.enrich.assert_not_called()

    def test_convert_rejects_package_without_geo_series_accession(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = self.write_json(tmpdir, {"series": {"title": "Missing accession"}})
            with self.assertRaisesRegex(ValueError, "package 1 has no GEO Series accession"):
                json2ae().convert(path)

    def test_convert_logs_stages_without_metadata_payload(self):
        constructor = MagicMock()
        constructor.miniml2magetab.return_value = "magetab"
        payload = package()
        payload["secret"] = "do-not-log"

        with tempfile.TemporaryDirectory() as tmpdir:
            path = self.write_json(tmpdir, payload)
            with self.assertLogs("meta_standards_converter.converters.json2ae", level="DEBUG") as logs:
                result = json2ae(
                    enricher=MagicMock(enrich=lambda data: data),
                    ae_constructor=constructor,
                ).convert(path)

        self.assertEqual(["magetab"], result)
        output = "\n".join(logs.output)
        self.assertIn("loaded 1 parsed package(s)", output)
        self.assertIn("enriching parsed package 1", output)
        self.assertIn("building MAGE-TAB package 1", output)
        self.assertNotIn("do-not-log", output)

    def test_fixture_parsed_package_matches_direct_ae_construction(self):
        fixture_path = os.path.join(ROOT, "tests", "GSE328265_family.xml")
        with open(fixture_path, encoding="utf-8") as handle:
            packages = GEOParser().parse(handle.read())
        insdc_fetcher = MagicMock()
        insdc_fetcher.fetch_sra_runs.return_value = []
        pubmed_fetcher = MagicMock()
        pubmed_fetcher.pubmed_summary.return_value = (None, None, None, None, None, None)

        with tempfile.TemporaryDirectory() as tmpdir:
            path = self.write_json(tmpdir, packages)
            converter_constructor = AEConstructor(
                idf_constructor=IDFConstructor(pubmed_fetcher=pubmed_fetcher),
                sdrf_constructor=SDRFConstructor(insdc_fetcher=insdc_fetcher)
            )
            actual = json2ae(ae_constructor=converter_constructor).convert(path, enrich=False)

        direct_constructor = AEConstructor(
            idf_constructor=IDFConstructor(pubmed_fetcher=pubmed_fetcher),
            sdrf_constructor=SDRFConstructor(insdc_fetcher=insdc_fetcher)
        )
        expected = [
            direct_constructor.miniml2magetab(data=copy.deepcopy(item))
            for item in packages
        ]
        self.assertEqual(expected, actual)


if __name__ == "__main__":
    unittest.main()
