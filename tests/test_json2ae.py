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
from unittest.mock import MagicMock, call


ROOT = os.path.dirname(os.path.dirname(__file__))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from meta_standards_converter.magetab.constructor import AEConstructor  # noqa: E402
from meta_standards_converter.magetab.idf import IDFConstructor  # noqa: E402
from meta_standards_converter.magetab.sdrf.constructor import SDRFConstructor
from meta_standards_converter.converters.json2ae import JSON2AEConverter  # noqa: E402
from meta_standards_converter.sources.json import (  # noqa: E402
    DatasetPackageGroup,
    SourceLoadResult,
)
from meta_standards_converter.miniml.geo_parser import GEOParser  # noqa: E402
from meta_standards_converter.miniml import MINiMLCodec, MINiMLPackage  # noqa: E402


def package(accession="GSE1"):
    return {
        "miniml_schema_version": "2.0",
        "source": {"format": "test"},
        "series": {
            "accession": [{"value": accession, "database": "GEO"}],
            "title": f"Study {accession}",
        },
        "sample": [{"iid": f"GSM-{accession}"}],
        "platform": [],
    }


def typed(payload):
    return MINiMLCodec().decode(payload).package


def atlas_v1(datasets):
    return {
        "schema_version": "1.0",
        "atlas": {"atlas_id": "atlas-test", "title": "Test", "theme": "test"},
        "run": {
            "run_id": "run-test",
            "created_at": "2026-07-31T00:00:00Z",
            "config": {"queries": [], "metadata_repositories": [], "options": {}},
            "status": "partial",
        },
        "datasets": datasets,
        "publications": [],
        "summary": {
            "dataset_count": len(datasets),
            "publication_count": 0,
            "completed_dataset_count": sum(
                item["status"] == "harmonized" for item in datasets
            ),
            "failed_dataset_count": sum(
                item["status"] == "failed" for item in datasets
            ),
        },
    }


def atlas_dataset(dataset_id, status, metadata, diagnostics=None):
    return {
        "dataset_id": dataset_id,
        "source_repository": "geo",
        "source_ordinal": 0,
        "status": status,
        "metadata": metadata,
        "publication_ids": [],
        "review": None,
        "harmonization": None,
        "diagnostics": diagnostics or [],
    }


class TestJSON2AEConverter(unittest.TestCase):
    def setUp(self):
        from unittest.mock import patch
        patcher = patch("meta_standards_converter.converters.json2ae.MAGETabWriter")
        self.writer = patcher.start().return_value
        self.addCleanup(patcher.stop)

    def write_json(self, directory, payload, name="input.json"):
        path = os.path.join(directory, name)
        if isinstance(payload, MINiMLPackage):
            payload = MINiMLCodec().encode(payload)
        elif isinstance(payload, list) and payload and all(
            isinstance(item, MINiMLPackage) for item in payload
        ):
            payload = MINiMLCodec().encode_many(payload)
        with open(path, "w", encoding="utf-8") as handle:
            json.dump(payload, handle)
        return path

    def test_convert_loads_list_enriches_and_builds_each_package_in_order(self):
        enricher = MagicMock()
        constructor = MagicMock()
        first = package("GSE1")
        second = package("GSE2")
        enriched_first = typed({**first, "extensions": {"enriched": True}})
        enriched_second = typed({**second, "extensions": {"enriched": True}})
        enricher.enrich.side_effect = [enriched_first, enriched_second]
        constructor.miniml2magetab.side_effect = ["first-magetab", "second-magetab"]

        with tempfile.TemporaryDirectory() as tmpdir:
            path = self.write_json(tmpdir, [first, second])
            result = JSON2AEConverter(enricher=enricher, ae_constructor=constructor).convert(path)

        self.assertEqual(["first-magetab", "second-magetab"], result)
        self.assertEqual([call(data=typed(first)), call(data=typed(second))], enricher.enrich.call_args_list)
        self.assertEqual(
            [call(data=enriched_first), call(data=enriched_second)],
            constructor.miniml2magetab.call_args_list,
        )
        self.writer.write.assert_not_called()

    def test_convert_accepts_one_package_object_and_can_skip_enrichment(self):
        enricher = MagicMock()
        constructor = MagicMock()
        constructor.miniml2magetab.return_value = "magetab"
        payload = package()

        with tempfile.TemporaryDirectory() as tmpdir:
            path = self.write_json(tmpdir, payload)
            result = JSON2AEConverter(enricher=enricher, ae_constructor=constructor).convert(
                path,
                enrich=False,
            )

        self.assertEqual(["magetab"], result)
        enricher.enrich.assert_not_called()
        constructor.miniml2magetab.assert_called_once_with(data=typed(payload))

    def test_convert_forwards_forced_platform_handler(self):
        constructor = MagicMock()
        constructor.miniml2magetab.return_value = "magetab"
        payload = package()

        with tempfile.TemporaryDirectory() as tmpdir:
            path = self.write_json(tmpdir, payload)
            result = JSON2AEConverter(ae_constructor=constructor).convert(
                path,
                enrich=False,
                platform_handler="bulk_sequencing",
            )

        self.assertEqual(["magetab"], result)
        constructor.miniml2magetab.assert_called_once_with(
            data=typed(payload),
            platform_handler="bulk_sequencing",
        )

    def test_convert_writes_each_magetab_when_out_is_supplied(self):
        constructor = MagicMock()
        constructor.miniml2magetab.side_effect = ["first", "second"]

        with tempfile.TemporaryDirectory() as tmpdir:
            path = self.write_json(tmpdir, [package("GSE1"), package("GSE2")])
            result = JSON2AEConverter(
                enricher=MagicMock(enrich=lambda data: data),
                ae_constructor=constructor,
            ).convert(path, out=tmpdir)

        self.assertEqual(["first", "second"], result)
        self.assertEqual(
            [call(magetab="first", out=tmpdir), call(magetab="second", out=tmpdir)],
            self.writer.write.call_args_list,
        )

    def test_convert_rejects_missing_file(self):
        with self.assertRaisesRegex(FileNotFoundError, "MINiML JSON file not found"):
            JSON2AEConverter().convert("missing.json")

    def test_convert_rejects_empty_package_list(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = self.write_json(tmpdir, [])
            with self.assertRaisesRegex(
                ValueError, "JSON source contains no convertible package groups"
            ):
                JSON2AEConverter().convert(path)

    def test_convert_rejects_non_object_package_before_enrichment(self):
        enricher = MagicMock()
        with tempfile.TemporaryDirectory() as tmpdir:
            path = self.write_json(tmpdir, [package(), "invalid"])
            with self.assertRaisesRegex(ValueError, "package 2 must be a JSON object"):
                JSON2AEConverter(enricher=enricher).convert(path)
        enricher.enrich.assert_not_called()

    def test_convert_rejects_package_without_geo_series_accession(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = self.write_json(
                tmpdir,
                {
                    "miniml_schema_version": "2.0",
                    "source": {"format": "test"},
                    "series": {"title": "Missing accession"},
                    "sample": [{"iid": "GSM1"}],
                },
            )
            with self.assertRaisesRegex(ValueError, "series requires iid or accession"):
                JSON2AEConverter().convert(path)

    def test_convert_accepts_non_geo_study_accession(self):
        constructor = MagicMock()
        constructor.miniml2magetab.return_value = "magetab"
        payload = package("E-MTAB-1")

        with tempfile.TemporaryDirectory() as tmpdir:
            path = self.write_json(tmpdir, payload)
            result = JSON2AEConverter(ae_constructor=constructor).convert(path, enrich=False)

        self.assertEqual(["magetab"], result)

    def test_convert_accepts_iid_only_study_and_uses_it_as_constructor_identity(self):
        constructor = MagicMock()
        constructor.miniml2magetab.return_value = "magetab"
        payload = package("E-MTAB-unused")
        payload["series"] = {"iid": "E-MTAB-ONLY", "title": "IID-only study"}

        with tempfile.TemporaryDirectory() as tmpdir:
            path = self.write_json(tmpdir, payload)
            result = JSON2AEConverter(ae_constructor=constructor).convert(path, enrich=False)

        decoded = typed(payload)
        self.assertEqual(["magetab"], result)
        constructor.miniml2magetab.assert_called_once_with(data=decoded)
        self.assertEqual(
            "E-MTAB-ONLY",
            AEConstructor()._series_accession(decoded.to_mapping()),
        )

    def test_convert_accepts_harmonized_v2_datasets_and_warns_for_skipped_states(self):
        constructor = MagicMock()
        constructor.miniml2magetab.side_effect = ["first", "second"]
        first = package("GSE1")
        second = package("E-MTAB-2")
        payload = atlas_v1(
            [
                atlas_dataset("GSE1", "harmonized", first),
                atlas_dataset("E-MTAB-2", "harmonized", second),
                atlas_dataset(
                    "GSE3",
                    "failed",
                    {},
                    diagnostics=[
                        {
                            "code": "harmonization_failed",
                            "message": "lookup failed",
                            "severity": "error",
                            "path": None,
                        }
                    ],
                ),
            ]
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            path = self.write_json(tmpdir, payload)
            with self.assertLogs(
                "meta_standards_converter.converters.json2ae",
                level="WARNING",
            ) as logs:
                result = JSON2AEConverter(ae_constructor=constructor).convert(
                    path,
                    enrich=False,
                )

        self.assertEqual(["first", "second"], result)
        self.assertEqual(
            [call(data=typed(first)), call(data=typed(second))],
            constructor.miniml2magetab.call_args_list,
        )
        self.assertIn(
            "GSE3: dataset status failed; no convertible metadata; "
            "harmonization_failed: lookup failed",
            "\n".join(logs.output),
        )

    def test_convert_rejects_atlas_without_harmonized_metadata(self):
        payload = atlas_v1([atlas_dataset("GSE1", "collected", {})])

        with tempfile.TemporaryDirectory() as tmpdir:
            path = self.write_json(tmpdir, payload)
            with self.assertRaisesRegex(
                ValueError,
                "no convertible package groups",
            ):
                JSON2AEConverter().convert(path)

    def test_convert_uses_injected_package_source(self):
        constructor = MagicMock()
        constructor.miniml2magetab.return_value = "magetab"
        source = MagicMock()
        payload = package("GSE1")
        source.load.return_value = SourceLoadResult(
            groups=(
                DatasetPackageGroup(
                    dataset_id="GSE1",
                    packages=(typed(payload),),
                    source_accession="GSE1",
                ),
            ),
        )

        result = JSON2AEConverter(
            ae_constructor=constructor,
            package_source=source,
        ).convert("virtual-atlas.json", enrich=False)

        self.assertEqual(["magetab"], result)
        source.load.assert_called_once_with("virtual-atlas.json")
        constructor.miniml2magetab.assert_called_once_with(data=typed(payload))

    def test_convert_still_rejects_malformed_geo_accession(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = self.write_json(tmpdir, package("GSE-not-numeric"))
            with self.assertRaisesRegex(ValueError, "no usable study accession"):
                JSON2AEConverter().convert(path)

    def test_convert_rejects_malformed_geo_iid(self):
        payload = package("E-MTAB-unused")
        payload["series"] = {"iid": "GSE-not-numeric", "title": "Invalid"}

        with tempfile.TemporaryDirectory() as tmpdir:
            path = self.write_json(tmpdir, payload)
            with self.assertRaisesRegex(ValueError, "no usable study accession"):
                JSON2AEConverter().convert(path)

    def test_convert_logs_stages_without_metadata_payload(self):
        constructor = MagicMock()
        constructor.miniml2magetab.return_value = "magetab"
        payload = package()
        payload["extensions"] = {"secret": "do-not-log"}

        with tempfile.TemporaryDirectory() as tmpdir:
            path = self.write_json(tmpdir, payload)
            with self.assertLogs("meta_standards_converter.converters.json2ae", level="DEBUG") as logs:
                result = JSON2AEConverter(
                    enricher=MagicMock(enrich=lambda data: data),
                    ae_constructor=constructor,
                ).convert(path)

        self.assertEqual(["magetab"], result)
        output = "\n".join(logs.output)
        self.assertIn("loaded 1 parsed package(s)", output)
        self.assertIn("enriching parsed package 1", output)
        self.assertIn("building MAGE-TAB package 1", output)
        self.assertNotIn("do-not-log", output)

    def test_fixture_parsed_package_matches_independent_magetab_expectations(self):
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
            actual = JSON2AEConverter(ae_constructor=converter_constructor).convert(path, enrich=False)

        self.assertEqual(1, len(actual))
        rows = {row[0]: row[1:] for row in actual[0] if row}
        self.assertEqual(["1.1"], rows["MAGE-TAB Version"])
        self.assertEqual(
            [
                "A CSF Disease-Associated Macrophage Signature defines "
                "Progressive Multiple Sclerosis"
            ],
            rows["Investigation Title"],
        )
        self.assertEqual(["2026-05-17"], rows["Public Release Date"])
        self.assertEqual(["GSE328265"], rows["Comment[SecondaryAccession]"])


if __name__ == "__main__":
    unittest.main()
