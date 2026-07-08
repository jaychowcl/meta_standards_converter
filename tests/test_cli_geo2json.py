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
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from io import StringIO
from unittest.mock import call, patch


ROOT = os.path.dirname(os.path.dirname(__file__))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from meta_standards_converter.cli.geo2json import main  # noqa: E402


class TestGeo2JSONCLI(unittest.TestCase):
    @patch("meta_standards_converter.cli.geo2json.geo2json")
    def test_one_accession_uses_defaults(self, geo2json_mock):
        converter = geo2json_mock.return_value
        converter.convert.return_value = [{"series": {"accession": "GSE1"}}]

        stdout = StringIO()
        with redirect_stdout(stdout):
            exit_code = main(["GSE234602"])

        self.assertEqual(0, exit_code)
        converter.convert.assert_called_once_with(
            gse="GSE234602",
            related_series=False,
            remove_empty=True,
            enrich=True,
            out=".",
        )
        self.assertEqual("", stdout.getvalue())

    @patch("meta_standards_converter.cli.geo2json.geo2json")
    def test_multiple_accessions_are_converted_in_order(self, geo2json_mock):
        converter = geo2json_mock.return_value
        converter.convert.return_value = [{"series": {"accession": "GSE1"}}]

        with redirect_stdout(StringIO()):
            exit_code = main(["GSE234602", "GSE34779"])

        self.assertEqual(0, exit_code)
        self.assertEqual(
            [
                call(gse="GSE234602", related_series=False, remove_empty=True, enrich=True, out="."),
                call(gse="GSE34779", related_series=False, remove_empty=True, enrich=True, out="."),
            ],
            converter.convert.call_args_list,
        )

    @patch("meta_standards_converter.cli.geo2json.geo2json")
    def test_flags_are_passed_to_converter(self, geo2json_mock):
        converter = geo2json_mock.return_value
        converter.convert.return_value = [{"series": {"accession": "GSE1"}}]

        with redirect_stdout(StringIO()):
            exit_code = main(["GSE234602", "--related", "--keep-empty", "--no-enrich", "--out", ".dev"])

        self.assertEqual(0, exit_code)
        converter.convert.assert_called_once_with(
            gse="GSE234602",
            related_series=True,
            remove_empty=False,
            enrich=False,
            out=".dev",
        )

    @patch("meta_standards_converter.cli.geo2json.geo2json")
    def test_related_aliases_are_supported(self, geo2json_mock):
        converter = geo2json_mock.return_value
        converter.convert.return_value = [{"series": {"accession": "GSE1"}}]

        with redirect_stdout(StringIO()):
            first_exit = main(["GSE234602", "--related-series"])
        with redirect_stdout(StringIO()):
            second_exit = main(["GSE34779", "--get-related-series"])

        self.assertEqual(0, first_exit)
        self.assertEqual(0, second_exit)
        self.assertEqual(
            [
                call(gse="GSE234602", related_series=True, remove_empty=True, enrich=True, out="."),
                call(gse="GSE34779", related_series=True, remove_empty=True, enrich=True, out="."),
            ],
            converter.convert.call_args_list,
        )

    @patch("meta_standards_converter.cli.geo2json.geo2json")
    def test_failed_accession_returns_one_and_continues(self, geo2json_mock):
        converter = geo2json_mock.return_value
        converter.convert.side_effect = [
            RuntimeError("network unavailable"),
            [{"series": {"accession": "GSE2"}}],
        ]

        stdout = StringIO()
        stderr = StringIO()
        with redirect_stdout(stdout), redirect_stderr(stderr):
            exit_code = main(["GSE1", "GSE2"])

        self.assertEqual(1, exit_code)
        self.assertEqual(
            [
                call(gse="GSE1", related_series=False, remove_empty=True, enrich=True, out="."),
                call(gse="GSE2", related_series=False, remove_empty=True, enrich=True, out="."),
            ],
            converter.convert.call_args_list,
        )
        self.assertIn(
            "ERROR meta_standards_converter.cli.geo2json: GSE1: JSON conversion failed",
            stdout.getvalue(),
        )
        self.assertIn("Traceback (most recent call last):", stdout.getvalue())
        self.assertIn("RuntimeError: network unavailable", stdout.getvalue())
        self.assertEqual("", stderr.getvalue())

    @patch("meta_standards_converter.cli.geo2json.geo2json")
    def test_verbose_emits_success_logs_to_stdout(self, geo2json_mock):
        converter = geo2json_mock.return_value
        converter.convert.return_value = [{"series": {"accession": "GSE1"}}]

        stdout = StringIO()
        with redirect_stdout(stdout):
            exit_code = main(["GSE234602", "-v"])

        self.assertEqual(0, exit_code)
        self.assertIn("INFO meta_standards_converter.cli.geo2json: GSE234602: JSON conversion started", stdout.getvalue())
        self.assertIn(
            "INFO meta_standards_converter.cli.geo2json: GSE234602: converted 1 JSON package(s) to .",
            stdout.getvalue(),
        )

    @patch("meta_standards_converter.cli.geo2json.geo2json")
    def test_double_verbose_emits_debug_logs_to_stdout(self, geo2json_mock):
        converter = geo2json_mock.return_value
        converter.convert.return_value = [{"series": {"accession": "GSE1"}}]

        stdout = StringIO()
        with redirect_stdout(stdout):
            exit_code = main(["GSE234602", "-vv"])

        self.assertEqual(0, exit_code)
        self.assertIn("DEBUG meta_standards_converter.cli.geo2json: Starting geo2json CLI", stdout.getvalue())

    @patch("meta_standards_converter.cli.geo2json.geo2json")
    def test_quiet_emits_only_errors(self, geo2json_mock):
        converter = geo2json_mock.return_value
        converter.convert.return_value = [{"series": {"accession": "GSE1"}}]

        stdout = StringIO()
        with redirect_stdout(stdout):
            exit_code = main(["GSE234602", "--quiet"])

        self.assertEqual(0, exit_code)
        self.assertEqual("", stdout.getvalue())

    @patch("meta_standards_converter.cli.geo2json.geo2json")
    def test_log_file_writes_configured_logs(self, geo2json_mock):
        converter = geo2json_mock.return_value
        converter.convert.return_value = [{"series": {"accession": "GSE1"}}]

        with tempfile.TemporaryDirectory() as tmpdir:
            log_path = os.path.join(tmpdir, "geo2json.log")
            with redirect_stdout(StringIO()):
                exit_code = main(["GSE234602", "-v", "--log-file", log_path])

            with open(log_path, encoding="utf-8") as handle:
                log_content = handle.read()

        self.assertEqual(0, exit_code)
        self.assertIn("INFO meta_standards_converter.cli.geo2json: GSE234602: JSON conversion started", log_content)
        self.assertIn("GSE234602: converted 1 JSON package(s) to .", log_content)


if __name__ == "__main__":
    unittest.main()
