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

from meta_standards_converter.cli.json2ae import main  # noqa: E402


PLATFORM_HANDLERS = (
    "plate_single_cell_sequencing",
    "droplet_single_cell_sequencing",
    "tenx_v2_droplet_single_cell_sequencing",
    "tenx_v3_droplet_single_cell_sequencing",
    "single_cell_sequencing",
    "spatial_sequencing",
    "bulk_sequencing",
    "sequencing",
    "array",
    "generic",
)


class TestJSON2AECLI(unittest.TestCase):
    @patch("meta_standards_converter.cli.json2ae.json2ae")
    def test_list_platform_handlers_requires_no_json_or_converter(self, json2ae_mock):
        stdout = StringIO()

        with redirect_stdout(stdout):
            exit_code = main(["--list-platform-handlers"])

        self.assertEqual(0, exit_code)
        self.assertEqual("".join(f"{name}\n" for name in PLATFORM_HANDLERS), stdout.getvalue())
        json2ae_mock.assert_not_called()

    @patch("meta_standards_converter.cli.json2ae.json2ae")
    def test_platform_handler_is_forwarded_to_every_json(self, json2ae_mock):
        json2ae_mock.return_value.convert.return_value = ["magetab"]

        with redirect_stdout(StringIO()):
            exit_code = main(
                ["GSE1.json", "GSE2.json", "--platform-handler", "bulk_sequencing"]
            )

        self.assertEqual(0, exit_code)
        self.assertEqual(
            [
                call(
                    json_path="GSE1.json",
                    enrich=True,
                    out=".",
                    platform_handler="bulk_sequencing",
                ),
                call(
                    json_path="GSE2.json",
                    enrich=True,
                    out=".",
                    platform_handler="bulk_sequencing",
                ),
            ],
            json2ae_mock.return_value.convert.call_args_list,
        )

    def test_invalid_platform_handler_is_rejected(self):
        with redirect_stderr(StringIO()), self.assertRaises(SystemExit) as raised:
            main(["GSE1.json", "--platform-handler", "invalid"])

        self.assertEqual(2, raised.exception.code)

    @patch("meta_standards_converter.cli.json2ae.json2ae")
    def test_one_json_uses_defaults(self, json2ae_mock):
        json2ae_mock.return_value.convert.return_value = ["magetab"]

        with redirect_stdout(StringIO()):
            exit_code = main(["GSE1.json"])

        self.assertEqual(0, exit_code)
        json2ae_mock.return_value.convert.assert_called_once_with(
            json_path="GSE1.json",
            enrich=True,
            out=".",
        )

    @patch("meta_standards_converter.cli.json2ae.json2ae")
    def test_multiple_json_files_are_converted_in_order(self, json2ae_mock):
        json2ae_mock.return_value.convert.return_value = ["magetab"]

        with redirect_stdout(StringIO()):
            exit_code = main(["GSE1.json", "GSE2.json"])

        self.assertEqual(0, exit_code)
        self.assertEqual(
            [
                call(json_path="GSE1.json", enrich=True, out="."),
                call(json_path="GSE2.json", enrich=True, out="."),
            ],
            json2ae_mock.return_value.convert.call_args_list,
        )

    @patch("meta_standards_converter.cli.json2ae.json2ae")
    def test_no_enrich_and_out_are_forwarded(self, json2ae_mock):
        json2ae_mock.return_value.convert.return_value = ["magetab"]

        with redirect_stdout(StringIO()):
            exit_code = main(["GSE1.json", "--no-enrich", "--out", "output"])

        self.assertEqual(0, exit_code)
        json2ae_mock.return_value.convert.assert_called_once_with(
            json_path="GSE1.json",
            enrich=False,
            out="output",
        )

    @patch("meta_standards_converter.cli.json2ae.json2ae")
    def test_failed_json_returns_one_and_continues(self, json2ae_mock):
        json2ae_mock.return_value.convert.side_effect = [RuntimeError("invalid JSON"), ["magetab"]]
        stdout = StringIO()
        stderr = StringIO()

        with redirect_stdout(stdout), redirect_stderr(stderr):
            exit_code = main(["bad.json", "good.json"])

        self.assertEqual(1, exit_code)
        self.assertEqual(
            [
                call(json_path="bad.json", enrich=True, out="."),
                call(json_path="good.json", enrich=True, out="."),
            ],
            json2ae_mock.return_value.convert.call_args_list,
        )
        self.assertIn("ERROR meta_standards_converter.cli.json2ae: bad.json: conversion failed", stdout.getvalue())
        self.assertIn("RuntimeError: invalid JSON", stdout.getvalue())
        self.assertEqual("", stderr.getvalue())

    @patch("meta_standards_converter.cli.json2ae.json2ae")
    def test_verbose_and_quiet_logging(self, json2ae_mock):
        json2ae_mock.return_value.convert.return_value = ["magetab"]

        verbose = StringIO()
        with redirect_stdout(verbose):
            verbose_exit = main(["GSE1.json", "-v"])
        quiet = StringIO()
        with redirect_stdout(quiet):
            quiet_exit = main(["GSE2.json", "--quiet"])

        self.assertEqual(0, verbose_exit)
        self.assertEqual(0, quiet_exit)
        self.assertIn("GSE1.json: conversion started", verbose.getvalue())
        self.assertIn("converted 1 MAGE-TAB output(s) to .", verbose.getvalue())
        self.assertEqual("", quiet.getvalue())

    @patch("meta_standards_converter.cli.json2ae.json2ae")
    def test_log_file_writes_configured_logs(self, json2ae_mock):
        json2ae_mock.return_value.convert.return_value = ["magetab"]

        with tempfile.TemporaryDirectory() as tmpdir:
            log_path = os.path.join(tmpdir, "json2ae.log")
            with redirect_stdout(StringIO()):
                exit_code = main(["GSE1.json", "-v", "--log-file", log_path])
            with open(log_path, encoding="utf-8") as handle:
                log_content = handle.read()

        self.assertEqual(0, exit_code)
        self.assertIn("GSE1.json: conversion started", log_content)


if __name__ == "__main__":
    unittest.main()
