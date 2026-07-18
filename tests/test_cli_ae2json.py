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
from contextlib import redirect_stderr, redirect_stdout
from io import StringIO
from unittest.mock import call, patch


ROOT = os.path.dirname(os.path.dirname(__file__))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from meta_standards_converter.cli.ae2json import main  # noqa: E402


class TestAE2JSONCLI(unittest.TestCase):
    @patch("meta_standards_converter.cli.ae2json.ae2json")
    def test_one_source_uses_defaults(self, converter_mock):
        converter_mock.return_value.convert.return_value = [{"series": {}}]

        with redirect_stdout(StringIO()):
            exit_code = main(["E-MTAB-1"])

        self.assertEqual(0, exit_code)
        converter_mock.return_value.convert.assert_called_once_with(
            source="E-MTAB-1", out=".", sdrf_sources=None
        )

    @patch("meta_standards_converter.cli.ae2json.ae2json")
    def test_multiple_sources_are_converted_in_order(self, converter_mock):
        converter_mock.return_value.convert.return_value = [{"series": {}}]

        with redirect_stdout(StringIO()):
            exit_code = main(["one.idf.txt", "https://example/two.idf.txt"])

        self.assertEqual(0, exit_code)
        self.assertEqual(
            [
                call(source="one.idf.txt", out=".", sdrf_sources=None),
                call(source="https://example/two.idf.txt", out=".", sdrf_sources=None),
            ],
            converter_mock.return_value.convert.call_args_list,
        )

    @patch("meta_standards_converter.cli.ae2json.ae2json")
    def test_sdrf_overrides_and_out_are_forwarded(self, converter_mock):
        converter_mock.return_value.convert.return_value = [{"series": {}}]

        with redirect_stdout(StringIO()):
            exit_code = main([
                "study.idf.txt", "--sdrf", "one.sdrf.txt", "--sdrf", "two.sdrf.txt", "--out", "output"
            ])

        self.assertEqual(0, exit_code)
        converter_mock.return_value.convert.assert_called_once_with(
            source="study.idf.txt",
            out="output",
            sdrf_sources=["one.sdrf.txt", "two.sdrf.txt"],
        )

    def test_sdrf_override_rejects_multiple_sources(self):
        with redirect_stderr(StringIO()):
            with self.assertRaises(SystemExit):
                main(["one.idf.txt", "two.idf.txt", "--sdrf", "override.sdrf.txt"])

    @patch("meta_standards_converter.cli.ae2json.ae2json")
    def test_failure_returns_one_and_continues(self, converter_mock):
        converter_mock.return_value.convert.side_effect = [RuntimeError("bad MAGE-TAB"), [{"series": {}}]]
        stdout = StringIO()

        with redirect_stdout(stdout):
            exit_code = main(["bad.idf.txt", "good.idf.txt"])

        self.assertEqual(1, exit_code)
        self.assertEqual(2, converter_mock.return_value.convert.call_count)
        self.assertIn("bad.idf.txt: conversion failed", stdout.getvalue())
        self.assertIn("RuntimeError: bad MAGE-TAB", stdout.getvalue())

    @patch("meta_standards_converter.cli.ae2json.ae2json")
    def test_verbose_logs_success(self, converter_mock):
        converter_mock.return_value.convert.return_value = [{"series": {}}]
        stdout = StringIO()

        with redirect_stdout(stdout):
            exit_code = main(["E-MTAB-1", "-v"])

        self.assertEqual(0, exit_code)
        self.assertIn("E-MTAB-1: conversion started", stdout.getvalue())
        self.assertIn("converted 1 package(s) to .", stdout.getvalue())


if __name__ == "__main__":
    unittest.main()
