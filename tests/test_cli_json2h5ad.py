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

from meta_standards_converter.cli.json2h5ad import main  # noqa: E402
from meta_standards_converter.converters.json2h5ad import ConversionResult  # noqa: E402


class TestJSON2H5ADCLI(unittest.TestCase):
    @patch("meta_standards_converter.cli.json2h5ad.json2h5ad")
    def test_workflow_and_asset_options_are_forwarded(self, json2h5ad_mock):
        converter = json2h5ad_mock.return_value
        converter.convert.return_value = "GSE1.h5ad"

        with redirect_stdout(StringIO()):
            exit_code = main([
                "GSE1.json",
                "--asset-manifest", "assets.csv",
                "--asset", "GSM1=local.h5ad",
                "--force-reprocess",
                "--pipeline", "scrnaseq",
                "--genome", "GRCh38",
                "--profile", "apptainer",
                "--revision", "4.1.0",
                "--params-file", "params.json",
                "--nextflow-config", "nextflow.config",
                "--work-dir", "work",
                "--resume",
                "--overwrite",
                "--matrix-orientation", "genes-by-observations",
            ])

        self.assertEqual(0, exit_code)
        converter.convert.assert_called_once_with(
            json_path="GSE1.json",
            out=".",
            asset_manifest="assets.csv",
            asset_specs=["GSM1=local.h5ad"],
            force_reprocess=True,
            pipeline="scrnaseq",
            genome="GRCh38",
            profile="apptainer",
            revision="4.1.0",
            params_file="params.json",
            nextflow_config="nextflow.config",
            work_dir="work",
            resume=True,
            overwrite=True,
            matrix_orientation="genes-by-observations",
        )

    @patch("meta_standards_converter.cli.json2h5ad.json2h5ad")
    def test_partial_conversion_returns_one(self, json2h5ad_mock):
        converter = json2h5ad_mock.return_value
        converter.convert.return_value = ConversionResult(
            study_accession="GSE1",
            sample_h5ads={"GSM1": "GSM1.h5ad"},
            failures=["combined output incompatible"],
        )

        with redirect_stdout(StringIO()):
            exit_code = main(["GSE1.json"])

        self.assertEqual(1, exit_code)

    @patch("meta_standards_converter.cli.json2h5ad.json2h5ad")
    def test_one_json_uses_defaults(self, json2h5ad_mock):
        converter = json2h5ad_mock.return_value
        converter.convert.return_value = "GSE1.h5ad"

        stdout = StringIO()
        with redirect_stdout(stdout):
            exit_code = main(["GSE1.json"])

        self.assertEqual(0, exit_code)
        converter.convert.assert_called_once_with(json_path="GSE1.json", out=".")
        self.assertEqual("", stdout.getvalue())

    @patch("meta_standards_converter.cli.json2h5ad.json2h5ad")
    def test_multiple_json_files_are_converted_in_order(self, json2h5ad_mock):
        converter = json2h5ad_mock.return_value
        converter.convert.return_value = "out.h5ad"

        with redirect_stdout(StringIO()):
            exit_code = main(["GSE1.json", "GSE2.json"])

        self.assertEqual(0, exit_code)
        self.assertEqual(
            [
                call(json_path="GSE1.json", out="."),
                call(json_path="GSE2.json", out="."),
            ],
            converter.convert.call_args_list,
        )

    @patch("meta_standards_converter.cli.json2h5ad.json2h5ad")
    def test_out_is_passed_to_converter(self, json2h5ad_mock):
        converter = json2h5ad_mock.return_value
        converter.convert.return_value = "out/GSE1.h5ad"

        with redirect_stdout(StringIO()):
            exit_code = main(["GSE1.json", "--out", "out"])

        self.assertEqual(0, exit_code)
        converter.convert.assert_called_once_with(json_path="GSE1.json", out="out")

    @patch("meta_standards_converter.cli.json2h5ad.json2h5ad")
    def test_failed_json_returns_one_and_continues(self, json2h5ad_mock):
        converter = json2h5ad_mock.return_value
        converter.convert.side_effect = [
            NotImplementedError("json2h5ad is not implemented yet"),
            "GSE2.h5ad",
        ]

        stdout = StringIO()
        stderr = StringIO()
        with redirect_stdout(stdout), redirect_stderr(stderr):
            exit_code = main(["GSE1.json", "GSE2.json"])

        self.assertEqual(1, exit_code)
        self.assertEqual(
            [
                call(json_path="GSE1.json", out="."),
                call(json_path="GSE2.json", out="."),
            ],
            converter.convert.call_args_list,
        )
        self.assertIn(
            "ERROR meta_standards_converter.cli.json2h5ad: GSE1.json: H5AD conversion failed",
            stdout.getvalue(),
        )
        self.assertIn("Traceback (most recent call last):", stdout.getvalue())
        self.assertIn("NotImplementedError: json2h5ad is not implemented yet", stdout.getvalue())
        self.assertEqual("", stderr.getvalue())

    @patch("meta_standards_converter.cli.json2h5ad.json2h5ad")
    def test_verbose_emits_success_logs_to_stdout(self, json2h5ad_mock):
        converter = json2h5ad_mock.return_value
        converter.convert.return_value = "GSE1.h5ad"

        stdout = StringIO()
        with redirect_stdout(stdout):
            exit_code = main(["GSE1.json", "-v"])

        self.assertEqual(0, exit_code)
        self.assertIn("INFO meta_standards_converter.cli.json2h5ad: GSE1.json: H5AD conversion started", stdout.getvalue())
        self.assertIn(
            "INFO meta_standards_converter.cli.json2h5ad: GSE1.json: converted to GSE1.h5ad",
            stdout.getvalue(),
        )

    @patch("meta_standards_converter.cli.json2h5ad.json2h5ad")
    def test_quiet_emits_only_errors(self, json2h5ad_mock):
        converter = json2h5ad_mock.return_value
        converter.convert.return_value = "GSE1.h5ad"

        stdout = StringIO()
        with redirect_stdout(stdout):
            exit_code = main(["GSE1.json", "--quiet"])

        self.assertEqual(0, exit_code)
        self.assertEqual("", stdout.getvalue())

    @patch("meta_standards_converter.cli.json2h5ad.json2h5ad")
    def test_log_file_writes_configured_logs(self, json2h5ad_mock):
        converter = json2h5ad_mock.return_value
        converter.convert.return_value = "GSE1.h5ad"

        with tempfile.TemporaryDirectory() as tmpdir:
            log_path = os.path.join(tmpdir, "json2h5ad.log")
            with redirect_stdout(StringIO()):
                exit_code = main(["GSE1.json", "-v", "--log-file", log_path])

            with open(log_path, encoding="utf-8") as handle:
                log_content = handle.read()

        self.assertEqual(0, exit_code)
        self.assertIn("INFO meta_standards_converter.cli.json2h5ad: GSE1.json: H5AD conversion started", log_content)
        self.assertIn("GSE1.json: converted to GSE1.h5ad", log_content)


if __name__ == "__main__":
    unittest.main()
