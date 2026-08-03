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
    @patch("meta_standards_converter.cli.json2h5ad.JSONDataOutputOrchestrator")
    def test_workflow_and_asset_options_are_forwarded(self, orchestrator_mock):
        orchestrator = orchestrator_mock.return_value
        orchestrator.export_h5ad.return_value = ConversionResult("GSE1", combined_h5ad="GSE1.h5ad")

        with redirect_stdout(StringIO()):
            exit_code = main([
                "GSE1.json",
                "--asset-manifest", "assets.csv",
                "--asset", "GSM1=local.h5ad",
                "--force-reprocess",
                "--pipeline", "scrnaseq",
                "--genome", "GRCh38",
                "--gff", "genes.gff3",
                "--profile", "apptainer",
                "--revision", "4.1.0",
                "--params-file", "params.json",
                "--nextflow-config", "nextflow.config",
                "--work-dir", "work",
                "--resume",
                "--processed-checkpoint-dir", "checkpoints",
                "--overwrite",
                "--allow-invalid",
                "--matrix-orientation", "genes-by-observations",
            ])

        self.assertEqual(0, exit_code)
        orchestrator.export_h5ad.assert_called_once_with(
            "GSE1.json",
            outdir=".",
            asset_manifest="assets.csv",
            asset_specs=["GSM1=local.h5ad"],
            force_reprocess=True,
            pipeline="scrnaseq",
            genome="GRCh38",
            gff="genes.gff3",
            profile="apptainer",
            revision="4.1.0",
            params_file="params.json",
            nextflow_config="nextflow.config",
            work_dir="work",
            resume=True,
            processed_checkpoint_dir="checkpoints",
            overwrite=True,
            allow_invalid=True,
            matrix_orientation="genes-by-observations",
        )

    @patch("meta_standards_converter.cli.json2h5ad.JSONDataOutputOrchestrator")
    def test_partial_conversion_returns_one(self, orchestrator_mock):
        orchestrator_mock.return_value.export_h5ad.return_value = ConversionResult(
            study_accession="GSE1",
            sample_h5ads={"GSM1": "GSM1.h5ad"},
            failures=["combined output incompatible"],
        )

        with redirect_stdout(StringIO()):
            exit_code = main(["GSE1.json"])

        self.assertEqual(1, exit_code)

    @patch("meta_standards_converter.cli.json2h5ad.JSONDataOutputOrchestrator")
    def test_one_json_uses_defaults(self, orchestrator_mock):
        orchestrator = orchestrator_mock.return_value
        orchestrator.export_h5ad.return_value = ConversionResult("GSE1", combined_h5ad="GSE1.h5ad")

        stdout = StringIO()
        with redirect_stdout(stdout):
            exit_code = main(["GSE1.json"])

        self.assertEqual(0, exit_code)
        orchestrator.export_h5ad.assert_called_once_with("GSE1.json", outdir=".")
        self.assertEqual("complete", __import__("json").loads(stdout.getvalue())["status"])

    @patch("meta_standards_converter.cli.json2h5ad.JSONDataOutputOrchestrator")
    def test_multiple_json_files_are_converted_in_order(self, orchestrator_mock):
        orchestrator = orchestrator_mock.return_value
        orchestrator.export_h5ad.return_value = ConversionResult("GSE1", combined_h5ad="out.h5ad")

        with redirect_stdout(StringIO()):
            exit_code = main(["GSE1.json", "GSE2.json"])

        self.assertEqual(0, exit_code)
        self.assertEqual(
            [
                call("GSE1.json", outdir="."),
                call("GSE2.json", outdir="."),
            ],
            orchestrator.export_h5ad.call_args_list,
        )

    @patch("meta_standards_converter.cli.json2h5ad.JSONDataOutputOrchestrator")
    def test_out_is_passed_to_converter(self, orchestrator_mock):
        orchestrator = orchestrator_mock.return_value
        orchestrator.export_h5ad.return_value = ConversionResult("GSE1", combined_h5ad="out/GSE1.h5ad")

        with redirect_stdout(StringIO()):
            exit_code = main(["GSE1.json", "--out", "out"])

        self.assertEqual(0, exit_code)
        orchestrator.export_h5ad.assert_called_once_with("GSE1.json", outdir="out")

    @patch("meta_standards_converter.cli.json2h5ad.JSONDataOutputOrchestrator")
    def test_failed_json_returns_one_and_continues(self, orchestrator_mock):
        orchestrator = orchestrator_mock.return_value
        orchestrator.export_h5ad.side_effect = [
            NotImplementedError("json2h5ad is not implemented yet"),
            ConversionResult("GSE2", combined_h5ad="GSE2.h5ad"),
        ]

        stdout = StringIO()
        stderr = StringIO()
        with redirect_stdout(stdout), redirect_stderr(stderr):
            exit_code = main(["GSE1.json", "GSE2.json"])

        self.assertEqual(1, exit_code)
        self.assertEqual(
            [
                call("GSE1.json", outdir="."),
                call("GSE2.json", outdir="."),
            ],
            orchestrator.export_h5ad.call_args_list,
        )
        self.assertIn(
            "ERROR meta_standards_converter.cli.json2h5ad: GSE1.json: H5AD conversion failed",
            stderr.getvalue(),
        )
        self.assertIn("Traceback (most recent call last):", stderr.getvalue())
        self.assertIn("NotImplementedError: json2h5ad is not implemented yet", stderr.getvalue())
        self.assertEqual("partial", __import__("json").loads(stdout.getvalue())["status"])

    @patch("meta_standards_converter.cli.json2h5ad.JSONDataOutputOrchestrator")
    def test_verbose_emits_success_logs_to_stderr(self, orchestrator_mock):
        orchestrator_mock.return_value.export_h5ad.return_value = ConversionResult("GSE1", combined_h5ad="GSE1.h5ad")

        stdout = StringIO()
        stderr = StringIO()
        with redirect_stdout(stdout), redirect_stderr(stderr):
            exit_code = main(["GSE1.json", "-v"])

        self.assertEqual(0, exit_code)
        self.assertIn("INFO meta_standards_converter.cli.json2h5ad: GSE1.json: H5AD conversion started", stderr.getvalue())
        self.assertIn(
            "INFO meta_standards_converter.cli.json2h5ad: GSE1.json: converted to GSE1.h5ad",
            stderr.getvalue(),
        )
        self.assertEqual("complete", __import__("json").loads(stdout.getvalue())["status"])

    @patch("meta_standards_converter.cli.json2h5ad.JSONDataOutputOrchestrator")
    def test_quiet_emits_only_errors(self, orchestrator_mock):
        orchestrator_mock.return_value.export_h5ad.return_value = ConversionResult("GSE1", combined_h5ad="GSE1.h5ad")

        stdout = StringIO()
        with redirect_stdout(stdout):
            exit_code = main(["GSE1.json", "--quiet"])

        self.assertEqual(0, exit_code)
        self.assertEqual("complete", __import__("json").loads(stdout.getvalue())["status"])

    @patch("meta_standards_converter.cli.json2h5ad.JSONDataOutputOrchestrator")
    def test_log_file_writes_configured_logs(self, orchestrator_mock):
        orchestrator_mock.return_value.export_h5ad.return_value = ConversionResult("GSE1", combined_h5ad="GSE1.h5ad")

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
