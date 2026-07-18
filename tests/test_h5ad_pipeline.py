# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = os.path.dirname(os.path.dirname(__file__))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from meta_standards_converter.converters.json2h5ad import (  # noqa: E402
    Asset,
    NFCoreRunner,
    ReferenceResolver,
)


def package(sample_ids=("GSM1",), organism="Homo sapiens", taxid="9606"):
    return {
        "series": {"accession": [{"value": "GSE1"}]},
        "sample": [
            {
                "accession": [{"value": sample_id}],
                "library_source": "TRANSCRIPTOMIC",
                "channel": [{"organism": [{"taxid": taxid, "value": organism}]}],
            }
            for sample_id in sample_ids
        ],
    }


def raw_asset(sample_id="GSM1"):
    return Asset(
        scope_id=sample_id,
        path=f"https://example/{sample_id}_R1.fastq.gz",
        kind="raw",
        members=(
            {"uri": f"https://example/{sample_id}_R1.fastq.gz", "run": "SRR1"},
            {"uri": f"https://example/{sample_id}_R2.fastq.gz", "run": "SRR1"},
        ),
    )


class TestReferenceResolver(unittest.TestCase):
    def test_explicit_genome_wins_without_confirmation(self):
        reference = ReferenceResolver().resolve(
            [package()], genome="GRCh37", accept_inferred=False
        )

        self.assertEqual({"genome": "GRCh37"}, reference)

    def test_inferred_reference_requires_confirmation(self):
        with self.assertRaisesRegex(ValueError, "GRCh38.*accept-inferred-reference"):
            ReferenceResolver().resolve([package()], accept_inferred=False)

    def test_confirmed_human_reference_is_grch38(self):
        reference = ReferenceResolver().resolve([package()], accept_inferred=True)

        self.assertEqual({"genome": "GRCh38", "inferred": True}, reference)

    def test_custom_reference_requires_fasta_and_gtf(self):
        with self.assertRaisesRegex(ValueError, "fasta.*gtf"):
            ReferenceResolver().resolve([package()], fasta="genome.fa")


class TestNFCoreRunner(unittest.TestCase):
    def test_auto_pipeline_groups_bulk_and_single_cell_samples(self):
        calls = []

        def command_runner(command, **kwargs):
            calls.append(command)
            return subprocess.CompletedProcess(command, 0, stdout="completed\n", stderr="")

        with tempfile.TemporaryDirectory() as tmpdir:
            sc_path = (
                Path(tmpdir)
                / "nfcore" / "GSE1" / "scrnaseq" / "results"
                / "simpleaf" / "mtx_conversions" / "GSM1_filtered_matrix.h5ad"
            )
            sc_path.parent.mkdir(parents=True)
            sc_path.touch()
            bulk_dir = Path(tmpdir) / "nfcore" / "GSE1" / "rnaseq" / "results" / "star_salmon"
            bulk_dir.mkdir(parents=True)
            (bulk_dir / "salmon.merged.gene_counts.tsv").write_text("gene\tGSM2\nENSG1\t1\n")
            data = package(("GSM1", "GSM2"))
            data["sample"][0]["library_source"] = "single cell transcriptomic"
            runner = NFCoreRunner(command_runner=command_runner, which=lambda name: f"/usr/bin/{name}")

            result = runner.process(
                {"GSM1": raw_asset("GSM1"), "GSM2": raw_asset("GSM2")},
                packages=[data],
                out=tmpdir,
                study_accession="GSE1",
                pipeline="auto",
                genome="GRCh38",
            )

            self.assertEqual({"GSM1", "GSM2"}, set(result.assets))
            self.assertEqual({"scrnaseq", "rnaseq"}, {run.pipeline for run in result.runs})
            self.assertEqual(2, len(calls))

    def test_scrnaseq_run_writes_samplesheet_and_discovers_filtered_h5ad(self):
        calls = []

        def command_runner(command, **kwargs):
            calls.append((command, kwargs))
            return subprocess.CompletedProcess(command, 0, stdout="completed\n", stderr="")

        with tempfile.TemporaryDirectory() as tmpdir:
            result_dir = Path(tmpdir) / "nfcore" / "GSE1" / "scrnaseq" / "results"
            h5ad_path = result_dir / "simpleaf" / "mtx_conversions" / "GSM1_filtered_matrix.h5ad"
            h5ad_path.parent.mkdir(parents=True)
            h5ad_path.touch()
            runner = NFCoreRunner(command_runner=command_runner, which=lambda name: f"/usr/bin/{name}")

            result = runner.process(
                {"GSM1": raw_asset()},
                packages=[package()],
                out=tmpdir,
                study_accession="GSE1",
                pipeline="scrnaseq",
                genome="GRCh38",
                profile="docker",
            )

            samplesheet = Path(tmpdir) / "nfcore" / "GSE1" / "scrnaseq" / "samplesheet.csv"
            self.assertEqual(
                "sample,fastq_1,fastq_2\n"
                "GSM1,https://example/GSM1_R1.fastq.gz,https://example/GSM1_R2.fastq.gz\n",
                samplesheet.read_text(),
            )
            command = calls[0][0]
            self.assertEqual("nextflow", command[0])
            self.assertIn("nf-core/scrnaseq", command)
            self.assertIn("4.1.0", command)
            self.assertIn("docker", command)
            self.assertEqual(str(h5ad_path), result.assets["GSM1"].path)
            self.assertEqual("h5ad", result.assets["GSM1"].kind)
            self.assertEqual(0, result.runs[0].returncode)

    def test_rnaseq_run_exposes_merged_counts_for_each_sample(self):
        def command_runner(command, **kwargs):
            return subprocess.CompletedProcess(command, 0, stdout="completed\n", stderr="")

        with tempfile.TemporaryDirectory() as tmpdir:
            result_dir = Path(tmpdir) / "nfcore" / "GSE1" / "rnaseq" / "results" / "star_salmon"
            result_dir.mkdir(parents=True)
            counts = result_dir / "salmon.merged.gene_counts.tsv"
            counts.write_text("gene\tGSM1\tGSM2\nENSG1\t1\t2\n")
            tpm = result_dir / "salmon.merged.gene_tpm.tsv"
            tpm.write_text("gene\tGSM1\tGSM2\nENSG1\t3\t4\n")
            runner = NFCoreRunner(command_runner=command_runner, which=lambda name: f"/usr/bin/{name}")

            result = runner.process(
                {"GSM1": raw_asset("GSM1"), "GSM2": raw_asset("GSM2")},
                packages=[package(("GSM1", "GSM2"))],
                out=tmpdir,
                study_accession="GSE1",
                pipeline="rnaseq",
                genome="GRCh38",
            )

            self.assertEqual({"GSM1", "GSM2"}, set(result.assets))
            self.assertEqual("rnaseq_counts", result.assets["GSM1"].role)
            self.assertEqual(str(tpm), result.assets["GSM1"].features_path)

    def test_failed_nextflow_run_raises_with_log_path(self):
        def command_runner(command, **kwargs):
            return subprocess.CompletedProcess(command, 1, stdout="", stderr="pipeline failed")

        with tempfile.TemporaryDirectory() as tmpdir:
            runner = NFCoreRunner(command_runner=command_runner, which=lambda name: f"/usr/bin/{name}")

            with self.assertRaisesRegex(RuntimeError, "Nextflow.*nextflow.log"):
                runner.process(
                    {"GSM1": raw_asset()},
                    packages=[package()],
                    out=tmpdir,
                    study_accession="GSE1",
                    pipeline="scrnaseq",
                    genome="GRCh38",
                )


if __name__ == "__main__":
    unittest.main()
