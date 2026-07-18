# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
import os
import gzip
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


ROOT = os.path.dirname(os.path.dirname(__file__))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from meta_standards_converter.converters.json2h5ad import (  # noqa: E402
    AnnotationConverter,
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
        with self.assertRaisesRegex(ValueError, "fasta.*annotation"):
            ReferenceResolver().resolve([package()], fasta="genome.fa")

    def test_catalogue_genome_accepts_gtf_override(self):
        reference = ReferenceResolver().resolve(
            [package()], genome="GRCh38", gtf="genes.gtf"
        )

        self.assertEqual({"genome": "GRCh38", "gtf": "genes.gtf"}, reference)

    def test_catalogue_genome_accepts_gff_override(self):
        reference = ReferenceResolver().resolve(
            [package()], genome="GRCh38", gff="genes.gff3"
        )

        self.assertEqual({"genome": "GRCh38", "gff": "genes.gff3"}, reference)

    def test_custom_fasta_accepts_gff_annotation(self):
        reference = ReferenceResolver().resolve(
            [package()], fasta="genome.fa", gff="genes.gff3"
        )

        self.assertEqual({"fasta": "genome.fa", "gff": "genes.gff3"}, reference)

    def test_gtf_and_gff_are_mutually_exclusive(self):
        with self.assertRaisesRegex(ValueError, "gtf.*gff"):
            ReferenceResolver().resolve(
                [package()], genome="GRCh38", gtf="genes.gtf", gff="genes.gff3"
            )

    def test_annotation_without_genome_or_fasta_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "annotation.*genome.*fasta"):
            ReferenceResolver().resolve([package()], gtf="genes.gtf")


class TestAnnotationConverter(unittest.TestCase):
    def test_gtf_is_retained_with_path_format_and_checksum(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            gtf = Path(tmpdir) / "genes.gtf"
            gtf.write_text("chr1\ttest\texon\t1\t2\t.\t+\t.\tgene_id \"g1\";\n")

            reference = AnnotationConverter().prepare(
                {"genome": "GRCh38", "gtf": str(gtf)}, Path(tmpdir) / "reference"
            )

            self.assertEqual(str(gtf), reference["gtf"])
            self.assertEqual(str(gtf), reference["annotation_source"])
            self.assertEqual("gtf", reference["annotation_format"])
            self.assertEqual(hashlib.sha256(gtf.read_bytes()).hexdigest(), reference["annotation_sha256"])

    def test_gff3_is_converted_to_checksum_addressed_gtf_and_reused(self):
        calls = []

        def runner(command, **kwargs):
            calls.append(command)
            Path(command[-1]).write_text("converted gtf\n")
            return subprocess.CompletedProcess(command, 0, stdout="", stderr="")

        with tempfile.TemporaryDirectory() as tmpdir:
            gff = Path(tmpdir) / "genes.gff3.gz"
            with gzip.open(gff, "wt", encoding="utf-8") as handle:
                handle.write("##gff-version 3\nchr1\ttest\tgene\t1\t2\t.\t+\t.\tID=g1\n")
            fasta = Path(tmpdir) / "genome.fa"
            fasta.write_text(">chr1\nAC\n")
            converter = AnnotationConverter(
                command_runner=runner,
                which=lambda name: f"/usr/bin/{name}",
            )
            reference_dir = Path(tmpdir) / "reference"

            first = converter.prepare({"fasta": str(fasta), "gff": str(gff)}, reference_dir)
            second = converter.prepare({"fasta": str(fasta), "gff": str(gff)}, reference_dir)

            digest = hashlib.sha256(gff.read_bytes()).hexdigest()
            self.assertEqual(str(reference_dir / f"{digest}.gtf"), first["gtf"])
            self.assertEqual(first, second)
            self.assertEqual("gff3", first["annotation_format"])
            self.assertEqual(digest, first["annotation_sha256"])
            self.assertNotIn("gff", first)
            self.assertEqual(1, len(calls))
            self.assertEqual(["gffread", str(gff), "-T", "-o", first["gtf"] + ".tmp"], calls[0])

    def test_gff3_requires_gffread(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            gff = Path(tmpdir) / "genes.gff3"
            gff.write_text("##gff-version 3\n")

            with self.assertRaisesRegex(RuntimeError, "gffread"):
                AnnotationConverter(which=lambda _name: None).prepare(
                    {"genome": "GRCh38", "gff": str(gff)}, Path(tmpdir) / "reference"
                )

    def test_failed_gff3_conversion_removes_temporary_output(self):
        def runner(command, **kwargs):
            Path(command[-1]).write_text("partial\n")
            return subprocess.CompletedProcess(command, 1, stdout="", stderr="bad annotation")

        with tempfile.TemporaryDirectory() as tmpdir:
            gff = Path(tmpdir) / "genes.gff3"
            gff.write_text("##gff-version 3\n")
            reference_dir = Path(tmpdir) / "reference"

            with self.assertRaisesRegex(RuntimeError, "bad annotation"):
                AnnotationConverter(command_runner=runner, which=lambda name: name).prepare(
                    {"genome": "GRCh38", "gff": str(gff)}, reference_dir
                )

            self.assertEqual([], list(reference_dir.glob("*.tmp")))

    def test_empty_gff3_conversion_is_rejected(self):
        def runner(command, **kwargs):
            return subprocess.CompletedProcess(command, 0, stdout="", stderr="")

        with tempfile.TemporaryDirectory() as tmpdir:
            gff = Path(tmpdir) / "genes.gff3"
            gff.write_text("##gff-version 3\n")

            with self.assertRaisesRegex(RuntimeError, "produced no GTF"):
                AnnotationConverter(command_runner=runner, which=lambda name: name).prepare(
                    {"genome": "GRCh38", "gff": str(gff)}, Path(tmpdir) / "reference"
                )


class TestNFCoreRunner(unittest.TestCase):
    def test_generated_params_keep_catalogue_genome_and_override_gtf(self):
        captured_params = []

        def command_runner(command, **kwargs):
            params_path = Path(command[command.index("-params-file") + 1])
            captured_params.append(json.loads(params_path.read_text()))
            result_dir = Path(captured_params[-1]["outdir"])
            output = result_dir / "simpleaf" / "mtx_conversions" / "GSM1_filtered_matrix.h5ad"
            output.parent.mkdir(parents=True)
            output.touch()
            return subprocess.CompletedProcess(command, 0, stdout="", stderr="")

        with tempfile.TemporaryDirectory() as tmpdir:
            gtf = Path(tmpdir) / "genes.gtf"
            gtf.write_text("chr1\ttest\texon\t1\t2\t.\t+\t.\tgene_id \"g1\";\n")
            supplied = Path(tmpdir) / "params.json"
            supplied.write_text(json.dumps({"genome": "wrong", "gtf": "wrong", "protocol": "10XV3"}))
            runner = NFCoreRunner(command_runner=command_runner, which=lambda name: f"/usr/bin/{name}")

            runner.process(
                {"GSM1": raw_asset()}, packages=[package()], out=tmpdir,
                study_accession="GSE1", pipeline="scrnaseq", genome="GRCh38",
                gtf=str(gtf), params_file=str(supplied),
            )

            self.assertEqual("GRCh38", captured_params[0]["genome"])
            self.assertEqual(str(gtf), captured_params[0]["gtf"])
            self.assertEqual("10XV3", captured_params[0]["protocol"])

    def test_samplesheet_upgrades_known_archive_ftp_urls_to_https(self):
        asset = Asset(
            scope_id="GSM1",
            path="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/GSM1_1.fastq.gz",
            kind="raw",
            members=(
                {
                    "uri": "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/GSM1_1.fastq.gz",
                    "run": "SRR1",
                },
                {
                    "uri": "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/GSM1_2.fastq.gz",
                    "run": "SRR1",
                },
            ),
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            samplesheet = Path(tmpdir) / "samplesheet.csv"
            NFCoreRunner()._write_samplesheet(samplesheet, {"GSM1": asset}, "scrnaseq")

            self.assertEqual(
                "sample,fastq_1,fastq_2\n"
                "GSM1,https://ftp.sra.ebi.ac.uk/vol1/fastq/GSM1_1.fastq.gz,"
                "https://ftp.sra.ebi.ac.uk/vol1/fastq/GSM1_2.fastq.gz\n",
                samplesheet.read_text(),
            )

    def test_required_rootless_docker_is_accepted(self):
        probes = []

        def runtime_runner(command, **kwargs):
            probes.append(command)
            return subprocess.CompletedProcess(
                command,
                0,
                stdout='["name=seccomp","name=rootless","name=cgroupns"]\n',
                stderr="",
            )

        runner = NFCoreRunner(
            which=lambda name: f"/usr/bin/{name}",
            runtime_runner=runtime_runner,
        )
        with patch.dict(
            os.environ,
            {"META_STANDARDS_REQUIRE_ROOTLESS_DOCKER": "1"},
        ):
            runner._preflight("docker")

        self.assertEqual(
            [["docker", "info", "--format", "{{json .SecurityOptions}}"]],
            probes,
        )

    def test_required_rootless_docker_rejects_rootful_daemon(self):
        def runtime_runner(command, **kwargs):
            return subprocess.CompletedProcess(
                command,
                0,
                stdout='["name=seccomp","name=cgroupns"]\n',
                stderr="",
            )

        runner = NFCoreRunner(
            which=lambda name: f"/usr/bin/{name}",
            runtime_runner=runtime_runner,
        )
        with patch.dict(
            os.environ,
            {"META_STANDARDS_REQUIRE_ROOTLESS_DOCKER": "true"},
        ):
            with self.assertRaisesRegex(RuntimeError, "rootless Docker daemon"):
                runner._preflight("docker")

    def test_required_rootless_docker_reports_unreachable_daemon(self):
        def runtime_runner(command, **kwargs):
            return subprocess.CompletedProcess(
                command,
                1,
                stdout="",
                stderr="Cannot connect to the Docker daemon",
            )

        runner = NFCoreRunner(
            which=lambda name: f"/usr/bin/{name}",
            runtime_runner=runtime_runner,
        )
        with patch.dict(
            os.environ,
            {"META_STANDARDS_REQUIRE_ROOTLESS_DOCKER": "yes"},
        ):
            with self.assertRaisesRegex(RuntimeError, "Cannot connect.*Docker daemon"):
                runner._preflight("docker")

    def test_relative_output_still_writes_absolute_nextflow_parameters(self):
        captured = []

        def command_runner(command, **kwargs):
            captured.append(command)
            return subprocess.CompletedProcess(command, 0, stdout="", stderr="")

        with tempfile.TemporaryDirectory(dir=os.getcwd()) as tmpdir:
            relative_out = os.path.relpath(tmpdir, os.getcwd())
            result_path = (
                Path(tmpdir).resolve()
                / "nfcore" / "GSE1" / "scrnaseq" / "results"
                / "simpleaf" / "mtx_conversions" / "GSM1_filtered_matrix.h5ad"
            )
            result_path.parent.mkdir(parents=True)
            result_path.touch()
            runner = NFCoreRunner(command_runner=command_runner, which=lambda name: f"/usr/bin/{name}")

            runner.process(
                {"GSM1": raw_asset()},
                packages=[package()],
                out=relative_out,
                study_accession="GSE1",
                pipeline="scrnaseq",
                genome="GRCh38",
            )

            params_index = captured[0].index("-params-file") + 1
            work_index = captured[0].index("-work-dir") + 1
            self.assertTrue(os.path.isabs(captured[0][params_index]))
            self.assertTrue(os.path.isabs(captured[0][work_index]))

    def test_auto_pipeline_groups_bulk_and_single_cell_samples(self):
        calls = []

        def command_runner(command, **kwargs):
            calls.append(command)
            return subprocess.CompletedProcess(command, 0, stdout="completed\n", stderr="")

        with tempfile.TemporaryDirectory() as tmpdir:
            gtf = Path(tmpdir) / "genes.gtf"
            gtf.write_text("chr1\ttest\texon\t1\t2\t.\t+\t.\tgene_id \"g1\";\n")
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
                gtf=str(gtf),
            )

            self.assertEqual({"GSM1", "GSM2"}, set(result.assets))
            self.assertEqual({"scrnaseq", "rnaseq"}, {run.pipeline for run in result.runs})
            self.assertEqual(2, len(calls))
            self.assertEqual(
                {str(gtf)},
                {asset.effective_annotation for asset in result.assets.values()},
            )

    def test_scrnaseq_run_writes_samplesheet_and_discovers_filtered_h5ad(self):
        calls = []

        def command_runner(command, **kwargs):
            calls.append((command, kwargs))
            return subprocess.CompletedProcess(
                command,
                0,
                stdout=(
                    "WARN: Unrecognized config option 'validation.example'\n"
                    "WARN: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
                    "  A multiline pipeline warning.\n"
                    "  More detail.\n"
                    "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
                    "completed\n"
                ),
                stderr="WARN: Unrecognized config option 'validation.example'\n",
            )

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
            self.assertIn("4.2.0", command)
            self.assertNotIn("4.1.0", command)
            self.assertIn("docker", command)
            self.assertEqual(str(h5ad_path), result.assets["GSM1"].path)
            self.assertEqual("h5ad", result.assets["GSM1"].kind)
            self.assertEqual(0, result.runs[0].returncode)
            self.assertEqual(
                [
                    "Unrecognized config option 'validation.example'",
                    "A multiline pipeline warning. More detail.",
                ],
                result.runs[0].warnings,
            )

    def test_scrnaseq_prefers_qcatch_filtered_quants_over_raw_matrix(self):
        def command_runner(command, **kwargs):
            return subprocess.CompletedProcess(command, 0, stdout="completed\n", stderr="")

        with tempfile.TemporaryDirectory() as tmpdir:
            result_dir = Path(tmpdir) / "nfcore" / "GSE1" / "scrnaseq" / "results"
            raw_path = (
                result_dir
                / "simpleaf" / "mtx_conversions" / "GSM1" / "GSM1_raw_matrix.h5ad"
            )
            filtered_path = (
                result_dir
                / "simpleaf" / "GSM1" / "qcatch" / "GSM1_filtered_quants.h5ad"
            )
            raw_path.parent.mkdir(parents=True)
            filtered_path.parent.mkdir(parents=True)
            raw_path.touch()
            filtered_path.touch()
            runner = NFCoreRunner(
                command_runner=command_runner,
                which=lambda name: f"/usr/bin/{name}",
            )

            result = runner.process(
                {"GSM1": raw_asset()},
                packages=[package()],
                out=tmpdir,
                study_accession="GSE1",
                pipeline="scrnaseq",
                genome="GRCh38",
            )

            self.assertEqual(str(filtered_path), result.assets["GSM1"].path)
            self.assertEqual(
                {str(raw_path), str(filtered_path)},
                set(result.retained_h5ads),
            )

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
