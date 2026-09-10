# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
import subprocess
from pathlib import Path
import anndata
from meta_standards_converter.converters import JSON2H5ADConverter
from meta_standards_converter.expression.nfcore import NFCoreRunner
from tests.support.contracts import PBMC, prepare


def test_raw_dispatch_reads_bounded_fake_pipeline_output(workspace):
    prepare(PBMC, workspace)
    calls = []
    def command_runner(command, **kwargs):
        calls.append(command)
        assert command[:3] == ["nextflow", "run", "nf-core/rnaseq"]
        folder = workspace / "out/nfcore/PBMC3K/rnaseq/results/star_salmon"
        folder.mkdir(parents=True)
        (folder / "salmon.merged.gene_counts.tsv").write_text("gene\tPBMC3K_SAMPLE\nENSG1\t7\nENSG2\t3\n")
        return subprocess.CompletedProcess(command, 0, stdout="completed\n", stderr="")
    runner = NFCoreRunner(command_runner=command_runner, which=lambda name: "/bounded-fake/" + name)
    result = JSON2H5ADConverter(pipeline_runner=runner).convert("miniml.json", out="out",
        asset_specs=["PBMC3K_SAMPLE=R1.fastq.gz", "PBMC3K_SAMPLE=R2.fastq.gz"],
        force_reprocess=True, pipeline="rnaseq", genome="GRCh38")
    assert len(calls) == 1
    assert result.pipeline_runs[0].returncode == 0
    assert not result.partial
    a = anndata.read_h5ad(result.sample_h5ads["PBMC3K_SAMPLE"])
    assert a.X.toarray().tolist() == [[7, 3]]
    assert list(a.var_names) == ["ENSG1", "ENSG2"]
    assert (workspace / "out/nfcore/PBMC3K/rnaseq/samplesheet.csv").read_text().splitlines()[-1].startswith("PBMC3K_SAMPLE,")
