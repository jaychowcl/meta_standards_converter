# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""CLI dispatch runs real converters; only HTTP and machine memory are controlled."""
import importlib
import json
import pytest
from tests.support.contracts import GEO, PBMC, AE, prepare, assert_expected

CASES = [
    ("geo2json", GEO, ["GSE328265", "--out", "out"]),
    ("geo2ae", GEO, ["GSE328265", "--out", "out", "--platform-handler", "single_cell_sequencing"]),
    ("json2ae", GEO, ["miniml.json", "--out", "out", "--no-enrich", "--platform-handler", "single_cell_sequencing"]),
    ("ae2json", AE, ["E-MTAB-6486.idf.txt", "--out", "out"]),
    ("json2tsv", GEO, ["miniml.json", "--out", "out"]),
    ("json2h5ad", PBMC, ["miniml.json", "--out", "out", "--asset", "PBMC3K_SAMPLE=counts.tsv", "--matrix-orientation", "genes-by-observations"]),
    ("json2obs", PBMC, ["miniml.json", "--outdir", "out", "--asset", "PBMC3K_SAMPLE=counts.tsv", "--matrix-orientation", "genes-by-observations", "--include-var", "--include-uns"]),
]

@pytest.mark.parametrize("name,case,argv", CASES, ids=[item[0] for item in CASES])
def test_cli_dispatch_publishes_reviewed_outputs(name, case, argv, workspace, replay, capsys):
    prepare(case, workspace)
    assert importlib.import_module("meta_standards_converter.cli." + name).main(argv) == 0
    assert_expected(case, name, workspace / "out", workspace)
    captured = capsys.readouterr()
    if name in {"json2tsv", "json2h5ad", "json2obs"}:
        assert json.loads(captured.out)["status"] == "complete"
    elif name == "ae2json":
        expected_warnings = json.loads((AE / "diagnostics.json").read_text())
        assert captured.out.splitlines() == ["WARNING meta_standards_converter.magetab.parser: " + text for text in expected_warnings]
        assert captured.err == ""
    else:
        assert captured.out == ""
    if name.startswith("geo2"):
        assert replay == ["geo:GSE328265", "pubmed:42129775", "sra:SRX32831930", "ena:SRX32831930"]
    else:
        assert replay == []
