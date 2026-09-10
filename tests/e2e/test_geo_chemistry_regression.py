# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Regression for MSC-TEST-001: explicit 5-prime evidence must survive rendering."""
import csv
from meta_standards_converter.converters import GEO2AEConverter
from tests.support.contracts import GEO

def test_auto_handler_must_not_contradict_explicit_five_prime_protocol(workspace, replay):
    assert "Single Cell 5′" in (GEO / "inputs/geo.xml").read_text()
    GEO2AEConverter().convert("GSE328265", out="out")
    with (workspace / "out/E-GEOD-328265.sdrf.txt").open() as handle:
        header, row = list(csv.reader(handle, delimiter="\t"))
    assert row[header.index("Comment[end bias]")] == "5 prime tag"
    assert "Comment[cdna read size]" not in header
    assert "Comment[umi barcode size]" not in header
