# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Scientific discrepancies are explicit strict xfails, not approved goldens."""
import csv
import pytest
from meta_standards_converter.converters import GEO2AEConverter
from tests.support.contracts import GEO

class ChemistryMismatch(AssertionError):
    """Only the verified scientific contradiction is an expected failure."""


@pytest.mark.xfail(strict=True, raises=ChemistryMismatch, reason="MSC-TEST-001: auto detection interprets explicit 5-prime kit as 10x v2 3-prime")
def test_auto_handler_must_not_contradict_explicit_five_prime_protocol(workspace, replay):
    assert "Single Cell 5′" in (GEO / "inputs/geo.xml").read_text()
    GEO2AEConverter().convert("GSE328265", out="out")
    with (workspace / "out/E-GEOD-328265.sdrf.txt").open() as handle:
        header, row = list(csv.reader(handle, delimiter="\t"))
    if row[header.index("Comment[end bias]")] == "3 prime tag":
        raise ChemistryMismatch("5-prime source was rendered as 3-prime")
