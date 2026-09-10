# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from meta_standards_converter.converters import JSON2OBSConverter
from tests.support.contracts import PBMC, prepare, assert_expected


def test_json2obs_public_counts_to_complete_sidecars(workspace):
    prepare(PBMC, workspace)
    result = JSON2OBSConverter().convert("miniml.json", outdir="out",
        asset_specs=["PBMC3K_SAMPLE=counts.tsv"], matrix_orientation="genes-by-observations",
        include_var=True, include_uns=True)
    assert not result.partial
    assert_expected(PBMC, "json2obs", workspace / "out", workspace)
