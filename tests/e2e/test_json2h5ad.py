# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from meta_standards_converter.converters import JSON2H5ADConverter
from tests.support.contracts import PBMC, prepare, assert_expected


def test_json2h5ad_public_counts_to_complete_catalogue(workspace):
    prepare(PBMC, workspace)
    result = JSON2H5ADConverter().convert("miniml.json", out="out",
        asset_specs=["PBMC3K_SAMPLE=counts.tsv"], matrix_orientation="genes-by-observations")
    assert not result.partial
    assert result.combined_h5ad is None
    assert_expected(PBMC, "json2h5ad", workspace / "out", workspace)
