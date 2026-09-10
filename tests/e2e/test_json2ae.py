# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from meta_standards_converter.converters import JSON2AEConverter
from tests.support.contracts import GEO, prepare, assert_expected


def test_json2ae_stored_miniml_to_complete_magetab(workspace, replay):
    prepare(GEO, workspace)
    JSON2AEConverter().convert("miniml.json", out="out", enrich=False, platform_handler="single_cell_sequencing")
    assert replay == []  # Stored enriched run evidence must suppress retrieval.
    assert_expected(GEO, "json2ae", workspace / "out", workspace)
