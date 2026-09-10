# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from meta_standards_converter.converters import JSON2TSVConverter
from tests.support.contracts import GEO, prepare, assert_expected


def test_json2tsv_stored_miniml_to_complete_table(workspace):
    prepare(GEO, workspace)
    result = JSON2TSVConverter().export_manifest("miniml.json", outdir="out")
    assert not result.partial
    assert_expected(GEO, "json2tsv", workspace / "out", workspace)


def test_atlas_groups_have_complete_ordered_csv_and_diagnostics(workspace):
    from tests.support.contracts import FIXTURES
    case = FIXTURES / "edge_cases/atlas-groups"
    prepare(case, workspace)
    JSON2TSVConverter().export_manifest("atlas.json", outdir="out", output_format="csv")
    assert_expected(case, "json2tsv", workspace / "out", workspace)
