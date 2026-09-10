# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
import json
from meta_standards_converter.converters import AE2JSONConverter
from tests.support.contracts import AE, prepare, assert_expected


def test_ae2json_local_resources_to_complete_miniml(workspace, caplog):
    prepare(AE, workspace)
    AE2JSONConverter().convert("E-MTAB-6486.idf.txt", out="out")
    assert [r.getMessage() for r in caplog.records if r.levelname == "WARNING"] == json.loads((AE / "diagnostics.json").read_text())
    assert_expected(AE, "ae2json", workspace / "out", workspace)


def test_constructed_material_type_edge_retains_full_output(workspace):
    from tests.support.contracts import FIXTURES
    case = FIXTURES / "edge_cases/magetab-normalization"
    prepare(case, workspace)
    AE2JSONConverter().convert("E-MTAB-6486.idf.txt", out="out")
    assert_expected(case, "ae2json", workspace / "out", workspace)
