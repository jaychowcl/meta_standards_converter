# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from meta_standards_converter.converters import GEO2JSONConverter
from tests.support.contracts import GEO, assert_expected


def test_geo2json_public_evidence_to_complete_artifacts(workspace, replay):
    GEO2JSONConverter().convert("GSE328265", out="out")
    assert replay == ["geo:GSE328265", "pubmed:42129775", "sra:SRX32831930", "ena:SRX32831930"]
    assert_expected(GEO, "geo2json", workspace / "out", workspace)
