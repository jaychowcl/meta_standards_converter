# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from tests.support.conversion import converted
from meta_standards_converter.converters import Converter, PreparationPolicy
from meta_standards_converter.magetab.exporter import MAGETabExporter
from meta_standards_converter.sources.metadata import GEOMetadataSource, MAGETabMetadataSource

from tests.support.contracts import GEO, assert_expected


def test_geo2json_public_evidence_to_complete_artifacts(workspace, replay):
    converted(Converter(), 'GSE328265', out_type='json', preparation=PreparationPolicy(enrichment='standard', expand_studies=False, linked_metadata=False, peer_provider=False), outdir='out')
    assert replay == ["geo:GSE328265", "pubmed:42129775", "sra:SRX32831930", "ena:SRX32831930"]
    assert_expected(GEO, "geo2json", workspace / "out", workspace)
