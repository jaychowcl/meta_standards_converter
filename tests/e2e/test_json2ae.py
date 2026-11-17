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

from tests.support.contracts import GEO, prepare, assert_expected


def test_json2ae_stored_miniml_to_complete_magetab(workspace, replay):
    prepare(GEO, workspace)
    converted(Converter(), 'miniml.json', out_type='magetab', preparation=PreparationPolicy(enrichment='off', expand_studies=False, linked_metadata=False, peer_provider=False), outdir='out', output_options={'platform_handler': 'single_cell_sequencing'})
    assert replay == []  # Stored enriched run evidence must suppress retrieval.
    assert_expected(GEO, "json2ae", workspace / "out", workspace)
