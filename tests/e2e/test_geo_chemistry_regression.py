# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Regression for MSC-TEST-001: explicit 5-prime evidence must survive rendering."""
from tests.support.conversion import converted
from meta_standards_converter.converters import Converter, PreparationPolicy
from meta_standards_converter.magetab.exporter import MAGETabExporter
from meta_standards_converter.sources.metadata import GEOMetadataSource, MAGETabMetadataSource
import csv

from tests.support.contracts import GEO

def test_auto_handler_must_not_contradict_explicit_five_prime_protocol(workspace, replay):
    assert "Single Cell 5′" in (GEO / "inputs/geo.xml").read_text()
    converted(Converter(), 'GSE328265', out_type='magetab', preparation=PreparationPolicy(enrichment='standard', expand_studies=False, linked_metadata=False, peer_provider=False), outdir='out')
    with (workspace / "out/E-GEOD-328265.sdrf.txt").open() as handle:
        header, row = list(csv.reader(handle, delimiter="\t"))
    assert "Comment[end bias]" not in header
    assert "Comment[library construction]" not in header
    assert 'Single Cell 5′' in (workspace / 'out/E-GEOD-328265.idf.txt').read_text()
    assert "Comment[cdna read size]" not in header
    assert "Comment[umi barcode size]" not in header
