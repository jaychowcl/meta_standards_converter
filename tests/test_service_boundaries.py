# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""MSC 6 service boundaries preserve metadata while exposing explicit collaborators."""
import importlib
from pathlib import Path


def test_geo_parser_is_network_free_and_preserves_fixture_packages():
    from meta_standards_converter.miniml.geo_parser import GEOParser
    parser = GEOParser()
    assert not hasattr(parser, 'geo_fetcher')
    packages = parser.parse((Path(__file__).parent / 'GSE328265_family.xml').read_text())
    assert packages
    assert all(package.miniml_schema_version == '3.0' for package in packages)


def test_services_have_direct_imports():
    for module, symbol in [
        ('sources.json', 'JSONPackageSource'),
        ('sources.insdc', 'INSDCWebfetcher'),
        ('sources.magetab', 'AEWebFetcher'),
        ('metadata.enrichment', 'MINiMLEnricher'),
        ('converters.geo2json', 'GEO2JSONConverter'),
        ('converters.geo2ae', 'GEO2AEConverter'),
        ('converters.ae2json', 'AE2JSONConverter'),
        ('converters.json2ae', 'JSON2AEConverter'),
    ]:
        assert getattr(importlib.import_module('meta_standards_converter.' + module), symbol)


def test_insdc_client_interception_preserves_run_parsing():
    from meta_standards_converter.sources.insdc import INSDCWebfetcher
    from xml.etree.ElementTree import fromstring
    calls = []
    class Client:
        def fetch_sra_xml(self, nrx):
            calls.append(('sra', nrx))
            return fromstring('<EXPERIMENT_PACKAGE_SET/>')
        def fetch_ena_file_report(self, accession):
            calls.append(('ena', accession))
            return []
    service = INSDCWebfetcher(client=Client())
    assert service.fetch_sra_runs('SRP1') == []
    assert calls == [('sra', 'SRP1'), ('ena', 'SRP1')]
