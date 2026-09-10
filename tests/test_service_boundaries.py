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


def test_magetab_writer_and_obs_converter_are_independent_services():
    from meta_standards_converter.magetab.writer import MAGETabWriter
    from meta_standards_converter.converters.json2obs import JSON2OBSConverter
    assert callable(MAGETabWriter().write)
    assert callable(JSON2OBSConverter().convert)


def test_injected_asset_reader_is_used_for_catalogue(tmp_path):
    from tests.test_json_outputs_orchestrator import _source
    from meta_standards_converter.converters.json2h5ad import JSON2H5ADConverter
    from meta_standards_converter.expression.readers import ProcessedAssetReader
    source, asset = _source(tmp_path)
    seen = []
    class Reader:
        def read(self, asset, *, orientation, localize):
            seen.append(asset.scope_id)
            return ProcessedAssetReader().read(asset, orientation=orientation, localize=localize)
    result = JSON2H5ADConverter(reader=Reader()).convert(
        str(source), out=str(tmp_path / 'catalogue'), asset_specs=[f'GSM1={asset}'])
    assert seen == ['GSM1']
    assert list(result.sample_h5ads) == ['GSM1']


def test_processed_checkpoint_version_change_preserves_original_destination(tmp_path):
    from meta_standards_converter.expression.assets import Asset
    from meta_standards_converter.expression.checkpoints import ProcessedCheckpointStore

    arguments = dict(sample_id="GSM1", source_json_sha256="fixed", sample={},
                     asset=Asset("GSM1", "matrix.h5ad", "h5ad"), orientation="auto")
    old = ProcessedCheckpointStore(lambda: "5.0.0").key(tmp_path, **arguments)
    new = ProcessedCheckpointStore(lambda: "6.0.0").key(tmp_path, **arguments)
    assert old[2] != new[2]
    assert old[0] != new[0]
    assert old[1] != new[1]
    assert new == ProcessedCheckpointStore(lambda: "6.0.0").key(tmp_path, **arguments)


def test_owning_packages_export_injectable_interfaces():
    from meta_standards_converter.expression import AssetDiscovery, AssetReader, SourcePlanner
    from meta_standards_converter.sources import GEOXMLParser, MAGETabSourceResolver, PackageLoader
    from meta_standards_converter.metadata import MetadataEnrichment
    from meta_standards_converter.converters import GEO2JSONConverter, JSON2OBSConverter
    assert all(symbol is not None for symbol in (
        AssetDiscovery, AssetReader, SourcePlanner, GEOXMLParser, MAGETabSourceResolver,
        PackageLoader, MetadataEnrichment, GEO2JSONConverter, JSON2OBSConverter,
    ))
