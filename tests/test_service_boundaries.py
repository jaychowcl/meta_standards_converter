# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""MSC 6 service boundaries preserve metadata while exposing explicit collaborators."""
from pathlib import Path


def test_geo_parser_is_network_free_and_preserves_fixture_packages():
    from meta_standards_converter.miniml.geo_parser import GEOParser
    parser = GEOParser()
    packages = parser.parse((Path(__file__).parent / 'GSE328265_family.xml').read_text())
    assert len(packages) == 1
    assert packages[0].series.iid == "GSE328265"
    assert packages[0]["sample"][0]["iid"] == "GSM9651991"
    assert all(package.miniml_schema_version == '3.1' for package in packages)


def test_insdc_client_interception_preserves_empty_response_and_call_order():
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


def test_injected_asset_reader_is_used_for_catalogue(tmp_path):
    from tests.support.expression import make_expression_source as _source
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


def test_preparation_owns_provider_calls_before_network_free_export():
    from meta_standards_converter.magetab.constructor import AEConstructor
    from meta_standards_converter.metadata.enrichment import MINiMLEnricher
    from tests.support.json2ae import package
    from copy import deepcopy
    from unittest.mock import Mock
    value = package()
    value['series']['pubmed_id'] = ['123']
    value['sample'][0]['relation'] = [{'type': 'SRA', 'target': 'SRX1'}]
    original = deepcopy(value)
    pubmed, insdc = Mock(), Mock()
    pubmed.pubmed_summary.return_value = (None,) * 6
    insdc.extract_sra_accessions.return_value = ['SRX1']
    insdc.fetch_sra_runs.return_value = []
    prepared = MINiMLEnricher(pubmed_fetcher=pubmed, insdc_fetcher=insdc).enrich(value)
    pubmed.pubmed_summary.assert_called_once_with(pubmed_id='123')
    insdc.fetch_sra_runs.assert_called_once_with(accession='SRX1')
    pubmed.reset_mock(); insdc.reset_mock()
    for technology in ['bulk_sequencing', 'array']:
        AEConstructor().miniml2magetab(prepared, platform_handler=technology)
    pubmed.pubmed_summary.assert_not_called()
    insdc.fetch_sra_runs.assert_not_called()
    assert value == original


def test_rendering_retained_evidence_never_constructs_provider_clients(monkeypatch):
    from meta_standards_converter.magetab.constructor import AEConstructor
    from meta_standards_converter.sources import pubmed, insdc
    from tests.support.json2ae import package
    from unittest.mock import Mock
    from copy import deepcopy
    forbidden = Mock(side_effect=AssertionError('rendering must not construct a client'))
    monkeypatch.setattr(pubmed, 'PubmedWebFetcher', forbidden)
    monkeypatch.setattr(insdc, 'INSDCWebfetcher', forbidden)
    value = package()
    value['series']['pubmed_id'] = ['123']
    value['sample'][0]['sra_run'] = []
    value['sample'][0]['relation'] = [{'type': 'SRA', 'target': 'SRX1'}]
    original = deepcopy(value)
    constructor = AEConstructor()
    plan = constructor.prepare(value)
    assert constructor.render(plan)
    forbidden.assert_not_called()
    assert value == original
