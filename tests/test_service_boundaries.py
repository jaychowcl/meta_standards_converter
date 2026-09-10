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
    assert all(package.miniml_schema_version == '3.0' for package in packages)




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






def test_magetab_explicit_evidence_clients_preserve_sra_then_publication_order():
    from meta_standards_converter.magetab.constructor import AEConstructor
    from meta_standards_converter.miniml import MINiMLCodec
    trace = []
    class INSDC:
        def extract_sra_accessions(self, value):
            return [value]
        def fetch_sra_runs(self, accession):
            trace.append(('sra', accession))
            return []
    class PubMed:
        def pubmed_summary(self, pubmed_id):
            trace.append(('pubmed', pubmed_id))
            return (None,) * 6
    package = MINiMLCodec().decode({
        'miniml_schema_version': '3.0', 'source': {'format': 'test'},
        'series': {'iid': 'GSE1', 'pubmed_id': ['123'], 'sample_ref': [{'ref': 'GSM1'}]},
        'sample': [{'iid': 'GSM1', 'relation': [{'type': 'SRA', 'target': 'SRX1'}],
                    'channel': [{'source': {'value': 'sample'}}]}],
    }).package
    constructor = AEConstructor(pubmed_client=PubMed(), insdc_client=INSDC())
    constructor.miniml2magetab(package, platform_handler='bulk_sequencing')
    assert trace == [('sra', 'SRX1'), ('pubmed', '123')]
    constructor.miniml2magetab(package, platform_handler='array')
    assert trace == [('sra', 'SRX1'), ('pubmed', '123'), ('pubmed', '123')]


def test_magetab_retained_evidence_suppresses_lookup_without_mutation():
    from unittest.mock import Mock
    from meta_standards_converter.metadata.enrichment import MAGETabEvidenceResolver
    pubmed, insdc = Mock(), Mock()
    evidence = MAGETabEvidenceResolver(pubmed_client=pubmed, insdc_client=insdc)
    handler = Mock()
    sample = {'iid': 'GSM1', 'sra_run': [], 'relation': [{'type': 'SRA', 'target': 'SRX1'}]}
    handler.ordered_samples.return_value = [sample]
    assert evidence.sample_runs(handler, 'bulk_sequencing') == {}
    assert evidence.publications({'series': {'pubmed_id': ['123'], 'pubmed_publication': [{}]}}) == []
    pubmed.pubmed_summary.assert_not_called()
    insdc.extract_sra_accessions.assert_not_called()
    insdc.fetch_sra_runs.assert_not_called()
    assert sample['sra_run'] == []
