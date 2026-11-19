# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Presets govern supplementary retrieval independently of the destination."""
from meta_standards_converter.converters import PreparationPolicy
from copy import deepcopy
import json

import pytest

from meta_standards_converter.converters import Converter
from tests.support.conversion import package

from tests.support.conversion import native, services


@pytest.mark.parametrize("preset,expected", [
    ("standard", ["peer", "linked", "standard"]),
    ("curators", ["standard"]),
    ("off", []),
])
@pytest.mark.parametrize("target", ["json", "csv"])
def test_preset_controls_all_enrichment(preset, expected, target):
    calls = []
    value = native()
    original = deepcopy(value)
    result = Converter(services=services(calls)).convert(value, out_type=target, preparation=PreparationPolicy(enrichment=preset, expand_studies=False))
    assert result.status == "complete", result.to_dict()
    assert [call[0] for call in calls] == expected
    assert value == original
    assert any(d.code.startswith("enrichment_") for d in result.items[0].diagnostics)


def test_file_and_object_share_preparation(tmp_path):
    source = tmp_path / "native.json"
    source.write_text(json.dumps(native()))
    a, b = [], []
    one = Converter(services=services(a)).convert(source, out_type='json', preparation=PreparationPolicy(enrichment='curators'))
    two = Converter(services=services(b)).convert(native(), out_type='json', preparation=PreparationPolicy(enrichment='curators'))
    assert one.status == two.status == "complete"
    assert a == b == [("standard", True, False)]
    assert one.items[0].payload == two.items[0].payload


def test_failed_enrichment_retains_valid_partial_payload():
    result = Converter(services=services([], failure=True)).convert(native(), out_type='json', preparation=PreparationPolicy(enrichment='curators'))
    assert result.status == "partial", result.to_dict()
    assert result.items[0].execution == "succeeded"
    assert result.items[0].payload[0].samples[0].iid == "SRS1"
    assert any(d.code == "source_partial" for d in result.items[0].diagnostics)


def test_off_magetab_does_not_resolve_hidden_evidence(monkeypatch):
    from meta_standards_converter.sources.pubmed import PubmedWebFetcher
    from unittest.mock import Mock
    call = Mock(side_effect=AssertionError("off must not retrieve publications"))
    monkeypatch.setattr(PubmedWebFetcher, "pubmed_summary", call)
    value = package()
    value["series"]["pubmed_id"] = ["123"]
    result = Converter().convert(value, out_type="magetab",
        preparation=PreparationPolicy(enrichment="off", expand_studies=False))
    assert result.status == "complete", result.to_dict()
    call.assert_not_called()


def test_curators_geo_excludes_sra_run_enrichment():
    calls = []
    value = package()
    value["source"]["format"] = "GEO"
    result = Converter(services=services(calls)).convert(value, out_type='json', preparation=PreparationPolicy(enrichment='curators'))
    assert result.status == "complete", result.to_dict()
    assert calls == [("standard", True, False)]


def test_publication_retrieval_is_deduplicated_per_convert_only(caplog):
    import logging
    caplog.set_level(logging.DEBUG, logger="meta_standards_converter")
    from meta_standards_converter.metadata.enrichment import MINiMLEnricher
    calls = []
    class Pubmed:
        def pubmed_summary(self, *, pubmed_id):
            calls.append(pubmed_id)
            return ('10.1/example', 'Authors', 'Title', 'published', 'EFO', 'EFO:1')
    enricher = MINiMLEnricher(pubmed_fetcher=Pubmed())
    converter = Converter(services={'enricher': enricher})
    sources = [package('GSE1', 'GSM1'), package('GSE2', 'GSM2')]
    for p in sources:
        p['series']['pubmed_id'] = ['12345']
    for _ in range(2):
        result = converter.convert(sources, out_type='json', preparation=PreparationPolicy(enrichment='curators'))
        assert result.status == 'complete', result.to_dict()
        assert all(p.series.to_mapping()['pubmed_publication'][0]['title'] == 'Title' for p in result.items[0].payload)
    assert calls == ['12345', '12345']
    completions = [record.msc_event for record in caplog.records
                   if hasattr(record, "msc_event") and record.msc_event["event_type"] == "conversion_finished"]
    assert len({event["run_id"] for event in completions}) == 2
    assert all(event["attributes"]["preparation_cache_hits"] >= 1 for event in completions)


def test_injected_services_are_used_in_stage_order():
    from meta_standards_converter.sources.metadata import GEOMetadataSource
    from meta_standards_converter.magetab.exporter import MAGETabExporter
    from meta_standards_converter.magetab.constructor import AEConstructor
    from meta_standards_converter.miniml.geo_parser import GEOParser
    from pathlib import Path
    calls = []
    class Fetcher:
        def fetch_gse_miniml(self, **kwargs):
            calls.append('fetch')
            return (Path(__file__).parents[1] / 'fixtures/studies/GSE100/inputs/geo.xml').read_text()
    class Parser:
        def parse(self, miniml, **kwargs):
            calls.append('parse')
            return GEOParser().parse(miniml)
    class Enricher:
        def enrich_selected(self, data, **kwargs):
            calls.append('enrich')
            return data
    class Constructor(AEConstructor):
        def miniml2magetab(self, **kwargs):
            calls.append('construct')
            return super().miniml2magetab(**kwargs)
    services = {'geo_source': GEOMetadataSource(fetcher=Fetcher(), parser=Parser()),
                'enricher': Enricher(), 'magetab_exporter': MAGETabExporter(Constructor())}
    result = Converter(services=services).convert('GSE100', out_type='magetab', preparation=PreparationPolicy(enrichment='curators', expand_studies=False))
    assert result.status == 'complete', result.to_dict()
    assert calls == ['fetch', 'parse', 'enrich', 'construct']


def test_saved_geo_provenance_is_reported():
    value = package()
    value['source']['format'] = 'GEO MINiML'
    result = Converter().convert(value, out_type='json', preparation=PreparationPolicy(enrichment='off'))
    assert result.items[0].provider == 'geo'
