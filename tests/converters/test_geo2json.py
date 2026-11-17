# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""GEO source loading and shared parent-publication preparation."""
from copy import deepcopy
from unittest.mock import Mock

import pytest

from meta_standards_converter.converters import Converter, PreparationPolicy, RuntimeOptions
from meta_standards_converter.miniml import MINiMLCodec
from meta_standards_converter.sources.metadata import GEOMetadataSource
from tests.converters.test_converter_expansion import geo


def family(*, child_ids=(), parents=('GSE2',), parent_ids=('12345',), reciprocal=True):
    child = geo('GSE1')
    child['series']['pubmed_id'] = list(child_ids)
    child['series']['relation'] = [{'type': 'SubSeries of', 'target': p} for p in parents]
    parent = geo('GSE2')
    parent['series']['pubmed_id'] = list(parent_ids)
    parent['series']['relation'] = [{'type': 'SuperSeries of', 'target': 'GSE1' if reciprocal else 'GSE99'}]
    return child, parent


def convert_family(child, parent, *, enrichment='standard', expand=False, failure=None):
    source, enricher = Mock(), Mock()
    def load(accession, **options):
        if accession == 'GSE2' and failure:
            raise failure
        return [MINiMLCodec().decode(child if accession == 'GSE1' else parent).package]
    source.load.side_effect = load
    seen = []
    def enrich(data, **options):
        seen.append(data.to_mapping())
        return data
    enricher.enrich_selected.side_effect = enrich
    enricher.publication_issues = []
    enricher._sra_failures = 0
    batch = Converter(services={'geo_source': source, 'enricher': enricher}).convert(
        'GSE1', out_type='json', preparation=PreparationPolicy(enrichment=enrichment, expand_studies=expand))
    return batch, source, seen


def test_source_fetches_and_parses_without_preparing():
    fetcher, parser = Mock(), Mock()
    fetcher.fetch_gse_miniml.return_value = '<MINiML />'
    parser.parse.return_value = []
    source = GEOMetadataSource(fetcher=fetcher, parser=parser)
    assert source.load('GSE1', remove_empty=False) == []
    fetcher.fetch_gse_miniml.assert_called_once_with(gse='GSE1')
    parser.parse.assert_called_once_with('<MINiML />', remove_empty=False)


@pytest.mark.parametrize('expand', [False, True])
def test_reciprocal_parent_is_inherited_before_enrichment_without_extra_output(expand):
    child, parent = family()
    before = deepcopy((child, parent))
    batch, source, seen = convert_family(child, parent, expand=expand)
    assert batch.status == 'complete', batch.to_dict()
    payload = batch.items[0].payload
    assert [p.series.iid for p in payload] == (['GSE1', 'GSE2'] if expand else ['GSE1'])
    assert payload[0]['series']['pubmed_id'] == ['12345']
    assert payload[0]['extensions']['publication_inheritance'] == {
        'relation_type': 'SubSeries of', 'source_series': 'GSE2', 'pubmed_ids': ['12345']}
    assert seen[0]['series']['pubmed_id'] == ['12345']
    assert [c.args[0] for c in source.load.call_args_list] == ['GSE1', 'GSE2']
    assert (child, parent) == before


@pytest.mark.parametrize('changes', [
    {'parents': ('GSE2','GSE3')}, {'parent_ids': ()},
    {'parent_ids': ('123','456')}, {'reciprocal': False},
])
def test_ambiguous_or_nonreciprocal_evidence_does_not_invent_publication(changes):
    batch, source, seen = convert_family(*family(**changes))
    assert batch.status == 'complete', batch.to_dict()
    package = batch.items[0].payload[0]
    assert not package['series'].get('pubmed_id')
    assert 'publication_inheritance' not in package.get('extensions', {})


@pytest.mark.parametrize('preset', ['standard','curators'])
def test_direct_publication_prevents_parent_retrieval(preset):
    batch, source, seen = convert_family(*family(child_ids=('999',)), enrichment=preset)
    assert batch.status == 'complete', batch.to_dict()
    assert batch.items[0].payload[0]['series']['pubmed_id'] == ['999']
    assert source.load.call_count == 1


def test_off_and_no_expansion_disable_parent_retrieval():
    batch, source, seen = convert_family(*family(), enrichment='off')
    assert batch.status == 'complete'
    assert source.load.call_count == 1 and seen == []


def test_failed_parent_retrieval_retains_valid_partial_metadata():
    batch, source, seen = convert_family(*family(), failure=OSError('secret'))
    item = batch.items[0]
    assert item.execution == 'succeeded' and item.validation == 'valid'
    assert item.completeness == 'partial'
    assert item.dataset_ids == ('GSE1',)
    assert any('parent GSE2 publication retrieval failed' in d.message for d in item.diagnostics)
    assert 'secret' not in str(batch.to_dict())


def test_geo_loader_declares_known_repository_when_xml_omits_database():
    fetcher = Mock()
    fetcher.fetch_gse_miniml.return_value = '''<MINiML><Series iid="GSE1">
        <Accession database="GEO">GSE1</Accession><Sample-Ref ref="GSM1"/>
        </Series><Sample iid="GSM1"><Accession database="GEO">GSM1</Accession></Sample></MINiML>'''
    package = GEOMetadataSource(fetcher=fetcher).load('GSE1')[0]
    MINiMLCodec().decode(package, strict=True)
    assert package['database'] == [{'iid': 'GEO', 'name': 'Gene Expression Omnibus'}]


def test_geo_loader_does_not_invent_unknown_database_declarations():
    from meta_standards_converter.miniml.codec import MINiMLCompatibilityError
    fetcher = Mock()
    fetcher.fetch_gse_miniml.return_value = '<MINiML><Series iid="GSE1"><Accession database="unknown">GSE1</Accession></Series></MINiML>'
    package = GEOMetadataSource(fetcher=fetcher).load("GSE1")[0]
    with pytest.raises(MINiMLCompatibilityError):
        MINiMLCodec().decode(package, strict=True)
