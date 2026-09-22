# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Facade presets select values without changing scientific identity joins."""
from copy import deepcopy
from meta_standards_converter.converters import Converter
from meta_standards_converter.miniml import MINiMLCodec
from meta_standards_converter.metadata.archive_enrichment import LinkedArchiveEnricher


def package(provider='SRA', layout='SINGLE'):
    accession = {'SRA': 'SRP1', 'GEO': 'GSE1', 'MAGE-TAB': 'E-MTAB-1'}[provider]
    data = {'miniml_schema_version': '3.0', 'source': {'format': provider},
        'database': [{'iid': p, 'name': p} for p in ['SRA','GEO','ArrayExpress']],
        'series': {'iid': accession, 'title': provider, 'accession': [{'value': accession}],
                   'sample_ref': [{'ref': 'SRS1'}], 'relation': [{'type': 'SRA', 'target': 'SRP1'}],
                   'assay_paths': [{'steps': [{'kind': 'source', 'name': 'SRS1', 'sample_ref': 'SRS1'},
                       {'kind': 'assay', 'name': 'SRX1'}, {'kind': 'scan', 'name': 'SRR1',
                        'comments': [{'name': 'LIBRARY_LAYOUT', 'value': layout}]}]}]},
        'sample': [{'iid': 'SRS1', 'accession': [{'value': 'SRS1', 'database': 'SRA'}],
                    'channel': [{'organism': [{'value': 'Homo sapiens', 'taxid': '9606'}]}],
                    'sra_run': [{'run': 'SRR1', 'experiment': 'SRX1', 'sample': 'SRS1',
                                 'library_strategy': 'RNA-Seq', 'library_layout': layout}]}]}
    return data


class Frozen:
    def __init__(self, data): self.data = data; self.calls = []
    def convert(self, accession, **kwargs):
        self.calls.append(accession)
        return [MINiMLCodec().decode(self.data, strict=True).package]


class NoPublication:
    publication_issues = []
    def enrich_selected(self, package, **kwargs): return package


def convert(data, incoming, *, mode='standard'):
    converter = Frozen(incoming)
    linked = LinkedArchiveEnricher(geo_converter=converter, ae_converter=converter)
    from meta_standards_converter.converters.archive_results import ArchiveImportResult, StudyImportOutcome
    class Peer:
        def convert(self, accession, **kwargs):
            return ArchiveImportResult(accession, [StudyImportOutcome(accession, accession, 'complete', MINiMLCodec().decode(data, strict=True).package)])
    facade = Converter(services={'linked_enricher': linked, 'enricher': NoPublication(), 'ena2json': Peer()})
    kwargs = {'enrichment': mode}
    result = facade.convert(MINiMLCodec().decode(data, strict=True).package, out_type='json',
        options={'expand_studies': False}, **kwargs)
    return result, converter.calls


def test_standard_selects_ae_layout_with_conflict_evidence_and_partial_valid():
    data, incoming = package(), package('MAGE-TAB', 'PAIRED')
    data['series']['relation'].append({'type': 'ArrayExpress', 'target': 'E-MTAB-1'})
    result, calls = convert(data, incoming)
    assert calls == ['E-MTAB-1']
    item = result.items[0]
    assert item.completeness == 'partial' and item.validation == 'valid', result.to_dict()
    output = item.payload[0]
    assert output['sample'][0]['sra_run'][0]['library_layout'] == 'PAIRED'
    assert any('conflict' in d.message.lower() and 'SINGLE' in d.message and 'PAIRED' in d.message for d in item.diagnostics)
    records = output['extensions']['insdc']['records']
    conflict = next(r['metadata'] for r in records if r['kind'] == 'source_conflict' and r['metadata']['field'] == 'library_layout')
    assert {s['value'] for s in conflict['statements']} == {'SINGLE', 'PAIRED'}
    assert conflict['selected']['provider'] == 'ArrayExpress'
    assert conflict['path'].endswith('/SRR1/library_layout')


def test_standard_direct_ae_does_not_lose_to_linked_geo():
    data, incoming = package('MAGE-TAB', 'PAIRED'), package('GEO', 'SINGLE')
    data['series']['relation'].append({'type': 'GEO', 'target': 'GSE1'})
    result, calls = convert(data, incoming)
    assert calls
    output = result.items[0].payload[0]
    assert output['series']['title'] == 'MAGE-TAB'
    assert output['sample'][0]['sra_run'][0]['library_layout'] == 'PAIRED'
    assert result.items[0].validation == 'valid'
    assert any('conflict' in d.message.lower() for d in result.items[0].diagnostics)


def test_curators_warns_from_supplied_evidence_without_repository_calls():
    data, incoming = package(), package('MAGE-TAB', 'PAIRED')
    data['extensions'] = {'insdc': {'version': '1.0', 'records': [
        {'provider': 'ArrayExpress', 'kind': 'MINiML', 'accession': 'E-MTAB-1', 'metadata': deepcopy(incoming)}]}}
    result, calls = convert(data, incoming, mode='curators')
    assert calls == []
    item = result.items[0]
    assert item.payload[0]['sample'][0]['sra_run'][0]['library_layout'] == 'SINGLE'
    assert item.completeness == 'partial' and item.validation == 'valid', result.to_dict()
    assert any('conflict' in d.message.lower() for d in item.diagnostics)


def test_off_preserves_recorded_conflict_without_reconciliation():
    data, incoming = package(), package('MAGE-TAB', 'PAIRED')
    data['extensions'] = {'insdc': {'version': '1.0', 'records': [
        {'provider': 'ArrayExpress', 'kind': 'MINiML', 'accession': 'E-MTAB-1', 'metadata': deepcopy(incoming)}]}}
    result, calls = convert(data, incoming, mode='off')
    assert calls == []
    assert result.status == 'complete'
    assert result.items[0].payload[0]['sample'][0]['sra_run'][0]['library_layout'] == 'SINGLE'
    assert not any('conflict' in d.message for d in result.items[0].diagnostics)


def test_equal_rank_peer_keeps_native_and_fills_missing_fact():
    from meta_standards_converter.metadata.archive_enrichment import merge_archive_metadata
    from meta_standards_converter.metadata.reconciliation import selection_policy
    a, b = package(), package(layout='PAIRED')
    b['source']['format'] = 'ENA'
    b['sample'][0]['sra_run'][0]['library_selection'] = 'cDNA'
    with selection_policy('standard'):
        merged, issues = merge_archive_metadata(MINiMLCodec().decode(a).package, MINiMLCodec().decode(b).package)
    run = merged.to_mapping()['sample'][0]['sra_run'][0]
    assert run['library_layout'] == 'SINGLE'
    assert run['library_selection'] == 'cDNA'
    assert any('conflict' in i for i in issues)


def test_standard_keeps_distinct_biological_occurrences():
    data, incoming = package(), package('MAGE-TAB')
    data['sample'][0]['channel'][0]['characteristics'] = [{'name': 'tissue', 'value': 'lung'}]
    incoming['sample'][0]['channel'][0]['characteristics'] = [{'name': 'tissue', 'value': 'heart'}]
    data['series']['relation'].append({'type': 'ArrayExpress', 'target': 'E-MTAB-1'})
    result, _ = convert(data, incoming)
    values = result.items[0].payload[0]['sample'][0]['channel'][0]['characteristics']
    assert [v['value'] for v in values if v['name'] == 'tissue'] == ['lung', 'heart']


def test_selected_source_provenance_prevents_later_lower_rank_overwrite():
    from meta_standards_converter.metadata.archive_enrichment import merge_archive_metadata
    from meta_standards_converter.metadata.reconciliation import selection_policy
    codec = MINiMLCodec()
    with selection_policy('standard'):
        first, _ = merge_archive_metadata(codec.decode(package()).package, codec.decode(package('MAGE-TAB', 'PAIRED')).package, prefer=True)
        second, issues = merge_archive_metadata(first, codec.decode(package('GEO', 'SINGLE')).package, prefer=True)
    data = second.to_mapping()
    assert data['series']['title'] == 'MAGE-TAB'
    assert data['sample'][0]['sra_run'][0]['library_layout'] == 'PAIRED'
    assert any('conflict' in i for i in issues)


def test_lower_rank_linked_workflow_does_not_replace_authored_ae_material():
    data, incoming = package('MAGE-TAB', 'PAIRED'), package('GEO', 'SINGLE')
    data['series']['relation'].append({'type': 'GEO', 'target': 'GSE1'})
    for record, material in [(data, 'total RNA'), (incoming, 'genomic DNA')]:
        record['series']['assay_paths'][0]['steps'].insert(1,
            {'kind': 'extract', 'name': 'extract1', 'material_type': {'value': material}})
    result, _ = convert(data, incoming)
    paths = result.items[0].payload[0]['series']['assay_paths']
    assert {s['material_type']['value'] for p in paths for s in p['steps'] if s['kind'] == 'extract'} == {'total RNA'}


def test_curators_compares_recorded_library_path_facts_without_fetching():
    data, incoming = package(), package('MAGE-TAB', 'PAIRED')
    # Deposited AE imports hold library facts on the authored assay path.
    incoming['sample'][0]['sra_run'][0].pop('library_layout')
    data['extensions'] = {'insdc': {'version': '1.0', 'records': [
        {'provider': 'ArrayExpress', 'kind': 'MINiML', 'accession': 'E-MTAB-1', 'metadata': incoming}]}}
    result, calls = convert(data, incoming, mode='curators')
    item = result.items[0]
    assert not calls
    assert item.payload[0]['sample'][0]['sra_run'][0]['library_layout'] == 'SINGLE'
    assert item.completeness == 'partial' and item.validation == 'valid'
    assert any('library_layout' in d.message and 'PAIRED' in d.message for d in item.diagnostics)


def test_curators_does_not_propagate_unbound_or_other_library_path_facts():
    for field, value in [('name', 'SRR2'), ('sample_ref', 'SRS9')]:
        data, incoming = package(), package('MAGE-TAB', 'PAIRED')
        incoming['sample'][0]['sra_run'][0].pop('library_layout')
        steps = incoming['series']['assay_paths'][0]['steps']
        if field == 'name': steps[-1]['name'] = value
        else: steps[0]['sample_ref'] = value
        data['extensions'] = {'insdc': {'version': '1.0', 'records': [
            {'provider': 'ArrayExpress', 'kind': 'MINiML', 'accession': 'E-MTAB-1', 'metadata': incoming}]}}
        result, calls = convert(data, incoming, mode='curators')
        assert not calls and result.status == 'complete'
        assert not any('library_layout' in d.message for d in result.items[0].diagnostics)
