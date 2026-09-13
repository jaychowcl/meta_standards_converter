# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from copy import deepcopy
from meta_standards_converter.miniml import MINiMLCodec
from meta_standards_converter.miniml.sra_parser import SRAParser
from meta_standards_converter.metadata.archive_enrichment import merge_archive_metadata, LinkedArchiveEnricher
from tests.test_native_archive_parsers import fixture_records


def native():
    return SRAParser().parse(fixture_records('sra'))


def linked(package, accession, title):
    data = package.to_mapping()
    data['source'] = {'format': 'MAGE-TAB' if accession.startswith('E-') else 'GEO'}
    data['series']['iid'] = accession
    data['series']['accession'].append({'value': accession, 'database': 'ArrayExpress' if accession.startswith('E-') else 'GEO'})
    data['series']['title'] = title
    data['sample'][0]['title'] = title
    old_iid = data['sample'][0]['iid']
    data['sample'][0]['iid'] = 'GSM1'
    for path in data['series'].get('assay_paths', []):
        for step in path['steps']:
            if step.get('sample_ref') == old_iid: step['sample_ref'] = 'GSM1'
    data['sample'][0]['channel'][0]['characteristics'] = [{'name': 'host', 'value': title, 'term_accession_number': 'ONT:1', 'term_source_ref': 'ONT'}]
    return MINiMLCodec().decode(data).package


def test_coherent_priority_identity_and_membership():
    package = native()
    other = linked(package, 'GSE1', 'enriched')
    merged, issues = merge_archive_metadata(package, other, prefer=True)
    data = merged.to_mapping()
    assert data['series']['title'] == 'enriched'
    assert data['series']['iid'] == package.series.iid
    assert data['sample'][0]['iid'] == package.samples[0].iid
    assert data['sample'][0]['channel'][0]['characteristics'][0]['value'] != 'enriched' or any(a['name'] == 'host' for a in data['sample'][0]['channel'][0]['characteristics'])
    host = next(a for a in data['sample'][0]['channel'][0]['characteristics'] if a['name'] == 'host')
    assert host['value'] == 'enriched'
    assert any(s.get('characteristics') and any(a.get('value') == 'enriched' for a in s['characteristics']) for p in data['series']['assay_paths'] for s in p['steps'])
    assert data['extensions']['insdc']['records'][:1] == package.to_mapping()['extensions']['insdc']['records'][:1]
    assert not issues


def test_missing_values_do_not_displace_informative_and_titles_do_not_join():
    package = native()
    other = linked(package, 'GSE1', 'not provided').to_mapping()
    other['sample'][0]['accession'] = [{'value': 'GSM9', 'database': 'GEO'}]
    other['sample'][0]['sra_accession'] = []
    other['sample'][0]['relation'] = []
    other['sample'][0]['sra_run'] = []
    merged, issues = merge_archive_metadata(package, MINiMLCodec().decode(other).package, prefer=True)
    assert merged.series.title == package.series.title
    assert merged.samples[0].title == package.samples[0].title
    assert len(merged.samples) == 1
    assert any('unmatched' in i for i in issues)


def test_ae_over_geo_and_unavailable_enrichment():
    package = native().to_mapping()
    package['series']['relation'] = [{'type': 'GEO', 'target': 'GSE1'}, {'type': 'ArrayExpress', 'target': 'E-MTAB-1'}]
    package = MINiMLCodec().decode(package).package
    class Converter:
        def convert(self, accession, **kwargs):
            return [linked(package, accession, 'AE' if accession.startswith('E-') else 'GEO')]
    enriched, issues = LinkedArchiveEnricher(geo_converter=Converter(), ae_converter=Converter()).enrich(package)
    assert enriched.series.title == 'AE'
    class Broken:
        def convert(self, *args, **kwargs):
            raise OSError('private detail')
    enriched, issues = LinkedArchiveEnricher(geo_converter=Broken(), ae_converter=Broken()).enrich(package)
    assert enriched.series.iid == package.series.iid
    assert issues and 'private detail' not in str(issues)


def test_ambiguous_sample_join_and_unrelated_peer_are_rejected():
    package = native()
    other = linked(package, 'GSE1', 'enriched').to_mapping()
    other['sample'].append(deepcopy(other['sample'][0]))
    other['sample'][1]['iid'] = 'GSM2'
    merged, issues = merge_archive_metadata(package, MINiMLCodec().decode(other).package, prefer=True)
    assert merged.samples[0].title == package.samples[0].title
    assert any('ambiguous' in i for i in issues)
    other['series']['accession'] = [{'value': 'SRP999', 'database': 'SRA'}]
    other['series']['iid'] = 'SRP999'
    merged, issues = merge_archive_metadata(package, MINiMLCodec().decode(other).package)
    assert merged == package
    assert issues


def test_enrichment_files_and_sample_protocols_reach_native_assay_paths():
    package = native()
    other = linked(package, 'GSE1', 'GEO').to_mapping()
    other['sample'][0]['supplementary_data'] = [{'value': 'https://example.org/counts.tsv', 'type': 'tsv'}]
    other['sample'][0]['channel'][0]['extract_protocol'] = 'Explicit extraction method'
    merged, issues = merge_archive_metadata(package, MINiMLCodec().decode(other).package, prefer=True)
    data = merged.to_mapping()
    assert any(s.get('link', {}).get('value') == 'https://example.org/counts.tsv' for p in data['series']['assay_paths'] for s in p['steps'])
    assert any(p.get('description') == 'Explicit extraction method' for p in data['series']['protocols'])
    assert any(s.get('protocol_ref', '').endswith('extract_protocol') for p in data['series']['assay_paths'] for s in p['steps'])


def test_peer_only_runs_are_added_on_exact_study_and_sample_match():
    package = native()
    other = package.to_mapping()
    other['sample'][0]['sra_run'][0]['run'] = 'SRR99'
    for path in other['series']['assay_paths']:
        for step in path['steps']:
            if step['kind'] == 'scan':
                step['name'] = 'SRR99'
                for comment in step.get('comments', []):
                    if comment['name'] in ('ENA_RUN', 'SRA_RUN'): comment['value'] = 'SRR99'
    merged, issues = merge_archive_metadata(package, MINiMLCodec().decode(other).package)
    assert {r['run'] for r in merged.to_mapping()['sample'][0]['sra_run']} == {'SRR11192680', 'SRR99'}
    assert any(s.get('name') == 'SRR99' for p in merged.to_mapping()['series']['assay_paths'] for s in p['steps'])


def test_processing_protocol_priority_is_coherent_and_scoped():
    package = native()
    geo = linked(package, 'GSE1', 'GEO').to_mapping()
    geo['sample'][0]['data_processing'] = 'GEO processing'
    merged, _ = merge_archive_metadata(package, MINiMLCodec().decode(geo).package, prefer=True)
    ae = linked(package, 'E-MTAB-1', 'AE').to_mapping()
    ae['sample'][0]['data_processing'] = 'AE processing'
    merged, _ = merge_archive_metadata(merged, MINiMLCodec().decode(ae).package, prefer=True)
    data = merged.to_mapping()
    names = {p['name']: p for p in data['series']['protocols']}
    for path in data['series']['assay_paths']:
        applied = [names[s['protocol_ref']].get('description') for s in path['steps'] if s.get('protocol_ref') in names]
        assert 'AE processing' in applied
        assert 'GEO processing' not in applied


def test_experiment_protocol_enrichment_does_not_leak_to_other_experiments():
    package = native().to_mapping()
    original = deepcopy(package['series']['assay_paths'][0])
    for step in original['steps']:
        if step['kind'] == 'assay': step['name'] = 'SRX99'
        if step['kind'] == 'scan': step['name'] = 'SRR99'; step['comments'] = []
    package['series']['assay_paths'].append(original)
    package = MINiMLCodec().decode(package).package
    extra = linked(package, 'E-MTAB-1', 'AE').to_mapping()
    path = deepcopy(extra['series']['assay_paths'][0])
    for step in path['steps']:
        if step.get('sample_ref'): step['sample_ref'] = 'GSM1'
    path['steps'].insert(2, {'kind': 'protocol_application', 'protocol_ref': 'specific'})
    extra['series']['assay_paths'] = [path]
    extra['series']['protocols'] = [{'name': 'specific', 'description': 'Specific experiment preparation', 'type': {'value': 'library construction protocol'}}]
    merged, issues = merge_archive_metadata(package, MINiMLCodec().decode(extra).package, prefer=True)
    for path in merged.to_mapping()['series']['assay_paths']:
        specific = any(s.get('protocol_ref', '').endswith(':specific') for s in path['steps'])
        if any(s.get('name') == 'SRX99' for s in path['steps']): assert not specific
        else: assert specific
