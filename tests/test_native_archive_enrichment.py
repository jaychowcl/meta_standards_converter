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
    data['sample'][0]['iid'] = 'GSM1'
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
