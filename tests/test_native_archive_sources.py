# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Discovery must expand identities, not silently import unrelated scopes."""
import xml.etree.ElementTree as ET
from meta_standards_converter.sources.sra import SRASource
from meta_standards_converter.sources.ena import ENASource
from meta_standards_converter.sources.archive_support import StudySeed


class HTTP:
    def __init__(self, handler):
        self.handler = handler
        self.calls = []
    def get(self, url, params=None, fmt='xml'):
        self.calls.append((url, params or {}, fmt))
        return self.handler(url, params or {}, fmt)


def test_sra_run_resolves_study_and_fetch_uses_study_query():
    xml = ET.fromstring('<EXPERIMENT_PACKAGE_SET><EXPERIMENT_PACKAGE><STUDY accession="SRP1"/><EXPERIMENT accession="SRX1"><STUDY_REF accession="SRP1"/></EXPERIMENT><RUN_SET><RUN accession="SRR1"/></RUN_SET></EXPERIMENT_PACKAGE></EXPERIMENT_PACKAGE_SET>')
    def handler(url, p, fmt):
        if 'esearch' in url:
            return {'esearchresult': {'count': '1', 'idlist': ['10'], 'querykey': '1', 'webenv': 'test'}}
        if 'efetch' in url:
            return xml
        return ET.fromstring('<eLinkResult/>')
    http = HTTP(handler)
    source = SRASource(http=http)
    result = source.resolve('SRR1')
    assert result.studies == [StudySeed('SRP1', 'SRP1')]
    source.fetch(result.studies[0])
    terms = [p['term'] for _, p, _ in http.calls if 'term' in p]
    assert terms == ['SRR1[ACCN]', 'SRP1[ACCN]']
    assert all('eutils.ncbi.nlm.nih.gov' in u for u, _, _ in http.calls)


def test_sra_search_walks_all_uid_pages():
    def handler(url, p, fmt):
        start = int(p.get('retstart', 0))
        return {'esearchresult': {'count': '3', 'idlist': ['1', '2'] if start == 0 else ['3']}}
    source = SRASource(http=HTTP(handler))
    source.page_size = 2
    assert source.search('SRP1[ACCN]')[0] == ['1', '2', '3']


def test_ena_sample_can_resolve_multiple_read_studies():
    http = HTTP(lambda u, p, f: [
        {'study_accession': 'PRJEB1', 'secondary_study_accession': 'ERP1'},
        {'study_accession': 'PRJEB2', 'secondary_study_accession': 'ERP2'}])
    resolution = ENASource(http=http).resolve('SAMEA123')
    assert [s.primary for s in resolution.studies] == ['PRJEB1', 'PRJEB2']
    assert 'sample_accession="SAMEA123"' in http.calls[0][1]['query']
    assert http.calls[0][1]['limit'] == 0


def test_ena_umbrella_follows_children_not_peer_projects():
    def handler(url, p, fmt):
        if '/xml/PRJEB1' in url:
            return ET.fromstring('<PROJECT_SET><PROJECT accession="PRJEB1"><UMBRELLA_PROJECT/><RELATED_PROJECTS><RELATED_PROJECT><CHILD_PROJECT accession="PRJEB2"/></RELATED_PROJECT><RELATED_PROJECT><RELATED_PROJECT accession="PRJEB9"/></RELATED_PROJECT></RELATED_PROJECTS></PROJECT></PROJECT_SET>')
        if '/xml/PRJEB2' in url:
            return ET.fromstring('<PROJECT_SET><PROJECT accession="PRJEB2"/></PROJECT_SET>')
        return [{'study_accession': 'PRJEB2', 'secondary_study_accession': 'ERP2'}]
    http = HTTP(handler)
    result = ENASource(http=http).resolve('PRJEB1')
    assert result.studies == [StudySeed('ERP2', 'PRJEB2')]
    assert all('PRJEB9' not in str(call) for call in http.calls)


def test_ena_partial_xml_retains_identified_study():
    def handler(url, p, fmt):
        if '/xml/' in url:
            raise OSError('unavailable')
        if 'returnFields' in url:
            return 'columnId\tdescription\ttype\nrun_accession\tRun\ttext\n'
        if '/count' in url:
            return [{'count': '1'}]
        return []
    records = ENASource(http=HTTP(handler)).fetch(StudySeed('ERP1', 'PRJEB1'))
    assert records.seed.primary == 'PRJEB1'
    assert records.issues


def test_sra_history_paging_uses_frozen_query():
    def handler(url, p, fmt):
        return {'esearchresult': {'count': '2', 'idlist': ['1' if p['retstart'] == 0 else '2'], 'webenv': 'frozen', 'querykey': '7'}}
    http = HTTP(handler); source = SRASource(http=http); source.page_size = 1
    source.search('SRP1[ACCN]')
    assert http.calls[1][1]['WebEnv'] == 'frozen'
    assert http.calls[1][1]['term'] == '#7'


def test_ena_xml_batch_missing_identifiers_are_reconciled():
    def handler(url, p, fmt):
        if '/xml/' in url:
            return ET.fromstring('<RUN_SET><RUN accession="ERR9"/></RUN_SET>')
        if 'returnFields' in url:
            return 'columnId\nrun_accession\n'
        if '/count' in url:
            return {'count': '1'}
        if p.get('result') == 'read_run':
            return [{'run_accession': 'ERR1'}]
        return []
    records = ENASource(http=HTTP(handler)).fetch(StudySeed('ERP1', 'PRJEB1'))
    assert any('ERR1' in issue for issue in records.issues)


def test_sra_sample_expansion_and_umbrella_children_are_explicit():
    def handler(url, p, fmt):
        if 'esearch' in url:
            return {'esearchresult': {'count': '1', 'idlist': ['1']}}
        if 'elink' in url:
            assert p['linkname'] == 'bioproject_bioproject_u2d'
            return ET.fromstring('<eLinkResult><LinkSet><LinkSetDb><Link><Id>2</Id></Link></LinkSetDb></LinkSet></eLinkResult>')
        if p['db'] == 'bioproject':
            if p['id'] == '1':
                return ET.fromstring('<RecordSet><Project><ProjectID><ArchiveID id="1" accession="PRJNA1"/></ProjectID><ProjectTypeTopAdmin/></Project></RecordSet>')
            return ET.fromstring('<RecordSet><Project><ProjectID><ArchiveID id="2" accession="PRJNA2"/></ProjectID></Project></RecordSet>')
        return ET.fromstring('<EXPERIMENT_PACKAGE_SET><EXPERIMENT_PACKAGE><EXPERIMENT><STUDY_REF accession="SRP2"/></EXPERIMENT></EXPERIMENT_PACKAGE><EXPERIMENT_PACKAGE><EXPERIMENT><STUDY_REF accession="SRP3"/></EXPERIMENT></EXPERIMENT_PACKAGE></EXPERIMENT_PACKAGE_SET>')
    source = SRASource(http=HTTP(handler))
    assert [s.study for s in source.resolve('SAMN1').studies] == ['SRP2', 'SRP3']
    source.resolve('PRJNA1')
    assert any(p.get('term') == 'PRJNA2[GPRJ]' for _, p, _ in source.http.calls)


def test_transport_enforces_aggregate_response_budget():
    import pytest
    from meta_standards_converter.sources.archive_support import ArchiveHTTP
    from meta_standards_converter.runtime_contracts import get_resource_profile
    class Response:
        headers = {}
        def raise_for_status(self): pass
        def close(self): pass
        def iter_content(self, **kwargs): yield b'<X/>'
    class Requester:
        def get(self, *args, **kwargs): return Response()
    profile = get_resource_profile('standard', overrides={'max_aggregate_download_bytes': 6})
    http = ArchiveHTTP('ena_portal', requester=Requester(), resource_profile=profile)
    http.get('https://www.ebi.ac.uk/ena/browser/api/xml/ERP1')
    with pytest.raises(Exception, match='aggregate|budget|limit'):
        http.get('https://www.ebi.ac.uk/ena/browser/api/xml/ERP2')


def test_sra_duplicate_packages_cannot_hide_missing_experiments():
    def handler(url, p, fmt):
        if 'esearch' in url:
            return {'esearchresult': {'count': '2', 'idlist': ['1', '2']}}
        if 'efetch' in url:
            return ET.fromstring('<EXPERIMENT_PACKAGE_SET><EXPERIMENT_PACKAGE><STUDY accession="SRP1"/><EXPERIMENT accession="SRX1"><STUDY_REF accession="SRP1"/></EXPERIMENT></EXPERIMENT_PACKAGE><EXPERIMENT_PACKAGE><EXPERIMENT accession="SRX1"><STUDY_REF accession="SRP1"/></EXPERIMENT></EXPERIMENT_PACKAGE></EXPERIMENT_PACKAGE_SET>')
        return ET.fromstring('<eLinkResult/>')
    records = SRASource(http=HTTP(handler)).fetch(StudySeed('SRP1', 'SRP1'))
    assert any('inventory' in issue for issue in records.issues)


def test_sra_assembly_links_retain_full_versioned_summaries():
    from meta_standards_converter.miniml.sra_parser import SRAParser
    def handler(url, p, fmt):
        if 'esearch' in url:
            return {'esearchresult': {'count': '1', 'idlist': ['1']}}
        if 'elink' in url:
            return ET.fromstring('<eLinkResult><LinkSet><LinkSetDb><Link><Id>900</Id></Link></LinkSetDb></LinkSet></eLinkResult>')
        if 'esummary' in url:
            assert p['db'] == 'assembly' and p['report'] == 'full'
            return {'result': {'900': {'assemblyaccession': 'GCA_000001405.29', 'synonym': {'refseq': 'GCF_000001405.40'}, 'ftppath_genbank': 'ftp://example.org/assembly-directory', 'meta': {'nested': ['retained']}}}}
        return ET.fromstring('<EXPERIMENT_PACKAGE_SET><EXPERIMENT_PACKAGE><STUDY accession="SRP1"/><EXPERIMENT accession="SRX1"><STUDY_REF accession="SRP1"/></EXPERIMENT></EXPERIMENT_PACKAGE></EXPERIMENT_PACKAGE_SET>')
    source = SRASource(http=HTTP(handler))
    records = source.fetch(StudySeed('SRP1', 'SRP1'))
    data = SRAParser().parse(records).to_mapping()
    assembly = next(r for r in data['extensions']['insdc']['records'] if r['kind'] == 'assembly')
    assert assembly['accession'] == 'GCA_000001405.29'
    assert assembly['metadata']['synonym']['refseq'] == 'GCF_000001405.40'
    assert assembly['metadata']['meta'] == {'nested': ['retained']}
    assert not data['series'].get('supplementary_data')
    assert any(r['type'] == 'assembly directory' for r in data['series']['relation'])
    assert all('example.org' not in url for url, _, _ in source.http.calls)
