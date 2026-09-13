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
