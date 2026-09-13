"""Regressions from the six-study archive fidelity evaluation."""
import xml.etree.ElementTree as ET
import pytest
from meta_standards_converter.sources.sra import SRASource
from meta_standards_converter.sources.ena import ENASource
from meta_standards_converter.sources.archive_support import StudySeed, Resolution, accession_kind
from meta_standards_converter.miniml.insdc_support import files_from_ena, accessions
from tests.test_native_archive_sources import HTTP

@pytest.mark.parametrize('accession,uid', [('PRJEB43529','713744'), ('PRJEB2208','203921'), ('PRJDA38027','38027')])
def test_project_fetch_resolves_exact_uid_and_validates(accession, uid):
    def handler(url,p,fmt):
        if 'esearch' in url:
            assert p['db']=='bioproject' and p['term']==accession+'[PRJA]'
            return {'esearchresult': {'count':'1','idlist':[uid]}}
        assert p['id']==uid
        return ET.fromstring(f'<RecordSet><DocumentSummary uid="{uid}"><Project><ProjectID><ArchiveID accession="{accession}" id="{uid}"/></ProjectID></Project></DocumentSummary></RecordSet>')
    result=Resolution()
    root=SRASource(http=HTTP(handler)).project_xml(accession,result)
    assert root is not None and not result.issues

@pytest.mark.parametrize('body', ['<DocumentSummary uid="713744"><error>not public</error></DocumentSummary>', '<DocumentSummary uid="713744"><Project><ProjectID><ArchiveID accession="PRJNA43529" id="713744"/></ProjectID></Project></DocumentSummary>', '<DocumentSummary uid="713744"/>'])
def test_bad_projects_rejected(body):
    http=HTTP(lambda u,p,f: {'esearchresult':{'count':'1','idlist':['713744']}} if 'esearch' in u else ET.fromstring('<RecordSet>'+body+'</RecordSet>'))
    result=Resolution()
    assert SRASource(http=http).project_xml('PRJEB43529',result) is None
    assert result.issues


def test_legacy_project_is_recognized_and_retained():
    assert accession_kind('PRJDA38027')==('PRJDA38027','project')
    node=ET.fromstring('<STUDY><IDENTIFIERS><EXTERNAL_ID>PRJDA38027</EXTERNAL_ID></IDENTIFIERS></STUDY>')
    assert {'value':'PRJDA38027','database':'BioProject'} in accessions(node,'DRP000001')

@pytest.mark.parametrize('path', ['ftp.example.org/33308_4#4.cram','ftp://ftp.example.org/33308_4%234.cram'])
def test_portal_ftp_paths_preserve_reserved_filename(path):
    file=files_from_ena({'submitted_ftp':path,'submitted_format':'CRAM'})[0]
    assert file['filename']=='33308_4#4.cram'
    assert file['uri']=='ftp://ftp.example.org/33308_4%234.cram'

@pytest.mark.parametrize('requested,version,missing', [('GCA_000209795','GCA_000209795.2',False),('GCA_000209795.1','GCA_000209795.2',True)])
def test_assembly_version_reconciliation(requested,version,missing):
    def handler(url,p,fmt):
        if '/xml/' in url:
            if 'GCA_' in url: return ET.fromstring(f'<ASSEMBLY_SET><ASSEMBLY accession="{version}"/></ASSEMBLY_SET>')
            acc=url.rsplit('/',1)[-1]
            return ET.fromstring(f'<STUDY_SET><STUDY accession="{acc}"/></STUDY_SET>')
        if 'returnFields' in url: return 'columnId\nassembly_accession\n'
        if '/count' in url: return {'count':'0'}
        if p.get('result')=='assembly': return [{'assembly_accession':requested}]
        return []
    records=ENASource(http=HTTP(handler)).fetch(StudySeed('DRP000001','PRJDA38027'))
    assert any('missing XML record '+requested in x for x in records.issues)==missing
