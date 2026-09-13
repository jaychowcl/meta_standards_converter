# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Source references are the only authority for native citation retrieval."""
import xml.etree.ElementTree as ET
import pytest
from meta_standards_converter.sources.archive_support import StudyRecords, StudySeed, publication_ids
from meta_standards_converter.sources.ena import ENASource
from meta_standards_converter.sources.sra import SRASource
from meta_standards_converter.miniml.sra_parser import SRAParser
from tests.test_native_archive_sources import HTTP
from tests.test_native_archive_parsers import fixture_records


def xref(target, pmid='123'):
    return {'Source':'EuropePMC', 'Source Primary Accession':'PMC456',
            'Source Secondary Accession':pmid, 'Source URL':'https://europepmc.org/articles/PMC456',
            'Target':'study', 'Target Primary Accession':target, 'Target Secondary Accession':''}


def test_ena_xrefs_json_first_row_pagination_cache_and_wrong_target():
    def handler(url,p,fmt):
        assert '/json/search' in url and fmt == 'json'
        return [xref('ERP1'), xref('ERP999','999')] if p['offset']==0 else [xref('ERP1','124')]
    http=HTTP(handler);s=ENASource(http=http);s.xref_page_size=2
    records=StudyRecords(StudySeed('ERP1','PRJEB1'))
    s.cross_references('ERP1', records)
    s.cross_references('ERP1', records)
    assert len(http.calls)==2
    assert publication_ids(records)=={'123','124'}
    assert any('target' in issue for issue in records.issues)
    assert all('999' not in str(r) for r in records.linked)


@pytest.mark.parametrize('change', ['dbfrom','input','dbto','linkname','error','missing'])
def test_sra_publication_links_validate_response_identity(change):
    xml='<eLinkResult><LinkSet><DbFrom>sra</DbFrom><IdList><Id>10</Id></IdList><LinkSetDb><DbTo>pubmed</DbTo><LinkName>sra_pubmed</LinkName><Link><Id>123</Id></Link></LinkSetDb></LinkSet></eLinkResult>'
    root=ET.fromstring(xml)
    if change=='dbfrom':root.find('.//DbFrom').text='biosample'
    if change=='input':root.find('.//IdList/Id').text='99'
    if change=='dbto':root.find('.//DbTo').text='pmc'
    if change=='linkname':root.find('.//LinkName').text='sra_pmc'
    if change=='error':ET.SubElement(root,'ERROR').text='bad response'
    if change=='missing':root.remove(root.find('LinkSet'))
    s=SRASource(http=HTTP(lambda *args:root))
    with pytest.raises(ValueError):s.publication_links('sra','pubmed',['10'],'sra_pubmed')


def test_sra_empty_publication_link_response_is_valid():
    root=ET.fromstring('<eLinkResult><LinkSet><DbFrom>sra</DbFrom><IdList><Id>10</Id></IdList></LinkSet></eLinkResult>')
    assert SRASource(http=HTTP(lambda *a:root)).publication_links('sra','pubmed',['10'],'sra_pubmed')=={'10':[]}


def test_sample_citation_is_not_a_study_publication():
    records=fixture_records('sra')
    sample=records.xml[1].find('.//BioSample')
    ET.SubElement(ET.SubElement(sample,'Links'),'Link',type='pubmed',label='paper').text='123'
    records.xml.append(ET.fromstring('<PubmedArticleSet><PubmedArticle><MedlineCitation><PMID>123</PMID><Article><ArticleTitle>Sample paper</ArticleTitle></Article></MedlineCitation></PubmedArticle></PubmedArticleSet>'))
    data=SRAParser().parse(records).to_mapping()
    assert '123' not in data['series'].get('pubmed_id',[])
    assert data['sample'][0]['pubmed_publication'][0]['title']=='Sample paper'
    assert {'type':'PubMed','target':'123'} in data['sample'][0]['relation']
    assert 'Sample paper' not in str(data['extensions'])


def test_doi_only_source_citation_retains_title_without_guessing():
    records=fixture_records('sra')
    project=records.xml[2].find('.//ProjectDescr')
    pub=ET.SubElement(project,'Publication',id='10.1234/explicit')
    ET.SubElement(pub,'DbType').text='eDOI'
    ET.SubElement(pub,'Title').text='Supplied citation'
    data=SRAParser().parse(records).to_mapping()
    assert data['series']['pubmed_publication'][0]['pubmed_id']==''
    assert data['series']['pubmed_publication'][0]['title']=='Supplied citation'
    assert data['series']['pubmed_publication'][0]['doi']=='10.1234/explicit'
    assert 'Supplied citation' not in str(data['extensions'])


def test_identifier_resolution_verifies_doi_and_pmc_namespace():
    from meta_standards_converter.sources.archive_publications import resolve_identifiers
    def handler(url,p,fmt):
        if 'idconv' in url:
            assert p['ids']=='PMC456'
            return {'records':[{'requested-id':'PMC456','pmcid':'PMC456','pmid':'123'}]}
        if 'esearch' in url:
            assert p['db']=='pubmed' and p['term']=='"10.1234/explicit"[AID]'
            return {'esearchresult':{'idlist':['999','124']}}
        return ET.fromstring('<PubmedArticleSet><PubmedArticle><MedlineCitation><PMID>999</PMID></MedlineCitation><PubmedData><ArticleIdList><ArticleId IdType="doi">10.1234/wrong</ArticleId></ArticleIdList></PubmedData></PubmedArticle><PubmedArticle><MedlineCitation><PMID>124</PMID></MedlineCitation><PubmedData><ArticleIdList><ArticleId IdType="doi">10.1234/explicit</ArticleId></ArticleIdList></PubmedData></PubmedArticle></PubmedArticleSet>')
    records=StudyRecords(StudySeed('SRP1','SRP1'))
    records.linked=[{'provider':'sra','kind':'publication_reference','accession':'SRP1','metadata':{'pmcid':'PMC456'}},
                    {'provider':'sra','kind':'publication_reference','accession':'SRP1','metadata':{'doi':'10.1234/explicit'}}]
    resolve_identifiers(records,HTTP(handler),'sra')
    assert publication_ids(records)=={'123','124'}


def test_failed_or_mismatched_identifier_resolution_retains_supplied_reference():
    from meta_standards_converter.sources.archive_publications import resolve_identifiers
    records=StudyRecords(StudySeed('SRP1','SRP1'))
    records.linked=[{'provider':'sra','kind':'publication_reference','accession':'SRP1','metadata':{'pmcid':'PMC456'}}]
    resolve_identifiers(records,HTTP(lambda *a:{'records':[{'requested-id':'PMC999','pmcid':'PMC999','pmid':'123'}]}),'sra')
    assert not publication_ids(records)
    assert records.issues and records.linked[0]['metadata']['pmcid']=='PMC456'


def test_cross_reference_enrichment_link_moves_into_core_before_pruning():
    from meta_standards_converter.miniml.ena_parser import ENAParser
    from meta_standards_converter.metadata.archive_enrichment import linked_accessions
    records=fixture_records('ena')
    records.linked.append({'provider':'ena','kind':'cross_references','accession':records.seed.study,
        'metadata':[{'Source':'ArrayExpress','Source Primary Accession':'E-MTAB-308',
                     'Source URL':'https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-308',
                     'Target':'study','Target Primary Accession':records.seed.primary}]})
    data=ENAParser().parse(records).to_mapping()
    assert {'value':'E-MTAB-308','database':'ArrayExpress'} in data['series']['accession']
    assert 'E-MTAB-308' in linked_accessions(data)
    assert not any(r['kind']=='cross_references' for r in data['extensions']['insdc']['records'])


def test_assembly_and_experiment_publications_keep_scope_and_details():
    records=fixture_records('sra')
    records.linked.extend([
        {'provider':'sra','kind':'assembly','accession':'GCA_123.2','metadata':{'assemblyaccession':'GCA_123.2'}},
        {'provider':'sra','kind':'publication_reference','accession':'GCA_123.2','metadata':{'pubmed_id':'123','title':'Assembly paper'}},
        {'provider':'sra','kind':'publication_reference','accession':'SRX7812918','metadata':{'pubmed_id':'456','title':'Experiment paper'}}])
    data=SRAParser().parse(records).to_mapping()
    assert not data['series'].get('pubmed_id')
    pubs=[r for r in data['series']['relation'] if r.get('publication')]
    assert any(r.get('assembly_ref')=='GCA_123.2' and r['publication']['title']=='Assembly paper' for r in pubs)
    assert any(r.get('experiment_ref')=='SRX7812918' for r in pubs)
    assert not any(r['kind']=='publication_reference' for r in data['extensions']['insdc']['records'])


def test_scoped_native_hydration_preserves_scope_and_uses_one_lookup_per_pmid():
    from meta_standards_converter.metadata.enrichment import MINiMLEnricher
    from meta_standards_converter.miniml import MINiMLCodec
    from tests.test_archive_publications import PubMed, NoRuns
    data=SRAParser().parse(fixture_records('sra')).to_mapping()
    data['sample'][0]['pubmed_publication']=[{'pubmed_id':'123'}]
    data['series'].setdefault('relation',[]).append({'type':'assembly publication','target':'123','assembly_ref':'GCA_123.2','publication':{'pubmed_id':'123'}})
    client=PubMed();service=MINiMLEnricher(pubmed_fetcher=client,insdc_fetcher=NoRuns())
    result=service.enrich(MINiMLCodec().decode(data).package).to_mapping()
    assert client.calls==['123']
    assert result['sample'][0]['pubmed_publication'][0]['title']=='Fetched title'
    assert not result['series'].get('pubmed_id')
    assert next(r['publication'] for r in result['series']['relation'] if r.get('publication'))['title']=='Fetched title'


def test_sra_link_retrieval_binds_uid_to_verified_experiment_project_and_sample():
    records=fixture_records('sra')
    biosample=records.xml[1].find('.//BioSample')
    project=records.xml[2].find('.//ProjectID/ArchiveID')
    def handler(url,p,fmt):
        if 'esummary' in url:
            return {'result':{'10':{'uid':'10','expxml':'<Experiment acc="SRX7812918"/>'}}}
        root=ET.Element('eLinkResult')
        for uid in p['id']:
            group=ET.SubElement(root,'LinkSet');ET.SubElement(group,'DbFrom').text=p['dbfrom']
            ET.SubElement(ET.SubElement(group,'IdList'),'Id').text=uid
            if p['db']=='pubmed':
                link=ET.SubElement(group,'LinkSetDb');ET.SubElement(link,'DbTo').text='pubmed'
                ET.SubElement(link,'LinkName').text=p['linkname']
                ET.SubElement(ET.SubElement(link,'Link'),'Id').text={'sra':'123','bioproject':'124','biosample':'125'}[p['dbfrom']]
        return root
    http=HTTP(handler);source=SRASource(http=http)
    source.linked_publications(records,['10'])
    refs={r['metadata']['pubmed_id']:r['accession'] for r in records.linked if r['kind']=='publication_reference'}
    assert refs=={'123':'SRX7812918','124':project.get('accession'),'125':biosample.get('accession')}
    data=SRAParser().parse(records).to_mapping()
    assert data['series']['pubmed_id']==['124']
    assert data['sample'][0]['pubmed_id']==['125']
    assert any(r.get('experiment_ref')=='SRX7812918' and r['target']=='123' for r in data['series']['relation'])
