import xml.etree.ElementTree as ET
from copy import deepcopy
from meta_standards_converter.miniml import MINiMLCodec
from meta_standards_converter.miniml.ena_parser import ENAParser
from meta_standards_converter.miniml.sra_parser import SRAParser
from meta_standards_converter.magetab.idf import IDFConstructor
from meta_standards_converter.metadata.archive_enrichment import merge_archive_metadata
from tests.test_native_archive_parsers import fixture_records
from tests.test_archive_fidelity_enrichment import workflow, enriched_native


def test_dates_keep_provider_and_submission_is_not_experiment():
    data = enriched_native().to_mapping()
    data['series']['status'] = [{'database':'SRA', 'release_date':'2020-01-01'}]
    extra = workflow().to_mapping()
    extra['series']['status'] = [{'release_date':'2019-01-01','submission_date':'2018-01-01'}]
    merged,_ = merge_archive_metadata(MINiMLCodec().decode(data).package,MINiMLCodec().decode(extra).package,prefer=True)
    statuses = merged.to_mapping()['series']['status']
    assert {s['database'] for s in statuses} == {'SRA','ArrayExpress'}
    rows = {r[0]:r[1:] for r in IDFConstructor()._idf_dates(merged.to_mapping())}
    assert rows['Public Release Date'] == ['2020-01-01']
    assert rows['Comment[ArrayExpressReleaseDate]'] == ['2019-01-01']
    assert rows['Date of Experiment'] == [] or rows['Date of Experiment'] == ['']


def test_downstream_factors_do_not_flatten_and_source_occurrences_survive():
    data = workflow().to_mapping()
    assert not any(c['name']=='AGE' for c in data['sample'][0]['channel'][0]['characteristics'])
    records=fixture_records('ena')
    node=records.xml[1].find('SAMPLE')
    attrs=node.find('SAMPLE_ATTRIBUTES')
    for _ in range(2):
        attr=ET.SubElement(attrs,'SAMPLE_ATTRIBUTE')
        ET.SubElement(attr,'TAG').text='organism'
        ET.SubElement(attr,'VALUE').text='human gut metagenome'
    data=ENAParser().parse(records).to_mapping()
    chars=data['sample'][0]['channel'][0]['characteristics']
    expected=sum(c['name']=='organism' for c in chars)+1
    assert all(sum(c['name']=='organism' for c in p['steps'][0]['characteristics'])==expected for p in data['series']['assay_paths'])


def test_assembly_description_is_not_a_protocol_and_sample_links_are_projected():
    records=fixture_records('ena');sample=records.xml[1].find('SAMPLE');acc=sample.get('accession')
    records.xml.append(ET.fromstring(f'<ASSEMBLY_SET><ASSEMBLY accession="GCA_1.2"><DESCRIPTION>Organism isolated from soil</DESCRIPTION><SAMPLE_REF accession="{acc}"/></ASSEMBLY></ASSEMBLY_SET>'))
    data=ENAParser().parse(records).to_mapping()
    assert not any(p['name']=='GCA_1.2:analysis' for p in data['series'].get('protocols', []))
    assert {'type':'assembly','target':'GCA_1.2'} in data['sample'][0]['relation']


def test_ncbi_assembly_reports_and_synonyms_are_sample_scoped():
    records=fixture_records('sra');sample=records.xml[0].find('.//SAMPLE')
    bio=next(n.text for n in sample.findall('IDENTIFIERS/*') if (n.text or '').startswith('SAM'))
    records.linked.append({'provider':'sra','kind':'assembly','accession':'GCF_1.2','metadata':{'assemblyaccession':'GCF_1.2','biosampleaccn':bio,'synonym':{'genbank':'GCA_1.2','refseq':'GCF_1.2','similarity':'identical'},'ftppath_stats_rpt':'ftp://example.org/stats.txt'}})
    data=SRAParser().parse(records).to_mapping()
    assert {'type':'assembly','target':'GCA_1.2'} in data['series']['relation']
    assert {'type':'assembly','target':'GCF_1.2'} in data['sample'][0]['relation']
    paths=[p for p in data['series']['assay_paths'] if any(s.get('link',{}).get('value')=='ftp://example.org/stats.txt' for s in p['steps'])]
    assert paths and all(not any(s['kind'] in ('assay','scan') for s in p['steps']) for p in paths)


def test_analysis_reference_fetch_is_bounded_and_validates_sample_binding():
    from meta_standards_converter.sources.ena import ENASource
    from tests.test_native_archive_sources import HTTP
    records=fixture_records('ena')
    sample=records.xml[1].find('SAMPLE')
    links=ET.SubElement(sample,'SAMPLE_LINKS');link=ET.SubElement(ET.SubElement(links,'SAMPLE_LINK'),'XREF_LINK')
    ET.SubElement(link,'DB').text='ENA-ANALYSIS';ET.SubElement(link,'ID').text='ERZ1'
    calls=[]
    def handler(url,p,fmt):
        calls.append(url)
        return ET.fromstring(f'<ANALYSIS_SET><ANALYSIS accession="ERZ1"><STUDY_REF accession="ERP999"/><SAMPLE_REF accession="{sample.get("accession")}"/></ANALYSIS></ANALYSIS_SET>')
    source=ENASource(http=HTTP(handler));source.fetch_associated_analyses(records)
    assert len(calls)==1 and calls[0].endswith('/ERZ1')
    assert any(n.tag=='ANALYSIS' for r in records.xml for n in r)
    assert not records.issues


def test_cross_reference_lists_split_only_verified_accession_literals():
    from meta_standards_converter.miniml.insdc_support import relations
    node=ET.fromstring('<SAMPLE><XREF_LINK><DB>ENA-STUDY</DB><ID>DRP1,ERP2</ID></XREF_LINK><XREF_LINK><DB>description</DB><ID>one,two</ID></XREF_LINK></SAMPLE>')
    assert relations(node)==[{'type':'ENA-STUDY','target':'DRP1'},{'type':'ENA-STUDY','target':'ERP2'},{'type':'description','target':'one,two'}]
