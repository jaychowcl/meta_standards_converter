# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
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
    assert not paths  # Assembly documentation is a supporting link, not a measurement node.
    link=next(f for f in data['sample'][0]['supplementary_data'] if f['value']=='ftp://example.org/stats.txt')
    assert link['assembly_accession']=='GCF_1.2' and link['type']=='assembly report'
    from tests.test_protocol_export import render
    table=next(r[1] for r in render(data) if r[0]=='SDRF File')
    assert 'ftp://example.org/stats.txt' in str(table)


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
        if p.get('result')=='analysis': return []
        return ET.fromstring(f'<ANALYSIS_SET><ANALYSIS accession="ERZ1"><STUDY_REF accession="ERP999"/><SAMPLE_REF accession="{sample.get("accession")}"/></ANALYSIS></ANALYSIS_SET>')
    source=ENASource(http=HTTP(handler));source.fetch_associated_analyses(records)
    assert len(calls)==2 and calls[0].endswith('/ERZ1') and calls[1].endswith('/search')
    assert any(n.tag=='ANALYSIS' for r in records.xml for n in r)
    assert not records.issues


def test_cross_reference_lists_split_only_verified_accession_literals():
    from meta_standards_converter.miniml.insdc_support import relations
    node=ET.fromstring('<SAMPLE><XREF_LINK><DB>ENA-STUDY</DB><ID>DRP1,ERP2</ID></XREF_LINK><XREF_LINK><DB>description</DB><ID>one,two</ID></XREF_LINK></SAMPLE>')
    assert relations(node)==[{'type':'ENA-STUDY','target':'DRP1'},{'type':'ENA-STUDY','target':'ERP2'},{'type':'description','target':'one,two'}]


def test_ena_range_links_are_resolved_before_individual_relations():
    from meta_standards_converter.sources.ena import ENASource
    from tests.test_native_archive_sources import HTTP
    records=fixture_records('ena');sample=records.xml[1].find('SAMPLE')
    link=ET.SubElement(sample,'XREF_LINK');ET.SubElement(link,'DB').text='ENA-RUN';ET.SubElement(link,'ID').text='ERR000001-ERR000002,ERR000004'
    def handler(url,p,fmt):
        ids=url.rsplit('/',1)[-1].split(',')
        return ET.fromstring('<RUN_SET>'+''.join(f'<RUN accession="{a}"/>' for a in ids)+'</RUN_SET>')
    source=ENASource(http=HTTP(handler));source.fetch_reference_ranges(records)
    data=ENAParser().parse(records).to_mapping()
    refs=data['sample'][0]['relation']
    assert all({'type':'ENA-RUN','target':a} in refs for a in ['ERR000001','ERR000002','ERR000004'])
    assert not any('-ERR' in r.get('target','') for r in refs)
    assert len(data['sample'][0]['sra_run'])==1


def test_associated_analysis_portal_files_project_without_inventing_urls():
    from meta_standards_converter.sources.ena import ENASource
    from tests.test_native_archive_sources import HTTP
    records=fixture_records('ena');sample=records.xml[1].find('SAMPLE');acc=sample.get('accession')
    link=ET.SubElement(sample,'XREF_LINK');ET.SubElement(link,'DB').text='ENA-ANALYSIS';ET.SubElement(link,'ID').text='ERZ1'
    def handler(url,p,fmt):
        if p.get('result')=='analysis':
            assert p['query']=='analysis_accession="ERZ1"'
            return [{'analysis_accession':'ERZ1','sample_accession':acc,'run_accession':'SRR11192680','submitted_ftp':'host/vol1/analysis/ERZ1/result.fa.gz','submitted_format':'FASTA','submitted_md5':'a'*32}]
        return ET.fromstring(f'<ANALYSIS_SET><ANALYSIS accession="ERZ1"><STUDY_REF accession="ERP999"/><SAMPLE_REF accession="{acc}"/><RUN_REF accession="SRR11192680"/><FILES><FILE filename="ERZ1/result.fa.gz" filetype="fasta" checksum_method="MD5" checksum="'+ 'a'*32 +'"/></FILES></ANALYSIS></ANALYSIS_SET>')
    ENASource(http=HTTP(handler)).fetch_associated_analyses(records)
    data=ENAParser().parse(records).to_mapping()
    assert any(f['value']=='ftp://host/vol1/analysis/ERZ1/result.fa.gz' for f in data['sample'][0].get('supplementary_data',[]))
    from tests.test_archive_residuals import nodes
    assert not any(n.get('tag')=='FILE' for n in nodes(data['extensions']))
    row=next((r['metadata'] for r in data['extensions']['insdc']['records'] if r['kind']=='analysis'),{})
    assert 'submitted_ftp' not in row
