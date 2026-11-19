# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from meta_standards_converter.miniml.archive_residuals import Projection
from meta_standards_converter.miniml.insdc_support import relations, tree
import xml.etree.ElementTree as ET


from meta_standards_converter.miniml.ena_parser import ENAParser

from meta_standards_converter.miniml.sra_parser import SRAParser


from tests.support.archive import fixture_records


def projection(xml, accepted):
    data = {'series': {'iid':'ERP1', 'accession':[{'value':'ERP1','database':'ENA'}], 'relation':accepted}}
    return Projection(data).xml(tree(ET.fromstring(xml)), 'STUDY', 'ERP1', 'ena')


def test_url_query_commas_remain_opaque_and_mapped_reference_disappears():
    literal = 'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=ERP1&amp;fields=run_accession,fastq_ftp'
    xml = f'<XREF_LINK><DB>ENA-FASTQ-FILES</DB><ID>{literal}</ID></XREF_LINK>'
    accepted = relations(ET.fromstring('<STUDY>'+xml+'</STUDY>'))
    assert len(accepted) == 1
    assert projection(xml, accepted) is None


def test_mapped_target_keeps_unknown_reference_siblings_and_attributes():
    xml = '<XREF_LINK qualifier="special"><DB version="1">ENA-RUN</DB><ID label="curated">ERR1</ID><NOTE>retain</NOTE></XREF_LINK>'
    residual = projection(xml, [{'type':'ENA-RUN','target':'ERR1'}])
    assert residual['attributes'] == {'qualifier':'special'}
    assert residual['children'] == [
        {'tag':'DB','attributes':{'version':'1'}}, {'tag':'ID','attributes':{'label':'curated'}},
        {'tag':'NOTE','text':'retain'}]


def test_partial_list_prunes_only_accepted_members_and_unknown_names_are_opaque():
    xml = '<XREF_LINK><DB>ENA-RUN</DB><ID>ERR1,ERR2</ID></XREF_LINK>'
    residual = projection(xml, [{'type':'ENA-RUN','target':'ERR1'}])
    assert residual['children'][1]['text'] == 'ERR2'
    xml = '<XREF_LINK><DB>ENA-CUSTOM</DB><ID>a,b</ID></XREF_LINK>'
    assert projection(xml, [{'type':'ENA-CUSTOM','target':'a'}])['children'][1]['text'] == 'a,b'


def test_malformed_url_and_unverified_ranges_remain_literal():
    from meta_standards_converter.miniml.reference_targets import parse_reference_targets
    assert parse_reference_targets('external','http://[unclosed') == ['http://[unclosed']
    assert parse_reference_targets('ENA-RUN','ERR1-ERR3') == ['ERR1-ERR3']
    assert parse_reference_targets('ENA-RUN','ERR1-ERR3',
        [{'database':'ENA-RUN','literal':'ERR1-ERR3','accessions':['ERR1','ERR2','ERR3']}]) == ['ERR1','ERR2','ERR3']


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
    from tests.support.magetab import render
    table=next(r[1] for r in render(data) if r[0]=='SDRF File')
    assert 'ftp://example.org/stats.txt' not in str(table)
    assert link in data['sample'][0]['supplementary_data']


def test_cross_reference_lists_split_only_verified_accession_literals():
    from meta_standards_converter.miniml.insdc_support import relations
    node=ET.fromstring('<SAMPLE><XREF_LINK><DB>ENA-STUDY</DB><ID>DRP1,ERP2</ID></XREF_LINK><XREF_LINK><DB>description</DB><ID>one,two</ID></XREF_LINK></SAMPLE>')
    assert relations(node)==[{'type':'ENA-STUDY','target':'DRP1'},{'type':'ENA-STUDY','target':'ERP2'},{'type':'description','target':'one,two'}]


def test_ena_range_links_are_resolved_before_individual_relations():
    from meta_standards_converter.sources.ena import ENASource
    from tests.support.transport import HTTP
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
