from meta_standards_converter.miniml.archive_residuals import Projection
from meta_standards_converter.miniml.insdc_support import relations, tree
import xml.etree.ElementTree as ET


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
