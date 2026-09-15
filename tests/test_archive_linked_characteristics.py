from copy import deepcopy
import xml.etree.ElementTree as ET
import pytest
from tests.test_native_archive_parsers import fixture_records
from meta_standards_converter.miniml.ena_parser import ENAParser


def record(items, native=None):
    records=fixture_records('ena'); node=records.xml[1].find('SAMPLE'); acc=node.get('accession')
    attrs=node.find('SAMPLE_ATTRIBUTES')
    for value in native or [{'name':'tissue','value':'brain'}]:
        attr=ET.SubElement(attrs,'SAMPLE_ATTRIBUTE')
        for tag, key in [('TAG','name'),('VALUE','value'),('UNITS','unit')]:
            if value.get(key): ET.SubElement(attr,tag).text=value[key]
    records.linked.append({'provider':'biosamples','kind':'sample','accession':acc,
        'metadata':{'accession':acc,'characteristics':{'tissue':items}}})
    return records


def test_linked_ontology_completes_unique_value_and_keeps_unknown_sibling():
    term='http://purl.obolibrary.org/obo/UBERON_0000955'
    data=ENAParser().parse(record([{'text':'brain','ontologyTerms':[term],'custom':'keep'}])).to_mapping()
    attrs=[v for v in data['sample'][0]['channel'][0]['characteristics'] if v['name']=='tissue']
    assert len(attrs)==1 and attrs[0]['term_accession_number']==term
    assert attrs[0]['term_source_ref']=='UBERON'
    assert all(any(v.get('term_accession_number')==term for v in p['steps'][0]['characteristics']) for p in data['series']['assay_paths'])
    residual=[r for r in data['extensions']['insdc']['records'] if r['provider']=='biosamples'][0]
    assert residual['metadata']['characteristics']['tissue']==[{'text':'brain','custom':'keep'}]


@pytest.mark.parametrize('items,native', [
    ([{'text':'brain','ontologyTerms':['UBERON:0000955']}], [{'name':'tissue','value':'brain'},{'name':'tissue','value':'brain'}]),
    ([{'text':'brain','ontologyTerms':['UBERON:0000955']},{'text':'brain','ontologyTerms':['UBERON:0000956']}], [{'name':'tissue','value':'brain'}]),
    ([{'text':'brain','unit':'g','ontologyTerms':['UBERON:0000955']}], [{'name':'tissue','value':'brain','unit':'mg'}]),
])
def test_ambiguous_or_conflicting_linked_groups_remain_residual(items,native):
    data=ENAParser().parse(record(items,native)).to_mapping()
    attrs=[v for v in data['sample'][0]['channel'][0]['characteristics'] if v['name']=='tissue']
    assert not any(v.get('term_accession_number') for v in attrs)
    residual=[r for r in data['extensions']['insdc']['records'] if r['provider']=='biosamples'][0]
    assert residual['metadata']['characteristics']['tissue']==items


def test_new_linked_characteristic_keeps_repeated_occurrences_and_units():
    items=[{'text':'brain','unit':'mg'}]*2
    data=ENAParser().parse(record(items,[{'name':'other','value':'x'}])).to_mapping()
    attrs=[v for v in data['sample'][0]['channel'][0]['characteristics'] if v['name']=='tissue']
    assert len(attrs)==2 and all(v['unit']['value']=='mg' for v in attrs)


def test_ebi_ontology_namespace_and_unit_are_completed_together():
    term='http://www.ebi.ac.uk/efo/EFO_0004472'
    data=ENAParser().parse(record([{'text':'brain','unit':'mg','ontologyTerms':[term]}])).to_mapping()
    value=next(v for v in data['sample'][0]['channel'][0]['characteristics'] if v['name']=='tissue')
    assert value['term_source_ref']=='EFO'
    assert value['term_accession_number']==term and value['unit']['value']=='mg'


def test_residual_annotations_retain_their_identifying_literal():
    data=ENAParser().parse(record([{'text':'brain','ontologyTerms':['A','B'],'custom':'x'}],
        [{'name':'tissue','value':'brain'},{'name':'tissue','value':'heart'}])).to_mapping()
    residual=next(r for r in data['extensions']['insdc']['records'] if r['provider']=='biosamples')
    assert residual['metadata']['characteristics']['tissue']==[{'text':'brain','ontologyTerms':['A','B'],'custom':'x'}]


def test_fully_mapped_biosamples_record_has_no_empty_residual():
    data=ENAParser().parse(record([{'text':'brain','ontologyTerms':['http://purl.obolibrary.org/obo/UBERON_0000955']}])).to_mapping()
    assert not any(r['provider']=='biosamples' for r in data['extensions']['insdc']['records'])


def test_organism_residual_keeps_conflicting_annotations():
    records=record([])
    name=records.xml[1].findtext('SAMPLE/SAMPLE_NAME/SCIENTIFIC_NAME')
    item={'text':name,'ontologyTerms':['different:identifier'],'custom':'retain me'}
    records.linked[-1]['metadata']['characteristics']={'organism':[item]}
    data=ENAParser().parse(records).to_mapping()
    residual=next(r for r in data['extensions']['insdc']['records'] if r['provider']=='biosamples')
    assert residual['metadata']['characteristics']['organism']==[item]


def test_organism_residual_preserves_unrepresented_repetitions():
    records=record([])
    name=records.xml[1].findtext('SAMPLE/SAMPLE_NAME/SCIENTIFIC_NAME')
    records.linked[-1]['metadata']['characteristics']={'organism':[{'text':name}]*3}
    data=ENAParser().parse(records).to_mapping()
    residual=next(r for r in data['extensions']['insdc']['records'] if r['provider']=='biosamples')
    channel=data['sample'][0]['channel'][0]
    mapped=len([v for v in channel['characteristics'] if v['name']=='organism' and v['value']==name])+len(channel['organism'])
    assert residual['metadata']['characteristics']['organism']==[{'text':name}]*(3-mapped)
