# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from copy import deepcopy
import json
import xml.etree.ElementTree as ET
from meta_standards_converter.miniml import MINiMLCodec
from meta_standards_converter.miniml.sra_parser import SRAParser
from meta_standards_converter.miniml.ena_parser import ENAParser
from meta_standards_converter.metadata.archive_enrichment import merge_archive_metadata, linked_accessions
from tests.test_native_archive_parsers import fixture_records
from tests.test_archive_fidelity_enrichment import enriched_native, workflow


def nodes(value):
    if isinstance(value,dict):
        yield value
        for v in value.values(): yield from nodes(v)
    elif isinstance(value,list):
        for v in value: yield from nodes(v)


def test_residual_only_contract_removes_mapped_leaves_not_unmapped_siblings():
    records=fixture_records('sra')
    study=records.xml[0].find('.//STUDY')
    ET.SubElement(study,'UNMAPPED_NOTE').text='preserve this note'
    data=SRAParser().parse(records).to_mapping();extra=data['extensions']['insdc']
    assert extra['version']=='2.0'
    assert not any(n.get('tag') in ('STUDY_TITLE','STUDY_ABSTRACT','LIBRARY_STRATEGY','SCIENTIFIC_NAME') for n in nodes(extra))
    assert 'preserve this note' in str(extra)
    assert not any(n.get('tag') in ('STUDY','SAMPLE','EXPERIMENT','RUN_SET') for r in extra['records'] if r['kind']=='EXPERIMENT_PACKAGE' for n in nodes(r['metadata']))
    MINiMLCodec().decode(data,strict=True)


def test_equal_text_in_unrelated_field_is_not_removed():
    records=fixture_records('ena');sample=records.xml[1].find('SAMPLE')
    text=sample.findtext('SAMPLE_NAME/SCIENTIFIC_NAME')
    ET.SubElement(sample,'CUSTOM_UNMAPPED').text=text
    data=ENAParser().parse(records).to_mapping()
    assert any(n.get('tag')=='CUSTOM_UNMAPPED' and n.get('text')==text for n in nodes(data['extensions']))


def test_unbound_magetab_annotations_survive_native_enrichment():
    extra = workflow().to_mapping()
    occurrence = {'sdrf':'source.sdrf','row_index':1,'column_index':4,'header':'Characteristics[note]',
                  'value':'ORPHAN-SENTINEL', 'sample_ref':extra['sample'][0]['iid']}
    extra['extensions'] = {'magetab':{'unbound_annotations':[occurrence]}}
    merged, issues = merge_archive_metadata(enriched_native(), MINiMLCodec().decode(extra).package, prefer=True)
    assert 'ORPHAN-SENTINEL' in str(merged.to_mapping()['extensions'])
    assert any('Unbound' in i for i in issues)


def test_enrichment_residual_keeps_displaced_values_not_complete_packages():
    original=enriched_native();extra=workflow().to_mapping()
    old=original.to_mapping()['series']['title']
    extra['series']['title']='replacement title'
    extra['sample'][0]['channel'][0]['characteristics'].append({'name':'unrelated note','value':'retain'})
    data,_=merge_archive_metadata(original,MINiMLCodec().decode(extra).package,prefer=True)
    data=data.to_mapping()
    assert data['series']['title']=='replacement title'
    assert old in str(data['extensions'])
    assert not any('miniml_schema_version' in r['metadata'] for r in data['extensions']['insdc']['records'])
    assert not any(n.get('tag')=='LIBRARY_STRATEGY' for n in nodes(data['extensions']))


def test_parallel_files_and_statistics_are_not_duplicated_in_residual():
    records=fixture_records('ena')
    records.indexed['read_run'][0].update(read_count='5',base_count='100')
    data=ENAParser().parse(records).to_mapping()
    row=next((r['metadata'] for r in data['extensions']['insdc']['records'] if r['kind']=='read_run'),{})
    assert not any(k in row for k in ('fastq_ftp','fastq_md5','fastq_bytes','read_count','base_count'))


def test_saved_v1_is_readable_and_saved_v2_enrichment_preserves_old_title():
    data=enriched_native().to_mapping();data['extensions']['insdc']['version']='1.0'
    package=MINiMLCodec().decode(json.loads(json.dumps(data))).package
    extra=workflow().to_mapping();old=data['series']['title'];extra['series']['title']='changed'
    merged,_=merge_archive_metadata(package,MINiMLCodec().decode(extra).package,prefer=True)
    assert old in str(merged.to_mapping()['extensions'])
    assert merged.to_mapping()['extensions']['insdc']['version']=='2.0'


def test_fully_mapped_record_is_omitted_and_repetitions_remain_in_core():
    from meta_standards_converter.miniml.archive_residuals import finalize
    data={'miniml_schema_version':'3.0','source':{'format':'ENA'},'series':{'iid':'ERP1','accession':[{'value':'ERP1','database':'ENA'}],'title':'title'},'database':[{'iid':'ENA','name':'ENA'}]}
    record={'provider':'ena','kind':'STUDY','accession':'ERP1','metadata':{'tag':'STUDY','attributes':{'accession':'ERP1'},'children':[{'tag':'DESCRIPTOR','children':[{'tag':'STUDY_TITLE','text':'title'}]}]}}
    assert finalize(data,[record]).to_mapping()['extensions']['insdc']['records']==[]


def test_private_evidence_does_not_serialize_and_unmapped_contact_detail_remains():
    package=SRAParser().parse(fixture_records('sra'))
    data=package.to_mapping()
    assert '_archive_source_records' not in json.dumps(data)
    assert 'sec_email' not in json.dumps(data['extensions'])
    assert any(c.get('extensions', {}).get('secondary_email') for c in data['contributor'])
    # The mapped primary contact email is absent from its native XML residual.
    assert not any(n.get('attributes',{}).get('email')=='bdaisley@uwo.ca' for n in nodes(data['extensions']))


def test_saved_v2_preserves_unmapped_siblings_across_repeated_enrichment():
    package=SRAParser().parse(fixture_records('sra'))
    extra=workflow().to_mapping();extra['series']['title']='new'
    first,_=merge_archive_metadata(package,MINiMLCodec().decode(extra).package,prefer=True)
    saved=MINiMLCodec().decode(json.loads(json.dumps(first.to_mapping()))).package
    second,_=merge_archive_metadata(saved,MINiMLCodec().decode(extra).package,prefer=True)
    assert 'sec_email' not in str(second.to_mapping()['extensions'])
    assert any(c.get('extensions', {}).get('secondary_email') for c in second.to_mapping()['contributor'])
    assert not any('miniml_schema_version' in r['metadata'] for r in second.to_mapping()['extensions']['insdc']['records'])


def test_study_residual_recognizes_description_dates_and_verified_attribute_links():
    from meta_standards_converter.miniml.archive_residuals import Projection
    from meta_standards_converter.miniml.insdc_support import tree
    data = {'series': {'iid': 'PRJEB1', 'accession': [{'value': 'ERP1'}, {'value': 'PRJEB1'}],
                       'summary': 'Submitted description with whitespace.',
                       'status': [{'database': 'ENA', 'release_date': '2021-03-09', 'last_update_date': '2021-03-20'}],
                       'relation': [{'type': 'ArrayExpress', 'target': 'E-MTAB-1'}]}}
    node = tree(ET.fromstring('''<STUDY accession="ERP1"><DESCRIPTOR>
      <STUDY_DESCRIPTION>Submitted description\n with whitespace.</STUDY_DESCRIPTION>
      </DESCRIPTOR><STUDY_ATTRIBUTES>
      <STUDY_ATTRIBUTE><TAG>ENA-FIRST-PUBLIC</TAG><VALUE>2021-03-09</VALUE></STUDY_ATTRIBUTE>
      <STUDY_ATTRIBUTE><TAG>ENA-LAST-UPDATE</TAG><VALUE>2021-03-09</VALUE></STUDY_ATTRIBUTE>
      <STUDY_ATTRIBUTE><TAG>ArrayExpress</TAG><VALUE>E-MTAB-1</VALUE><CUSTOM>keep</CUSTOM></STUDY_ATTRIBUTE>
      <STUDY_ATTRIBUTE><TAG>UNMAPPED</TAG><VALUE>2021-03-09</VALUE></STUDY_ATTRIBUTE>
      </STUDY_ATTRIBUTES></STUDY>'''))
    residual = Projection(data).xml(node, 'STUDY', 'ERP1', 'ena')
    rendered = str(residual)
    assert 'STUDY_DESCRIPTION' not in rendered
    assert 'ENA-FIRST-PUBLIC' not in rendered
    assert 'E-MTAB-1' not in rendered
    assert 'ENA-LAST-UPDATE' in rendered  # Differing browser and index dates survive.
    assert 'UNMAPPED' in rendered and 'CUSTOM' in rendered


def test_residual_list_matching_counts_separate_source_occurrences():
    from meta_standards_converter.miniml.archive_residuals import diff
    occurrence = {'name': 'organism', 'value': 'Danio rerio'}
    assert diff([occurrence, occurrence], [occurrence]) == [occurrence]


def test_residual_character_matching_does_not_reuse_one_destination_occurrence():
    from meta_standards_converter.miniml.archive_residuals import Projection
    sample = {'iid': 'S1', 'channel': [{'characteristics': [{'name': 'organism', 'value': 'Danio rerio'}]}]}
    projection = Projection({'series': {}, 'sample': [sample]})
    assert projection.character(sample, 'organism', 'Danio rerio')
    assert not projection.character(sample, 'organism', 'Danio rerio')


def test_single_channel_residual_prunes_mapped_siblings_after_one_conflict():
    from meta_standards_converter.miniml.archive_residuals import diff
    source={'sample':[{'iid':'SRS1','channel':[{'organism':[{'value':'Mus musculus','taxid':'10090'}],
        'molecule':{'value':'polyA RNA'},'extract_protocol':'supplied extraction',
        'characteristics':[{'name':'strain','value':'FVN/Swiss'},{'name':'cell type','value':'Astrocyte'}]}]}]}
    target=deepcopy(source)
    target['sample'][0]['channel'][0]['characteristics'][1]={'name':'cell type','value':'astrocyte','term_source_ref':'CL','term_accession_number':'CL:0000127'}
    residual=diff(source,target)
    channel=residual['sample'][0]['channel'][0]
    assert channel=={'characteristics':[{'name':'cell type','value':'Astrocyte'}]}
    # Multiple channels without explicit identities remain ambiguous.
    target['sample'][0]['channel'].append(deepcopy(target['sample'][0]['channel'][0]))
    residual=diff(source,target)
    assert residual['sample'][0]['channel'][0]['extract_protocol']=='supplied extraction'


def test_residual_keeps_displaced_biological_source_but_not_package_source_marker():
    from meta_standards_converter.miniml.archive_residuals import diff
    source={'source':{'format':'GEO'},'sample':[{'iid':'SRS1','channel':[{'source':{'value':'liver'},'molecule':{'value':'RNA'}}]}]}
    target=deepcopy(source);target['source']['format']='SRA';target['sample'][0]['channel'][0]['source']={'value':'brain'}
    assert diff(source,target)=={'sample':[{'iid':'SRS1','channel':[{'source':{'value':'liver'}}]}]}


def test_filename_uri_transform_prunes_only_represented_source_attributes():
    records = fixture_records('sra')
    run = records.xml[0].find('.//RUN')
    existing = run.find('SRAFiles')
    if existing is not None: run.remove(existing)
    group = ET.SubElement(run, 'SRAFiles')
    file = ET.SubElement(group, 'SRAFile', filename='a#1.fastq.gz', url='https://files.example.org/a#1.fastq.gz', semantic_name='fastq')
    ET.SubElement(file, 'Alternatives', url='https://mirror.example.org/a#1.fastq.gz')
    ET.SubElement(file, 'UNMAPPED_NOTE').text='source-only sibling'
    data = SRAParser().parse(records).to_mapping()
    files = [f for s in data['sample'] for r in s.get('sra_run', []) for f in r.get('files', []) if f.get('filename') == 'a#1.fastq.gz']
    assert len(files) == 2 and all('%231.fastq.gz' in f['uri'] for f in files)
    residual = data['extensions']['insdc']
    assert 'source-only sibling' in str(residual)
    assert not any(n.get('attributes', {}).get('url') in {'https://files.example.org/a#1.fastq.gz','https://mirror.example.org/a#1.fastq.gz'} for n in nodes(residual))
