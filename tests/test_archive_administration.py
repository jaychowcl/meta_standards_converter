# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from copy import deepcopy
import xml.etree.ElementTree as ET
from tests.test_archive_export_cleanup import legacy
from tests.test_native_archive_parsers import fixture_records
from meta_standards_converter.miniml.ena_parser import ENAParser


def test_administration_preserves_occurrences_and_cleans_copies():
    from meta_standards_converter.miniml.archive_administration import normalize_administration
    data=legacy();sample=data['sample'][0]
    items=[{'name':'INSDC status','value':'public'},{'name':'INSDC status','value':'public'},
           {'name':'INSDC center name','value':'SC'},{'name':'INSDC center alias','value':'SC'}]
    sample['channel'][0]['characteristics']+=deepcopy(items)
    for path in data['series']['assay_paths']:path['steps'][0]['characteristics']+=deepcopy(items)
    normalize_administration(data)
    statuses=[s for s in sample['status'] if s.get('comment')==[items[0]]]
    assert len(statuses)==2
    orgs={o['iid']:o for o in data['organization']}
    links=[r for r in sample['relation'] if r['type']=='archive center' and orgs[r['target']].get('role','').startswith('INSDC')]
    assert len(links)==2
    assert {orgs[r['target']]['role'] for r in links}=={'INSDC center name','INSDC center alias'}
    assert all(c not in items for c in sample['channel'][0]['characteristics'])
    assert all(c not in items for p in data['series']['assay_paths'] for st in p['steps'] for c in st.get('characteristics',[]))
    before=deepcopy(data);normalize_administration(data);assert data==before


def test_fresh_independent_center_and_status_sources_remain_scoped():
    records=fixture_records('ena');node=records.xml[1].find('SAMPLE');acc=node.get('accession');node.set('center_name','SC')
    attrs=node.find('SAMPLE_ATTRIBUTES');a=ET.SubElement(attrs,'SAMPLE_ATTRIBUTE')
    ET.SubElement(a,'TAG').text='INSDC status';ET.SubElement(a,'VALUE').text='public';ET.SubElement(a,'CUSTOM').text='keep'
    records.linked.append({'provider':'biosamples','kind':'sample','accession':acc,'metadata':{'accession':acc,'status':'PUBLIC','characteristics':{'INSDC center name':[{'text':'SC'}]}}})
    data=ENAParser().parse(records).to_mapping();sample=data['sample'][0]
    assert len([r for r in sample['relation'] if r['type']=='archive center'])==2
    assert any(s.get('database')=='INSDC' and s.get('comment')==[{'name':'INSDC status','value':'public'}] for s in sample['status'])
    assert any(s.get('database')=='BioSamples' and s.get('comment')==[{'name':'status','value':'PUBLIC'}] for s in sample['status'])
    assert 'keep' in str(data['extensions'])
    assert not any(c['name']=='INSDC status' for c in sample['channel'][0]['characteristics'])
