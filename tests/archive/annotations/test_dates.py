# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from tests.support.conversion import converted
from meta_standards_converter.converters import Converter, PreparationPolicy
from meta_standards_converter.magetab.exporter import MAGETabExporter
from copy import deepcopy
import json
import xml.etree.ElementTree as ET
import pytest
from meta_standards_converter.miniml import MINiMLCodec
from meta_standards_converter.miniml.ena_parser import ENAParser
from meta_standards_converter.magetab.constructor import AEConstructor
from meta_standards_converter.metadata.archive_enrichment import merge_archive_metadata
from tests.support.archive import fixture_records
from tests.support.archive import NoRuns, PubMed
from tests.support.archive import nodes
from tests.support.archive import enriched_native, workflow

from tests.support.archive import legacy


from meta_standards_converter.magetab.idf import IDFConstructor

from tests.support.archive import workflow, enriched_native


def test_administrative_dates_are_scoped_counted_and_idempotent():
    from meta_standards_converter.miniml.archive_dates import normalize_archive_dates
    data=legacy();series_status=deepcopy(data['series'].get('status'))
    normalize_archive_dates(data)
    sample=data['sample'][0];statuses=sample['status']
    assert statuses.count({'database':'ENA','release_date':'2020-01-01'})==2
    assert {'database':'ENA','last_update_date':'2020-02-03T12:01:02Z'} in statuses
    assert {'database':'INSDC','last_update_date':'2020-02-04'} in statuses
    assert {'database':'BioSample','last_update_date':'2020-04-05'} in statuses
    assert data['series'].get('status')==series_status
    assert all(c['name'] not in ('ENA-FIRST-PUBLIC','ENA-LAST-UPDATE','INSDC last update') for c in sample['channel'][0]['characteristics'])
    assert {'name':'collection_date','value':'2019-12'} in sample['channel'][0]['characteristics']
    assert 'INSDC' in {d['iid'] for d in data['database']}
    once=deepcopy(data);normalize_archive_dates(data);assert data==once


def test_idf_partition_happens_after_overlay_and_ae_accession_is_not_secondary():
    data=legacy();data['series']['accession'] += [{'value':'E-MTAB-308','database':'ArrayExpress'},{'value':'GSE1','database':'GEO'}]
    data['series']['comments']=[{'name':'UserNote','value':'keep'}]
    rows=AEConstructor().miniml2magetab(MINiMLCodec().decode(data).package)
    first=next(i for i,r in enumerate(rows) if r[0].startswith('Comment['))
    assert all(r[0].startswith('Comment[') for r in rows[first:])
    values={r[0]:r[1:] for r in rows}
    assert 'E-MTAB-308' not in values['Comment[SecondaryAccession]']
    assert values['Comment[ArrayExpressAccession]']==['E-MTAB-308']
    assert ('GSE1','GEO') in list(zip(values['Comment[SecondaryAccession]'],values['Comment[SecondaryAccessionTermSourceRef]']))
    assert values['Comment[UserNote]']==['keep']


@pytest.mark.parametrize('enrich',[False,True])
def test_existing_json_exports_without_admin_characteristics_or_mutation(tmp_path,enrich):

    from meta_standards_converter.metadata.enrichment import MINiMLEnricher
    path=tmp_path/'legacy.json';path.write_text(json.dumps(legacy()));before=path.read_bytes()
    client=PubMed();service=MINiMLEnricher(pubmed_fetcher=client,insdc_fetcher=NoRuns())
    rows=converted(Converter(services={'enricher': service, 'magetab_exporter': MAGETabExporter(AEConstructor())}), str(path), out_type='magetab', preparation=PreparationPolicy(enrichment='standard' if enrich else 'off', expand_studies=False, linked_metadata=False, peer_provider=False))[0]
    table=next(r[1] for r in rows if r[0]=='SDRF File')
    assert not any(any(key in label for key in ('ENA-FIRST-PUBLIC','ENA-LAST-UPDATE','INSDC last update')) for label in table[0])
    assert 'Characteristics[collection_date]' in table[0]
    assert path.read_bytes()==before and not client.calls


def test_fresh_dates_project_to_status_with_unmapped_xml_siblings_retained():
    records=fixture_records('ena');sample=records.xml[1].find('SAMPLE')
    attributes=sample.find('SAMPLE_ATTRIBUTES')
    attr=ET.SubElement(attributes,'SAMPLE_ATTRIBUTE')
    ET.SubElement(attr,'TAG').text='ENA-LAST-UPDATE';ET.SubElement(attr,'VALUE').text='2021-02-03T00:01:02Z'
    ET.SubElement(attr,'CUSTOM').text='keep sibling'
    data=ENAParser().parse(records).to_mapping();s=data['sample'][0]
    assert {'database':'ENA','last_update_date':'2021-02-03T00:01:02Z'} in s['status']
    assert not any(c['name']=='ENA-LAST-UPDATE' for c in s['channel'][0]['characteristics'])
    assert 'keep sibling' in str(data['extensions'])
    assert not any(n.get('tag')=='VALUE' and n.get('text')=='2021-02-03T00:01:02Z' for n in nodes(data['extensions']))


def test_enrichment_keeps_conflicting_dates_in_the_same_archive():
    data=enriched_native().to_mapping();data['series']['status']=[{'database':'ENA','release_date':'2020-01-01'}]
    extra=workflow().to_mapping();extra['series']['status']=[{'database':'ENA','release_date':'2021-01-01'}]
    merged,_=merge_archive_metadata(MINiMLCodec().decode(data).package,MINiMLCodec().decode(extra).package,prefer=True)
    assert {s['release_date'] for s in merged.to_mapping()['series']['status'] if s['database']=='ENA'}=={'2020-01-01','2021-01-01'}


@pytest.mark.parametrize('name', ['ENA first public','ENA last update'])
def test_space_separated_dates_preserve_precision_and_residual_siblings(name):
    records=fixture_records('ena');attributes=records.xml[1].find('SAMPLE/SAMPLE_ATTRIBUTES')
    attr=ET.SubElement(attributes,'SAMPLE_ATTRIBUTE')
    ET.SubElement(attr,'TAG').text=name;ET.SubElement(attr,'VALUE').text='2020-03'
    ET.SubElement(attr,'CUSTOM').text='preserve me'
    data=ENAParser().parse(records).to_mapping()
    field='release_date' if name.endswith('public') else 'last_update_date'
    assert {'database':'ENA',field:'2020-03'} in data['sample'][0]['status']
    assert not any(c['name']==name for c in data['sample'][0]['channel'][0]['characteristics'])
    assert 'preserve me' in str(data['extensions'])
    assert not any(n.get('tag')=='VALUE' and n.get('text')=='2020-03' for n in nodes(data['extensions']))
    from tests.support.magetab import render
    table=next(r[1] for r in render(data) if r[0]=='SDRF File')
    assert not any(name in header for header in table[0])


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


def test_analysis_reference_fetch_is_bounded_and_validates_sample_binding():
    from meta_standards_converter.sources.ena import ENASource
    from tests.support.transport import HTTP
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
