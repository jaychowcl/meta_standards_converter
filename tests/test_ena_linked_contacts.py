# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Native ENA imports use explicitly identified BioSample contacts without a peer import."""
from copy import deepcopy
import xml.etree.ElementTree as ET
from meta_standards_converter.sources.ena import ENASource
from meta_standards_converter.miniml.ena_parser import ENAParser
from tests.test_native_archive_parsers import fixture_records
from tests.test_native_archive_sources import HTTP


def test_linked_biosample_contact_has_sample_scope_and_pruned_residual():
    records = fixture_records('ena')
    before = ENAParser().parse(deepcopy(records)).to_mapping()
    owner = '<Owner><Name>Native laboratory</Name><Contacts><Contact sec_email="other@example.org"><Name><First>Named</First><Last>Scientist</Last></Name></Contact></Contacts></Owner>'
    def handler(url, params, fmt):
        assert params['db'] == 'biosample'
        if 'esearch' in url:
            assert params['term'] == 'SAMN14218700[Accession]'
            return {'esearchresult': {'count':'1', 'idlist':['14218700']}}
        return ET.fromstring('<BioSampleSet><BioSample id="14218700" accession="SAMN14218700">'+owner+'<Unknown>retain</Unknown></BioSample></BioSampleSet>')
    source = ENASource(http=HTTP(handler))
    source.fetch_biosamples(['SAMN14218700'], records)
    data = ENAParser().parse(records).to_mapping()
    contact = next(c for c in data['contributor'] if c.get('person', {}).get('first') == 'Named')
    assert {'ref':contact['iid']} in data['sample'][0]['contact_ref']
    assert contact['extensions']['secondary_email'] == 'other@example.org'
    assert data['series']['iid'] == before['series']['iid']
    assert data['sample'][0]['sra_run'] == before['sample'][0]['sra_run']
    assert 'retain' in str(data['extensions'])
    assert 'Named' not in str(data['extensions'])
    assert all(p['db']=='biosample' for _,p,_ in source.http.calls)


def test_ena_fetch_requests_only_its_inventory_biosamples(monkeypatch):
    source = ENASource(http=HTTP(lambda u,p,f: []))
    monkeypatch.setattr(source, 'verified_xml', lambda *args: None)
    monkeypatch.setattr(source, 'search', lambda kind,*args: [{'sample_accession':'SAMN1','run_accession':'SRR1'}] if kind=='read_run' else [])
    called=[]
    monkeypatch.setattr(source,'fetch_biosamples',lambda accessions,records: called.extend(accessions))
    from meta_standards_converter.sources.archive_support import StudySeed
    source.fetch(StudySeed('SRP1','PRJNA1'))
    assert called == ['SAMN1']
