# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Reduced real GEO records exercise both converters without provider access."""
import hashlib
import json
import csv
from dataclasses import asdict
import xml.etree.ElementTree as ET
from pathlib import Path
import pytest
from meta_standards_converter.converters import GEO2AEConverter, GEO2JSONConverter, JSON2AEConverter, AE2JSONConverter
from tests.support.contracts import install_replay

ROOT=Path(__file__).resolve().parents[1]/'fixtures/chemistry/routing'
CASES=json.loads((ROOT/'manifest.json').read_text())

@pytest.mark.parametrize('case',CASES,ids=lambda c:c['id'])
def test_real_routing_regressions(case,workspace,monkeypatch):
    raw=(ROOT/(case['id']+'.xml')).read_bytes()
    assert hashlib.sha256(raw).hexdigest()==case['sha256']
    calls=install_replay(monkeypatch,documents={case['gse']:raw})
    converter=GEO2AEConverter()
    converter.convert(case['gse'],out='geo')
    GEO2JSONConverter().convert(case['gse'],out='json')
    json_converter=JSON2AEConverter()
    json_converter.convert(str(next((workspace/'json').glob('*.json'))),out='json_ae',enrich=False)
    assert calls==['geo:'+case['gse']]*2
    for suffix in ('idf.txt','sdrf.txt'):
        geo=next((workspace/'geo').glob('*.'+suffix)).read_text()
        assert geo==next((workspace/'json_ae').glob('*.'+suffix)).read_text()
        expected=ROOT/(case['id']+'.'+suffix)
        assert geo==expected.read_text()
    rows=list(csv.reader(next((workspace/'geo').glob('*.sdrf.txt')).open(),delimiter='\t'))
    values=dict(zip(rows[0],rows[1]))
    assert values.get('Comment[library construction]')==('10x 3 prime v2' if case['id']=='GSM5388031' else '10x 3 prime v3.1')
    assert values['Comment[end bias]']=='3 prime tag'
    if case['id']=='mixed':
        col=rows[0].index('Comment[library construction]')
        assert [row[col] for row in rows[1:]]==['10x 3 prime v3.1','10x Visium']
        assert rows[0].count('Comment[sample barcode read]')==2
        assert [v for h,v in zip(rows[0],rows[1]) if h=='Comment[sample barcode read]']==['index1','index2']
        assert [v for h,v in zip(rows[0],rows[1]) if h=='Comment[sample barcode size]']==['10','10']
    else:
        for key in ('cell barcode size','umi barcode size','cdna read','sample barcode read'):
            assert 'Comment['+key+']' not in values
    retained=AE2JSONConverter().convert(str(next((workspace/'geo').glob('*.idf.txt'))))[0].to_mapping()
    assert retained['sample']
    original_protocols=[node.text for node in ET.fromstring(raw).iter() if node.tag.endswith('}Extract-Protocol')]
    retained_protocols=[c.get('extract_protocol') for sample in retained['sample'] for c in sample['channel']]
    for protocol in original_protocols:
        assert ' '.join(protocol.split()) in retained_protocols
    # The real MAGE-TAB parser retains every extract comment, including repeated
    # index attributes, independently for each source row.
    for row, path in zip(rows[1:],retained['series']['assay_paths']):
        extract=next(step for step in path['steps'] if step['kind']=='extract')
        start=rows[0].index('Extract Name');end=rows[0].index('Assay Name')
        expected=[(h[8:-1],v) for h,v in list(zip(rows[0],row))[start:end] if h.startswith('Comment[')]
        assert [(c['name'],c['value']) for c in extract.get('comments',[])]==expected
    assert len(retained['series']['assay_paths'])==len(rows)-1
    for conversion in (converter,json_converter):
        assert asdict(conversion.ae_constructor.sdrf_constructor.last_sdrf_audit)==case['audit']
