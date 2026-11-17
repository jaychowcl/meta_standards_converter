# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Reduced real GEO records exercise both converters without provider access."""
from tests.support.conversion import converted
from meta_standards_converter.converters import Converter, PreparationPolicy
from meta_standards_converter.magetab.exporter import MAGETabExporter
from meta_standards_converter.sources.metadata import GEOMetadataSource, MAGETabMetadataSource
import hashlib
import json
import re
import csv
from dataclasses import asdict
import xml.etree.ElementTree as ET
from pathlib import Path
import pytest

from tests.support.contracts import install_replay

ROOT=Path(__file__).resolve().parents[1]/'fixtures/chemistry/routing'
CASES=json.loads((ROOT/'manifest.json').read_text())

@pytest.mark.parametrize('case',CASES,ids=lambda c:c['id'])
def test_real_routing_regressions(case,workspace,monkeypatch):
    from meta_standards_converter.magetab.idf import IDFConstructor
    monkeypatch.setattr(IDFConstructor, "_current_idf_date", lambda self: "2026-09-10")
    raw=(ROOT/(case['id']+'.xml')).read_bytes()
    assert hashlib.sha256(raw).hexdigest()==case['sha256']
    calls=install_replay(monkeypatch,documents={case['gse']:raw})
    exporter = MAGETabExporter()
    converted(Converter(services={'magetab_exporter': exporter}), case['gse'], out_type='magetab', preparation=PreparationPolicy(enrichment='standard', expand_studies=False, linked_metadata=False, peer_provider=False), outdir='geo')
    converted(Converter(), case['gse'], out_type='json', preparation=PreparationPolicy(enrichment='standard', expand_studies=False, linked_metadata=False, peer_provider=False), outdir='json')
    json_exporter = MAGETabExporter()
    converted(Converter(services={'magetab_exporter': json_exporter}), str(next((workspace / 'json').glob('*.json'))), out_type='magetab', preparation=PreparationPolicy(enrichment='off', expand_studies=False, linked_metadata=False, peer_provider=False), outdir='json_ae')
    assert calls==['geo:'+case['gse']]*2
    for suffix in ('idf.txt','sdrf.txt'):
        geo=next((workspace/'geo').glob('*.'+suffix)).read_text()
        assert geo==next((workspace/'json_ae').glob('*.'+suffix)).read_text()
        expected=ROOT/(case['id']+'.'+suffix)
        assert geo==expected.read_text()
    rows=list(csv.reader(next((workspace/'geo').glob('*.sdrf.txt')).open(),delimiter='\t'))
    values=dict(zip(rows[0],rows[1]))
    for key in ('library construction', 'end bias', 'cell barcode size', 'umi barcode size',
                'cdna read', 'sample barcode read', 'sample barcode size'):
        assert 'Comment[' + key + ']' not in values
    retained=converted(Converter(), str(next((workspace / 'geo').glob('*.idf.txt'))), out_type='json', preparation=PreparationPolicy(enrichment='off', expand_studies=False, linked_metadata=False, peer_provider=False))[0].to_mapping()
    assert retained['sample']
    original_protocols=[node.text for node in ET.fromstring(raw).iter() if node.tag.endswith('}Extract-Protocol')]
    retained_protocols=[c.get('extract_protocol') for sample in retained['sample'] for c in sample['channel']]
    for protocol in original_protocols:
        assert re.sub(r'\r\n|[\t\r\n]', ' ', protocol) in retained_protocols
    # Material Type owns these populated comments, including repeated index
    # attributes, independently for each source row.
    assert len(retained['series']['assay_paths'])==len(rows)-1
    for row, path in zip(rows[1:],retained['series']['assay_paths']):
        extract=next(step for step in path['steps'] if step['kind']=='extract')
        start=rows[0].index('Extract Name');end=rows[0].index('Assay Name')
        expected=[(h[8:-1],v) for h,v in list(zip(rows[0],row))[start:end] if h.startswith('Comment[') and v.strip()]
        assert rows[0][start+1]=='Material Type'
        assert extract['material_type']['value']==row[start+1]
        assert not extract.get('comments')
        assert [(c['name'],c['value']) for c in extract['material_type'].get('comments',[])]==expected
    for rendering in (exporter, json_exporter):
        assert asdict(rendering.constructor.sdrf_constructor.last_sdrf_audit)==case['audit']
