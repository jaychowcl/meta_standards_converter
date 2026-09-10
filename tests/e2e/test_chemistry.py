# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Both real converters parse and publish recorded chemistry in synthetic envelopes."""
import csv
import io
import hashlib
import json
from pathlib import Path
import pytest
from meta_standards_converter.converters import GEO2AEConverter, GEO2JSONConverter, JSON2AEConverter
from tests.support.contracts import install_replay

ROOT=Path(__file__).resolve().parents[1]/'fixtures/chemistry'
CASES=json.loads((ROOT/'manifest.json').read_text())

@pytest.mark.parametrize('case',CASES,ids=lambda c:c['id'])
def test_chemistry_complete_converter_outputs(case,workspace,monkeypatch):
    xml=(ROOT/(case['id']+'.xml')).read_bytes()
    assert hashlib.sha256(xml).hexdigest()==case['envelope']['sha256']
    accession=case['envelope']['gse']
    calls=install_replay(monkeypatch,documents={accession:xml})
    geo=GEO2AEConverter()
    geo.convert(accession,out='geo')
    # The stored-source route passes through the real GEO parser and enrichment.
    GEO2JSONConverter().convert(accession,out='json')
    json_path=next((workspace/'json').glob('*.json'))
    converter=JSON2AEConverter()
    converter.convert(str(json_path),out='json_ae',enrich=False)
    assert calls==['geo:'+accession,'geo:'+accession]
    geo_sdrf=next((workspace/'geo').glob('*.sdrf.txt')).read_text()
    json_sdrf=next((workspace/'json_ae').glob('*.sdrf.txt')).read_text()
    assert geo_sdrf==json_sdrf
    assert geo_sdrf==(ROOT/(case['id']+'.sdrf.txt')).read_text()
    warnings=geo.ae_constructor.sdrf_constructor.last_sdrf_audit.warnings
    codes=[w.split(' chemistry ',1)[1].split(':',1)[0] for w in warnings if ' chemistry ' in w]
    assert codes==case.get('diagnostics',[])

    from meta_standards_converter.converters import AE2JSONConverter
    reparsed=AE2JSONConverter().convert(str(next((workspace/'geo').glob('*.idf.txt'))))
    retained=reparsed[0].to_mapping()
    # Original protocol and every repeated SDRF column survive the real parser.
    fields=json.loads((ROOT/(case['id']+'.json')).read_text())
    # IDF protocol registration has always flattened whitespace. Compare the
    # published description exactly, rather than changing that separate contract.
    idf_rows=list(csv.reader(next((workspace/'geo').glob('*.idf.txt')).open(),delimiter='\t'))
    descriptions=next(r[1:] for r in idf_rows if r[0]=='Protocol Description')
    assert ' '.join(fields['protocol'].split()) in descriptions
    assert ' '.join(fields['protocol'].split()) in [c.get('extract_protocol') for s in retained['sample'] for c in s['channel']]
    from meta_standards_converter.miniml import MINiMLCodec
    path=workspace/'roundtrip.json'
    path.write_text(json.dumps(MINiMLCodec().encode_many(reparsed)))
    JSON2AEConverter().convert(str(path),out='roundtrip',enrich=False)
    original=list(csv.reader(io.StringIO(geo_sdrf),delimiter='\t'))
    again=list(csv.reader(next((workspace/'roundtrip').glob('*.sdrf.txt')).open(),delimiter='\t'))
    # The existing semantic renderer groups characteristics before comments.
    # Every comment value and repeated occurrence must survive in exact order.
    assert [(h,v) for h,v in zip(*again) if h.startswith('Comment[')] == [(h,v) for h,v in zip(*original) if h.startswith('Comment[')]
    extract=next(step for step in retained['series']['assay_paths'][0]['steps'] if step['kind']=='extract')
    start=original[0].index('Extract Name')
    end=original[0].index('Assay Name')
    assert [(c['name'],c['value']) for c in extract['comments']] == [(h[8:-1],v) for h,v in list(zip(*original))[start:end] if h.startswith('Comment[')]
