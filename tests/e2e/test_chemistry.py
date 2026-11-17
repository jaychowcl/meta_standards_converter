# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Both real converters parse and publish recorded chemistry in synthetic envelopes."""
from tests.support.conversion import converted
from meta_standards_converter.converters import Converter, PreparationPolicy
from meta_standards_converter.magetab.exporter import MAGETabExporter
from meta_standards_converter.sources.metadata import GEOMetadataSource, MAGETabMetadataSource
import csv
import io
import hashlib
import json
import re
from pathlib import Path
import pytest

from tests.support.contracts import install_replay

ROOT=Path(__file__).resolve().parents[1]/'fixtures/chemistry'
CASES=json.loads((ROOT/'manifest.json').read_text())

@pytest.mark.parametrize('case',CASES,ids=lambda c:c['id'])
def test_chemistry_complete_converter_outputs(case,workspace,monkeypatch):
    xml=(ROOT/(case['id']+'.xml')).read_bytes()
    assert hashlib.sha256(xml).hexdigest()==case['envelope']['sha256']
    accession=case['envelope']['gse']
    calls=install_replay(monkeypatch,documents={accession:xml})
    exporter = MAGETabExporter()
    converted(Converter(services={'magetab_exporter': exporter}), accession, out_type='magetab', preparation=PreparationPolicy(enrichment='standard', expand_studies=False, linked_metadata=False, peer_provider=False), outdir='geo')
    # The stored-source route passes through the real GEO parser and enrichment.
    converted(Converter(), accession, out_type='json', preparation=PreparationPolicy(enrichment='standard', expand_studies=False, linked_metadata=False, peer_provider=False), outdir='json')
    json_path=next((workspace/'json').glob('*.json'))
    converted(Converter(), str(json_path), out_type='magetab', preparation=PreparationPolicy(enrichment='off', expand_studies=False, linked_metadata=False, peer_provider=False), outdir='json_ae')
    assert calls==['geo:'+accession,'geo:'+accession]
    geo_sdrf=next((workspace/'geo').glob('*.sdrf.txt')).read_text()
    json_sdrf=next((workspace/'json_ae').glob('*.sdrf.txt')).read_text()
    assert geo_sdrf==json_sdrf
    assert geo_sdrf==(ROOT/(case['id']+'.sdrf.txt')).read_text()
    warnings=exporter.constructor.sdrf_constructor.last_sdrf_audit.warnings
    codes=[w.split(' chemistry ',1)[1].split(':',1)[0] for w in warnings if ' chemistry ' in w]
    assert codes == []  # GEO export retains prose without deriving chemistry claims.


    reparsed=converted(Converter(), str(next((workspace / 'geo').glob('*.idf.txt'))), out_type='json', preparation=PreparationPolicy(enrichment='off', expand_studies=False, linked_metadata=False, peer_provider=False))
    retained=reparsed[0].to_mapping()
    # Original protocol and every repeated SDRF column survive the real parser.
    fields=json.loads((ROOT/(case['id']+'.json')).read_text())
    # Only TSV transport characters change; ordinary spaces remain authored.
    idf_rows=list(csv.reader(next((workspace/'geo').glob('*.idf.txt')).open(),delimiter='\t'))
    descriptions=next(r[1:] for r in idf_rows if r[0]=='Protocol Description')
    published = re.sub(r'\r\n|[\t\r\n]', ' ', fields['protocol'])
    assert published in descriptions
    assert published in [c.get('extract_protocol') for s in retained['sample'] for c in s['channel']]
    from meta_standards_converter.miniml import MINiMLCodec
    path=workspace/'roundtrip.json'
    path.write_text(json.dumps(MINiMLCodec().encode_many(reparsed)))
    converted(Converter(), str(path), out_type='magetab', preparation=PreparationPolicy(enrichment='off', expand_studies=False, linked_metadata=False, peer_provider=False), outdir='roundtrip')
    original=list(csv.reader(io.StringIO(geo_sdrf),delimiter='\t'))
    again=list(csv.reader(next((workspace/'roundtrip').glob('*.sdrf.txt')).open(),delimiter='\t'))
    # The existing semantic renderer groups characteristics before comments.
    # Every comment value and repeated occurrence must survive in exact order.
    assert [(h,v) for h,v in zip(*again) if h.startswith('Comment[')] == [(h,v) for h,v in zip(*original) if h.startswith('Comment[')]
    extract=next(step for step in retained['series']['assay_paths'][0]['steps'] if step['kind']=='extract')
    start=original[0].index('Extract Name')
    end=original[0].index('Assay Name')
    # Material Type precedes these comments, so it owns every occurrence.
    assert original[0][start+1] == 'Material Type'
    assert extract['material_type']['value'] == original[1][start+1]
    assert not extract.get('comments')
    assert [(c['name'],c['value']) for c in extract['material_type'].get('comments', [])] == [(h[8:-1],v) for h,v in list(zip(*original))[start:end] if h.startswith('Comment[')]
