# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Source-backed chemistry contracts; no provider calls or inferred read recipes."""
import hashlib
import json
from pathlib import Path
import pytest
from meta_standards_converter.magetab.chemistry import resolve_chemistry

FIXTURES = Path(__file__).resolve().parents[1] / 'fixtures/chemistry'
CASES = json.loads((FIXTURES / 'manifest.json').read_text())

def resolve(protocol, **kwargs):
    return resolve_chemistry({'iid':'GSM1', 'channel':[{'extract_protocol':protocol}]}, **kwargs)

@pytest.mark.parametrize('case', CASES, ids=lambda c:c['id'])
def test_recorded_evidence(case):
    raw=(FIXTURES / (case['id']+'.json')).read_bytes()
    assert hashlib.sha256(raw).hexdigest()==case['sha256']
    fields=json.loads(raw)
    result=resolve_chemistry({'iid':case['id'], 'data_processing':fields['processing'], 'channel':[{'extract_protocol':fields['protocol']}]})
    assert result.family == case['family']
    assert list(result.versions)==case['versions']
    if case['id']=='GSM8288494':
        assert 'ambiguous_family' in [d.code for d in result.diagnostics]
    if case['family']:
        assert result.manufacturer=='10x Genomics'
        assert result.evidence
    if case['id'] in ('pbmc_5_v2','pbmc_3_v31'):
        attrs=result.attributes
        assert ('cdna read size','90') in attrs
        assert ('sample barcode read','index1') in attrs
        assert ('sample barcode read','index2') in attrs
        assert attrs.count(('sample barcode size','10'))==2
        assert ('umi barcode size','10' if case['id']=='pbmc_5_v2' else '12') in attrs
        assert not any('offset' in k for k,v in attrs)

@pytest.mark.parametrize('prime', ["5'", '5′', '5’', '5ʹ', 'five prime', '5-prime'])
def test_prime_and_version_alternatives(prime):
    r=resolve(f'Chromium Single Cell {prime} Library & Gel Bead Kit v.1.1/v2')
    assert r.family=='5 prime'
    assert r.versions==('1.1','2')
    assert 'ambiguous_version' in [d.code for d in r.diagnostics]

@pytest.mark.parametrize('text', ["10x Chromium v3", "Cell Ranger v3 processes 10x data", "Not prepared with Chromium Single Cell 5' Kit v2", "Compatible with Chromium Single Cell 5' Kit v2", "https://example.org/Chromium_Single_Cell_5_v2.fastq.gz"])
def test_no_unsupported_chemistry(text):
    r=resolve(text)
    assert r.family is None
    assert not any(k != "single cell isolation" for k,v in r.attributes)

def test_scope_and_read_roles():
    r=resolve("Chromium Single Cell 5' Kit v2", series={'summary':"Chromium Single Cell 3' Kit v3"}, run={'read_lengths':['26','90']})
    assert r.family=='5 prime'
    assert not any('read size' in k for k,v in r.attributes)
    sample={'iid':'GSM1','channel':[{'extract_protocol':"Chromium Single Cell 3' Kit v2"},{'extract_protocol':"Chromium Single Cell 5' Kit v3"}]}
    assert resolve_chemistry(sample, channel=sample['channel'][1]).family=='5 prime'

def test_explicit_shared_scope_only():
    assert resolve('',series={'overall_design':"Chromium Single Cell 5' Kit v2"}).family is None
    assert resolve('',series={'overall_design':"All libraries were prepared using Chromium Single Cell 5' Kit v2"}).family=='5 prime'

def test_library_specific_recipe_not_shared():
    protocol="Chromium Single Cell 5' Kit v2. Gene Expression libraries: Read 2: 90 cycles (transcript). V(D)J libraries: Read 2: 150 cycles (transcript)."
    r=resolve(protocol)
    assert ('cdna read size','90') not in r.attributes
    assert ('cdna read size','150') not in r.attributes
    sample={'iid':'GSM1','channel':[{'extract_protocol':protocol}]}
    r=resolve_chemistry(sample, run={'library_name':'Gene Expression'})
    assert ('cdna read size','90') in r.attributes

@pytest.mark.parametrize('handler', ['single_cell_sequencing','droplet_single_cell_sequencing','tenx_v2_droplet_single_cell_sequencing','tenx_v3_droplet_single_cell_sequencing'])
def test_forced_handlers_use_source_and_never_presets(handler):
    from meta_standards_converter.magetab.sdrf.constructor import SDRFConstructor
    sample={'iid':'GSM1','accession':[{'value':'GSM1'}], 'channel':[{'source':'source','extract_protocol':"Chromium Single Cell 5' Kit v2"}]}
    constructor=SDRFConstructor()
    rows=constructor._miniml2sdrf({'sample':[sample]},technology_type=handler)
    assert rows[1][rows[0].index('Comment[end bias]')]=='5 prime tag'
    for name in ('cdna read size','cell barcode size','sample barcode size','umi barcode size','primer','LIBRARY_STRAND'):
        assert 'Comment['+name+']' not in rows[0]

def test_two_sample_chemistries_remain_separate():
    from meta_standards_converter.magetab.sdrf.constructor import SDRFConstructor
    samples=[{'iid':f'GSM{i}', 'accession':[{'value':f'GSM{i}'}], 'channel':[{'source':'source','extract_protocol':f"Chromium Single Cell {prime}' Kit v2"}]} for i,prime in enumerate([3,5],1)]
    rows=SDRFConstructor()._miniml2sdrf({'sample':samples},technology_type='droplet_single_cell_sequencing')
    assert [r[rows[0].index('Comment[end bias]')] for r in rows[1:]]==['3 prime tag','5 prime tag']

def test_spatial_preserves_existing_rendering():
    from meta_standards_converter.magetab.sdrf.constructor import SDRFConstructor
    sample={'iid':'GSM1','description':'10x Visium spatial transcriptomics','channel':[{'source':'source'}]}
    rows=SDRFConstructor()._miniml2sdrf({'sample':[sample]},technology_type='spatial_sequencing')
    assert rows[1][rows[0].index('Comment[library construction]')]=='10x Visium'

@pytest.mark.parametrize('protocol', ["Chromium Single Cell 5' Kit v2; compatible with Single Cell 3' Kit v3", "Chromium Single Cell 5' Kit v2; not Single Cell 3' Kit v3"])
def test_excluded_clauses_do_not_override_preparation(protocol):
    assert resolve(protocol).family=='5 prime'


def test_explicit_offsets_and_conflicting_lengths():
    r=resolve("Chromium Single Cell 5' Kit v2\nCell barcode offset: 0\nUMI barcode offset: 16\nRead 2: 90 cycles (transcript)\nRead 2: 150 cycles (transcript)")
    assert ('cell barcode offset','0') in r.attributes
    assert ('umi barcode offset','16') in r.attributes
    assert not any(k=='cdna read size' for k,v in r.attributes)
    assert any(d.code=='ambiguous_cdna_read_size' for d in r.diagnostics)

def test_vendor_in_compatibility_clause_does_not_claim_other_kit():
    r=resolve("OtherVendor Single Cell 3' Kit v2; compatible with Chromium data")
    assert r.manufacturer is None
    assert r.family is None


def test_unknown_kit_reports_uncertainty():
    r=resolve('Chromium Single Cell unsupported library kit v2')
    assert r.manufacturer=='10x Genomics'
    assert r.family is None
    assert r.versions==()
    assert 'unresolved_family' in [d.code for d in r.diagnostics]


def test_explicit_read_equals_and_dual_indices():
    r=resolve("Chromium Single Cell 5' Kit v2\nR2=90 bp (cDNA)\ni5=10 cycles\ni7=10 cycles")
    assert ('cdna read size','90') in r.attributes
    assert r.index_configuration=='dual'
