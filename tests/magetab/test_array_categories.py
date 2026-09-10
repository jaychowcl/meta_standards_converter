# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""GEO platform categories route arrays without interpreting sample prose."""
from copy import deepcopy
import json
from pathlib import Path
from types import SimpleNamespace

import pytest
from meta_standards_converter.magetab.technology import _detect_base_technology, resolve_technology
from meta_standards_converter.magetab.constructor import AEConstructor
from meta_standards_converter.converters import GEO2AEConverter, JSON2AEConverter
from meta_standards_converter.miniml.geo_parser import GEOParser
from tests.support.contracts import install_replay

CATEGORIES = ('in situ oligonucleotide', 'spotted oligonucleotide',
              'mixed spotted oligonucleotide', 'spotted DNA/cDNA',
              'spotted peptide or protein', 'antibody', 'tissue', 'oligonucleotide beads')

@pytest.mark.parametrize('technology', CATEGORIES)
@pytest.mark.parametrize('normalize', [False, True])
def test_exact_geo_array_categories(technology, normalize):
    if normalize:
        technology = '  ' + '\t  '.join(technology.upper().split()) + ' \n'
    assert _detect_base_technology({'platform': [{'technology': technology}]}) == 'array'

@pytest.mark.parametrize('technology', ['MS', 'RT-PCR', 'other', 'antibody sequencing', 'tissue profiling'])
def test_non_array_categories(technology):
    assert _detect_base_technology({'platform': [{'technology': technology}]}) == 'generic'

@pytest.mark.parametrize('field', ['title', 'description', 'extract_protocol', 'characteristics'])
def test_sample_prose_cannot_supply_array_category(field):
    sample = {field: 'antibody tissue spotted DNA/cDNA'}
    assert _detect_base_technology({'sample': [sample]}) == 'generic'
    sample['type'] = 'SRA'
    assert _detect_base_technology({'sample': [sample]}) == 'bulk_sequencing'

@pytest.mark.parametrize('signal', [
    {'type': 'SRA'}, {'library_strategy': 'RNA-Seq'},
    {'relation': [{'type': 'SRA', 'target': 'SRX1'}]},
])
def test_sequencing_precedence(signal):
    assert _detect_base_technology({'platform': [{'technology': 'antibody'}], 'sample': [signal]}) == 'bulk_sequencing'

def test_sample_platform_filtering_in_mixed_study():
    data = {'platform': [{'iid': 'GPL1', 'technology': 'spotted DNA/cDNA'},
                         {'iid': 'GPL2', 'technology': 'high-throughput sequencing'}]}
    array = {'iid': 'GSM1', 'platform_ref': [{'ref': 'GPL1'}]}
    seq = {'iid': 'GSM2', 'platform_ref': [{'ref': 'GPL2'}]}
    data['sample'] = [array, seq]
    original = deepcopy(data)
    assert _detect_base_technology(data) == 'bulk_sequencing'
    assert resolve_technology(array, data=data).handler == 'array'
    assert resolve_technology(seq, data=data).handler == 'bulk_sequencing'
    assert data == original

def test_gse100_both_converters_match_forced_array(tmp_path, monkeypatch):
    raw = (Path(__file__).parents[1] / 'fixtures/studies/GSE100/inputs/geo.xml').read_bytes()
    install_replay(monkeypatch, documents={'GSE100': raw})
    package = GEOParser().parse(raw.decode(), remove_empty=True)[0]
    assert len(package.to_mapping()['sample']) == 4
    assert _detect_base_technology(package.to_mapping()) == 'array'
    source = tmp_path / 'GSE100.json'
    source.write_text(json.dumps(package.to_mapping()))
    # Metadata only: source is real; external enrichment is deliberately disabled.
    def constructor():
        return AEConstructor(evidence_resolver=SimpleNamespace(sample_runs=lambda *args: {}, publications=lambda data: []))
    geo = GEO2AEConverter(enricher=SimpleNamespace(enrich=lambda data: data), ae_constructor=constructor())
    js = JSON2AEConverter(ae_constructor=constructor())
    outputs = [geo.convert('GSE100'), geo.convert('GSE100', platform_handler='array'),
               js.convert(str(source), enrich=False), js.convert(str(source), enrich=False, platform_handler='array')]
    assert all(output == outputs[0] for output in outputs)
    rows = next(row[1] for row in outputs[0][0] if row[0] == 'SDRF File')
    assert len(rows) == 9
    assert len({row[0] for row in rows[1:]}) == 4
    assert 'Assay Name' in rows[0]
    assert {row[rows[0].index('Array Design REF')] for row in rows[1:]} == {'GPL221'}
    generic = js.convert(str(source), enrich=False, platform_handler='generic')
    assert 'Array Design REF' not in next(row[1] for row in generic[0] if row[0] == 'SDRF File')[0]
