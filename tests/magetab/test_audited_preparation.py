from copy import deepcopy

import pytest

from meta_standards_converter.magetab.technology import resolve_technology
from meta_standards_converter.magetab.chemistry import resolve_chemistry
from tests.magetab.test_scoped_routing import sample


@pytest.mark.parametrize('identity', ['sci-RNA-seq', 'single cell transcriptional profiling by combinatorial indexing'])
def test_combinatorial_identity_is_single_cell_without_invented_format(identity):
    s = sample(identity)
    assert resolve_technology(s).handler == 'single_cell_sequencing'
    assert resolve_chemistry(s).manufacturer is None


def test_explicit_bulk_control_overrides_shared_single_cell_preparation():
    s = sample('E18.5, 200 cell bulk control', 'Single-cell RNA sequencing using Fluidigm C1.')
    result = resolve_technology(s)
    assert result.handler == 'bulk_sequencing'
    assert result.control_role == 'bulk control'


@pytest.mark.parametrize('title', ['no cell control', 'Negative control L139_EMPTY_93'])
def test_empty_control_retains_only_evidenced_preparation(title):
    s = sample(title)
    data = {'sample': [s], 'series': {'title': 'Single-cell RNA sequencing'}}
    before = deepcopy(data)
    result = resolve_technology(s, data=data)
    assert result.handler == 'sequencing'
    assert result.control_role == 'empty control'
    assert data == before
    s['channel'][0]['extract_protocol'] = 'Single-cell RNA sequencing libraries were prepared using Drop-seq.'
    result = resolve_technology(s, data=data)
    assert result.handler == 'droplet_single_cell_sequencing'
    assert result.control_role == 'empty control'


@pytest.mark.parametrize('protocol,expected', [
    ('Single mES cells were FACS sorted (BD Influx; BD Biosciences) to each well. '
     'Subsequent steps were performed as described in Smart-Seq2.', 'plate_single_cell_sequencing'),
    ('Smart-Seq2 library preparation on a Fluidigm C1 microfluidic device.', 'single_cell_sequencing'),
    ('Cells were cultured in 96-well plates before RNA isolation.', 'single_cell_sequencing'),
])
def test_plate_evidence_requires_cell_deposition_not_kit_or_culture(protocol, expected):
    assert resolve_technology(sample('scRNA-seq', protocol)).handler == expected


def test_rna_preparation_cannot_classify_genomic_chip_library():
    s = sample('H3K4me3 ChIP-seq', 'Single-cell RNA sequencing libraries were made with Chromium Single Cell 3\' Library Kit v3.1.')
    run = {'library_strategy': 'ChIP-Seq', 'library_source': 'GENOMIC'}
    original = deepcopy(s)
    decision = resolve_technology(s, run=run)
    assert decision.handler == 'sequencing'
    assert any(d.code == 'incompatible_preparation_evidence' for d in decision.diagnostics)
    assert resolve_chemistry(s, run=run).manufacturer is None
    assert s == original


def test_explicit_single_cell_chip_is_not_banned_by_genomic_library_source():
    s = sample('single-cell ChIP-seq', 'Individual cells were sorted into wells for ChIP library preparation.')
    run = {'library_strategy': 'ChIP-Seq', 'library_source': 'GENOMIC'}
    assert resolve_technology(s, run=run).handler == 'plate_single_cell_sequencing'
