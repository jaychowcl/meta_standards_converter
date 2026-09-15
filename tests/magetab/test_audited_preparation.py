# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from copy import deepcopy
import json
from pathlib import Path

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


def test_structured_bulk_control_identity_overrides_shared_preparation():
    s = sample('control 1', 'Drop-seq libraries were prepared.')
    s['channel'][0]['characteristics'] = [{'name': 'control type', 'value': '200 cell bulk control'}]
    decision = resolve_technology(s)
    assert decision.handler == 'bulk_sequencing' and decision.control_role == 'bulk control'


@pytest.mark.parametrize('case', json.loads((Path(__file__).parents[1]/'fixtures/audited_preparations.json').read_text()))
def test_reduced_actual_archive_preparation_scopes(case):
    data = case['data']; s = data['sample'][0]
    before = deepcopy(data)
    result = resolve_technology(s, run=s['sra_run'][0], data=data)
    assert result.handler == case['expected']
    assert data == before


def test_empty_control_uses_uniquely_bound_library_preparation():
    s = sample('no cell control'); run = {'run':'SRR1', 'experiment':'SRX1'}
    series = {'protocols':[{'name':'P1','type':{'value':'library construction protocol'},
                           'description':'Single-cell RNA sequencing libraries were prepared using Drop-seq.'}],
              'assay_paths':[{'steps':[{'kind':'source','name':'GSM1','sample_ref':'GSM1'},
                {'kind':'protocol_application','protocol_ref':'P1'},
                {'kind':'assay','name':'SRX1','sample_ref':'GSM1'}, {'kind':'scan','name':'SRR1'}]}]}
    result = resolve_technology(s, run=run, data={'series':series})
    assert result.handler == 'droplet_single_cell_sequencing' and result.control_role == 'empty control'
    series['protocols'][0]['description'] = ('Single-cell RNA sequencing used Drop-seq. '
        'A no cell control was separately prepared in a PCR tube without any cell input.')
    assert resolve_technology(s, run=run, data={'series':series}).handler == 'sequencing'
def test_bound_protocol_index_is_operation_local_and_shared(monkeypatch):
    from meta_standards_converter.magetab import preparation
    from copy import deepcopy
    series = {'protocols':[{'name':'p', 'type':{'value':'library preparation'}, 'description':'Dropseq'}],
              'assay_paths':[{'steps':[{'kind':'source','name':'s','sample_ref':'s'},
                {'kind':'assay','name':'e'}, {'kind':'scan','name':'r'}, {'kind':'protocol_application','protocol_ref':'p'}]}]}
    original = deepcopy(series)
    builds = []
    builder = preparation._bound_protocol_index
    def counted(value):
        builds.append(value); return builder(value)
    monkeypatch.setattr(preparation, '_bound_protocol_index', counted)
    @preparation.preparation_operation
    def operation():
        for _ in range(10):
            assert preparation.bound_protocols({'iid':'s'}, {'run':'r','experiment':'e'}, series)
    operation()
    assert len(builds) == 1
    operation()
    assert len(builds) == 2
    assert series == original


@pytest.mark.parametrize('protocol,expected', [
    ('Single mES cells were FACS sorted (BD Influx; BD Biosciences) to each well. Subsequent steps used Smart-Seq2.', 'plate_single_cell_sequencing'),
    ('Smart-Seq2 library preparation on a Fluidigm C1 microfluidic device.', 'single_cell_sequencing'),
    ('Single mES cells were sorted into wells. Other libraries used Drop-seq.', 'single_cell_sequencing'),
])
def test_scoped_preparation_specializes_single_cell_source_without_title_signal(protocol, expected):
    s = sample('sample 17', protocol)
    run = {'library_strategy':'RNA-Seq', 'library_source':'TRANSCRIPTOMIC SINGLE CELL'}
    original = deepcopy(s)
    assert resolve_technology(s, run=run).handler == expected
    assert s == original


def test_bound_plate_preparation_can_specialize_source_signal_without_shared_study_override():
    s = sample('sample 17'); run = {'run':'SRR1','experiment':'SRX1','library_source':'TRANSCRIPTOMIC SINGLE CELL'}
    data = {'series':{'title':'A study also discusses Drop-seq',
        'protocols':[{'name':'P1','type':{'value':'nucleic acid extraction protocol'},
                     'description':'Single mES cells were FACS sorted to each well.'}],
        'assay_paths':[{'steps':[{'kind':'source','name':s['iid'],'sample_ref':s['iid']},
            {'kind':'protocol_application','protocol_ref':'P1'}, {'kind':'assay','name':'SRX1'}, {'kind':'scan','name':'SRR1'}]}]}}
    assert resolve_technology(s, run=run, data=data).handler == 'plate_single_cell_sequencing'


def test_conflicting_bound_preparation_formats_report_ambiguity():
    s = sample('sample 17', 'Single cells were FACS sorted into wells. These libraries used Drop-seq.')
    result = resolve_technology(s, run={'library_source':'TRANSCRIPTOMIC SINGLE CELL'})
    assert result.handler == 'single_cell_sequencing'
    assert any(d.code == 'ambiguous_preparation' for d in result.diagnostics)
