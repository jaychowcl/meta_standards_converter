# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from copy import deepcopy

import pytest

from meta_standards_converter.magetab.chemistry import resolve_chemistry
from meta_standards_converter.magetab.technology import resolve_technology
from meta_standards_converter.magetab.sdrf.constructor import SDRFConstructor
from tests.magetab.test_scoped_routing import sample


CITE_PROTOCOL = ('Cord blood mononuclear cells (CBMCs) were isolated from cord blood as described in Breton et al. 2015. '
    'Library was constructed as described in the protocol of the manufacturer (10x Genomics Single Cell) '
    'or for Dropseq (Macosko et al., 2015) with the modifications described in Stoeckius et al. 2017')


@pytest.mark.parametrize('provider', ['ENA', 'SRA', 'GEO'])
@pytest.mark.parametrize('method', ['Dropseq', 'Drop-seq', '10x'])
@pytest.mark.parametrize('role', ['RNA', 'ADT'])
def test_cite_library_identity_controls_droplet_routing_and_manufacturer(provider, method, role):
    s = sample(f'{role}: CBMC_8K_13AB_{method}', CITE_PROTOCOL)
    data = {'source': {'format': provider}, 'sample': [s], 'series': {'title': 'CITE-seq'}}
    before = deepcopy(data)
    assert resolve_technology(s, data=data).handler == 'droplet_single_cell_sequencing'
    result = resolve_chemistry(s, series=data['series'])
    assert result.manufacturer == ('10x Genomics' if method == '10x' else None)
    assert result.versions == ()
    assert result.family is None
    rows = SDRFConstructor()._miniml2sdrf(data)
    if method != '10x':
        assert not any(value == '10x technology' for row in rows[1:] for value in row)
    assert data == before


@pytest.mark.parametrize('provider', ['ENA', 'SRA', 'GEO'])
@pytest.mark.parametrize('prose', ['Cells were purified by immunopanning plates into a single cell suspension.',
    'A single cell solution was sorted into purified cell populations.'])
def test_cortical_bulk_suspensions_do_not_establish_single_cell_sequencing(provider, prose):
    s = sample('Astrocyte 1', 'RNA was purified from 100ng total RNA using TruSeq RNA Sample Prep Kit.')
    data = {'source': {'format': provider}, 'sample': [s], 'series': {'summary': prose}}
    assert resolve_technology(s, data=data).handler == 'bulk_sequencing'


@pytest.mark.parametrize('title,protocol,expected', [
    ('scRNA-seq', '', 'single_cell_sequencing'),
    ('single-cell RNA sequencing', '', 'single_cell_sequencing'),
    ('snRNA-seq', '', 'single_cell_sequencing'),
    ('scRNA-seq', 'Single cells were sorted into individual 96-well plates for library preparation.', 'plate_single_cell_sequencing'),
    ('plate-based single-cell RNA-seq', '', 'plate_single_cell_sequencing'),
    ('scRNA-seq', 'Not a plate-based preparation.', 'single_cell_sequencing'),
    ('scRNA-seq', 'Cells were grown in 96-well culture plates.', 'single_cell_sequencing'),
    ('RNA: donor_Dropseq', '', 'droplet_single_cell_sequencing'),
    ('Chromium genomic library', '', 'bulk_sequencing'),
    ('10x RNA-seq', '', 'bulk_sequencing'),
])
def test_affirmative_preparation_format_required(title, protocol, expected):
    s = sample(title, protocol)
    assert resolve_technology(s).handler == expected


@pytest.mark.parametrize('title', ['Chromium WGS', 'Chromium genomic libraries', 'bulk RNA-seq', 'not single-cell sequencing'])
def test_explicit_bulk_identity_blocks_shared_alternative_methods(title):
    s = sample(title, CITE_PROTOCOL)
    data = {'sample': [s], 'series': {'title': 'Bulk and single-cell RNA sequencing'}}
    assert resolve_technology(s, data=data).handler == 'bulk_sequencing'
    assert resolve_chemistry(s).manufacturer is None


def test_run_method_identity_precedes_shared_sample_label():
    s = sample('CITE-seq_10x', CITE_PROTOCOL)
    drop = {'library_name': 'ADT_Dropseq'}
    tenx = {'library_name': 'RNA_10x'}
    assert resolve_chemistry(s, run=drop).manufacturer is None
    assert resolve_chemistry(s, run=tenx).manufacturer == '10x Genomics'
    assert resolve_technology(s, run=drop).handler == 'droplet_single_cell_sequencing'


def test_conflicting_method_identity_is_not_resolved_from_shared_protocol():
    s = sample('scRNA-seq_Dropseq_10x', CITE_PROTOCOL)
    decision = resolve_technology(s)
    assert decision.handler == 'single_cell_sequencing'
    assert decision.diagnostics
    chemistry = resolve_chemistry(s)
    assert chemistry.manufacturer is None
    assert chemistry.diagnostics


def test_modern_explicit_chemistry_is_unchanged():
    s = sample('scRNA-seq', "Chromium Single Cell 3' Library Kit v3.1", code='SC3Pv3HT')
    chemistry = resolve_chemistry(s)
    assert (chemistry.manufacturer, chemistry.family, chemistry.versions) == ('10x Genomics', '3 prime', ('3.1',))
    assert resolve_technology(s).handler == 'droplet_single_cell_sequencing'


def test_generic_vendor_mention_without_single_cell_context_has_no_isolation_attribute():
    result = resolve_chemistry(sample('Chromium genome', 'Chromium genomic DNA library preparation'))
    assert result.manufacturer is None
    assert ('single cell isolation', '10x technology') not in result.attributes


@pytest.mark.parametrize('text', ['Bulk RNA-seq. 10x single-cell RNA sequencing', '10x single-cell RNA sequencing. Bulk RNA-seq'])
def test_conflicting_identity_clauses_do_not_depend_on_sentence_order(text):
    s = sample(text)
    result = resolve_technology(s)
    assert result.handler in {'sequencing', 'single_cell_sequencing'}
    assert result.diagnostics
    assert resolve_chemistry(s).manufacturer is None


def test_run_assay_identity_beats_sample_bulk_even_without_named_format():
    s = sample('bulk RNA-seq')
    assert resolve_technology(s, run={'library_name': 'single-cell RNA-seq'}).handler == 'single_cell_sequencing'
    s['title'] = 'single-cell RNA-seq'
    assert resolve_technology(s, run={'library_name': 'bulk RNA-seq'}).handler == 'bulk_sequencing'


def test_explicit_generic_droplet_preparation_has_a_droplet_handler():
    s = sample('scRNA-seq', 'Droplet-based single-cell RNA sequencing was performed.')
    assert resolve_technology(s).handler == 'droplet_single_cell_sequencing'
    assert resolve_chemistry(s).manufacturer is None


def test_explicit_exclusion_blocks_lower_priority_chemistry():
    s = sample('scRNA-seq; not 10x', "Chromium Single Cell 3' Library Kit v3.1")
    assert resolve_technology(s).handler == 'single_cell_sequencing'
    result = resolve_chemistry(s)
    assert result.manufacturer is None
    assert result.versions == ()
    assert result.diagnostics


def test_saved_native_json_exports_offline_without_changing_input(tmp_path):
    import json
    from meta_standards_converter.magetab.constructor import AEConstructor
    from meta_standards_converter.converters.json2ae import JSON2AEConverter
    from tests.test_native_archive_enrichment import native

    class Forbidden:
        def __getattr__(self, name):
            raise AssertionError('No retrieval is allowed for saved native JSON with no enrichment')

    data = native().to_mapping()
    data['sample'][0]['title'] = 'RNA: CD8_high_Dropseq'
    data['sample'][0]['channel'][0]['extract_protocol'] = CITE_PROTOCOL
    path = tmp_path / 'saved.json'
    path.write_text(json.dumps(data))
    before = path.read_bytes()
    converter = JSON2AEConverter(enricher=Forbidden(), ae_constructor=AEConstructor(insdc_client=Forbidden()))
    first = converter.convert(str(path), enrich=False)
    assert converter.convert(str(path), enrich=False) == first
    assert path.read_bytes() == before
    assert resolve_technology(data['sample'][0], data=data).handler == 'droplet_single_cell_sequencing'
    table = next(r[1] for r in first[0] if r[0] == 'SDRF File')
    assert not any(cell == '10x technology' for row in table[1:] for cell in row)


@pytest.mark.parametrize('method', ['10x', 'Dropseq'])
def test_named_droplet_method_and_generic_format_are_compatible(method):
    s = sample(f'{method} droplet-based scRNA-seq')
    assert resolve_technology(s).handler == 'droplet_single_cell_sequencing'
    assert resolve_chemistry(s).manufacturer == ('10x Genomics' if method == '10x' else None)


def test_generic_droplet_format_keeps_explicit_vendor_exclusion():
    s = sample('scRNA-seq; droplet-based; not 10x', "Chromium Single Cell 3' Library Kit v3.1")
    assert resolve_technology(s).handler == 'droplet_single_cell_sequencing'
    assert resolve_chemistry(s).manufacturer is None


def test_shared_plate_and_tenx_alternative_is_unresolved():
    s = sample('scRNA-seq', "Plate-based single-cell RNA-seq or Chromium Single Cell 3' Library Kit v3.1")
    result = resolve_chemistry(s)
    assert result.manufacturer is None
    assert result.versions == ()
    assert result.diagnostics
    assert resolve_technology(s).handler == 'single_cell_sequencing'


@pytest.mark.parametrize('alternative', ['Plate-based single-cell RNA-seq', 'Dropseq libraries were prepared'])
def test_shared_alternative_preparations_are_scoped_across_sentences(alternative):
    s = sample('scRNA-seq', alternative + ". Chromium Single Cell 3' Library Kit v3.1")
    result = resolve_chemistry(s)
    assert result.manufacturer is None
    assert result.versions == ()
    assert result.diagnostics
    s['title'] = 'scRNA-seq_10x'
    assert resolve_chemistry(s).versions == ('3.1',)


def test_explicit_channel_chemistry_precedes_lower_bulk_title():
    s = sample('bulk RNA-seq', code='SC3Pv3HT')
    assert resolve_technology(s).handler == 'droplet_single_cell_sequencing'
    result = resolve_chemistry(s)
    assert (result.manufacturer, result.versions) == ('10x Genomics', ('3.1',))


def test_explicit_channel_chemistry_conflicts_with_same_scope_bulk_run():
    s = sample('RNA-seq', code='SC3Pv3HT')
    run = {'library_name': 'bulk RNA-seq'}
    result = resolve_technology(s, run=run)
    assert result.handler == 'sequencing'
    assert result.diagnostics
    chemistry = resolve_chemistry(s, run=run)
    assert chemistry.manufacturer is None
    assert chemistry.diagnostics


def test_unknown_chemistry_identifier_does_not_displace_bulk_evidence():
    s = sample('bulk RNA-seq', code='SC3Pv3HT-unknown')
    assert resolve_technology(s).handler == 'bulk_sequencing'
