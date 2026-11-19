# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Description ordering is local to a material and preserves complete fields."""
from copy import deepcopy

import pytest

from meta_standards_converter.magetab.sdrf.model import SDRFAttr, SDRFNode, SDRFPath
from meta_standards_converter.magetab.sdrf.renderer import SDRFRenderer
from meta_standards_converter.magetab.semantics import render_miniml_assay_documents
from tests.support.deposited_fidelity import assert_deposited_sdrf
from tests.support.magetab import parse


@pytest.mark.parametrize('preserve_order', [False, True])
def test_description_moves_with_companions_within_each_material(preserve_order):
    materials = ['Source Name', 'Sample Name', 'Extract Name', 'Labeled Extract Name']
    parts, header, values = [], [], []
    for index, kind in enumerate(materials):
        parts.append(SDRFNode(kind, kind, f'm{index}', [
            SDRFAttr('Comment[node]', f'node {index}'),
            SDRFAttr('Characteristics[dose]', '0'),
            SDRFAttr('Material Type', f'type {index}', [
                SDRFAttr('Term Source REF', 'EFO'),
                SDRFAttr('Term Accession Number', f'EFO:{index}'),
                SDRFAttr('Comment[material]', f'material {index}')]),
            SDRFAttr('Label', 'label'),
            SDRFAttr('Description', f' description {index} ', [
                SDRFAttr('Comment[description]', f'detail {index}')]),
        ]))
        header.extend([kind, 'Comment[node]', 'Characteristics[dose]',
                       'Description', 'Comment[description]', 'Material Type',
                       'Term Source REF', 'Term Accession Number', 'Comment[material]', 'Label'])
        values.extend([f'm{index}', f'node {index}', '0', f' description {index} ',
                       f'detail {index}', f'type {index}', 'EFO', f'EFO:{index}',
                       f'material {index}', 'label'])
    parts.extend([
        SDRFNode('Assay Name', 'a', 'a1', [
            SDRFAttr('Description', 'assay description'),
            SDRFAttr('Factor Value[dose]', '0', [SDRFAttr('Unit', 'hour')])]),
        SDRFNode('Array Data File', 'f', 'reads.fastq'),
    ])
    header.extend(['Assay Name', 'Description', 'Array Data File', 'Factor Value[dose]', 'Unit'])
    values.extend(['a1', 'assay description', 'reads.fastq', '0', 'hour'])
    paths = [SDRFPath(parts)]
    renderer = SDRFRenderer(preserve_order=preserve_order)
    columns = renderer.plan_columns(paths)
    original = deepcopy((paths, columns))

    assert renderer.render_paths(columns, paths) == [header, values]
    assert (paths, columns) == original
    assert renderer.render_paths(columns, paths) == [header, values]


@pytest.mark.parametrize('preserve_order', [False, True])
def test_description_discovered_in_later_row_precedes_merged_material(preserve_order):
    paths = [
        SDRFPath([SDRFNode('Source Name', 's1', 's1', [
            SDRFAttr('Material Type', 'tissue'), SDRFAttr('Label', 'first')])]),
        SDRFPath([SDRFNode('Source Name', 's2', 's2', [
            SDRFAttr('Material Type', 'RNA'), SDRFAttr('Label', 'second'),
            SDRFAttr('Description', ' later description ')])]),
        SDRFPath([SDRFNode('Source Name', 's3', 's3', [
            SDRFAttr('Description', 'description without material')])]),
    ]
    renderer = SDRFRenderer(preserve_order=preserve_order)
    columns = renderer.plan_columns(paths)
    original = deepcopy((paths, columns))
    assert renderer.render_paths(columns, paths) == [
        ['Source Name', 'Description', 'Material Type', 'Label'],
        ['s1', '', 'tissue', 'first'],
        ['s2', ' later description ', 'RNA', 'second'],
        ['s3', 'description without material', '', ''],
    ]
    assert (paths, columns) == original


@pytest.mark.parametrize('preserve_order', [False, True])
@pytest.mark.parametrize('fields', [[], ['Description'], ['Material Type'],
                                   ['Description', 'Characteristics[note]', 'Material Type']])
def test_missing_field_or_existing_order_is_unchanged(preserve_order, fields):
    attrs = [SDRFAttr('Comment[Sample_description]', 'GEO description')]
    attrs.extend(SDRFAttr(field, f'value {index}') for index, field in enumerate(fields))
    paths = [SDRFPath([SDRFNode('Source Name', 's', 's1', attrs)])]
    renderer = SDRFRenderer(preserve_order=preserve_order)
    assert renderer.render_paths(renderer.plan_columns(paths), paths) == [
        ['Source Name', 'Comment[Sample_description]', *fields],
        ['s1', 'GEO description', *(f'value {index}' for index in range(len(fields)))],
    ]


@pytest.mark.parametrize('preserve_order', [False, True])
def test_repeated_descriptions_keep_occurrences_before_first_material(preserve_order):
    paths = [SDRFPath([SDRFNode('Source Name', 's', 's1', [
        SDRFAttr('Description', 'first'), SDRFAttr('Provider', 'provider'),
        SDRFAttr('Material Type', 'RNA'), SDRFAttr('Label', 'label'),
        SDRFAttr('Description', 'second'), SDRFAttr('Material Type', 'DNA'),
        SDRFAttr('Description', 'third'),
    ])])]
    renderer = SDRFRenderer(preserve_order=preserve_order)
    assert renderer.render_paths(renderer.plan_columns(paths), paths) == [
        ['Source Name', 'Description', 'Provider', 'Description', 'Description',
         'Material Type', 'Label', 'Material Type'],
        ['s1', 'first', 'provider', 'second', 'third', 'RNA', 'label', 'DNA'],
    ]


@pytest.mark.parametrize('preserve_order', [False, True])
def test_semantic_export_and_reimport_preserve_material_owners(preserve_order):
    header = ['Source Name', 'Comment[source]', 'Material Type', 'Term Source REF',
              'Term Accession Number', 'Comment[material]', 'Description',
              'Sample Name', 'Material Type', 'Description',
              'Extract Name', 'Material Type', 'Comment[material]', 'Description', 'Assay Name']
    row = ['s1', 'source note', 'tissue', 'EFO', 'EFO:1', 'source material', 'source description',
           'sample1', 'cells', 'sample description',
           'extract1', 'RNA', 'extract material', 'extract description', 'a1']
    source = parse(header, [row]).to_mapping()['series']['assay_paths']
    original = deepcopy(source)
    table = render_miniml_assay_documents(source, preserve_order=preserve_order)['x.sdrf.txt']
    for owner in ['Source Name', 'Sample Name', 'Extract Name']:
        start = table[0].index(owner)
        assert table[0].index('Description', start) < table[0].index('Material Type', start)
    assert_deposited_sdrf([header, row], table, accession='synthetic')
    assert parse(table[0], table[1:]).to_mapping()['series']['assay_paths'] == original
    assert source == original


def test_plate_handler_uses_same_description_order():
    from meta_standards_converter.magetab.sdrf.constructor import SDRFConstructor
    data = {
        'source': {'format': 'GEO MINiML'},
        'series': {'iid': 'GSE1', 'sample_ref': [{'ref': 'GSM1'}]},
        'sample': [{'iid': 'GSM1', 'description': 'plate sample',
                    'channel': [{'source': 'cells', 'molecule': 'total RNA'}]}],
    }
    table = SDRFConstructor()._miniml2sdrf(data, technology_type='plate_single_cell_sequencing')
    header = table[0]
    assert header.index('Description') < header.index('Material Type')
    assert table[1][header.index('Description')] == 'plate sample'
    assert table[1][header.index('Comment[Sample_description]')] == 'plate sample'
