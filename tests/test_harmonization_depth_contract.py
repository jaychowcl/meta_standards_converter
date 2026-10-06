# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================

"""XSD-derived, compatibility-first model for MINiML-compatible JSON."""
"""Depth-aware output and immutable legacy patch compatibility."""
from copy import deepcopy
import json
from pathlib import Path

import pytest

from meta_standards_converter.miniml import (
    HarmonizedValue, MINiMLCodec, MINiMLHarmonizationPatch,
    append_harmonized_value, apply_miniml_harmonization_patch,
    canonical_miniml_document, miniml_source_fingerprint,
    iter_harmonized_values,
)
from meta_standards_converter.magetab.harmonized import columns, read_group


def legacy():
    return json.loads((Path(__file__).parent / 'fixtures/contracts/legacy-depth-patch.json').read_text())


def test_depths_have_independent_collision_groups_and_preserve_unknown():
    dest = {'value': 'raw'}
    for label, depth in [('T cell', 0), ('lymphocyte', 1), ('immune cell', 1), ('cell', 2)]:
        append_harmonized_value(dest, HarmonizedValue('cell_type', label, hierarchy_depth=depth))
    assert dest['hz_cell_type'] == 'T cell'
    assert dest['hz_cell_type__ancestor_depth_1'] == 'lymphocyte'
    assert dest['hz_cell_type__ancestor_depth_1(1)'] == 'immune cell'
    assert dest['hz_cell_type__ancestor_depth_2'] == 'cell'
    assert dest['hz_cell_type__ancestor_depth_1_hierarchy_depth(1)'] == 1
    assert HarmonizedValue('tissue', 'lung').to_mapping() == {'hz_tissue': 'lung'}
    assert [(v.field, v.hierarchy_depth, v.index) for v in iter_harmonized_values(dest)] == [
        ('cell_type', 0, 0), ('cell_type', 1, 0), ('cell_type', 1, 1), ('cell_type', 2, 0)]


def test_new_names_infer_depth_but_reject_conflicts_and_reserved_field():
    assert iter_harmonized_values({'hz_cell_type__ancestor_depth_2': 'cell'})[0].hierarchy_depth == 2
    with pytest.raises(ValueError, match='depth'):
        iter_harmonized_values({'hz_cell_type__ancestor_depth_2': 'cell', 'hz_cell_type__ancestor_depth_2_hierarchy_depth': 1})
    with pytest.raises(ValueError, match='reserved'):
        HarmonizedValue('cell_type__ancestor_depth_2', 'cell')


def test_legacy_hash_is_stable_and_new_writer_upgrades_without_mutation():
    fixture = legacy()
    source = deepcopy(fixture['source'])
    assert miniml_source_fingerprint(source, schema_version='3.0') == fixture['sha256']
    assert canonical_miniml_document(source, schema_version='3.0') == fixture['canonical']
    result = MINiMLCodec.encode(MINiMLCodec().decode(source).package)
    assert result['miniml_schema_version'] == '3.1'
    assert result['sample'][0]['channel'][0]['source']['hz_cell_type__ancestor_depth_1'] == 'lymphocyte'
    assert source == fixture['source']
    assert MINiMLCodec.encode(MINiMLCodec().decode(result).package) == result


def test_old_patch_verifies_before_upgrade_and_remains_idempotent():
    fixture = legacy()
    result = apply_miniml_harmonization_patch(fixture['source'], fixture['patch'])
    assert result['miniml_schema_version'] == '3.1'
    assert result['sample'][0]['channel'][0]['source']['hz_cell_type__ancestor_depth_2'] == 'leukocyte'
    assert apply_miniml_harmonization_patch(result, fixture['patch']) == result
    changed = deepcopy(fixture['source'])
    changed['sample'][0]['channel'][0]['source']['value'] = 'B cell'
    with pytest.raises(ValueError, match='fingerprint'):
        apply_miniml_harmonization_patch(changed, fixture['patch'])


def test_new_patch_binds_new_canonicalization_and_roundtrips_sdrf():
    fixture = legacy()
    patch = MINiMLHarmonizationPatch(base_sha256=miniml_source_fingerprint(fixture['source']), adds=tuple(fixture['patch']['adds']))
    assert (patch.schema_version, patch.miniml_schema_version) == ('3.2', '3.1')
    result = apply_miniml_harmonization_patch(fixture['source'], patch)
    values = result['sample'][0]['channel'][0]['source']
    pairs = columns(values)
    header, row = map(list, zip(*pairs))
    found = []
    index = 0
    while index < len(header):
        _, value, index = read_group(header, row, index)
        found.append(value)
    assert [(v.field, v.hierarchy_depth) for v in found] == [('cell_type', 0), ('cell_type', 1), ('cell_type', 2)]


def test_named_rows_upgrade_names_and_companions_together():
    source = legacy()['source']
    source['sample'][0]['channel'][0]['characteristics'] = [
        {'name': 'raw', 'value': 'unchanged'},
        {'name': 'hz_cell_type(1)', 'value': 'lymphocyte'},
        {'name': 'hz_cell_type_hierarchy_depth(1)', 'value': 1},
    ]
    result = canonical_miniml_document(source)
    rows = result['sample'][0]['channel'][0]['characteristics']
    assert rows[0] == {'name': 'raw', 'value': 'unchanged'}
    assert rows[1]['name'] == 'hz_cell_type__ancestor_depth_1'
    assert rows[2]['name'] == 'hz_cell_type__ancestor_depth_1_hierarchy_depth'
