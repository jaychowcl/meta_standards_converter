# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Reference-aware selection preserves experiments, roles and source positions."""
from copy import deepcopy
import json

import pytest

from meta_standards_converter.miniml import MINiMLCodec, MINiMLCompatibilityError
from meta_standards_converter.miniml.selection import (
    SAMPLE_SELECTION_EXTENSION_KEY, select_miniml_samples, sample_role,
    selection_content_fingerprint,
)


def mixed_package():
    def path(*ids):
        return {'steps': [{'kind': 'source', 'name': iid, 'sample_ref': iid} for iid in ids]}
    return {'miniml_schema_version': '3.1', 'source': {'format': 'test'},
        'series': {'iid': 'study', 'sample_ref': [{'ref': s} for s in 'ABCD'],
            'protocols': [{'name': 'authored', 'description': 'original'}],
            'assay_paths': [path('B', 'C'), path('A', 'B'), path('D'),
                            {'steps': [{'kind': 'assay', 'name': 'unbound'}]}],
            'variable': [{'factor': 'tissue', 'sample_ref': [{'ref': 'A'}, {'ref': 'D'}]},
                         {'factor': 'age', 'sample_ref': [{'ref': 'D'}]}],
            'repeats': [{'sample_ref': [{'ref': 'B'}, {'ref': 'D'}]}]},
        'sample': [{'iid': s, 'title': s} for s in 'ABCD']}


def resolve(value, pointer):
    for part in pointer.strip('/').split('/'):
        part = part.replace('~1', '/').replace('~0', '~')
        value = value[int(part)] if isinstance(value, list) else value[part]
    return value


def test_selection_closure_subsets_provenance_and_unchanged_source():
    source = mixed_package(); before = deepcopy(source)
    selected = select_miniml_samples(source, ['A'])
    result = selected.package.to_mapping()
    assert source == before
    assert selected.analysis_ids == ('A',)
    assert selected.context_ids == ('B', 'C')
    assert selected.excluded_ids == ('D',)
    assert [s['iid'] for s in result['sample']] == list('ABC')
    assert result['series']['assay_paths'] == [source['series']['assay_paths'][i] for i in (0, 1, 3)]
    assert result['series']['protocols'] == source['series']['protocols']
    assert result['series']['variable'] == [{'factor': 'tissue', 'sample_ref': [{'ref': 'A'}]}]
    assert result['series']['repeats'] == [{'sample_ref': [{'ref': 'B'}]}]
    MINiMLCodec().decode(result, strict=True)
    for output, original in selected.selected_to_source.items():
        assert resolve(result, output) == resolve(source, original)
    assert selected.selected_to_source['/series/assay_paths/2/steps/0/name'] == '/series/assay_paths/3/steps/0/name'
    marker = result['extensions'][SAMPLE_SELECTION_EXTENSION_KEY]
    assert marker['selected_sha256'] == selection_content_fingerprint(result)
    assert marker['source_sha256'] == selection_content_fingerprint(source)
    decoded = MINiMLCodec().decode(json.loads(json.dumps(result)), strict=True).package
    assert sample_role(decoded, 'A') == 'analysis'
    assert sample_role(decoded, 'B') == 'context'
    assert sample_role(source, 'D') == 'analysis'
    assert select_miniml_samples(decoded, ['A']).package == decoded
    with pytest.raises(ValueError, match='context'):
        select_miniml_samples(decoded, ['B'])


def test_empty_selection_keeps_only_unbound_study_paths():
    selected = select_miniml_samples(mixed_package(), [])
    result = selected.package.to_mapping()
    assert not selected.package.samples
    assert result['series']['assay_paths'] == [{'steps': [{'kind': 'assay', 'name': 'unbound'}]}]
    assert not result['series'].get('sample_ref')
    assert not result['series'].get('variable')
    assert not result['series'].get('repeats')
    assert selected.excluded_ids == tuple('ABCD')


def test_selection_does_not_hide_invalid_excluded_references():
    source = mixed_package()
    source['series']['assay_paths'][2]['steps'][0]['sample_ref'] = 'missing'
    with pytest.raises(MINiMLCompatibilityError):
        select_miniml_samples(source, ['A'])
    with pytest.raises(ValueError, match='unknown'):
        select_miniml_samples(mixed_package(), ['missing'])
