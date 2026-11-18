from copy import deepcopy

from meta_standards_converter.miniml import harmonization as hz


def test_companion_rewrite_groups_keys_once_and_preserves_values(monkeypatch):
    # Count parsing work, not machine-dependent timing. Companion lookup used
    # to scan all N fields for every value, making large annotations quadratic.
    document = {'raw': 'unchanged', 'extensions': {'hz_unparsed': 'preserved'}}
    for i in range(1, 101):
        document.update({f'hz_tissue({i})': f'value {i}',
                         f'hz_tissue_id({i})': f'UBERON:{i}',
                         f'hz_tissue_onto({i})': 'UBERON'})
    original = deepcopy(document)
    parse = hz.parse_harmonized_key_details
    calls = []
    def counted(key):
        calls.append(key)
        return parse(key)
    monkeypatch.setattr(hz, 'parse_harmonized_key_details', counted)
    assert hz.rewrite_harmonized_names(document) == original
    assert document == original
    assert len(calls) <= 6 * 300
