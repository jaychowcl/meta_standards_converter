# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""File and object APIs must share projection without temporary JSON files."""
import copy
import json

from meta_standards_converter.sources.json import JSONPackageSource
from meta_standards_converter.converters.json2ae import JSON2AEConverter
from meta_standards_converter.converters.json2tsv import JSON2TSVConverter
from tests.converters.test_json2tsv import package, atlas_v1_payload, atlas_v1_dataset


def test_decode_accepts_packages_and_envelopes_without_mutating(tmp_path):
    source = JSONPackageSource()
    for value in (package(), [package()], {'miniml_json': package()},
                  atlas_v1_payload([atlas_v1_dataset('GSE1', package())])):
        original = copy.deepcopy(value)
        path = tmp_path / 'source.json'
        path.write_text(json.dumps(value))
        loaded = source.decode(value, fallback='source')
        assert loaded == source.load(path)
        assert value == original
        assert source.decode(loaded.groups[0].packages).groups[0].packages == loaded.groups[0].packages


def test_loaded_table_matches_file_api_and_has_no_writes(tmp_path):
    source = tmp_path / 'source.json'
    source.write_text(json.dumps(package()))
    converter = JSON2TSVConverter()
    loaded = converter.package_source.load(source)
    table = converter.project_loaded(loaded)
    assert len(table.rows) == 1
    assert table.rows[0]['msc.sample.accession'] == 'GSM1'
    assert sorted(p.name for p in tmp_path.iterdir()) == ['source.json']
    a = converter.convert_source(source, tmp_path / 'a.tsv')
    b = converter.convert_loaded(loaded, tmp_path / 'b.tsv')
    assert (tmp_path / 'a.tsv').read_bytes() == (tmp_path / 'b.tsv').read_bytes()
    assert a.columns == b.columns == table.columns


def test_loaded_magetab_retains_injected_collaborators(tmp_path):
    calls = []
    class Constructor:
        def miniml2magetab(self, data, **kwargs):
            calls.append((data, kwargs))
            return [['Investigation Title', data.series.title]]
    converter = JSON2AEConverter(ae_constructor=Constructor())
    loaded = JSONPackageSource().decode(package())
    result = converter.convert_loaded(loaded, enrich=False, platform_handler='generic')
    assert result == [[['Investigation Title', loaded.groups[0].packages[0].series.title]]]
    assert calls[0][1] == {'platform_handler': 'generic'}
    assert not list(tmp_path.iterdir())
