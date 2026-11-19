# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Source-derived regressions for blank characteristic annotation ownership."""
import csv
import json
from pathlib import Path

import pytest

from meta_standards_converter.converters import Converter, PreparationPolicy
from scripts.audit_recorded_provider import primitive

FIXTURE = Path(__file__).parents[1] / 'fixtures/edge_cases/blank-characteristic-comments'
FIELDS = ('Sample_description', 'Sample_characteristics', 'Sample_source_name')


def convert(tmp_path, route, *, populated=False):
    source = list(csv.reader((FIXTURE / 'E-GEOD-13398.sdrf.txt').open(), delimiter='\t'))
    if populated: source[1][source[0].index('Characteristics [site]')] = 'reference material'
    idf = tmp_path / 'E-GEOD-13398.idf.txt'
    idf.write_bytes((FIXTURE / idf.name).read_bytes())
    with (tmp_path / 'E-GEOD-13398.sdrf.txt').open('w', newline='') as stream:
        csv.writer(stream, delimiter='\t', lineterminator='\n').writerows(source)
    policy = PreparationPolicy(enrichment='off', expand_studies=False, linked_metadata=False, peer_provider=False)
    value = idf
    if route == 'via-json':
        result = Converter().convert(idf, out_type='json', preparation=policy)
        value = tmp_path / 'saved.json'
        value.write_text(json.dumps(primitive(result.items[0].payload[0])))
    result = Converter().convert(value, out_type='magetab', preparation=policy)
    document = result.items[0].payload[0]
    actual = next(row[1] for row in document if row[0] == 'SDRF File')
    return source, actual


def cell(rows, field, source_name):
    header = rows[0]
    column = next(i for i,h in enumerate(header) if h.replace(' ', '') == f'Comment[{field}]')
    row = next(r for r in rows[1:] if r[header.index('Source Name')] == source_name)
    return row[column]


@pytest.mark.parametrize('route', ['direct', 'via-json'])
@pytest.mark.parametrize('field', FIELDS)
def test_blank_characteristic_does_not_drop_supplied_channel_comments(tmp_path, route, field):
    source, actual = convert(tmp_path, route)
    name = source[1][source[0].index('Source Name')]
    expected = cell(source, field, name)
    assert expected.strip()
    assert cell(actual, field, name) == expected


@pytest.mark.parametrize('route', ['direct', 'via-json'])
def test_populated_characteristic_retains_the_same_comment_qualifiers(tmp_path, route):
    source, actual = convert(tmp_path, route, populated=True)
    name = source[1][source[0].index('Source Name')]
    for field in FIELDS: assert cell(actual, field, name) == cell(source, field, name)


def test_blank_owner_keeps_qualifiers_and_does_not_steal_adjacent_comments():
    from tests.support.magetab import parse, render, import_tables
    package = parse(
        ['Source Name', 'Characteristics[site]', 'Term Source REF', 'Term Accession Number',
         'Comment[note]', 'Characteristics[site]', 'Comment[note]', 'Assay Name'],
        [['S1', '', 'EFO', 'EFO:1', 'blank owner', 'lung', 'populated owner', 'A1']])
    for value in (package, import_tables(render(package))[0]):
        node = value.to_mapping()['series']['assay_paths'][0]['steps'][0]
        attributes = [c for c in node['characteristics'] if c['name'] == 'site']
        assert len(attributes) == 2
        assert attributes[0]['value'] == ''
        assert attributes[0]['term_accession_number'] == 'EFO:1'
        assert attributes[0]['comments'] == [{'name': 'note', 'value': 'blank owner'}]
        assert attributes[1]['comments'] == [{'name': 'note', 'value': 'populated owner'}]
