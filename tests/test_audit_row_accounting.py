"""Occurrence-aware audit: distinct repository file evidence is not deduplicated."""
from collections import Counter
from copy import deepcopy
import json
from pathlib import Path

import pytest

from tests.support.recorded_replay import (
    ReplayTransport, assert_occurrences, occurrence_counter, table_occurrences,
)


ROOT = Path(__file__).parent / 'fixtures/recorded_provider'
MANIFEST = ROOT / 'row-evidence.manifest.json'


def ledger(study):
    replay = ReplayTransport(MANIFEST)
    value = json.loads(replay.read(f'{study}.rows'))
    assert len(value['rows']) == value['source_rows']
    headers = [name for _, name in value['columns']]
    return headers, value['rows']


@pytest.mark.parametrize('study,count', [
    ('E-MTAB-16253', 8), ('E-MTAB-14566', 779), ('E-MTAB-14416', 52),
])
def test_doubled_rows_are_distinct_original_and_ena_file_occurrences(study, count):
    headers, rows = ledger(study)
    at = headers.index
    roles = Counter(row[at('Comment[FASTQ_ROLE]')] for row in rows)
    assert roles == {'Original': count, 'GENERATED_FILE': count}
    assert max(occurrence_counter(table_occurrences([headers, *rows])).values()) == 1
    identities = [[row[at(key)] for key in (
        'Source Name', 'Assay Name', 'Comment[ENA_RUN]', 'Comment[FASTQ_URI]'
    )] for row in rows]
    assert max(occurrence_counter(identities).values()) == 1
    md5 = {
        role: {row[at('Comment[FASTQ_MD5]')] for row in rows
               if row[at('Comment[FASTQ_ROLE]')] == role}
        for role in roles
    }
    assert '' not in md5['Original'] | md5['GENERATED_FILE']
    assert md5['Original'].isdisjoint(md5['GENERATED_FILE'])
    # Each run keeps both repository representations; equal run is not equal file.
    original = [row[at('Comment[ENA_RUN]')] for row in rows
                if row[at('Comment[FASTQ_ROLE]')] == 'Original']
    generated = [row[at('Comment[ENA_RUN]')] for row in rows
                 if row[at('Comment[FASTQ_ROLE]')] == 'GENERATED_FILE']
    assert_occurrences(original, generated)


@pytest.mark.parametrize('mutation', ['drop', 'duplicate', 'md5', 'uri', 'role', 'factor', 'zero'])
def test_audit_comparison_rejects_occurrence_and_scientific_cell_loss(mutation):
    headers, rows = ledger('E-GEOD-60450')
    expected = table_occurrences([headers, *rows])
    changed = deepcopy(expected)
    if mutation == 'drop':
        changed.pop()
    elif mutation == 'duplicate':
        changed.append(deepcopy(changed[0]))
    else:
        labels = {'md5': 'Comment[FASTQ_MD5]', 'uri': 'Comment[FASTQ_URI]',
                  'role': 'Comment[FASTQ_ROLE]', 'factor': 'Factor Value[developmental stage]',
                  'zero': 'Comment[FASTQ_BYTES]'}
        cell = next(cell for cell in changed[0] if cell[0] == labels[mutation])
        cell[2] = 0 if mutation == 'zero' else 'changed'
    with pytest.raises(AssertionError, match='lost occurrences'):
        assert_occurrences(expected, changed)


def test_audit_preserves_duplicate_columns_literal_blanks_and_duplicate_rows():
    table = [['Source Name', 'Factor Value[dose]', 'Factor Value[dose]', 'Comment[bytes]'],
             ['s1', ' 0 ', '', 0], ['s1', ' 0 ', '', 0]]
    expected = table_occurrences(table)
    assert expected[0][1:3] == [['Factor Value[dose]', 1, ' 0 '],
                              ['Factor Value[dose]', 2, '']]
    assert list(occurrence_counter(expected).values()) == [2]
    with pytest.raises(AssertionError):
        assert_occurrences(expected, expected[:1])
    changed = deepcopy(expected)
    changed[0][1][2] = '0'
    with pytest.raises(AssertionError):
        assert_occurrences(expected, changed)


@pytest.mark.parametrize('study,count,suffix,fastq_count', [
    ('E-GEOD-52564', 17, '.xls', 34), ('E-GEOD-60450', 12, '.txt', 12),
])
def test_geo_supplement_branches_are_distinct_from_authored_uncompressed_files(
    study, count, suffix, fastq_count,
):
    headers, rows = ledger(study)
    at = headers.index
    supplemental = [row for row in rows if row[at('Comment[SAMPLE_FILE_NAME]')]]
    authored = [row for row in rows if row[at('Derived Array Data File')]]
    assert len(rows) == fastq_count + count
    assert len(supplemental) == count
    assert len(authored) == fastq_count
    assert all(row[at('Comment[SAMPLE_FILE_NAME]')].endswith(suffix + '.gz')
               and not row[at('Comment[ENA_RUN]')] for row in supplemental)
    assert all(row[at('Derived Array Data File')].endswith(suffix) for row in authored)
    supplied_uris = {row[at('Comment[SAMPLE_FILE_URI]')] for row in supplemental}
    authored_uris = {row[at('Comment[File URI]')] for row in authored}
    assert supplied_uris.isdisjoint(authored_uris)
    assert len(supplied_uris) == len(authored_uris) == count
    assert {row[at('Source Name')] for row in supplemental} == {
        row[at('Source Name')] for row in authored}


@pytest.mark.parametrize('study', [
    'E-MTAB-16253', 'E-MTAB-14566', 'E-MTAB-14416', 'E-GEOD-52564', 'E-GEOD-60450',
])
def test_reduced_source_package_renders_exact_scoped_recorded_file_occurrences(study):
    from meta_standards_converter.miniml import MINiMLCodec
    from meta_standards_converter.magetab.constructor import AEConstructor

    replay = ReplayTransport(MANIFEST)
    metadata = replay.files[f'{study}.package']
    headers, rows = ledger(study)
    expected_rows = [row for row in rows
                     if row[headers.index('Source Name')] == metadata['selected_sample']
                     and row[headers.index('Comment[ENA_RUN]')] in ('', metadata['selected_run'])]
    expected = table_occurrences([headers, *expected_rows])
    with replay.offline():
        package = MINiMLCodec().decode(json.loads(replay.read(f'{study}.package')), strict=True).package
        rendered = AEConstructor().miniml2magetab(package)
    replay.assert_complete()
    actual_table = next(row[1] for row in rendered if row and row[0] == 'SDRF File')
    # The full study's wider table may have columns absent in the scoped sample.
    # Project by header occurrence, preserving repeated archive-file columns.
    columns = [(name, ordinal) for name, ordinal, _ in expected[0]]
    actual = []
    for row in table_occurrences(actual_table):
        cells = {(name, ordinal): value for name, ordinal, value in row}
        actual.append([[name, ordinal, cells.get((name, ordinal), '')]
                       for name, ordinal in columns])
    assert_occurrences(expected, actual)


@pytest.mark.parametrize('target', ['tsv', 'csv'])
def test_delimited_projection_materializes_each_typed_package_once(target, monkeypatch):
    from meta_standards_converter.converters import Converter
    from meta_standards_converter.converters.json2tsv import JSON2DelimitedConverter
    from meta_standards_converter.miniml import MINiMLPackage

    replay = ReplayTransport(MANIFEST)
    data = json.loads(replay.read('E-MTAB-14566.package'))
    project = JSON2DelimitedConverter.project_loaded
    materialize = MINiMLPackage.to_mapping
    active = False
    calls = []

    def counted(self):
        if active:
            calls.append(self)
        return materialize(self)

    def in_projection(self, *args, **kwargs):
        nonlocal active
        active = True
        try:
            return project(self, *args, **kwargs)
        finally:
            active = False

    monkeypatch.setattr(MINiMLPackage, 'to_mapping', counted)
    monkeypatch.setattr(JSON2DelimitedConverter, 'project_loaded', in_projection)
    with replay.offline():
        result = Converter().convert(data, out_type=target, enrichment='off',
                                     options={'expand_studies': False})
    replay.assert_complete()
    assert result.status == 'complete', result.to_dict()
    # Strict preparation reconstructs candidates; each distinct typed package
    # still materializes at most once within the projection operation.
    assert calls
    assert max(Counter(map(id, calls)).values()) == 1
    table = next(iter(result.items[0].payload.values()))
    assert len(table.rows) == 1
    assert table.rows[0]['msc.sample.accession'] == replay.files['E-MTAB-14566.package']['selected_sample']
