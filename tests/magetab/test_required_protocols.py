# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Required structural protocols survive explicit MINiML workflow export."""
from copy import deepcopy
import json

import pytest

from meta_standards_converter.converters import Converter, PreparationPolicy
from meta_standards_converter.magetab.constructor import AEConstructor
from meta_standards_converter.miniml import MINiMLCodec
from tests.support.conversion import converted
from tests.support.magetab import import_tables


COLLECTION = 'sample collection protocol'
SEQUENCING = 'nucleic acid sequencing protocol'


def definition(name, kind, description='', **extra):
    return dict(name=name, type={'value': kind}, description=description, **extra)


def package(technologies=('sequencing',)):
    data = {'miniml_schema_version': '3.1', 'source': {'format': 'MAGE-TAB'},
            'series': {'iid': 'E-MTAB-1', 'title': 'Protocol coverage',
                       'protocols': [definition('P-E-MTAB-1-1', 'nucleic acid extraction protocol', 'Extract RNA.')],
                       'assay_paths': [], 'sample_ref': []}, 'sample': [], 'platform': []}
    for i, technology in enumerate(technologies, 1):
        sid = f'S{i}'
        data['series']['sample_ref'].append({'ref': sid})
        sample = {'iid': sid, 'channel': [{'source': {'value': sid}}], 'platform_ref': {'ref': f'PL{i}'}}
        if technology == 'sequencing':
            sample['library_strategy'] = 'RNA-Seq'
        data['sample'].append(sample)
        data['platform'].append({'iid': f'PL{i}', 'technology': {'sequencing': 'high-throughput sequencing',
                                   'array': 'expression array', 'generic': 'other'}[technology]})
        data['series']['assay_paths'].append({'steps': [
            {'kind': 'source', 'name': sid, 'sample_ref': sid},
            {'kind': 'protocol_application', 'protocol_ref': 'P-E-MTAB-1-1',
             'comments': [{'name': 'authored', 'value': 'keep attached'}]},
            {'kind': 'extract', 'name': f'E{i}'},
            {'kind': 'assay', 'name': f'A{i}'},
            {'kind': 'scan', 'name': f'R{i}'},
        ]})
    return data


def export(data, **options):
    constructor = AEConstructor()
    rows = constructor.miniml2magetab(MINiMLCodec().decode(data).package, **options)
    return rows, constructor.sdrf_constructor.last_sdrf_audit


def unpack(rows):
    values = {r[0]: r[1:] for r in rows}
    names = values['Protocol Name']
    definitions = {name: {label: cells[i] for label, cells in values.items()
                          if label.startswith('Protocol ') and i < len(cells)} for i, name in enumerate(names)}
    table = values['SDRF File'][0]
    paths = {row[0]: [definitions[value] for label, value in zip(table[0], row)
                     if label == 'Protocol REF' and value] for row in table[1:]}
    return definitions, table, paths


def assert_blank(record, kind):
    assert record['Protocol Type'] == kind
    assert record['Protocol Term Source REF'] == 'EFO'
    assert record['Protocol Term Accession Number'] == {'sample collection protocol': 'EFO_0005518',
                                                       'nucleic acid sequencing protocol': 'EFO_0004170'}[kind]
    for label in ('Protocol Description', 'Protocol Hardware', 'Protocol Software', 'Protocol Parameters',
                  'Protocol Contact', 'Protocol Performer'):
        assert record.get(label) in (None, '')


@pytest.mark.parametrize('present', [(), (COLLECTION,), (SEQUENCING,), (COLLECTION, SEQUENCING)])
@pytest.mark.parametrize('description', ['Authored method.', ''])
def test_explicit_lists_complete_missing_types_and_references(present, description):
    data = package()
    for i, kind in enumerate(present):
        data['series']['protocols'].append(definition(f'P-existing-{i}', kind, description))
        data['series']['assay_paths'][0]['steps'].insert(1, {'kind': 'protocol_application', 'protocol_ref': f'P-existing-{i}'})
    before = deepcopy(data)
    rows, _ = export(data)
    definitions, table, paths = unpack(rows)
    for kind in (COLLECTION, SEQUENCING):
        applied = [p for p in paths['S1'] if p['Protocol Type'] == kind]
        assert len(applied) == 1
        if kind not in present:
            assert_blank(applied[0], kind)
        else:
            assert applied[0]['Protocol Description'] == description
    assert len(definitions) == 3
    assert table[1][table[0].index('Comment[authored]')] == 'keep attached'
    assert rows == export(data)[0]
    assert data == before


def test_uncovered_workflows_share_blanks_without_borrowing_authored_procedures():
    data = package(('sequencing',) * 3)
    data['series']['protocols'] += [definition('P-collect', COLLECTION, 'Only sample one.'),
                                  definition('P-sequence', SEQUENCING, 'Only sample one.', hardware=['Machine A'])]
    first = data['series']['assay_paths'][0]['steps']
    first.insert(1, {'kind': 'protocol_application', 'protocol_ref': 'P-collect'})
    first.insert(-1, {'kind': 'protocol_application', 'protocol_ref': 'P-sequence'})
    rows, _ = export(data)
    _, _, paths = unpack(rows)
    for kind in (COLLECTION, SEQUENCING):
        by_sample = {sid: next(p for p in ps if p['Protocol Type'] == kind) for sid, ps in paths.items()}
        assert by_sample['S1']['Protocol Description'] == 'Only sample one.'
        assert_blank(by_sample['S2'], kind)
        assert by_sample['S2']['Protocol Name'] == by_sample['S3']['Protocol Name']
        assert by_sample['S1']['Protocol Name'] != by_sample['S2']['Protocol Name']


@pytest.mark.parametrize('accession', ['EFO:0005518', 'EFO_0005518', 'http://www.ebi.ac.uk/efo/EFO_0005518'])
def test_existing_blank_ontology_definition_is_reused_and_ids_do_not_collide(accession):
    data = package()
    blank = definition('P-E-MTAB-1-2', 'deposited label')
    blank['type'].update(term_source_ref='EFO', term_accession_number=accession)
    data['series']['protocols'].append(blank)
    rows, _ = export(data)
    definitions, _, paths = unpack(rows)
    assert [p for p in paths['S1'] if p['Protocol Name'] == blank['name']]
    assert len(definitions) == 3
    assert definitions[blank['name']]['Protocol Type'] == 'deposited label'
    assert definitions[blank['name']]['Protocol Term Accession Number'] == accession


@pytest.mark.parametrize('technology', ['array', 'generic', 'sequencing'])
def test_selected_technology_limits_sequencing_completion(technology):
    rows, _ = export(package((technology,)))
    _, _, paths = unpack(rows)
    types = [p['Protocol Type'] for p in paths['S1']]
    assert COLLECTION in types
    assert (SEQUENCING in types) == (technology == 'sequencing')


def test_mixed_scopes_and_forced_handler():
    data = package(('sequencing', 'array', 'generic'))
    _, _, paths = unpack(export(data)[0])
    assert [any(p['Protocol Type'] == SEQUENCING for p in paths[sid]) for sid in ('S1', 'S2', 'S3')] == [True, False, False]
    _, _, forced = unpack(export(data, platform_handler='bulk_sequencing')[0])
    assert all(any(p['Protocol Type'] == SEQUENCING for p in ps) for ps in forced.values())
    _, _, forced = unpack(export(data, platform_handler='generic')[0])
    assert not any(p['Protocol Type'] == SEQUENCING for ps in forced.values() for p in ps)


@pytest.mark.parametrize('last_node', ['scan', 'array_data_file', 'array_data_matrix_file', 'assay', 'extract', None])
def test_sparse_boundaries_and_audit(last_node):
    data = package()
    steps = [{'kind': 'source', 'name': 'S1', 'sample_ref': 'S1'}]
    if last_node:
        steps.append({'kind': last_node, 'name': 'reads.fastq' if last_node == 'array_data_file' else 'node'})
    data['series']['assay_paths'][0]['steps'] = steps
    rows, audit = export(data)
    definitions, table, paths = unpack(rows)
    assert {p['Protocol Type'] for p in definitions.values()} >= {COLLECTION, SEQUENCING}
    assert any(p['Protocol Type'] == COLLECTION for p in paths['S1'])
    if last_node:
        sequence_ref = next(p['Protocol Name'] for p in paths['S1'] if p['Protocol Type'] == SEQUENCING)
        index = table[1].index(sequence_ref)
        boundary = table[1].index('reads.fastq' if last_node == 'array_data_file' else 'node')
        assert (index < boundary) == (last_node in {'scan', 'array_data_file', 'array_data_matrix_file'})
    else:
        assert not any(p['Protocol Type'] == SEQUENCING for p in paths['S1'])
        assert any('nucleic acid sequencing' in w and 'boundary' in w for w in audit.warnings)


def test_json2ae_writes_blanks_and_roundtrip_is_idempotent(tmp_path):
    data = package()
    source = tmp_path / 'input.json'
    source.write_text(json.dumps(data))
    before = source.read_bytes()
    tables = converted(Converter(), str(source), in_type='json', out_type='magetab', outdir=str(tmp_path / 'out'),
                       preparation=PreparationPolicy(enrichment='off', expand_studies=False, linked_metadata=False, peer_provider=False))
    rows = tables[0]
    definitions, _, paths = unpack(rows)
    assert {p['Protocol Type'] for p in paths['S1']} >= {COLLECTION, SEQUENCING}
    for p in definitions.values():
        if p['Protocol Type'] in {COLLECTION, SEQUENCING}:
            assert_blank(p, p['Protocol Type'])
    idf_path, = (tmp_path / 'out').glob('*.idf.txt')
    parsed, = converted(Converter(), str(idf_path), out_type='json',
                        preparation=PreparationPolicy(enrichment='off', expand_studies=False, linked_metadata=False, peer_provider=False))
    again, _ = export(parsed)
    after_definitions, _, after_paths = unpack(again)
    assert definitions == after_definitions
    assert paths == after_paths
    assert source.read_bytes() == before
    assert list((tmp_path / 'out').glob('*.idf.txt'))
    assert list((tmp_path / 'out').glob('*.sdrf.txt'))


def test_empty_explicit_list_is_completed_without_modifying_input():
    data = package()
    data['series']['protocols'] = []
    data['series']['assay_paths'][0]['steps'].pop(1)
    before = deepcopy(data)
    definitions, _, paths = unpack(export(data)[0])
    assert len(definitions) == 2
    assert {p['Protocol Type'] for p in paths['S1']} == {COLLECTION, SEQUENCING}
    assert data == before


def test_generated_paths_and_explicit_definitions_share_declared_identities():
    data = package()
    data['series']['assay_paths'] = []
    data['sample'][0]['channel'][0]['extract_protocol'] = 'Extract RNA.'
    before = deepcopy(data)
    definitions, _, paths = unpack(export(data)[0])
    assert {p['Protocol Type'] for p in paths['S1']} >= {COLLECTION, SEQUENCING}
    assert definitions['P-E-MTAB-1-1']['Protocol Description'] == 'Extract RNA.'
    assert len([p for p in definitions.values() if p['Protocol Description'] == 'Extract RNA.']) == 1
    assert data == before


def test_multiple_authored_methods_and_hardware_only_method_are_not_borrowed():
    data = package(('sequencing',) * 3)
    data['series']['protocols'] += [definition('P-author-one', COLLECTION, 'First collection'),
                                  definition('P-author-two', COLLECTION, 'Second collection'),
                                  definition('P-hardware', SEQUENCING, hardware=['Authored instrument'])]
    for i, name in enumerate(('P-author-one', 'P-author-two')):
        data['series']['assay_paths'][i]['steps'].insert(1, {'kind': 'protocol_application', 'protocol_ref': name})
    data['series']['assay_paths'][0]['steps'].insert(-1, {'kind': 'protocol_application', 'protocol_ref': 'P-hardware'})
    definitions, _, paths = unpack(export(data)[0])
    assert [p['Protocol Name'] for p in paths['S1'] if p['Protocol Type'] == COLLECTION] == ['P-author-one']
    assert [p['Protocol Name'] for p in paths['S2'] if p['Protocol Type'] == COLLECTION] == ['P-author-two']
    assert_blank(next(p for p in paths['S3'] if p['Protocol Type'] == COLLECTION), COLLECTION)
    assert_blank(next(p for p in paths['S2'] if p['Protocol Type'] == SEQUENCING), SEQUENCING)
    assert definitions['P-hardware']['Protocol Hardware'] == 'Authored instrument'


def test_collection_without_source_anchor_retains_definition_and_warns():
    data = package()
    data['series']['assay_paths'][0]['steps'] = [{'kind': 'assay', 'name': 'A1'}, {'kind': 'scan', 'name': 'R1'}]
    rows, audit = export(data)
    definitions, _, paths = unpack(rows)
    assert any(p['Protocol Type'] == COLLECTION for p in definitions.values())
    assert not any(p['Protocol Type'] == COLLECTION for p in paths['A1'])
    assert any('sample collection' in w and 'boundary' in w for w in audit.warnings)


def test_sequencing_uses_first_acquisition_boundary_when_raw_file_precedes_scan():
    data = package()
    data['series']['assay_paths'][0]['steps'].insert(-1, {'kind': 'array_data_file', 'name': 'reads.fastq'})
    definitions, table, _ = unpack(export(data)[0])
    sequence = next(name for name, p in definitions.items() if p['Protocol Type'] == SEQUENCING)
    assert table[1].index(sequence) < table[1].index('reads.fastq') < table[1].index('R1')


@pytest.mark.parametrize('scope', ['channel', 'run'])
def test_path_scoped_technology_does_not_leak_between_channels_or_runs(scope):
    data = package(('generic',))
    sample = data['sample'][0]
    first = data['series']['assay_paths'][0]
    second = deepcopy(first)
    second['steps'][-1]['name'] = 'R2'
    data['series']['assay_paths'].append(second)
    if scope == 'run':
        sample['sra_run'] = [{'run': 'R1', 'library_name': 'single-cell RNA-seq'},
                             {'run': 'R2', 'library_name': 'unclassified'}]
    else:
        sample['channel'] = [{'source': {'value': 'c1'}, 'extract_protocol': 'Single-cell RNA sequencing.'},
                             {'source': {'value': 'c2'}}]
        first['steps'][0]['name'] = 'c1'
        second['steps'][0]['name'] = 'c2'
    definitions, table, _ = unpack(export(data)[0])
    applied = {row[table[0].index('Scan Name')]: [definitions[v]['Protocol Type'] for h, v in zip(table[0], row)
                                               if h == 'Protocol REF' and v] for row in table[1:]}
    assert SEQUENCING in applied['R1']
    assert SEQUENCING not in applied['R2']


def test_multiple_documents_and_runs_keep_reference_identity():
    data = package()
    first = data['series']['assay_paths'][0]
    first['document'] = 'one.sdrf.txt'
    second = deepcopy(first)
    second['document'] = 'two.sdrf.txt'
    second['steps'][-1]['name'] = 'R2'
    data['series']['assay_paths'].append(second)
    definitions, table, _ = unpack(export(data)[0])
    assert len(table) == 3
    refs = [[v for h, v in zip(table[0], row) if h == 'Protocol REF' and v] for row in table[1:]]
    assert refs[0] == refs[1]
    assert all(ref in definitions for ref in refs[0])


def test_completion_preserves_different_authored_protocol_positions_across_rows():
    data = package(('sequencing', 'sequencing'))
    data['series']['protocols'].append(definition('P-authored-seq', SEQUENCING, 'Authored sequencing.'))
    data['series']['assay_paths'][0]['steps'].insert(1, {'kind': 'protocol_application', 'protocol_ref': 'P-authored-seq'})
    rows, _ = export(data)
    definitions, table, _ = unpack(rows)
    first, second = table[1:]
    assert first.index('P-authored-seq') < first.index('E1')
    sequence = next(name for name, p in definitions.items() if p['Protocol Type'] == SEQUENCING and not p['Protocol Description'])
    assert second.index('A2') < second.index(sequence) < second.index('R2')
    assert second.index('P-E-MTAB-1-1') < second.index('E2')
    parsed, _ = import_tables(rows)
    assert unpack(export(parsed.to_mapping())[0])[1] == table
