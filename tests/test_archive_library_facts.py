# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from copy import deepcopy

from meta_standards_converter.miniml import MINiMLCodec
from meta_standards_converter.metadata.archive_enrichment import merge_archive_metadata
from tests.test_archive_fidelity_enrichment import enriched_native, workflow


def linked_layout(layout='PAIRED'):
    extra = workflow().to_mapping()
    step = next(s for s in extra['series']['assay_paths'][0]['steps'] if s['kind'] == 'extract')
    step.setdefault('comments', []).append({'name': 'LIBRARY_LAYOUT', 'value': layout})
    return extra


def native_layout(layout='SINGLE'):
    data = enriched_native().to_mapping()
    run = data['sample'][0]['sra_run'][0]
    run['library_layout'] = layout
    run['read_count'] = '123'
    if layout == 'PAIRED':
        run['nominal_length'] = '350'
        run['nominal_sdev'] = '30'
    return MINiMLCodec().decode(data).package


def test_accepted_library_layout_replaces_native_run_and_path_projection():
    data, issues = merge_archive_metadata(native_layout(), MINiMLCodec().decode(linked_layout()).package, prefer=True)
    data = data.to_mapping()
    assert not issues
    assert data['sample'][0]['sra_run'][0]['library_layout'] == 'PAIRED'
    assert data['sample'][0]['sra_run'][0]['read_count'] == '123'
    values = [c['value'] for p in data['series']['assay_paths'] for s in p['steps']
              for c in s.get('comments', []) if c['name'] == 'LIBRARY_LAYOUT']
    assert values and set(values) == {'PAIRED'}
    assert 'SINGLE' in str(data['extensions']['insdc'])


def test_selecting_single_layout_removes_incompatible_paired_geometry():
    result, issues = merge_archive_metadata(native_layout('PAIRED'), MINiMLCodec().decode(linked_layout('SINGLE')).package, prefer=True)
    run = result.to_mapping()['sample'][0]['sra_run'][0]
    assert run['library_layout'] == 'SINGLE'
    assert 'nominal_length' not in run and 'nominal_sdev' not in run


def test_rejected_ambiguous_workflows_cannot_change_run_layout():
    extra = linked_layout()
    second = deepcopy(extra['series']['assay_paths'][0])
    next(s for s in second['steps'] if s['kind'] == 'extract')['name'] = 'different material'
    extra['series']['assay_paths'].append(second)
    extra['sample'][0]['sra_run'][0]['library_layout'] = 'PAIRED'
    data, issues = merge_archive_metadata(native_layout(), MINiMLCodec().decode(extra).package, prefer=True)
    assert issues
    assert data.to_mapping()['sample'][0]['sra_run'][0]['library_layout'] == 'SINGLE'


def test_conflicting_run_and_path_facts_are_reported_without_arbitrary_choice():
    extra = linked_layout('PAIRED')
    extra['sample'][0]['sra_run'][0]['library_layout'] = 'SINGLE'
    data, issues = merge_archive_metadata(native_layout(), MINiMLCodec().decode(extra).package, prefer=True)
    assert any('library_layout' in issue for issue in issues)
    assert data.to_mapping()['sample'][0]['sra_run'][0]['library_layout'] == 'SINGLE'


def test_missing_layout_does_not_displace_informative_native_layout():
    result, _ = merge_archive_metadata(native_layout('PAIRED'), MINiMLCodec().decode(linked_layout('not provided')).package, prefer=True)
    assert result.to_mapping()['sample'][0]['sra_run'][0]['library_layout'] == 'PAIRED'


def test_explicit_geo_run_facts_apply_without_an_ae_workflow():
    extra = workflow().to_mapping()
    extra['source']['format'] = 'GEO'
    extra['series']['assay_paths'] = []
    extra['sample'][0]['sra_run'][0]['library_layout'] = 'PAIRED'
    result, _ = merge_archive_metadata(native_layout(), MINiMLCodec().decode(extra).package, prefer=True)
    assert result.to_mapping()['sample'][0]['sra_run'][0]['library_layout'] == 'PAIRED'


def test_saved_bound_ae_extract_facts_are_synchronized_without_input_mutation():
    from meta_standards_converter.miniml.archive_libraries import synchronize_library_facts
    data, _ = merge_archive_metadata(native_layout(), MINiMLCodec().decode(linked_layout()).package, prefer=True)
    data = data.to_mapping()
    data['sample'][0]['sra_run'][0]['library_layout'] = 'SINGLE'
    before = deepcopy(data)
    copy = deepcopy(data)
    synchronize_library_facts(copy)
    assert copy['sample'][0]['sra_run'][0]['library_layout'] == 'PAIRED'
    first = deepcopy(copy)
    synchronize_library_facts(copy)
    assert first == copy and data == before


def test_conflicting_layout_cannot_supply_geometry_to_retained_native_group():
    from meta_standards_converter.miniml.archive_libraries import resolve_library_facts
    native = {'library_layout': 'PAIRED', 'nominal_length': '350'}
    selected = resolve_library_facts(native, [{'library_layout':'SINGLE','nominal_length':'999'},
        {'library_layout':'PAIRED'}], [], 'SRR1')
    assert selected == native


def test_source_run_with_conflicting_experiment_cannot_supply_facts():
    extra = linked_layout()
    extra['sample'][0]['sra_run'][0].update(experiment='SRX999999', library_source='TRANSCRIPTOMIC')
    data, issues = merge_archive_metadata(native_layout(), MINiMLCodec().decode(extra).package, prefer=True)
    assert issues
    assert data.to_mapping()['sample'][0]['sra_run'][0]['library_source'] == native_layout().to_mapping()['sample'][0]['sra_run'][0]['library_source']


def test_replaced_insert_measurement_does_not_inherit_previous_standard_deviation():
    from meta_standards_converter.miniml.archive_libraries import resolve_library_facts
    result = resolve_library_facts({'library_layout':'PAIRED', 'nominal_length':'350', 'nominal_sdev':'30'},
                                   [{'library_layout':'PAIRED', 'nominal_length':'500'}], [], 'SRR1')
    assert result == {'library_layout':'PAIRED', 'nominal_length':'500'}


def test_shared_experiment_path_does_not_take_last_runs_layout():
    from meta_standards_converter.miniml.archive_libraries import synchronize_library_facts
    data = native_layout().to_mapping()
    sample = data['sample'][0]
    run = deepcopy(sample['sra_run'][0]); run.update(run='SRR222', library_layout='PAIRED')
    sample['sra_run'].append(run)
    for path in data['series']['assay_paths']:
        path['steps'] = [s for s in path['steps'] if s['kind'] not in ('scan','array_data_file')]
    reverse = deepcopy(data); reverse['sample'][0]['sra_run'].reverse()
    issues = []
    synchronize_library_facts(data, issues); synchronize_library_facts(reverse, [])
    assert issues
    assert data['series']['assay_paths'] == reverse['series']['assay_paths']


def shared_library(native_values, field='library_layout', incoming='not provided'):
    data = native_layout().to_mapping()
    sample = data['sample'][0]
    sample['sra_run'] = [dict(deepcopy(sample['sra_run'][0]), run=f'SRR{111 + i}', **{field: value})
                         for i, value in enumerate(native_values)]
    data['series']['assay_paths'] = data['series']['assay_paths'][:1]
    for path in data['series']['assay_paths']:
        path['steps'] = [s for s in path['steps'] if s['kind'] not in ('scan', 'array_data_file')]
        for step in path['steps']:
            step['comments'] = [c for c in step.get('comments', []) if c['name'].lower() != field]
    extra = linked_layout(incoming)
    extra['sample'][0]['accession'] = deepcopy(sample['accession'])
    extra['sample'][0]['sra_run'] = []
    for path in extra['series']['assay_paths']:
        path['steps'] = [s for s in path['steps'] if s['kind'] != 'scan']
        for step in path['steps']:
            for comment in step.get('comments', []):
                if comment['name'] == 'LIBRARY_LAYOUT': comment['name'] = field.upper()
    return data, extra


def test_missing_shared_workflow_fact_preserves_each_native_run_and_reports_disagreement():
    from meta_standards_converter.miniml.archive_libraries import synchronize_library_facts
    for field, values in [('library_layout', ['SINGLE', 'PAIRED']),
                          ('library_strategy', ['RNA-Seq', 'WGS']),
                          ('library_source', ['TRANSCRIPTOMIC', 'GENOMIC']),
                          ('library_selection', ['cDNA', 'RANDOM'])]:
        data, extra = shared_library(values, field)
        data['sample'][0]['sra_run'][1]['nominal_length'] = '350'
        for reverse in (False, True):
            current = deepcopy(data)
            if reverse: current['sample'][0]['sra_run'].reverse()
            original = deepcopy(current)
            result, issues = merge_archive_metadata(MINiMLCodec().decode(current).package,
                MINiMLCodec().decode(extra).package, prefer=True)
            result = result.to_mapping()
            assert {r['run']: r[field] for r in result['sample'][0]['sra_run']} == dict(zip(['SRR111','SRR112'], values))
            assert any('shared workflow' in issue and field in issue for issue in issues)
            assert field not in result['sample'][0]
            extract_values = [c['value'] for p in result['series']['assay_paths'] for s in p['steps']
                              if s['kind'] == 'extract' for c in s.get('comments', []) if c['name'].lower() == field]
            assert not (set(extract_values) & set(values))
            if field == 'library_layout':
                assert next(r for r in result['sample'][0]['sra_run'] if r['run'] == 'SRR112')['nominal_length'] == '350'
            saved = deepcopy(result)
            synchronize_library_facts(result, [])
            assert result == saved
            assert current == original


def test_explicit_shared_workflow_fact_may_override_both_runs_without_stale_geometry():
    data, extra = shared_library(['SINGLE', 'PAIRED'], incoming='SINGLE')
    data['sample'][0]['sra_run'][1]['nominal_length'] = '350'
    result, issues = merge_archive_metadata(MINiMLCodec().decode(data).package,
        MINiMLCodec().decode(extra).package, prefer=True)
    runs = result.to_mapping()['sample'][0]['sra_run']
    assert {r['library_layout'] for r in runs} == {'SINGLE'}
    assert all('nominal_length' not in r for r in runs)
    assert not issues


def test_missing_authored_fact_does_not_become_native_consensus_in_disguise():
    data, extra = shared_library(['PAIRED', 'PAIRED'])
    result, _ = merge_archive_metadata(MINiMLCodec().decode(data).package,
        MINiMLCodec().decode(extra).package, prefer=True)
    result = result.to_mapping()
    assert {r['library_layout'] for r in result['sample'][0]['sra_run']} == {'PAIRED'}
    assert all(c['value'] != 'PAIRED' for p in result['series']['assay_paths'] for s in p['steps']
               if s['kind'] == 'extract' for c in s.get('comments', []) if c['name'] == 'LIBRARY_LAYOUT')
