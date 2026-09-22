# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Operation counts guard scalability without machine-dependent timing assertions."""
from copy import deepcopy
from meta_standards_converter.metadata.modality import resolve_modality
from meta_standards_converter.metadata.preparation import preparation_operation
from meta_standards_converter.miniml import MINiMLCodec
from meta_standards_converter.miniml.model import MINiMLPackage
from tests.test_native_archive_enrichment import native


def test_one_package_materialization_per_outer_operation(monkeypatch):
    data = native().to_mapping()
    sample = data['sample'][0]
    first = sample['sra_run'][0]
    sample['sra_run'] = [{**deepcopy(first), 'run': f'SRR{i}', 'description': 'Visium spatial sequencing'} for i in range(1, 10)]
    package = MINiMLCodec().decode(data).package
    original = MINiMLPackage.to_mapping
    calls = []
    def counted(self):
        calls.append(self)
        return original(self)
    monkeypatch.setattr(MINiMLPackage, 'to_mapping', counted)
    @preparation_operation
    def project():
        a = resolve_modality(sample, data=package)
        b = resolve_modality(sample, data=package)
        assert a == b
        return a
    assert project().value == 'spatial'
    assert len(calls) == 1
    project()
    assert len(calls) == 2  # A second conversion must own fresh state.


def test_nested_resolvers_share_protocol_index_and_keep_run_scope(monkeypatch):
    from meta_standards_converter.metadata import preparation
    original = preparation._bound_protocol_index
    calls = []
    def counted(series):
        calls.append(series)
        return original(series)
    monkeypatch.setattr(preparation, '_bound_protocol_index', counted)
    sample = {'iid': 's1', 'sra_run': [{'run': 'SRR1', 'experiment': 'SRX1'}, {'run': 'SRR2', 'experiment': 'SRX1'}]}
    data = {'series': {'protocols': [
        {'name': 'P1', 'type': {'value': 'library construction'}, 'description': 'Visium spatial sequencing'},
        {'name': 'P2', 'type': {'value': 'library construction'}, 'description': 'bulk RNA sequencing'}],
        'assay_paths': [ {'steps': [{'kind': 'source', 'name': 's1', 'sample_ref': 's1'},
                                   {'kind': 'assay', 'name': 'SRX1'}, {'kind': 'scan', 'name': run},
                                   {'kind': 'protocol_application', 'protocol_ref': proto}]}
                         for run,proto in [('SRR1','P1'), ('SRR2','P2')]]}}
    @preparation_operation
    def project():
        return [resolve_modality(sample, data=data) for _ in range(2)]
    decisions = project()
    assert all(d.value == 'unknown' for d in decisions)  # Conflicting runs of one library.
    assert all(d.diagnostics for d in decisions)
    assert len(calls) == 1


def test_workflow_scope_extraction_is_linear_in_number_of_paths(monkeypatch):
    from meta_standards_converter.metadata import archive_workflows as workflows
    sample = {'iid': 's1', 'channel': [{'characteristics': []}], 'sra_run': [
        {'run': f'SRR{i}', 'experiment': 'SRX1'} for i in range(1, 31)]}
    paths = [{'steps': [{'kind': 'source', 'name': 's1', 'sample_ref': 's1'},
                        {'kind': 'assay', 'name': 'SRX1'}, {'kind': 'scan', 'name': r['run']}]} for r in sample['sra_run']]
    data = {'source': {'format': 'SRA'}, 'series': {'iid': 'SRP1', 'assay_paths': paths}, 'sample': [sample]}
    extra = deepcopy(data)
    extra['source']['format'] = 'MAGE-TAB'
    extra['series']['iid'] = 'E-MTAB-1'
    original = workflows.path_ids
    calls = []
    def counted(path):
        calls.append(path)
        return original(path)
    monkeypatch.setattr(workflows, 'path_ids', counted)
    issues = []
    workflows.merge_workflows(data, extra, {'s1': 's1'}, {}, True, issues)
    assert not issues
    assert {s['name'] for p in data['series']['assay_paths'] for s in p['steps'] if s['kind'] == 'scan'} == {f'SRR{i}' for i in range(1,31)}
    assert len(calls) <= 12 * len(paths)


def test_cached_native_source_records_do_not_materialize_package(monkeypatch):
    from meta_standards_converter.miniml.archive_residuals import source_records
    package = native()
    expected = package._archive_source_records
    def forbidden(self): raise AssertionError('Unneeded full package serialization')
    monkeypatch.setattr(MINiMLPackage, 'to_mapping', forbidden)
    assert source_records(package) == list(expected)


def test_residual_workflow_candidates_are_indexed_before_copying(monkeypatch):
    from meta_standards_converter.miniml import archive_workflow_residuals as residuals
    count = 30
    paths = [{'steps': [
        {'kind': 'source', 'name': 's1', 'sample_ref': 's1'},
        {'kind': 'assay', 'name': 'SRX1', 'sample_ref': 's1'},
        {'kind': 'scan', 'name': f'SRR{i}', 'sample_ref': 's1'},
        {'kind': 'array_data_file', 'name': f'{i}.fq', 'link': {'value': f'ftp://example.org/{i}.fq'}}
    ]} for i in range(1, count + 1)]
    source = {'series': {'assay_paths': paths}}
    target = {'source': {'format': 'SRA'}, 'series': {'iid': 'SRP1', 'assay_paths': deepcopy(paths)},
              'sample': [{'iid': 's1', 'sra_run': [{'run': f'SRR{i}', 'experiment': 'SRX1'} for i in range(1, count + 1)]}]}
    from meta_standards_converter.metadata.archive_workflows import mark_file_origins
    mark_file_origins({'source': {'format': 'MAGE-TAB'}, 'series': {'iid': 'E-MTAB-1', 'assay_paths': target['series']['assay_paths']}})
    original = residuals.deepcopy
    calls = []
    def counted(value):
        if isinstance(value, dict) and 'steps' in value:
            calls.append(value)
        return original(value)
    monkeypatch.setattr(residuals, 'deepcopy', counted)
    result = residuals.prune_bound_workflows(source, target, deepcopy(source), 'MAGE-TAB', 'E-MTAB-1')
    assert result['series']['assay_paths'] == []
    assert len(calls) <= 3 * count


def test_anndata_groups_share_one_operation_materialization(monkeypatch, tmp_path):
    from types import SimpleNamespace
    from meta_standards_converter.converters.json2h5ad import JSON2H5ADConverter
    from meta_standards_converter.metadata.preparation import package_mapping
    package = native()
    original = MINiMLPackage.to_mapping
    calls = []
    def counted(self):
        calls.append(self)
        return original(self)
    monkeypatch.setattr(MINiMLPackage, 'to_mapping', counted)
    converter = JSON2H5ADConverter()
    def project(packages, **kwargs):
        # Two sample projections inside the same expression conversion.
        for _ in range(2):
            @preparation_operation
            def sample_projection(): return package_mapping(packages[0])
            sample_projection()
        return None
    monkeypatch.setattr(converter, '_convert_packages', project)
    loaded = SimpleNamespace(warnings=[], groups=[SimpleNamespace(dataset_id='SRP1', packages=[package],
                              source_packages=[package], harmonization_resolution=None)])
    converter._convert_groups(loaded, source_json='input.json', out=str(tmp_path))
    assert len(calls) == 1
