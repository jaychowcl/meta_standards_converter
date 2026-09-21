"""Library-scoped biological modality, independently of renderer defaults."""
from collections import Counter
from copy import deepcopy
from pathlib import Path
import pytest
from meta_standards_converter.metadata.modality import resolve_modality
from meta_standards_converter.magetab.parser import AEParser
from meta_standards_converter.sources.magetab import AEWebFetcher


@pytest.mark.parametrize('text,expected', [
    ('Visium spatial gene expression', 'spatial'), ('GeoMx DSP RNA profiling', 'spatial'),
    ('single-nucleus ATAC-seq', 'single_nucleus'), ('single-cell ATAC-seq', 'single_cell'),
    ('SPLiT-seq combinatorial indexing', 'single_cell'), ('Evercode whole transcriptome', 'single_cell'),
    ('single-cell combinatorial indexing RNA sequencing', 'single_cell'),
    ('bulk RNA-seq control for a single-cell study', 'bulk'),
    ('10x Genomics', 'unknown'), ('Illumina RNA-seq', 'unknown'),
    ('single-cell and spatial RNA sequencing', 'unknown'),
    ('without single-cell preparation; bulk RNA libraries', 'bulk'),
])
def test_positive_scoped_evidence_and_counterexamples(text, expected):
    sample = {'iid': 's', 'title': text, 'library_strategy': 'RNA-Seq'}
    before = deepcopy(sample)
    decision = resolve_modality(sample)
    assert decision.value == expected
    assert sample == before
    if expected != 'unknown': assert decision.evidence


def test_mixed_requires_different_established_libraries_not_conflicting_claims():
    sample = {'iid': 's', 'sra_run': [
        {'run': 'r1', 'experiment': 'e1', 'library_name': 'bulk RNA-seq'},
        {'run': 'r2', 'experiment': 'e2', 'library_name': 'single-nucleus RNA-seq'}]}
    assert resolve_modality(sample).value == 'mixed'
    sample['sra_run'][1]['experiment'] = 'e1'
    result = resolve_modality(sample)
    assert result.value == 'unknown'
    assert result.diagnostics and result.diagnostics[0].paths
    sample['sra_run'][1]['library_name'] = 'unannotated'
    assert resolve_modality(sample).value == 'unknown'


@pytest.mark.parametrize('accession,expected', [
    ('E-MTAB-14560', {'single_cell': 4, 'spatial': 16}),
    ('E-MTAB-14566', {'spatial': 2}), ('E-MTAB-14960', {'spatial': 4}),
    ('E-MTAB-13131', {'single_nucleus': 1}), ('E-MTAB-16856', {'single_cell': 2}),
])
def test_recorded_current_technologies(accession, expected):
    path = Path(__file__).parents[1] / 'fixtures/fidelity' / (accession + '.idf.txt')
    value = AEParser().parse(AEWebFetcher().resolve(str(path))).to_mapping()
    assert Counter(resolve_modality(s, data=value).value for s in value['sample']) == expected


def test_run_identity_precedes_shared_protocol_and_preserves_evidence_path():
    sample = {'iid': 's', 'channel': [{'extract_protocol': 'Visium spatial gene expression'}],
              'sra_run': [{'run': 'r', 'library_name': 'Chromium single-cell RNA-seq'}]}
    result = resolve_modality(sample)
    assert result.value == 'single_cell'
    assert any('r' in fact.path and 'library_name' in fact.path for fact in result.evidence)


def test_tabular_and_anndata_share_modality_detail_without_changing_legacy_service(tmp_path):
    import anndata
    import numpy as np
    from meta_standards_converter.converters.json2h5ad import JSON2H5ADConverter
    from meta_standards_converter.converters.json2tsv import JSON2TSVConverter
    from meta_standards_converter.sources.json import JSONPackageSource
    from meta_standards_converter.expression.assets import Asset
    from tests.converters.test_json2tsv import package
    value = package()
    sample = value['sample'][0]
    sample['title'] = 'single-nucleus ATAC-seq'
    sample['library_strategy'] = 'ATAC-seq'
    table = JSON2TSVConverter().project_loaded(JSONPackageSource().decode(value))
    assert table.rows[0]['msc.expression.modality_detail'] == 'single_nucleus'
    assert table.rows[0]['msc.library.strategy'] == 'ATAC-seq'
    source = tmp_path / 'input.tsv'
    source.write_text('gene\tcell\na\t1\n')
    adata = anndata.AnnData(np.ones((1, 1)))
    JSON2H5ADConverter().normalizer.normalize(adata, sample, value, 'GSE1',
        Asset('GSM1', str(source), 'matrix'), [], tmp_path)
    assert adata.obs['msc.expression.modality_detail'].tolist() == ['single_nucleus']
    assert 'modality_evidence' in adata.uns['meta_standards_converter']
    adata.write_h5ad(tmp_path / 'test.h5ad')
    assert anndata.read_h5ad(tmp_path / 'test.h5ad').obs['msc.expression.modality_detail'].tolist() == ['single_nucleus']
