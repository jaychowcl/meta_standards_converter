# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
import json
from pathlib import Path

import anndata
import pandas
from scipy import sparse

from meta_standards_converter.converters import JSONDataOutputOrchestrator
from meta_standards_converter.converters.json2h5ad import (
    BatchConversionResult,
    ConversionResult,
)


def _source(tmp_path: Path) -> tuple[Path, Path]:
    h5ad = tmp_path / "source.h5ad"
    adata = anndata.AnnData(
        X=sparse.csr_matrix([[1, 0], [0, 2]]),
        obs=pandas.DataFrame(
            {"author_cluster": ["alpha", "beta"]},
            index=["GSM1-cell-1", "GSM1-cell-2"],
        ),
        var=pandas.DataFrame(
            {"gene_symbol": ["A", "B"]}, index=["ENSG1", "ENSG2"]
        ),
    )
    adata.uns["source_note"] = "kept"
    adata.write_h5ad(h5ad)
    payload = {
        "miniml_schema_version": "2.0",
        "source": {"format": "test"},
        "series": {"accession": [{"value": "GSE1"}]},
        "sample": [
            {
                "iid": "GSM1",
                "accession": [{"value": "GSM1"}],
                "title": "sample one",
            }
        ],
    }
    source = tmp_path / "GSE1.json"
    source.write_text(json.dumps(payload), encoding="utf-8")
    return source, h5ad


def test_orchestrator_exports_combined_obs_and_optional_metadata_sidecars(tmp_path):
    source, h5ad = _source(tmp_path)
    outdir = tmp_path / "components"

    result = JSONDataOutputOrchestrator().export_anndata_metadata(
        str(source),
        outdir=str(outdir),
        asset_specs=[f"GSM1={h5ad}"],
        include_var=True,
        include_uns=True,
    )

    obs = pandas.read_csv(result.obs_path)
    var = pandas.read_csv(result.var_path)
    uns = json.loads(Path(result.uns_path).read_text(encoding="utf-8"))
    manifest = json.loads(Path(result.manifest_path).read_text(encoding="utf-8"))

    assert obs["cell_id"].tolist() == ["GSM1-cell-1", "GSM1-cell-2"]
    assert obs["author_cluster"].tolist() == ["alpha", "beta"]
    assert obs["msc.sample.accession"].tolist() == ["GSM1", "GSM1"]
    assert var["feature_id"].tolist() == ["ENSG1", "ENSG2"]
    assert uns["values"]["msc_metadata"]["sample_values"]["type"] == "dataframe"
    assert manifest["operation"] == "anndata_metadata"
    assert manifest["status"] == "complete"
    assert result.obs.shape[0] == 2
    assert result.var.shape[0] == 2
    assert result.uns["msc_metadata"]["schema_version"] == "1.0"


def test_orchestrator_obs_export_omits_unrequested_sidecars(tmp_path):
    source, h5ad = _source(tmp_path)

    result = JSONDataOutputOrchestrator().export_anndata_metadata(
        str(source),
        outdir=str(tmp_path / "obs-only"),
        asset_specs=[f"GSM1={h5ad}"],
    )

    assert Path(result.obs_path).is_file()
    assert result.var is None and result.var_path is None
    assert result.uns is None and result.uns_path is None


def test_orchestrator_aggregates_obs_metadata_without_combining_expression(tmp_path):
    sources = []
    samples = []
    for sample_id in ("GSM1", "GSM2"):
        h5ad = tmp_path / f"{sample_id}.h5ad"
        anndata.AnnData(
            X=sparse.csr_matrix([[1]]),
            obs=pandas.DataFrame(index=[f"{sample_id}-cell"]),
            var=pandas.DataFrame(index=[f"{sample_id}-feature"]),
        ).write_h5ad(h5ad)
        sources.append(f"{sample_id}={h5ad}")
        samples.append(
            {
                "iid": sample_id,
                "accession": [{"value": sample_id}],
                "title": sample_id,
            }
        )
    source = tmp_path / "GSE1.json"
    source.write_text(
        json.dumps(
            {
                "miniml_schema_version": "2.0",
                "source": {"format": "test"},
                "series": {"accession": [{"value": "GSE1"}]},
                "sample": samples,
            }
        ),
        encoding="utf-8",
    )

    result = JSONDataOutputOrchestrator().export_anndata_metadata(
        source,
        outdir=tmp_path / "obs-catalogue",
        asset_specs=sources,
    )

    assert result.obs.index.tolist() == ["GSM1-cell", "GSM2-cell"]
    assert result.obs["msc.sample.accession"].tolist() == ["GSM1", "GSM2"]
    assert result.failures == ()


def test_orchestrator_manifest_writes_selected_format_and_json_summary(tmp_path):
    source, _h5ad = _source(tmp_path)

    result = JSONDataOutputOrchestrator().export_manifest(
        source,
        outdir=tmp_path / "manifest",
        output_format="csv",
    )

    assert Path(result.output_path).suffix == ".csv"
    assert Path(result.manifest_path).is_file()
    assert json.loads(Path(result.manifest_path).read_text())["artifacts"][
        "table"
    ] == result.output_path


def test_batch_obs_keeps_dataset_directory_when_other_group_fails(
    tmp_path, monkeypatch
):
    class FakeH5ADConverter:
        def convert(self, *_args, **_kwargs):
            return BatchConversionResult(
                conversions={"GSE1": ConversionResult("GSE1")},
                failures=["GSE2: conversion failed"],
            )

    orchestrator = JSONDataOutputOrchestrator(h5ad_converter=FakeH5ADConverter())
    targets = []

    def export_components(_conversion, destination, **_options):
        targets.append(destination)
        return object()

    monkeypatch.setattr(orchestrator, "_export_components", export_components)

    result = orchestrator.export_anndata_metadata(
        tmp_path / "atlas.json", outdir=tmp_path / "outputs"
    )

    assert targets == [tmp_path / "outputs" / "GSE1"]
    assert result.partial
