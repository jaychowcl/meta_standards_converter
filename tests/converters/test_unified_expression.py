# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Expression routes are lazy, explicit about processing, and preserve caller data."""

import json
from pathlib import Path

import anndata
import scanpy
import numpy as np
import pandas as pd
import pytest

from meta_standards_converter.converters import Converter, InputSpec
from meta_standards_converter.converters.json2h5ad import JSON2H5ADConverter
from meta_standards_converter.sources.json import JSONPackageSource
from tests.converters.test_json2tsv import package


def adata():
    return anndata.AnnData(
        np.array([[1.0, 2.0], [3.0, 4.0]]),
        obs=pd.DataFrame({"sample": ["a", "a"]}, index=["c1", "c2"]),
        var=pd.DataFrame(index=["g1", "g2"]),
    )


@pytest.mark.parametrize("file", [False, True])
def test_standalone_h5ad_obs_and_metadata_requirement(file, tmp_path):
    value = adata()
    if file:
        path = tmp_path / "input.h5ad"
        value.write_h5ad(path)
        value = path
    converted = Converter().convert(value, out_type="obs")
    assert converted.status == "complete", converted.to_dict()
    assert converted.items[0].payload["obs"].index.tolist() == ["c1", "c2"]
    assert Converter().convert(value, out_type="json").status == "failed"
    published = Converter().convert(
        value,
        out_type="obs",
        outdir=tmp_path / "out",
        output_options={"include_var": True, "include_uns": True},
    )
    assert published.status == "complete", published.to_dict()
    assert len([p for p in (tmp_path / "out").iterdir() if p.is_file()]) >= 4


def test_anndata_not_mutated_and_embedded_metadata_is_validated(tmp_path):
    value = adata()
    value.uns["msc_miniml"] = {
        "schema_version": "1.0",
        "packages_json": json.dumps([package()]),
    }
    original = value.copy()
    result = Converter().convert(value, out_type="json")
    assert result.status == "complete", result.to_dict()
    assert value.uns == original.uns and value.obs.equals(original.obs)
    output = tmp_path / "copy.h5ad"
    assert (
        Converter().convert(value, out_type="h5ad", outfile=output).status == "complete"
    )
    np.testing.assert_array_equal(anndata.read_h5ad(output).X, value.X)
    value.uns["msc_miniml"]["packages_json"] = "{}"
    assert Converter().convert(value, out_type="json").status == "failed"


def test_matrix_requires_orientation_and_preserves_axes(tmp_path):
    path = tmp_path / "counts.tsv"
    path.write_text("gene\tc1\tc2\ng1\t1\t3\ng2\t2\t4\n")
    assert Converter().convert(path, out_type="h5ad").status == "failed"
    result = Converter().convert(
        path, out_type="h5ad", input_options={"orientation": "genes-by-observations"}
    )
    assert result.status == "complete", result.to_dict()
    converted = result.items[0].payload
    assert converted.obs_names.tolist() == ["c1", "c2"]
    assert converted.var_names.tolist() == ["g1", "g2"]
    np.testing.assert_array_equal(converted.X.toarray(), adata().X)


def test_loaded_h5ad_conversion_does_not_require_source_json(tmp_path):
    path = tmp_path / "counts.tsv"
    path.write_text("gene\tc1\tc2\ng1\t1\t3\ng2\t2\t4\n")
    value = package()
    value["sample"][0]["supplementary_data"] = [{"value": str(path)}]
    loaded = JSONPackageSource().decode(value)
    converter = JSON2H5ADConverter(
        available_memory=lambda: 10**12, memory_estimator=lambda *_: 100
    )
    converted = converter.convert_loaded(
        loaded, out=str(tmp_path / "out"), matrix_orientation="genes-by-observations"
    )
    assert not converted.partial
    matrix = anndata.read_h5ad(converted.sample_h5ads["GSM1"])
    assert (
        json.loads(matrix.uns["msc_miniml"]["packages_json"])[0]["sample"][0]["iid"]
        == "GSM1"
    )
    assert not list(tmp_path.glob("*.json"))


def test_metadata_h5ad_requires_destination_and_known_assets(tmp_path):
    path = tmp_path / "counts.tsv"
    path.write_text("gene\tc1\ng1\t1\n")
    value = package()
    value["sample"][0]["supplementary_data"] = [{"value": str(path)}]
    converter = Converter(
        services={
            "json2h5ad": JSON2H5ADConverter(
                available_memory=lambda: 10**12, memory_estimator=lambda *_: 100
            )
        }
    )
    assert converter.convert(value, out_type="h5ad").status == "failed"
    result = converter.convert(
        value,
        out_type="h5ad",
        outdir=tmp_path / "out",
        output_options={"matrix_orientation": "genes-by-observations"},
    )
    assert result.status == "complete", result.to_dict()
    assert any(p.endswith(".h5ad") for p in result.items[0].artifacts.values())


def test_raw_assets_in_metadata_do_not_launch_pipeline(tmp_path):
    value = package()
    value["sample"][0]["supplementary_data"] = [
        {"value": "https://ftp.ebi.ac.uk/sample.fastq.gz"}
    ]

    class Never:
        def process(self, *args, **kwargs):
            raise AssertionError("pipeline must not start")

    converter = Converter(
        services={"json2h5ad": JSON2H5ADConverter(pipeline_runner=Never())}
    )
    result = converter.convert(value, out_type="h5ad", outdir=tmp_path / "out")
    assert result.status == "failed"
    assert any(d.code == "processing_required" for d in result.items[0].diagnostics)


def test_standalone_rejects_settings_it_cannot_apply():
    with pytest.raises(ValueError):
        Converter().convert(
            adata(), out_type="h5ad", output_options={"force_reprocess": True}
        )


def test_explicit_asset_injection_is_retained(tmp_path):
    from meta_standards_converter.expression.assets import Asset

    path = tmp_path / "counts.tsv"
    path.write_text("gene\tc1\ng1\t1\n")
    converter = Converter(
        services={
            "json2h5ad": JSON2H5ADConverter(
                available_memory=lambda: 10**12, memory_estimator=lambda *_: 100
            )
        }
    )
    result = converter.convert(
        package(),
        out_type="h5ad",
        outdir=tmp_path / "out",
        output_options={
            "explicit_assets": [
                Asset(
                    "GSM1",
                    str(path),
                    "matrix",
                    source="cli",
                    orientation="genes-by-observations",
                )
            ]
        },
    )
    assert result.status == "complete", result.to_dict()


def test_discovery_groups_prefixed_10x_companions(tmp_path):
    from scipy import sparse
    from scipy.io import mmwrite

    mmwrite(str(tmp_path / "sample_matrix.mtx"), sparse.coo_matrix([[1, 2], [3, 4]]))
    (tmp_path / "sample_features.tsv").write_text(
        "g1\tG1\tGene Expression\ng2\tG2\tGene Expression\n"
    )
    (tmp_path / "sample_barcodes.tsv").write_text("c1\nc2\n")
    result = Converter().convert(tmp_path, out_type="h5ad")
    assert len(result.items) == 1, result.to_dict()
    assert result.status == "complete", result.to_dict()
    assert result.items[0].payload.shape == (2, 2)


def test_explicit_processing_option_runs_injected_pipeline(tmp_path):
    from meta_standards_converter.expression.assets import Asset
    from meta_standards_converter.expression.catalogue import RawProcessingResult

    reads = tmp_path / "input.fastq"
    reads.write_text("@read\nAC\n+\n!!\n")
    matrix = tmp_path / "counts.tsv"
    matrix.write_text("gene\tc1\ng1\t1\n")
    calls = []

    class Runner:
        def process(self, assets, **options):
            calls.append((assets, options))
            return RawProcessingResult(
                {
                    "GSM1": Asset(
                        "GSM1",
                        str(matrix),
                        "matrix",
                        orientation="genes-by-observations",
                    )
                },
                [],
                [],
            )

    h5ad = JSON2H5ADConverter(
        pipeline_runner=Runner(),
        available_memory=lambda: 10**12,
        memory_estimator=lambda *_: 100,
    )
    result = Converter(services={"json2h5ad": h5ad}).convert(
        InputSpec(reads, metadata=package()),
        out_type="h5ad",
        outdir=tmp_path / "out",
        runtime_options={"allow_processing": True},
        output_options={"genome": "GRCh38"},
    )
    assert result.status == "complete", result.to_dict()
    assert len(calls) == 1 and calls[0][1]["genome"] == "GRCh38"


def test_failed_component_dataset_does_not_erase_successful_sibling(tmp_path):
    from meta_standards_converter.expression.components import AnnDataComponentExporter

    path = tmp_path / "counts.tsv"
    path.write_text("gene\tc1\ng1\t1\n")
    first = package()
    second = package("GSE2", "GSM2")
    for value in (first, second):
        value["sample"][0]["supplementary_data"] = [{"value": str(path)}]

    class Components:
        def export(self, conversion, *args, **kwargs):
            if conversion.study_accession == "GSE2":
                raise OSError("fake publication failure")
            return AnnDataComponentExporter().export(conversion, *args, **kwargs)

    h5ad = JSON2H5ADConverter(
        available_memory=lambda: 10**12, memory_estimator=lambda *_: 100
    )
    result = Converter(
        services={"json2h5ad": h5ad, "components": Components()}
    ).convert(
        [first, second],
        out_type="obs",
        outdir=tmp_path / "out",
        output_options={"matrix_orientation": "genes-by-observations"},
    )
    assert result.status == "partial", result.to_dict()
    assert (tmp_path / "out/GSE1/GSE1.obs.csv").exists()
    assert not (tmp_path / "out/GSE2/GSE2.obs.csv").exists()


def test_raw_handler_rejects_malformed_fastq_before_pipeline(tmp_path):
    reads = tmp_path / "bad.fastq"
    reads.write_text("not fastq\n")
    result = Converter().convert(
        InputSpec(reads, metadata=package()),
        out_type="h5ad",
        outdir=tmp_path / "out",
        runtime_options={"allow_processing": True},
    )
    assert result.status == "failed"
    assert any(d.code == "invalid_fastq" for d in result.items[0].diagnostics)
    assert not (tmp_path / "out").exists()


def test_loaded_obs_api_matches_saved_metadata(tmp_path):
    from meta_standards_converter.converters.json2obs import JSON2OBSConverter

    path = tmp_path / "counts.tsv"
    path.write_text("gene\tc1\ng1\t1\n")
    value = package()
    value["sample"][0]["supplementary_data"] = [{"value": str(path)}]
    h5ad = JSON2H5ADConverter(
        available_memory=lambda: 10**12, memory_estimator=lambda *_: 100
    )
    result = JSON2OBSConverter(h5ad_converter=h5ad).convert_loaded(
        JSONPackageSource().decode(value),
        outdir=tmp_path / "out",
        matrix_orientation="genes-by-observations",
    )
    assert not result.partial and result.obs.shape[0] == 1
    assert Path(result.obs_path).exists()


def test_compressed_h5ad_is_bounded_and_read_backed_for_obs(tmp_path):
    import gzip

    path = tmp_path / "data.h5ad"
    adata().write_h5ad(path)
    compressed = tmp_path / "data.h5ad.gz"
    with gzip.open(compressed, "wb") as stream:
        stream.write(path.read_bytes())
    result = Converter().convert(compressed, out_type="obs")
    assert result.status == "complete", result.to_dict()
    rejected = Converter().convert(
        compressed,
        out_type="obs",
        runtime_options={"resource_overrides": {"max_expanded_archive_bytes": 100}},
    )
    assert rejected.status == "failed"
    assert any(d.code == "input_too_large" for d in rejected.items[0].diagnostics)


def test_mixed_tenx_directory_keeps_other_inputs_and_rejects_ambiguous_features(
    tmp_path,
):
    from scipy import sparse
    from scipy.io import mmwrite

    mmwrite(str(tmp_path / "matrix.mtx"), sparse.coo_matrix([[1, 2], [3, 4]]))
    (tmp_path / "features.tsv").write_text(
        "g1\tG1\tGene Expression\ng2\tG2\tGene Expression\n"
    )
    (tmp_path / "barcodes.tsv").write_text("c1\nc2\n")
    (tmp_path / "study.json").write_text(json.dumps(package()))
    result = Converter().convert(tmp_path, out_type="json")
    assert len(result.items) == 2, result.to_dict()
    assert any(item.dataset_ids == ("GSE1",) for item in result.items)
    (tmp_path / "genes.tsv").write_text("g1\tG1\ng2\tG2\n")
    rejected = Converter().convert(tmp_path / "matrix.mtx", out_type="h5ad")
    assert rejected.status == "failed"
    assert any(d.code == "ambiguous_companion" for d in rejected.items[0].diagnostics)


def test_failed_expression_requirements_preserve_independent_dataset(tmp_path):
    path = tmp_path / "counts.tsv"
    path.write_text("gene\tc1\ng1\t1\n")
    first = package()
    first["sample"][0]["supplementary_data"] = [{"value": str(path)}]
    second = package("GSE2", "GSM2")
    h5ad = JSON2H5ADConverter(
        available_memory=lambda: 10**12, memory_estimator=lambda *_: 100
    )
    result = Converter(services={"json2h5ad": h5ad}).convert(
        [first, second],
        out_type="h5ad",
        outdir=tmp_path / "out",
        output_options={"matrix_orientation": "genes-by-observations"},
    )
    assert result.status == "partial", result.to_dict()
    assert (tmp_path / "out/GSE1/GSM1.h5ad").is_file()
    assert any(d.code == "conversion_partial" for d in result.items[0].diagnostics)
