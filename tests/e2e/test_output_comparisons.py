# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Deliberate corruptions must be detected by the same full-output assertion."""
import json
import anndata
import numpy as np
import pytest
from meta_standards_converter.converters import JSON2H5ADConverter
from tests.support.contracts import PBMC, prepare, assert_expected

@pytest.mark.parametrize("change", ["count", "observation_order", "feature_order", "identity", "dtype", "retained_metadata", "diagnostic"])
def test_comparison_rejects_scientific_and_contract_drift(change, workspace):
    prepare(PBMC, workspace)
    result = JSON2H5ADConverter().convert("miniml.json", out="out", asset_specs=["PBMC3K_SAMPLE=counts.tsv"], matrix_orientation="genes-by-observations")
    assert_expected(PBMC, "json2h5ad", workspace / "out", workspace)
    path = result.sample_h5ads["PBMC3K_SAMPLE"]
    a = anndata.read_h5ad(path)
    if change == "count": a.X[0, 4] = 8
    elif change == "observation_order": a = a[::-1].copy()
    elif change == "feature_order": a = a[:, ::-1].copy()
    elif change == "identity": a.obs_names = ["changed", *list(a.obs_names)[1:]]
    elif change == "dtype": a.X = a.X.astype(np.float32)
    elif change == "retained_metadata": del a.uns["msc_miniml"]
    else:
        manifest = workspace / "out/PBMC3K.json2h5ad.json"
        data = json.loads(manifest.read_text()); data["warnings"].append("unexpected semantic diagnostic")
        manifest.write_text(json.dumps(data))
    if change != "diagnostic": a.write_h5ad(path)
    with pytest.raises(AssertionError):
        assert_expected(PBMC, "json2h5ad", workspace / "out", workspace)


@pytest.mark.parametrize("actual,expected", [(True, 1), (1, 1.0), ("", None)])
def test_comparison_preserves_json_scalar_type_and_missing_value_distinctions(actual, expected):
    from tests.support.contracts import assert_contract_equal
    with pytest.raises(AssertionError):
        assert_contract_equal({"field": actual}, {"field": expected})
