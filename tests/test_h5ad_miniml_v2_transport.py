# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from __future__ import annotations

import json

import pytest

from meta_standards_converter.converters import JSON2H5ADConverter


def _package() -> dict:
    return {
        "miniml_schema_version": "2.0",
        "source": {"format": "GEO MINiML"},
        "sample": [{"iid": "GSM1"}],
        "series": {"accession": [{"value": "GSE1"}]},
    }


def test_h5ad_miniml_envelope_losslessly_transports_typed_packages(tmp_path) -> None:
    anndata = pytest.importorskip("anndata")
    package = _package()
    adata = anndata.AnnData(shape=(1, 1))

    JSON2H5ADConverter()._attach_miniml(
        adata,
        [package],
        source_json=str(tmp_path / "input.json"),
        source_json_sha256=None,
        artifact_parent=tmp_path,
    )

    assert json.loads(adata.uns["msc_miniml"]["packages_json"]) == [package]
    destination = tmp_path / "transport.h5ad"
    adata.write_h5ad(destination)
    restored = anndata.read_h5ad(destination)
    assert json.loads(restored.uns["msc_miniml"]["packages_json"]) == [package]
