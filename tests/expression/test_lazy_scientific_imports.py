# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
import json
import builtins
from pathlib import Path
from meta_standards_converter.converters import JSON2H5ADConverter

def test_public_processed_h5ad_conversion_does_not_import_scanpy(
    tmp_path, monkeypatch
):
    import anndata
    import pandas
    from scipy import sparse

    source_h5ad = tmp_path / "source.h5ad"
    anndata.AnnData(
        X=sparse.csr_matrix([[1, 2]]),
        obs=pandas.DataFrame(index=["cell-1"]),
        var=pandas.DataFrame(index=["ENSG1", "ENSG2"]),
    ).write_h5ad(source_h5ad)
    package = {
        "miniml_schema_version": "3.0",
        "source": {"format": "test"},
        "series": {"accession": [{"value": "GSE1"}]},
        "sample": [
            {
                "iid": "GSM1",
                "accession": [{"value": "GSM1"}],
                "supplementary_data": [{"value": str(source_h5ad)}],
            }
        ],
    }
    source_json = tmp_path / "GSE1.json"
    source_json.write_text(json.dumps(package), encoding="utf-8")

    original_import = builtins.__import__

    def reject_scanpy(name, *args, **kwargs):
        if name == "scanpy" or name.startswith("scanpy."):
            raise ImportError("scanpy must stay lazy for processed H5AD")
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", reject_scanpy)

    result = JSON2H5ADConverter().convert(
        str(source_json), out=str(tmp_path / "h5ad")
    )

    assert Path(result.sample_h5ads["GSM1"]).is_file()
    assert result.combined_h5ad is None
    assert Path(result.manifest_path).is_file()
    assert not result.partial
