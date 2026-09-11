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

def make_expression_source(tmp_path: Path) -> tuple[Path, Path]:
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
        "miniml_schema_version": "3.0",
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
