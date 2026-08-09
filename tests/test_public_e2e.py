# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Offline end-to-end coverage through the supported converter interfaces."""

import json
import builtins
from pathlib import Path
from unittest.mock import Mock

from meta_standards_converter.ae_handlers.ae_constructor import AEConstructor
from meta_standards_converter.ae_handlers.ae_idf_handlers import IDFConstructor
from meta_standards_converter.ae_handlers.ae_sdrf_handlers import SDRFConstructor
from meta_standards_converter.converters.ae2json import ae2json
from meta_standards_converter.converters.geo2ae import geo2ae
from meta_standards_converter.converters.geo2json import geo2json
from meta_standards_converter.converters.json2ae import json2ae
from meta_standards_converter.converters.json2h5ad import JSON2H5ADConverter
from meta_standards_converter.geo_handlers.geo_parser import GEOParser
from meta_standards_converter.miniml import MINiMLCodec


ROOT = Path(__file__).resolve().parents[1]


class StaticGEOFetcher:
    def __init__(self, miniml):
        self.miniml = miniml

    def fetch_gse_miniml(self, gse):
        assert gse == "GSE328265"
        return self.miniml


class IdentityEnricher:
    def enrich(self, data):
        return data


def offline_constructor():
    pubmed = Mock()
    pubmed.pubmed_summary.return_value = (None, None, None, None, None, None)
    insdc = Mock()
    insdc._extract_sra.return_value = []
    insdc.fetch_sra_runs.return_value = []
    return AEConstructor(
        idf_constructor=IDFConstructor(pubmed_fetcher=pubmed),
        sdrf_constructor=SDRFConstructor(insdc_fetcher=insdc),
    )


def test_public_metadata_converters_write_interoperable_artifacts(tmp_path):
    miniml = (ROOT / "tests" / "GSE328265_family.xml").read_text(encoding="utf-8")
    fetcher = StaticGEOFetcher(miniml)
    parser = GEOParser(geo_fetcher=fetcher)
    json_out = tmp_path / "json"

    packages = geo2json(
        geo_fetcher=fetcher,
        parser=parser,
        enricher=IdentityEnricher(),
    ).convert("GSE328265", out=str(json_out))

    json_path = json_out / "GSE328265.json"
    assert json.loads(json_path.read_text(encoding="utf-8")) == MINiMLCodec().encode_many(packages)

    geo_magetab = tmp_path / "geo-magetab"
    geo2ae(
        geo_fetcher=fetcher,
        parser=parser,
        enricher=IdentityEnricher(),
        ae_constructor=offline_constructor(),
    ).convert("GSE328265", out=str(geo_magetab))
    assert (geo_magetab / "E-GEOD-328265.idf.txt").is_file()
    assert (geo_magetab / "E-GEOD-328265.sdrf.txt").is_file()

    json_magetab = tmp_path / "json-magetab"
    json2ae(ae_constructor=offline_constructor()).convert(
        str(json_path), out=str(json_magetab), enrich=False
    )
    idf = json_magetab / "E-GEOD-328265.idf.txt"
    assert idf.is_file()

    roundtrip_out = tmp_path / "roundtrip"
    roundtrip = ae2json().convert(str(idf), out=str(roundtrip_out))
    assert roundtrip[0]["series"]["title"].startswith(
        "A CSF Disease-Associated Macrophage Signature"
    )
    assert list(roundtrip_out.glob("*.json"))


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
        "miniml_schema_version": "2.0",
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
    assert Path(result.manifest_path).is_file()
    assert not result.partial
