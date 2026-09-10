# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Offline provider replay and lossless, reviewable output assertions.

Only the temporary workspace and JSON2OBS staging-directory token are normalized.
No ordering, scientific values, diagnostics, or release versions are discarded.
"""
import io
from datetime import date
import json
import re
import shutil
import tarfile
from pathlib import Path

import requests

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures"
GEO = FIXTURES / "studies/GSE328265"
PBMC = FIXTURES / "studies/pbmc3k"
AE = FIXTURES / "studies/E-MTAB-6486"


def prepare(case, workspace):
    for source in (case / "inputs").iterdir():
        shutil.copyfile(source, workspace / source.name)


def install_replay(monkeypatch, geo_xml=None, *, documents=None, fail_on=()):
    from meta_standards_converter.helpers.request_helper import RateLimitedRequester
    calls = []
    documents = documents or {"GSE328265": geo_xml or (GEO / "inputs/geo.xml").read_bytes()}
    def get(_requester, url, **kwargs):
        params = kwargs.get("params", {})
        accession = url.rsplit("/", 1)[-1].removesuffix("_family.xml.tgz")
        if url.endswith("_family.xml.tgz") and accession in documents:
            calls.append("geo:" + accession)
            xml = documents[accession]
            stream = io.BytesIO()
            with tarfile.open(fileobj=stream, mode="w:gz") as archive:
                entry = tarfile.TarInfo(accession + "_family.xml")
                entry.size = len(xml)
                archive.addfile(entry, io.BytesIO(xml))
            content = stream.getvalue()
        elif url.endswith("/esummary.fcgi") and params.get("db") == "pubmed" and params.get("id") == "42129775":
            calls.append("pubmed:42129775")
            content = (GEO / "provider_responses/pubmed.xml").read_bytes()
        elif url.endswith("/efetch.fcgi") and params.get("db") == "sra" and params.get("id") == "SRX32831930":
            calls.append("sra:SRX32831930")
            content = (GEO / "provider_responses/sra.xml").read_bytes()
        elif url.endswith("/filereport") and params.get("accession") == "SRX32831930":
            assert params["result"] == "read_run"
            assert params["fields"] == "run_accession,fastq_ftp,fastq_md5,fastq_bytes"
            calls.append("ena:SRX32831930")
            content = (GEO / "provider_responses/ena.json").read_bytes()
        else:
            raise AssertionError(f"Unexpected provider request: {url} {params}")
        if calls[-1] in fail_on:
            raise requests.RequestException("fixture provider unavailable")
        response = requests.Response()
        response.status_code = 200
        response.url = url
        response.headers["Content-Length"] = str(len(content))
        response._content = content
        response._content_consumed = True
        return response
    monkeypatch.setattr(RateLimitedRequester, "get", get)
    return calls


def primitive(value):
    import numpy as np
    import pandas as pd
    from scipy import sparse
    if isinstance(value, pd.DataFrame):
        return {"kind": "dataframe", "index": primitive(value.index),
                "columns": primitive(value.columns),
                "dtypes": [str(dtype) for dtype in value.dtypes],
                "data": [[primitive(x) for x in row] for row in value.itertuples(index=False, name=None)],
                "categories": {str(col): {"values": primitive(value[col].cat.categories), "ordered": value[col].cat.ordered}
                               for col in value if isinstance(value[col].dtype, pd.CategoricalDtype)}}
    if isinstance(value, pd.Index):
        return {"values": [primitive(x) for x in value], "name": primitive(value.name), "dtype": str(value.dtype)}
    if sparse.issparse(value):
        return {"kind": value.format, "dtype": str(value.dtype), "shape": list(value.shape), "values": value.toarray().tolist()}
    if isinstance(value, np.ndarray):
        return {"kind": "array", "dtype": str(value.dtype), "shape": list(value.shape), "values": primitive(value.tolist())}
    if isinstance(value, np.generic):
        return primitive(value.item())
    if isinstance(value, dict):
        return {str(k): primitive(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [primitive(x) for x in value]
    if isinstance(value, float) and not np.isfinite(value):
        return {"nonfinite": str(value)}
    return value


def h5ad_content(path):
    import anndata
    a = anndata.read_h5ad(path)
    return {"X": primitive(a.X), "obs": primitive(a.obs), "var": primitive(a.var),
            "uns": primitive(dict(a.uns)), "layers": primitive(dict(a.layers)),
            "obsm": primitive(dict(a.obsm)), "varm": primitive(dict(a.varm)),
            "obsp": primitive(dict(a.obsp)), "varp": primitive(dict(a.varp)),
            "raw": None if a.raw is None else {"X": primitive(a.raw.X), "var": primitive(a.raw.var)}}


def _submission_date(match):
    assert match.group(1) == date.today().isoformat(), "Submission date must be today's date"
    return "Comment[ArrayExpressSubmissionDate]\t<RUN_DATE>"


def normalize(value, workspace):
    if isinstance(value, dict):
        return {k: normalize(v, workspace) for k, v in value.items()}
    if isinstance(value, list):
        return [normalize(v, workspace) for v in value]
    if isinstance(value, str):
        if "Comment[ArrayExpressSubmissionDate]\t" in value:
            value = re.sub(r"(?m)^Comment\[ArrayExpressSubmissionDate\]\t([^\n]+)$",
                           _submission_date, value)
        value = value.replace(str(workspace), "<WORK>")
        return re.sub(r"\.out\.json2obs-[^/\\\s\"]+", ".out.json2obs-<STAGING>", value)
    return value


def artifacts(out, workspace):
    result = {}
    for path in sorted(out.rglob("*")):
        relative = path.relative_to(out)
        if not path.is_file() or any(p.startswith(".") for p in relative.parts):
            continue
        if path.suffix == ".h5ad":
            content = h5ad_content(path)
        elif path.suffix == ".json":
            content = json.loads(path.read_text())
        else:
            content = path.read_text()
        result[str(relative)] = normalize(content, workspace)
    assert result, "Conversion must publish artifacts"
    return result


def assert_expected(case, converter, out, workspace):
    expected_root = case / "expected" / converter
    assert expected_root.is_dir(), f"Missing reviewed expectation: {expected_root}"
    expected = {}
    for path in sorted(expected_root.iterdir()):
        name = path.name.removesuffix(".json") if path.name.endswith(".h5ad.json") else path.name
        expected[name] = json.loads(path.read_text()) if path.suffix == ".json" else path.read_text()
    assert artifacts(out, workspace) == expected
