# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from __future__ import annotations


import gzip


import os

import re

import shutil


import tempfile

from dataclasses import replace


from pathlib import Path
from typing import Protocol, Callable, Any


from urllib.parse import urlparse


from .assets import Asset

_TENX_MEMBER = re.compile(
    r"^(?P<prefix>.+?)[._](?P<role>"
    r"matrix\.mtx|barcodes\.tsv|genes\.tsv|features\.tsv)$",
    re.IGNORECASE,
)

def scientific_modules():
    try:
        import anndata
        import numpy
        import pandas
        from scipy import sparse
    except ImportError as exc:
        raise RuntimeError(
            "json2h5ad requires optional dependencies; install "
            "meta-standards-converter[h5ad]."
        ) from exc
    return anndata, numpy, pandas, sparse


def scanpy_module():
    try:
        import scanpy
    except ImportError as exc:
        raise RuntimeError(
            "10x matrix input requires Scanpy; install "
            "meta-standards-converter[h5ad]."
        ) from exc
    return scanpy


def read_h5ad(anndata, path: str):
    if not str(path).lower().endswith(".gz"):
        return anndata.read_h5ad(path)
    temporary_path = None
    try:
        with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as temporary:
            temporary_path = temporary.name
            with gzip.open(path, "rb") as compressed:
                shutil.copyfileobj(compressed, temporary, length=1024 * 1024)
        return anndata.read_h5ad(temporary_path)
    finally:
        if temporary_path:
            Path(temporary_path).unlink(missing_ok=True)


def underlying_suffix(path: str) -> str:
    name = Path(path).name.lower()
    for suffix in (".gz", ".bz2", ".xz", ".zip"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    return Path(name).suffix


def _tenx_member(value: str) -> tuple[str, str] | None:
    filename = os.path.basename(urlparse(str(value)).path)
    for suffix in (".gz", ".bz2", ".xz", ".zip"):
        if filename.casefold().endswith(suffix):
            filename = filename[: -len(suffix)]
            break
    match = _TENX_MEMBER.fullmatch(filename)
    if match is None:
        return None
    raw_role = match.group("role").casefold()
    role = (
        "matrix"
        if raw_role == "matrix.mtx"
        else "barcodes"
        if raw_role == "barcodes.tsv"
        else "features"
    )
    return match.group("prefix").casefold(), role

def _tenx_local_name(role: str, source: str, *, legacy: bool) -> str:
    if legacy:
        return {
            "matrix": "matrix.mtx",
            "barcodes": "barcodes.tsv",
            "features": "genes.tsv",
        }[role]
    filename = os.path.basename(urlparse(str(source)).path).casefold()
    compressed = ".gz" if filename.endswith(".gz") else ""
    if role == "matrix":
        return f"matrix.mtx{compressed}"
    if role == "barcodes":
        return f"barcodes.tsv{compressed}"
    stem = "features" if ".features.tsv" in filename else "genes"
    return f"{stem}.tsv{compressed}"

def _copy_10x_member(source: str, destination: Path, *, decompress: bool) -> None:
    if decompress and str(source).casefold().endswith(".gz"):
        with gzip.open(source, "rb") as input_stream:
            with destination.open("wb") as output_stream:
                shutil.copyfileobj(input_stream, output_stream)
        return
    shutil.copyfile(source, destination)


class AssetReader(Protocol):
    def read(self, asset: Asset, *, orientation: str = "auto", localize: Callable[[str], str]) -> Any: ...


class ProcessedAssetReader:
    def read(self, asset: Asset, *, orientation: str = "auto", localize):
        if (
            underlying_suffix(asset.path) == ".mtx"
            and asset.features_path
            and asset.barcodes_path
        ):
            with tempfile.TemporaryDirectory(prefix="msc-10x-") as directory:
                prepared = Path(directory)
                paths = {
                    "matrix": localize(asset.path, md5=asset.md5),
                    "barcodes": localize(asset.barcodes_path),
                    "features": localize(asset.features_path),
                }
                legacy = ".genes.tsv" in os.path.basename(
                    urlparse(str(asset.features_path)).path
                ).casefold()
                for role, source in paths.items():
                    destination = prepared / _tenx_local_name(
                        role,
                        asset.features_path if role == "features" else source,
                        legacy=legacy,
                    )
                    _copy_10x_member(source, destination, decompress=legacy)
                return self.read(
                    replace(
                        asset,
                        path=str(prepared),
                        features_path=None,
                        barcodes_path=None,
                    ),
                    orientation=orientation, localize=localize,
                )
        anndata, numpy, pandas, sparse = scientific_modules()
        path = localize(asset.path, md5=asset.md5)
        if asset.kind == "h5ad":
            adata = read_h5ad(anndata, path)
            if asset.study_scope:
                accession_column = next(
                    (
                        column
                        for column in (
                            "msc.sample.accession",
                            "geo_accession",
                            "sample_id",
                            "sample",
                            "gsm_accession",
                        )
                        if column in adata.obs
                    ),
                    None,
                )
                if not accession_column:
                    raise ValueError(
                        f"Study H5AD {asset.path} cannot be mapped to samples; "
                        "obs needs msc.sample.accession, geo_accession, "
                        "sample_id, sample, or gsm_accession."
                    )
                mask = adata.obs[accession_column].astype(str).str.upper() == asset.scope_id.upper()
                if not mask.any():
                    raise ValueError(f"Study H5AD {asset.path} contains no observations for {asset.scope_id}.")
                adata = adata[mask].copy()
        elif underlying_suffix(path) == ".h5":
            scanpy = scanpy_module()
            adata = scanpy.read_10x_h5(path, gex_only=True)
        elif underlying_suffix(path) == ".mtx" or Path(path).is_dir():
            scanpy = scanpy_module()
            matrix_dir = path if Path(path).is_dir() else str(Path(path).parent)
            adata = scanpy.read_10x_mtx(matrix_dir, var_names="gene_ids", make_unique=True)
        else:
            separator = "," if underlying_suffix(path) == ".csv" else "\t"
            frame = pandas.read_csv(path, sep=separator, index_col=0)
            feature_annotations = None
            if asset.role == "rnaseq_counts":
                if asset.scope_id not in frame.columns:
                    raise ValueError(
                        f"RNA-seq count matrix {asset.path} has no column for {asset.scope_id}."
                    )
                if "gene_name" in frame.columns:
                    feature_annotations = frame[["gene_name"]].copy()
                frame = frame[[asset.scope_id]]
            try:
                values = frame.apply(pandas.to_numeric, errors="raise")
            except (TypeError, ValueError) as exc:
                raise ValueError(f"Matrix {asset.path} contains nonnumeric values.") from exc
            raw = values.to_numpy()
            if raw.size == 0:
                raise ValueError(f"Matrix {asset.path} is empty.")
            if not numpy.isfinite(raw).all() or (raw < 0).any():
                raise ValueError(f"Matrix {asset.path} must contain finite nonnegative values.")
            if asset.study_scope and asset.role != "rnaseq_counts":
                if orientation == "genes-by-observations" and asset.scope_id in values.columns:
                    values = values[[asset.scope_id]]
                    raw = values.to_numpy()
                elif orientation == "observations-by-genes" and asset.scope_id in values.index:
                    values = values.loc[[asset.scope_id]]
                    raw = values.to_numpy()
                else:
                    raise ValueError(
                        f"Study matrix {asset.path} cannot be mapped to {asset.scope_id}."
                    )
            if orientation == "auto":
                raise ValueError(
                    f"Matrix orientation for {asset.path} is ambiguous; specify "
                    "genes-by-observations or observations-by-genes."
                )
            if orientation == "genes-by-observations":
                raw = raw.T
                obs_names = values.columns.astype(str)
                var_names = values.index.astype(str)
            elif orientation == "observations-by-genes":
                obs_names = values.index.astype(str)
                var_names = values.columns.astype(str)
            else:
                raise ValueError(f"Unsupported matrix orientation: {orientation}")
            adata = anndata.AnnData(
                X=sparse.csr_matrix(raw),
                obs=pandas.DataFrame(index=obs_names),
                var=pandas.DataFrame(index=var_names),
            )
            if feature_annotations is not None:
                adata.var["gene_name"] = (
                    feature_annotations.reindex(var_names)["gene_name"].astype(str).to_numpy()
                )
            if asset.role == "rnaseq_counts" and asset.features_path:
                tpm = pandas.read_csv(asset.features_path, sep="\t", index_col=0)
                if asset.scope_id not in tpm.columns:
                    raise ValueError(
                        f"RNA-seq TPM matrix {asset.features_path} has no column for {asset.scope_id}."
                    )
                tpm_values = tpm[[asset.scope_id]].reindex(values.index)
                tpm_values = tpm_values.apply(pandas.to_numeric, errors="raise")
                if tpm_values.isna().any().any():
                    raise ValueError("RNA-seq TPM features do not align with count features.")
                adata.layers["tpm"] = sparse.csr_matrix(tpm_values.to_numpy().T)
        if adata.n_obs == 0 or adata.n_vars == 0:
            raise ValueError(f"Processed asset {asset.path} contains an empty matrix.")
        if not sparse.issparse(adata.X):
            adata.X = sparse.csr_matrix(adata.X)
        else:
            adata.X = adata.X.tocsr()
        adata.var_names_make_unique()
        return adata
