# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from __future__ import annotations
import csv
import logging
import os
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import urlparse


logger = logging.getLogger(__name__)

@dataclass(frozen=True)
class Asset:
    """One processed or raw data source associated with a sample or study."""

    scope_id: str
    path: str
    kind: str
    role: str = "primary"
    source: str = "json"
    members: tuple[dict, ...] = ()
    features_path: str | None = None
    barcodes_path: str | None = None
    orientation: str = "auto"
    md5: str | None = None
    study_scope: str | None = None
    reference: str | None = None
    annotation_source: str | None = None
    annotation_format: str | None = None
    annotation_sha256: str | None = None
    effective_annotation: str | None = None


class AssetManifest:
    """Load explicit asset mappings from CSV/TSV or compact CLI specifications."""

    def load(self, path: str) -> list[Asset]:
        delimiter = "\t" if Path(path).suffix.lower() in {".tsv", ".tab"} else ","
        with open(path, encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter=delimiter))
        if not rows or not {"scope_id", "path"}.issubset(rows[0]):
            raise ValueError("Asset manifest requires scope_id and path columns.")
        processed = []
        raw_groups = {}
        for row in rows:
            scope_id = (row.get("scope_id") or "").strip()
            asset_path = (row.get("path") or "").strip()
            if not scope_id or not asset_path:
                raise ValueError("Asset manifest scope_id and path values cannot be blank.")
            kind = (row.get("kind") or classify_asset(asset_path) or "").strip().lower()
            if not kind and os.path.isdir(asset_path):
                kind = "matrix"
            if kind not in {"h5ad", "matrix", "raw"}:
                raise ValueError(f"Unsupported asset kind for {asset_path}: {kind or 'unknown'}")
            role = (row.get("role") or "primary").strip()
            if kind == "raw":
                key = (scope_id, role)
                raw_groups.setdefault(key, []).append({
                    "uri": asset_path,
                    "read": (row.get("read") or "").strip() or None,
                    "lane": (row.get("lane") or "").strip() or None,
                    "run": (row.get("run") or row.get("lane") or "").strip() or None,
                    "md5": (row.get("md5") or "").strip() or None,
                })
                continue
            processed.append(Asset(
                scope_id=scope_id,
                path=asset_path,
                kind=kind,
                role=role,
                source="manifest",
                features_path=(row.get("features_path") or "").strip() or None,
                barcodes_path=(row.get("barcodes_path") or "").strip() or None,
                orientation=(row.get("orientation") or "auto").strip(),
                md5=(row.get("md5") or "").strip() or None,
            ))
        for (scope_id, role), members in raw_groups.items():
            processed.append(Asset(
                scope_id=scope_id,
                path=members[0]["uri"],
                kind="raw",
                role=role,
                source="manifest",
                members=tuple(members),
            ))
        return processed

    def parse_spec(self, spec: str) -> Asset:
        if "=" not in spec:
            raise ValueError("--asset must use ACCESSION=PATH_OR_URL syntax.")
        scope_id, path = (value.strip() for value in spec.split("=", 1))
        if not scope_id or not path:
            raise ValueError("--asset accession and path cannot be blank.")
        kind = classify_asset(path)
        if not kind and os.path.isdir(path):
            kind = "matrix"
        if not kind:
            raise ValueError(f"Cannot infer asset kind from {path}; use an asset manifest.")
        if kind == "raw":
            member = {"uri": path, "read": None, "run": None}
            return Asset(scope_id, path, kind, source="cli", members=(member,))
        return Asset(scope_id, path, kind, source="cli")


def classify_asset( path: str | None) -> str | None:
    if not path:
        return None
    filename = os.path.basename(urlparse(str(path)).path).lower()
    for suffix in (".gz", ".bz2", ".xz", ".zip"):
        if filename.endswith(suffix):
            filename = filename[: -len(suffix)]
            break
    if filename.endswith(".h5ad"):
        return "h5ad"
    if filename.endswith((".h5", ".mtx", ".csv", ".tsv", ".txt")):
        return "matrix"
    if filename.endswith((".fastq", ".fq")):
        return "raw"
    return None
