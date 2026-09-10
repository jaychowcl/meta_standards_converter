# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from __future__ import annotations


import heapq

import logging

import os


from dataclasses import replace


from typing import Mapping, Protocol

from urllib.parse import urlparse


from meta_standards_converter.expression.assets import Asset

from .readers import _tenx_member
logger = logging.getLogger(__name__)

class AssetDiscovery(Protocol):
    def discover(self, packages: list[dict]) -> list[Asset]: ...


class SourcePlanner:
    """Discover assets and select the best available source for every sample."""

    SOURCE_RANK = {"json": 0, "cli": 1, "manifest": 2, "nfcore": 3}
    KIND_RANK = {"raw": 0, "matrix": 1, "h5ad": 2}

    def __init__(self, discovery: AssetDiscovery | None = None):
        self.discovery = discovery or DefaultAssetDiscovery()

    def plan(
        self,
        packages: list[dict],
        explicit_assets: list[Asset] | None = None,
        force_reprocess: bool = False,
    ) -> dict[str, Asset]:
        assets = self.discovery.discover(packages)
        assets.extend(explicit_assets or [])
        assets = self._group_10x_assets(assets)
        assets = self._coalesce_raw_assets(assets)
        samples = self.samples(packages)
        study_by_sample = self._study_by_sample(packages)
        assets_by_scope = self._index_assets_by_scope(assets)
        planned = {}

        for sample_id in samples:
            study_id = study_by_sample.get(sample_id)
            scoped = assets_by_scope.get(sample_id, ())
            study_scoped = assets_by_scope.get(study_id, ()) if study_id else ()
            candidates = [
                asset
                if asset.scope_id == sample_id
                else replace(asset, scope_id=sample_id, study_scope=study_id)
                for _order, asset in heapq.merge(
                    scoped,
                    study_scoped,
                    key=lambda item: item[0],
                )
            ]
            if force_reprocess:
                candidates = [asset for asset in candidates if asset.kind == "raw"]
                if not candidates:
                    raise ValueError(f"{sample_id} has no raw FASTQ files for forced reprocessing.")
            if not candidates:
                raise ValueError(f"{sample_id} has no supported H5AD, matrix, or raw FASTQ source.")
            planned[sample_id] = max(
                candidates,
                key=lambda asset: (
                    self.SOURCE_RANK.get(asset.source, -1),
                    self.KIND_RANK.get(asset.kind, -1),
                ),
            )
        return planned

    @staticmethod
    def _index_assets_by_scope(
        assets: list[Asset],
    ) -> dict[str, list[tuple[int, Asset]]]:
        by_scope: dict[str, list[tuple[int, Asset]]] = {}
        for order, asset in enumerate(assets):
            by_scope.setdefault(asset.scope_id, []).append((order, asset))
        return by_scope

    def _group_10x_assets(self, assets: list[Asset]) -> list[Asset]:
        groups: dict[
            tuple[str, str, str], dict[str, tuple[int, Asset]]
        ] = {}
        for index, asset in enumerate(assets):
            member = _tenx_member(asset.path)
            if member is None:
                continue
            prefix, role = member
            groups.setdefault(
                (asset.scope_id, prefix, asset.source), {}
            )[role] = (index, asset)
        complete = {
            key: members
            for key, members in groups.items()
            if {"matrix", "barcodes", "features"} <= set(members)
        }
        member_groups = {
            index: key
            for key, members in complete.items()
            for index, _asset in members.values()
        }
        emitted: set[tuple[str, str, str]] = set()
        result: list[Asset] = []
        for index, asset in enumerate(assets):
            key = member_groups.get(index)
            if key is None:
                result.append(asset)
                continue
            if key in emitted:
                continue
            emitted.add(key)
            members = complete[key]
            matrix = members["matrix"][1]
            result.append(
                replace(
                    matrix,
                    role="10x_mtx",
                    barcodes_path=members["barcodes"][1].path,
                    features_path=members["features"][1].path,
                )
            )
        return result

    def _coalesce_raw_assets(self, assets: list[Asset]) -> list[Asset]:
        retained = [asset for asset in assets if asset.kind != "raw"]
        groups = {}
        for asset in assets:
            if asset.kind != "raw":
                continue
            key = (asset.scope_id, asset.source, asset.role)
            members = list(asset.members) or [{"uri": asset.path, "md5": asset.md5}]
            groups.setdefault(key, []).extend(members)
        for (scope_id, source, role), members in groups.items():
            deduped = []
            seen = set()
            for member in members:
                path = member.get("uri") or member.get("filename")
                if not path or path in seen:
                    continue
                seen.add(path)
                deduped.append(member)
            retained.append(Asset(
                scope_id=scope_id,
                path=deduped[0].get("uri") or deduped[0].get("filename"),
                kind="raw",
                role=role,
                source=source,
                members=tuple(deduped),
            ))
        return retained

    def _study_by_sample(self, packages: list[dict]) -> dict[str, str]:
        result = {}
        for package in packages:
            series = package.get("series") if isinstance(package, Mapping) else None
            study_id = None
            if isinstance(series, dict):
                for accession in self._as_list(series.get("accession")):
                    value = self._value(accession)
                    if isinstance(value, str) and value.upper().startswith("GSE"):
                        study_id = value.upper()
                        break
            if not study_id:
                continue
            for sample in self._as_list(package.get("sample")):
                if isinstance(sample, dict):
                    sample_id = self.sample_accession(sample)
                    if sample_id:
                        result[sample_id] = study_id
        return result


    def samples(self, packages: list[dict]) -> list[str]:
        values = []
        for package in packages:
            for sample in self._as_list(package.get("sample")):
                if isinstance(sample, dict):
                    accession = self.sample_accession(sample)
                    if accession and accession not in values:
                        values.append(accession)
        return values

    def sample_accession(self, sample: dict) -> str | None:
        for accession in self._as_list(sample.get("accession")):
            value = self._value(accession)
            if isinstance(value, str) and value.upper().startswith("GSM"):
                return value.upper()
        value = sample.get("iid")
        return str(value) if value else None

    def classify(self, path: str | None) -> str | None:
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

    def _value(self, value):
        if isinstance(value, dict):
            return value.get("value") or value.get("uri") or value.get("filename")
        return value

    def _as_list(self, value) -> list:
        if value is None:
            return []
        return value if isinstance(value, list) else [value]


class DefaultAssetDiscovery:
    def discover(self, packages: list[dict]) -> list[Asset]:
        assets = []
        for package in packages:
            for sample in self._as_list(package.get("sample")):
                if not isinstance(sample, dict):
                    continue
                sample_id = self.sample_accession(sample)
                if not sample_id:
                    continue
                for key in ("supplementary_data", "raw_data"):
                    for entry in self._as_list(sample.get(key)):
                        path = self._value(entry)
                        kind = self.classify(path)
                        if path and kind:
                            assets.append(Asset(
                                sample_id,
                                path,
                                kind,
                                md5=entry.get("md5") if isinstance(entry, dict) else None,
                            ))

                fastqs = []
                for run in self._as_list(sample.get("sra_run")):
                    if not isinstance(run, dict):
                        continue
                    for fastq in self._as_list(run.get("fastq_files")):
                        if isinstance(fastq, dict) and (fastq.get("uri") or fastq.get("filename")):
                            fastqs.append(dict(fastq, run=run.get("run")))
                if fastqs:
                    assets.append(
                        Asset(
                            sample_id,
                            fastqs[0].get("uri") or fastqs[0].get("filename"),
                            "raw",
                            members=tuple(fastqs),
                        )
                    )
        return assets

    def samples(self, packages: list[dict]) -> list[str]:
        values = []
        for package in packages:
            for sample in self._as_list(package.get("sample")):
                if isinstance(sample, dict):
                    accession = self.sample_accession(sample)
                    if accession and accession not in values:
                        values.append(accession)
        return values

    def sample_accession(self, sample: dict) -> str | None:
        for accession in self._as_list(sample.get("accession")):
            value = self._value(accession)
            if isinstance(value, str) and value.upper().startswith("GSM"):
                return value.upper()
        value = sample.get("iid")
        return str(value) if value else None

    def classify(self, path: str | None) -> str | None:
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

    def _value(self, value):
        if isinstance(value, dict):
            return value.get("value") or value.get("uri") or value.get("filename")
        return value

    def _as_list(self, value) -> list:
        if value is None:
            return []
        return value if isinstance(value, list) else [value]
