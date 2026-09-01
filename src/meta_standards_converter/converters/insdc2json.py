# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Shared study-scoped SRA/ENA conversion and explicit origin enrichment."""

from __future__ import annotations

from collections.abc import Mapping
from copy import deepcopy
import logging
from pathlib import Path
from typing import Any

from meta_standards_converter.miniml import MINiMLCodec, MINiMLPackage
from meta_standards_converter.runtime_contracts import ResourceProfile, get_resource_profile


logger = logging.getLogger(__name__)


class StudyConversionFailures(ValueError):
    """One or more independently resolved studies failed conversion."""

    def __init__(self, failures: Mapping[str, BaseException]) -> None:
        self.failures = dict(failures)
        rendered = "; ".join(
            f"{study}: {error}" for study, error in self.failures.items()
        )
        super().__init__(rendered)


def _items(value) -> list:
    if value is None:
        return []
    return list(value) if isinstance(value, (list, tuple)) else [value]


def _row_key(value: Any) -> str:
    if isinstance(value, Mapping):
        for key in ("value", "pubmed_id", "iid", "name", "ref"):
            if value.get(key) is not None:
                return f"{key}:{value[key]}"
        return repr(sorted(value.items()))
    return repr(value)


def _union(first, second) -> list:
    values = []
    seen = set()
    for item in [*_items(first), *_items(second)]:
        key = _row_key(item)
        if key not in seen:
            seen.add(key)
            values.append(deepcopy(item))
    return values


def _sample_identifiers(sample: Mapping) -> set[str]:
    values: set[str] = set()
    for candidate in [
        sample.get("iid"),
        *_items(sample.get("accession")),
        *_items(sample.get("sra_accession")),
        *_items(sample.get("ena_accession")),
    ]:
        if isinstance(candidate, Mapping):
            candidate = candidate.get("value") or candidate.get("target")
        if candidate:
            values.add(str(candidate).strip().upper())
    for run in _items(sample.get("sra_run")):
        if not isinstance(run, Mapping):
            continue
        for key in ("run", "experiment", "sample", "biosample", "geo_sample"):
            if run.get(key):
                values.add(str(run[key]).strip().upper())
    return values


def _merge_sample(base: dict, origin: Mapping) -> dict:
    merged = deepcopy(base)
    structural = {
        "iid",
        "status",
        "sra_run",
        "sra_accession",
        "ena_accession",
        "platform_ref",
        "library_layout",
        "library_selection",
        "library_source",
        "library_strategy",
        "instrument_model",
    }
    for key, value in origin.items():
        if key in structural or value in (None, [], {}, ""):
            continue
        merged[key] = deepcopy(value)
    merged["accession"] = _union(origin.get("accession"), base.get("accession"))
    for key in structural - {"iid"}:
        if key in base:
            merged[key] = deepcopy(base[key])
    merged["iid"] = base.get("iid") or origin.get("iid")
    return merged


def merge_origin_package(
    base: MINiMLPackage,
    origin: MINiMLPackage,
    *,
    source: str,
) -> MINiMLPackage:
    """Overlay only explicit one-to-one origin matches and keep INSDC structure."""

    codec = MINiMLCodec()
    base_data = codec.encode(base)
    origin_data = codec.encode(origin)
    insdc = base_data["extensions"]["insdc"]
    diagnostics = insdc["enrichment"].setdefault("diagnostics", [])

    series = base_data["series"]
    origin_series = origin_data["series"]
    for key in (
        "title",
        "summary",
        "overall_design",
        "type",
        "contributor_ref",
        "contributor",
        "contact_ref",
        "contact",
        "variable",
        "repeats",
        "quality_controls",
        "replicate_types",
        "normalization_types",
    ):
        if origin_series.get(key) not in (None, [], {}, ""):
            series[key] = deepcopy(origin_series[key])
    series["accession"] = _union(origin_series.get("accession"), series.get("accession"))
    series["pubmed_id"] = _union(origin_series.get("pubmed_id"), series.get("pubmed_id"))
    series["pubmed_publication"] = _union(
        origin_series.get("pubmed_publication"), series.get("pubmed_publication")
    )
    series["web_link"] = _union(origin_series.get("web_link"), series.get("web_link"))
    series["protocols"] = _union(
        origin_series.get("protocols"), series.get("protocols")
    )
    series["assay_paths"] = _union(
        origin_series.get("assay_paths"), series.get("assay_paths")
    )

    base_samples = list(base_data.get("sample", []))
    origin_samples = list(origin_data.get("sample", []))
    base_ids = [_sample_identifiers(sample) for sample in base_samples]
    origin_ids = [_sample_identifiers(sample) for sample in origin_samples]
    merged_origin: set[int] = set()
    for base_index, identifiers in enumerate(base_ids):
        matches = [
            index
            for index, candidate in enumerate(origin_ids)
            if identifiers & candidate
        ]
        if len(matches) == 1:
            origin_index = matches[0]
            reverse = [
                index
                for index, candidate in enumerate(base_ids)
                if candidate & origin_ids[origin_index]
            ]
            if len(reverse) == 1:
                base_samples[base_index] = _merge_sample(
                    base_samples[base_index], origin_samples[origin_index]
                )
                merged_origin.add(origin_index)
                continue
        if matches:
            diagnostics.append(
                {
                    "code": "ambiguous_sample_alignment",
                    "base_sample": base_samples[base_index].get("iid"),
                    "origin_samples": [origin_samples[index].get("iid") for index in matches],
                }
            )
    for index, sample in enumerate(origin_samples):
        if index not in merged_origin:
            base_samples.append(deepcopy(sample))
            diagnostics.append(
                {
                    "code": "unmatched_origin_sample",
                    "origin_sample": sample.get("iid"),
                }
            )
    base_data["sample"] = base_samples

    for collection in ("database", "organization", "contributor"):
        base_data[collection] = _union(
            origin_data.get(collection), base_data.get(collection)
        )
    # INSDC sequencer definitions remain authoritative; origin array platforms
    # are additive when their identifiers are distinct.
    base_data["platform"] = _union(
        base_data.get("platform"), origin_data.get("platform")
    )

    documents = list(base_data["source"].get("documents", []))
    existing_names = {item["name"] for item in documents}
    for document in origin_data.get("source", {}).get("documents", []):
        copied = deepcopy(document)
        original_name = copied["name"]
        copied["name"] = f"{source}-{original_name}"
        suffix = 2
        while copied["name"] in existing_names:
            copied["name"] = f"{source}-{suffix}-{original_name}"
            suffix += 1
        existing_names.add(copied["name"])
        copied["kind"] = f"{source}_{copied['kind']}"
        documents.append(copied)
    base_data["source"]["documents"] = documents
    insdc["enrichment"]["status"] = "applied"
    insdc["enrichment"]["source"] = source
    return codec.decode(base_data).package


class INSDC2JSONConverter:
    provider = "provider"
    base_suffix = ".json"

    def __init__(
        self,
        *,
        fetcher,
        parser,
        geo_converter=None,
        ae_converter=None,
        resource_profile: str | ResourceProfile = "standard",
        resource_overrides=None,
    ) -> None:
        self.resource_profile = get_resource_profile(
            resource_profile, overrides=resource_overrides
        )
        self.fetcher = fetcher
        self.parser = parser
        self.geo_converter = geo_converter
        self.ae_converter = ae_converter
        self._written_studies: set[tuple[str, str]] = set()

    def _geo(self):
        if self.geo_converter is None:
            from meta_standards_converter.converters.geo2json import geo2json

            self.geo_converter = geo2json(resource_profile=self.resource_profile)
        return self.geo_converter

    def _ae(self):
        if self.ae_converter is None:
            from meta_standards_converter.converters.ae2json import ae2json

            self.ae_converter = ae2json(resource_profile=self.resource_profile)
        return self.ae_converter

    @staticmethod
    def _origin_package(packages: list[MINiMLPackage], accession: str) -> MINiMLPackage:
        candidates = []
        for package in packages:
            data = MINiMLCodec().encode(package)
            values = {
                str(item.get("value") if isinstance(item, Mapping) else item).upper()
                for item in _items(data["series"].get("accession"))
            }
            if (
                accession.upper() in values
                or data["series"].get("iid", "").upper() == accession.upper()
            ):
                candidates.append(package)
        if len(candidates) != 1:
            raise ValueError(
                f"origin enrichment for {accession} did not resolve exactly one study"
            )
        return candidates[0]

    def convert(
        self,
        accession: str,
        *,
        enrich_geo: bool = False,
        enrich_ae: bool = False,
        out: str | Path | None = None,
    ) -> list[MINiMLPackage]:
        if enrich_geo and enrich_ae:
            raise ValueError("GEO and ArrayExpress enrichment are mutually exclusive")
        mode = "geo" if enrich_geo else "ae" if enrich_ae else "none"
        unique_results = {}
        for result in self.fetcher.fetch(accession):
            unique_results.setdefault(result.study_accession, result)
        results = [unique_results[key] for key in sorted(unique_results)]
        packages: list[MINiMLPackage] = []
        failures: dict[str, BaseException] = {}
        for result in results:
            try:
                package = self._convert_result(
                    result,
                    mode=mode,
                    enrich_geo=enrich_geo,
                    enrich_ae=enrich_ae,
                )
            except Exception as error:
                failures[result.study_accession] = error
                continue
            packages.append(package)
            if out is not None:
                suffix = self.base_suffix
                if enrich_geo:
                    suffix = suffix.removesuffix(".json") + ".geo.json"
                elif enrich_ae:
                    suffix = suffix.removesuffix(".json") + ".ae.json"
                key = (result.study_accession, suffix)
                if key not in self._written_studies:
                    MINiMLCodec().dump(
                        [package], Path(out) / f"{result.study_accession}{suffix}"
                    )
                    self._written_studies.add(key)
        if failures:
            raise StudyConversionFailures(failures)
        return packages

    def _convert_result(
        self,
        result,
        *,
        mode: str,
        enrich_geo: bool,
        enrich_ae: bool,
    ) -> MINiMLPackage:
        package = self.parser.parse(result, enrichment=mode)
        data = MINiMLCodec().encode(package)
        extension = data["extensions"]["insdc"]
        links = extension["enrichment"]["linked_accessions"]
        if enrich_geo:
            if not links["geo"]:
                raise ValueError(
                    f"study {result.study_accession} has no compatible GEO link"
                )
            origin_accession = links["geo"][0]
            origin_packages = self._geo().convert(
                gse=origin_accession,
                enrich=False,
            )
            package = merge_origin_package(
                package,
                self._origin_package(origin_packages, origin_accession),
                source="geo",
            )
        elif enrich_ae:
            if not links["arrayexpress"]:
                raise ValueError(
                    f"study {result.study_accession} has no compatible ArrayExpress link"
                )
            origin_accession = links["arrayexpress"][0]
            origin_packages = self._ae().convert(origin_accession)
            package = merge_origin_package(
                package,
                self._origin_package(origin_packages, origin_accession),
                source="arrayexpress",
            )
        else:
            for warning in extension["warnings"]:
                if warning["code"].startswith("linked_"):
                    logger.warning(
                        "%s: %s",
                        result.study_accession,
                        warning["message"],
                    )
        return package


__all__ = [
    "INSDC2JSONConverter",
    "StudyConversionFailures",
    "merge_origin_package",
]
