# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Project native MSC MINiML assay parameters without replacing raw values."""

from __future__ import annotations

import re
from typing import Any, Mapping

from meta_standards_converter.miniml import harmonized_value_mappings


PARAMETER_FIELDS = (
    "value", "unit", "harmonized_value", "harmonized_value_id",
    "harmonized_value_ontology", "harmonized_unit", "harmonized_unit_id",
    "harmonized_unit_ontology",
)


def _parameter_summary(
    package: Mapping[str, Any], sample: Mapping[str, Any]
) -> dict[str, tuple[Any, ...]]:
    """Return deterministic per-parameter values for one bound sample."""
    result: dict[str, list[Any]] = {}
    for row in _parameter_rows(package, sample=sample):
        slug = _slug(row["name"])
        if not slug:
            continue
        for field in PARAMETER_FIELDS:
            value = row.get(field)
            if value in (None, ""):
                continue
            key = f"msc.assay.parameter.{slug}.{field}"
            values = result.setdefault(key, [])
            if value not in values:
                values.append(value)
    return {key: tuple(values) for key, values in result.items()}


def _protocols_for_sample(
    package: Mapping[str, Any], sample: Mapping[str, Any]
) -> list[Mapping[str, Any]]:
    """Return declared protocols in sample-path order, or study order if unbound."""
    series = package.get("series")
    if not isinstance(series, Mapping):
        return []
    protocols = [
        item
        for item in _as_list(series.get("protocols"))
        if isinstance(item, Mapping) and item.get("name")
    ]
    if not protocols:
        return []
    paths = _bound_assay_paths(series, sample)
    if not paths:
        return protocols
    by_name = {str(item["name"]): item for item in protocols}
    selected: list[Mapping[str, Any]] = []
    seen: set[str] = set()
    for path in paths:
        for step in _as_list(path.get("steps")):
            if not isinstance(step, Mapping) or step.get("kind") != "protocol_application":
                continue
            reference = str(step.get("protocol_ref") or "")
            protocol = by_name.get(reference)
            if protocol is not None and reference not in seen:
                selected.append(protocol)
                seen.add(reference)
    return selected


def _material_types_for_sample(
    package: Mapping[str, Any], sample: Mapping[str, Any]
) -> list[Any]:
    """Return material types from assay nodes bound to the selected sample."""
    series = package.get("series")
    if not isinstance(series, Mapping):
        return []
    values: list[Any] = []
    for path in _bound_assay_paths(series, sample):
        for step in _as_list(path.get("steps")):
            if not isinstance(step, Mapping) or step.get("kind") == "protocol_application":
                continue
            value = step.get("material_type")
            if value not in (None, "") and value not in values:
                values.append(value)
    return values


def _bound_assay_paths(
    series: Mapping[str, Any], sample: Mapping[str, Any]
) -> list[Mapping[str, Any]]:
    identities = _sample_identities(sample)
    if not identities:
        return []
    result = []
    for path in _as_list(series.get("assay_paths")):
        if not isinstance(path, Mapping):
            continue
        bound = {
            str(step.get("sample_ref") or step.get("name"))
            for step in _as_list(path.get("steps"))
            if isinstance(step, Mapping)
            and step.get("kind") == "sample"
            and (step.get("sample_ref") or step.get("name")) not in (None, "")
        }
        if identities.intersection(bound):
            result.append(path)
    return result


def _parameter_rows(
    package: Mapping[str, Any], sample: Mapping[str, Any] | None = None
) -> list[dict[str, Any]]:
    """Return typed parameter occurrences, optionally restricted to one sample."""
    series = package.get("series")
    if not isinstance(series, Mapping):
        return []
    identities = _sample_identities(sample) if sample is not None else set()
    rows = []
    for path_index, assay in enumerate(_as_list(series.get("assay_paths"))):
        if not isinstance(assay, Mapping):
            continue
        steps = [step for step in _as_list(assay.get("steps")) if isinstance(step, Mapping)]
        bound = {
            str(step.get("sample_ref") or step.get("name"))
            for step in steps
            if step.get("kind") == "sample"
            and (step.get("sample_ref") or step.get("name")) not in (None, "")
        }
        if identities and not identities.intersection(bound):
            continue
        for step_index, step in enumerate(steps):
            if step.get("kind") != "protocol_application":
                continue
            for parameter in _as_list(step.get("parameter_values")):
                if not isinstance(parameter, Mapping):
                    continue
                unit = parameter.get("unit")
                unit = unit if isinstance(unit, Mapping) else {"value": unit}
                value_annotation = _first_annotation(parameter)
                unit_annotation = _first_annotation(unit, field="unit")
                rows.append({
                    "sample_accession": next(iter(sorted(identities or bound)), ""),
                    "document": assay.get("document", ""),
                    "path_index": path_index,
                    "step_index": step_index,
                    "protocol_ref": step.get("protocol_ref", ""),
                    "name": parameter.get("name", ""),
                    "value": parameter.get("value", ""),
                    "unit": unit.get("value", ""),
                    "term_source_ref": parameter.get("term_source_ref", ""),
                    "term_accession_number": parameter.get("term_accession_number", ""),
                    "harmonized_value": value_annotation.get("value", ""),
                    "harmonized_value_id": value_annotation.get("term_accession_number", ""),
                    "harmonized_value_ontology": value_annotation.get("term_source_ref", ""),
                    "harmonized_unit": unit_annotation.get("value", ""),
                    "harmonized_unit_id": unit_annotation.get("term_accession_number", ""),
                    "harmonized_unit_ontology": unit_annotation.get("term_source_ref", ""),
                })
    return rows


def _first_annotation(value: Any, *, field: str | None = None) -> Mapping[str, Any]:
    for annotation in harmonized_value_mappings(value):
        if isinstance(annotation, Mapping) and (
            field is None or annotation.get("field") == field
        ):
            return annotation
    return {}


def _as_list(value: Any) -> list[Any]:
    if value is None:
        return []
    return value if isinstance(value, list) else [value]


def _sample_identities(sample: Mapping[str, Any] | None) -> set[str]:
    if not isinstance(sample, Mapping):
        return set()
    identities = {str(sample.get("iid") or "")}
    accession = sample.get("accession")
    for item in accession if isinstance(accession, list) else [accession]:
        value = item.get("value") if isinstance(item, Mapping) else item
        if value not in (None, ""):
            identities.add(str(value))
    identities.discard("")
    return identities


def _slug(value: Any) -> str:
    return re.sub(r"_+", "_", re.sub(r"[^a-z0-9]+", "_", str(value).lower())).strip("_")
