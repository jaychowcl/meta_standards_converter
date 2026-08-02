# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Project enriched MAGE-TAB assay attributes without replacing raw values."""

from __future__ import annotations

import re
from typing import Any, Mapping


PARAMETER_FIELDS = (
    "value", "unit", "hz_value", "hz_value_id", "hz_value_onto",
    "hz_unit", "hz_unit_id", "hz_unit_onto",
)


def _parameter_summary(
    package: Mapping[str, Any], sample: Mapping[str, Any]
) -> dict[str, tuple[Any, ...]]:
    """Return deterministic per-parameter values for one bound sample."""
    result: dict[str, list[Any]] = {}
    for row in _parameter_rows(package, sample=sample):
        if row["attribute_type"] != "parameter value":
            continue
        slug = _slug(row["name"])
        if not slug:
            continue
        for field in PARAMETER_FIELDS:
            value = row.get(field)
            if value in (None, ""):
                continue
            key = f"msc.mage_tab.parameter.{slug}.{_column_field(field)}"
            values = result.setdefault(key, [])
            if value not in values:
                values.append(value)
    return {key: tuple(values) for key, values in result.items()}


def _parameter_rows(
    package: Mapping[str, Any], sample: Mapping[str, Any] | None = None
) -> list[dict[str, Any]]:
    """Return lossless typed-attribute occurrence rows, optionally sample-bound."""
    model = _model(package)
    if model is None:
        return []
    identities = _sample_identities(sample) if sample is not None else set()
    rows = []
    for assay in model.get("assay_paths", []):
        if not isinstance(assay, Mapping):
            continue
        binding = assay.get("binding") if isinstance(assay.get("binding"), Mapping) else {}
        if identities and not identities.intersection(
            str(value) for value in binding.values() if value not in (None, "")
        ):
            continue
        for step in assay.get("steps", []):
            if not isinstance(step, Mapping) or step.get("kind") != "attribute":
                continue
            rows.append({
                "sample_accession": next(iter(sorted(identities)), ""),
                "assay_path_id": assay.get("id", ""),
                "sdrf": assay.get("sdrf", ""),
                "row_index": assay.get("row_index", -1),
                "column_index": step.get("column_index", -1),
                "attribute_type": step.get("attribute_type", ""),
                "name": step.get("name", ""),
                "value": step.get("value", ""),
                "unit": step.get("unit", ""),
                "term_source_ref": step.get("term_source_ref", ""),
                "term_accession_number": step.get("term_accession_number", ""),
                "hz_field": step.get("hz_field", ""),
                "hz_value": step.get("hz_value", ""),
                "hz_value_id": step.get("hz_value_id", ""),
                "hz_value_onto": step.get("hz_value_onto", ""),
                "hz_unit": step.get("hz_unit", ""),
                "hz_unit_id": step.get("hz_unit_id", ""),
                "hz_unit_onto": step.get("hz_unit_onto", ""),
            })
    return rows


def _model(package: Mapping[str, Any]) -> Mapping[str, Any] | None:
    mage_tab = package.get("mage_tab")
    model = mage_tab.get("model") if isinstance(mage_tab, Mapping) else None
    return model if isinstance(model, Mapping) and model.get("schema_version") == 1 else None


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


def _column_field(field: str) -> str:
    return field[:-5] + "ontology" if field.endswith("_onto") else field
