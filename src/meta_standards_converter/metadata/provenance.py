# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Deterministic converter-facing views of retained MSC patch provenance."""

from __future__ import annotations

from typing import Any, Mapping

from meta_standards_converter.miniml import iter_harmonization_operations


PROVENANCE_FIELDS = (
    "value",
    "id",
    "ontology",
    "hierarchy_depth",
    "source_field",
    "source_label",
    "source_path",
    "match_kind",
    "patch_id",
)


def patch_provenance_columns(
    package: Mapping[str, Any],
    sample: Mapping[str, Any] | str,
    *,
    occupied: set[str] | None = None,
    operations: tuple[dict, ...] | None = None,
) -> dict[str, Any]:
    """Return stable indexed ``msc.harmonization.*`` columns for one sample."""

    sample_id = sample if isinstance(sample, str) else str(sample.get("iid") or "")
    result: dict[str, Any] = {}
    counts: dict[str, int] = {}
    used_prefixes = {
        key.rsplit(".", 1)[0]
        for key in (occupied or set())
        if key.startswith("msc.harmonization.") and "." in key
    }
    if operations is None:
        if not package.get("extensions", {}).get("msc_harmonization"):
            return {}
        operations = iter_harmonization_operations(package, sample=sample_id)
    for operation in operations:
        harmonized = operation["harmonized_value"]
        field = str(harmonized["field"])
        index = counts.get(field, 0)
        while True:
            suffix = "" if index == 0 else f"({index})"
            prefix = f"msc.harmonization.{field}{suffix}"
            if prefix not in used_prefixes:
                break
            index += 1
        counts[field] = index + 1
        used_prefixes.add(prefix)
        evidence = operation.get("source_evidence")
        evidence = evidence if isinstance(evidence, Mapping) else {}
        values = {
            "value": harmonized.get("value"),
            "id": harmonized.get("term_accession_number"),
            "ontology": harmonized.get("term_source_ref"),
            "hierarchy_depth": harmonized.get("hierarchy_depth"),
            "source_field": evidence.get("source_field"),
            "source_label": evidence.get("source_label"),
            "source_path": evidence.get("source_value_path"),
            "match_kind": evidence.get("match_kind"),
            "patch_id": operation.get("patch_id"),
        }
        for name in PROVENANCE_FIELDS:
            result[f"{prefix}.{name}"] = values[name]
    return result


__all__ = ["PROVENANCE_FIELDS", "patch_provenance_columns"]
