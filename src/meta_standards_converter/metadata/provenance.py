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


"""Repository identity from explicit provenance, without accession translation."""
def repository_name(data):
    fmt = str(data.get("source", {}).get("format", "")).casefold()
    if fmt in {"sra", "ena"}:
        return fmt.upper()
    if fmt in {"geo", "geo miniml", "miniml"}:
        return "GEO"
    if fmt in {"arrayexpress", "biostudies"}:
        return "ArrayExpress"
    series = data.get("series", {})
    series = series if isinstance(series, list) else [series]
    identifiers = [v for item in series for v in [str(item.get("iid", "")),
                   *[str(a.get("value", "")) for a in item.get("accession", [])]]]
    if fmt == "mage-tab" and any(v.startswith(("E-MTAB-", "E-GEOD-")) for v in identifiers):
        return "ArrayExpress"
    if not fmt and any(v.startswith("GSE") and v[3:].isdigit() for v in identifiers):
        return "GEO"
    return None
