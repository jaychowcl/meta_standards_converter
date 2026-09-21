# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Validate and apply Agentic Curator harmonization override profiles."""

from __future__ import annotations

import copy
import re
from dataclasses import dataclass
from typing import Any, Mapping, Sequence

from meta_standards_converter.miniml import MINiMLCodec, harmonized_value_mappings


_FIXED_DESTINATIONS = {
    "organism",
    "organism_part",
    "developmental_stage",
    "disease",
    "genotype",
    "source",
    "biomaterial_provider",
    "material_type",
    "molecule",
}
_CHARACTERISTIC_TAGS = {
    "organism_part": "organism part",
    "developmental_stage": "developmental stage",
    "disease": "disease",
    "genotype": "genotype",
}
_DYNAMIC_DESTINATION = re.compile(r"^characteristics\.([a-z][a-z0-9_]*[a-z0-9]|[a-z])$")
_SOURCE_FIELD = re.compile(r"^[a-z][a-z0-9_]*$")


@dataclass(frozen=True)
class HarmonizationSelection:
    sample_accession: str
    destination: str
    value: str
    identifier: str | None
    ontology: str | None
    source_field: str
    hierarchy_depth: int | None
    status: str = "selected"


@dataclass(frozen=True)
class HarmonizationResolution:
    packages: tuple[dict[str, Any], ...]
    profile: Mapping[str, Any] | None = None
    selections: tuple[HarmonizationSelection, ...] = ()
    warnings: tuple[str, ...] = ()
    enabled: bool = False
    applied: bool = False


def resolve_harmonization_overrides(
    packages: Sequence[Mapping[str, Any]],
    profile: Mapping[str, Any] | None,
    *,
    enabled: bool,
) -> HarmonizationResolution:
    """Return a harmonization-aware deep copy derived from typed annotations."""
    if profile is not None and not isinstance(profile, Mapping):
        raise ValueError("replacement_profile must be an object")
    original = tuple(copy.deepcopy(dict(package)) for package in packages)
    if not enabled or profile is None:
        return HarmonizationResolution(original, profile=profile, enabled=enabled)
    try:
        normalized = validate_harmonization_overrides(profile)
    except ValueError as error:
        return HarmonizationResolution(
            original,
            profile=copy.deepcopy(profile),
            warnings=(f"Harmonization overrides disabled: {error}",),
            enabled=True,
        )

    raw = tuple(_canonical_package(package) for package in packages)

    selections: list[HarmonizationSelection] = []
    for package in raw:
        for sample in _as_list(package.get("sample")):
            if not isinstance(sample, dict):
                continue
            sample_accession = str(sample.get("iid") or _value(sample.get("accession")) or "")
            for channel in _as_list(sample.get("channel")):
                if not isinstance(channel, dict):
                    continue
                for destination, sources in normalized["replacements"].items():
                    selected_source = None
                    selected_values: list[dict[str, Any]] = []
                    for source in sources:
                        values = _harmonized_values(channel, source)
                        if values:
                            selected_source = source
                            selected_values = values
                            break
                    if selected_source is None:
                        continue
                    _apply_destination(channel, destination, selected_values)
                    selections.extend(
                        HarmonizationSelection(
                            sample_accession=sample_accession,
                            destination=destination,
                            value=str(item["value"]),
                            identifier=_optional_string(item.get("id")),
                            ontology=_optional_string(item.get("onto")),
                            source_field=selected_source,
                            hierarchy_depth=_depth(item.get("hierarchy_depth")),
                        )
                        for item in selected_values
                    )
        # Apply the same policy to explicitly bound path-local source evidence.
        # Resolving each occurrence independently preserves differing path metadata.
        from meta_standards_converter.magetab.harmonized import _path_binding
        for path in package.get("series", {}).get("assay_paths", []):
            binding = _path_binding(package, path)
            if binding is None:
                continue
            node, channel = binding
            for destination, sources in normalized["replacements"].items():
                for source in sources:
                    values = _harmonized_values(node, source) or _harmonized_values(channel, source)
                    if values:
                        _apply_path_destination(node, destination, values)
                        break
    return HarmonizationResolution(
        raw,
        profile=normalized,
        selections=tuple(selections),
        enabled=True,
        applied=True,
    )


def _apply_path_destination(node, destination, values):
    if destination == "organism":
        # Assay nodes express organisms as characteristics, not Channel.organism.
        _apply_destination(node, "characteristics.organism", values)
    elif destination == "source":
        comments = node.setdefault("comments", [])
        comments[:] = [c for c in comments if c.get("name") != "Sample_source_name"]
        comments.extend({"name": "Sample_source_name", "value": item["value"]} for item in values)
    elif destination in {"molecule", "material_type"}:
        node["material_type"] = _ontology_container(values[0])
        if len(values) > 1:
            node["material_type"]["value"] = " | ".join(str(item["value"]) for item in values)
    elif destination == "biomaterial_provider":
        node["provider"] = " | ".join(str(item["value"]) for item in values)
    else:
        _apply_destination(node, destination, values)


def _canonical_package(package: Mapping[str, Any]) -> dict[str, Any]:
    """Give every consumer one canonical 3.0 view, including readable 2.0 input."""
    version = package.get("miniml_schema_version")
    if version in {"2.0", "3.0"}:
        return MINiMLCodec().encode(MINiMLCodec().decode(package).package)
    return copy.deepcopy(dict(package))


def validate_harmonization_overrides(profile: Mapping[str, Any]) -> dict[str, Any]:
    if not isinstance(profile, Mapping):
        raise ValueError("profile must be an object")
    if profile.get("schema_version") != "1.0":
        raise ValueError("schema_version must be '1.0'")
    replacements = profile.get("replacements")
    if not isinstance(replacements, Mapping) or not replacements:
        raise ValueError("replacements must be a non-empty object")
    normalized: dict[str, list[str]] = {}
    for destination, sources in replacements.items():
        destination = str(destination)
        if destination not in _FIXED_DESTINATIONS and not _DYNAMIC_DESTINATION.fullmatch(destination):
            raise ValueError(f"unsupported destination {destination!r}")
        if not isinstance(sources, list) or not sources:
            raise ValueError(f"destination {destination!r} requires a non-empty source list")
        rendered = [str(source) for source in sources]
        if any(not _SOURCE_FIELD.fullmatch(source) for source in rendered):
            raise ValueError(f"destination {destination!r} contains an invalid source field")
        if len(set(rendered)) != len(rendered):
            raise ValueError(f"destination {destination!r} contains duplicate source fields")
        normalized[destination] = rendered
    return {"schema_version": "1.0", "replacements": normalized}


def _harmonized_values(channel: Mapping[str, Any], source: str) -> list[dict[str, Any]]:
    values: list[dict[str, Any]] = []
    containers = [channel]
    for field in ("source", "molecule", "material_type"):
        item = channel.get(field)
        if isinstance(item, Mapping):
            containers.append(item)
    containers.extend(item for item in _as_list(channel.get("organism")) if isinstance(item, Mapping))
    containers.extend(item for item in _as_list(channel.get("characteristics")) if isinstance(item, Mapping))
    for container in containers:
        for annotation in harmonized_value_mappings(container):
            if not isinstance(annotation, Mapping) or annotation.get("field") != source:
                continue
            _append_value(values, {
                "value": annotation.get("value"),
                "id": annotation.get("term_accession_number"),
                "onto": annotation.get("term_source_ref"),
                "hierarchy_depth": annotation.get("hierarchy_depth"),
            })
    for annotation in harmonized_value_mappings(
        [
            item
            for item in _as_list(channel.get("characteristics"))
            if isinstance(item, Mapping)
        ]
    ):
        if annotation.get("field") != source:
            continue
        _append_value(values, {
            "value": annotation.get("value"),
            "id": annotation.get("term_accession_number"),
            "onto": annotation.get("term_source_ref"),
            "hierarchy_depth": annotation.get("hierarchy_depth"),
        })
    return values


def _append_value(values: list[dict[str, Any]], item: Mapping[str, Any]) -> None:
    value = item.get("value")
    if value is None or not str(value).strip():
        return
    normalized = {
        "value": value,
        "id": item.get("id"),
        "onto": item.get("onto"),
        "hierarchy_depth": item.get("hierarchy_depth"),
    }
    identity = tuple(str(normalized.get(key) or "") for key in normalized)
    if not any(tuple(str(existing.get(key) or "") for key in existing) == identity for existing in values):
        values.append(normalized)


def _apply_destination(channel: dict[str, Any], destination: str, values: list[dict[str, Any]]) -> None:
    if destination in {"organism", "source", "molecule", "material_type"}:
        from meta_standards_converter.miniml.harmonization import iter_harmonized_values, named_harmonized_rows
        from dataclasses import replace
        rows = channel.setdefault("characteristics", [])
        used = {(v.field, v.index) for v in iter_harmonized_values(rows)}
        for container in _as_list(channel.get(destination)):
            if not isinstance(container, Mapping):
                continue
            for value in iter_harmonized_values(container):
                while (value.field, value.index) in used:
                    value = replace(value, index=value.index + 1)
                used.add((value.field, value.index))
                rows.extend(named_harmonized_rows([value]))
    if destination == "organism":
        channel["organism"] = [_ontology_container(item, taxon=True) for item in values]
        return
    if destination in {"source", "molecule", "biomaterial_provider", "material_type"}:
        rendered = [item["value"] for item in values]
        channel[destination] = rendered[0] if len(rendered) == 1 else rendered
        return
    tag = _CHARACTERISTIC_TAGS.get(destination)
    if tag is None:
        match = _DYNAMIC_DESTINATION.fullmatch(destination)
        tag = match.group(1) if match else destination
    rows = [item for item in _as_list(channel.get("characteristics")) if isinstance(item, dict)]
    from .fields import BIOLOGICAL_ALIASES
    aliases = BIOLOGICAL_ALIASES.get(destination, {tag.casefold().replace(" ", "_")})
    rows = [item for item in rows if re.sub(r"[^a-z0-9]+", "_",
            str(item.get("name") or item.get("tag") or "").casefold()).strip("_") not in aliases]
    rows.extend({"name": tag, **_ontology_container(item)} for item in values)
    channel["characteristics"] = rows


def _ontology_container(item: Mapping[str, Any], *, taxon: bool = False) -> dict[str, Any]:
    result = {"value": item["value"]}
    identifier = _optional_string(item.get("id"))
    ontology = _optional_string(item.get("onto"))
    if taxon and identifier and ontology and ontology.casefold() == "ncbitaxon":
        result["taxid"] = identifier.rsplit(":", 1)[-1].rsplit("_", 1)[-1]
    if ontology:
        result["term_source_ref"] = ontology
    if identifier:
        result["term_accession_number"] = identifier
    return result


def _as_list(value: Any) -> list[Any]:
    if value is None:
        return []
    return value if isinstance(value, list) else [value]


def _value(value: Any) -> Any:
    for item in _as_list(value):
        if isinstance(item, Mapping) and item.get("value") is not None:
            return item["value"]
        if item is not None:
            return item
    return None


def _ordinal(values: list[Any], ordinal: int) -> Any:
    return values[ordinal] if ordinal < len(values) else None


def _optional_string(value: Any) -> str | None:
    return str(value) if value is not None and str(value).strip() else None


def _depth(value: Any) -> int | None:
    if value is None or value == "":
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None
