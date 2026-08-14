# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Canonical MINiML 3.0 harmonized-value wire helpers."""

from __future__ import annotations

from dataclasses import dataclass
import re
from typing import Any, Iterable, Mapping, MutableMapping, MutableSequence, Sequence


_HZ_PATTERN = re.compile(
    r"^hz_(?P<body>[a-z][a-z0-9_]*?)(?P<qualifier>_hierarchy_depth|_id|_onto)?"
    r"(?:\((?P<index>[1-9][0-9]*)\))?$"
)
_LEGACY_COLLISION_PATTERN = re.compile(r"^hz_[a-z][a-z0-9_]*_[1-9][0-9]*$")


def _model_error(message: str) -> ValueError:
    # Import lazily so this module remains the leaf dependency of model.py.
    from .model import MINiMLModelError

    return MINiMLModelError(message)


@dataclass(frozen=True)
class HarmonizedValue:
    """One typed harmonized value represented by a flat ``hz_*`` wire group."""

    field: str
    value: Any
    term_source_ref: str | None = None
    term_accession_number: str | None = None
    hierarchy_depth: int | None = None
    index: int = 0

    def __post_init__(self) -> None:
        if not isinstance(self.field, str) or re.fullmatch(
            r"[a-z][a-z0-9_]*", self.field
        ) is None:
            raise _model_error(f"invalid harmonized field {self.field!r}")
        if self.value is None or (isinstance(self.value, str) and not self.value.strip()):
            raise _model_error("harmonized value must be nonblank")
        if isinstance(self.index, bool) or not isinstance(self.index, int) or self.index < 0:
            raise _model_error("harmonized value index must be a non-negative integer")
        if self.hierarchy_depth is not None and (
            isinstance(self.hierarchy_depth, bool)
            or not isinstance(self.hierarchy_depth, int)
            or self.hierarchy_depth < 0
        ):
            raise _model_error(
                "harmonized hierarchy depth must be a non-negative integer"
            )
        for name, value in (
            ("term_source_ref", self.term_source_ref),
            ("term_accession_number", self.term_accession_number),
        ):
            if value is not None and not isinstance(value, str):
                raise _model_error(f"harmonized {name} must be a string")

    @property
    def suffix(self) -> str:
        return "" if self.index == 0 else f"({self.index})"

    def to_mapping(self) -> dict[str, Any]:
        prefix = f"hz_{self.field}"
        result: dict[str, Any] = {f"{prefix}{self.suffix}": self.value}
        if self.term_accession_number is not None:
            result[f"{prefix}_id{self.suffix}"] = self.term_accession_number
        if self.term_source_ref is not None:
            result[f"{prefix}_onto{self.suffix}"] = self.term_source_ref
        if self.hierarchy_depth is not None:
            result[f"{prefix}_hierarchy_depth{self.suffix}"] = self.hierarchy_depth
        return result

    def to_annotation_mapping(self) -> dict[str, Any]:
        """Return a consumer-neutral semantic mapping for projection code."""
        result = {"field": self.field, "value": self.value, "index": self.index}
        if self.term_accession_number is not None:
            result["term_accession_number"] = self.term_accession_number
        if self.term_source_ref is not None:
            result["term_source_ref"] = self.term_source_ref
        if self.hierarchy_depth is not None:
            result["hierarchy_depth"] = self.hierarchy_depth
        return result


def is_harmonized_key(value: Any) -> bool:
    return isinstance(value, str) and value.startswith("hz_")


def parse_harmonized_key(key: str) -> tuple[str, str, int]:
    if _LEGACY_COLLISION_PATTERN.fullmatch(key):
        raise _model_error(
            f"invalid harmonized field {key!r}; collisions use parentheses"
        )
    match = _HZ_PATTERN.fullmatch(key)
    if match is None:
        raise _model_error(f"invalid harmonized field {key!r}")
    qualifier = match.group("qualifier") or ""
    field = match.group("body")
    role = {
        "": "value",
        "_id": "id",
        "_onto": "onto",
        "_hierarchy_depth": "hierarchy_depth",
    }[qualifier]
    index = int(match.group("index") or 0)
    return field, role, index


def parse_harmonized_mapping(value: Mapping[str, Any]) -> tuple[HarmonizedValue, ...]:
    groups: dict[tuple[str, int], dict[str, Any]] = {}
    for raw_key, item in value.items():
        if not is_harmonized_key(raw_key):
            continue
        field, role, index = parse_harmonized_key(str(raw_key))
        group = groups.setdefault((field, index), {})
        if role in group:
            raise _model_error(f"duplicate harmonized {role} for hz_{field}")
        group[role] = item

    result: list[HarmonizedValue] = []
    for (field, index), group in sorted(
        groups.items(), key=lambda item: (item[0][0], item[0][1])
    ):
        if "value" not in group:
            raise _model_error(
                f"harmonized companions for hz_{field}({index}) exist without a "
                "corresponding value"
            )
        result.append(
            HarmonizedValue(
                field=field,
                value=group["value"],
                term_source_ref=(
                    None if group.get("onto") is None else str(group["onto"])
                ),
                term_accession_number=(
                    None if group.get("id") is None else str(group["id"])
                ),
                hierarchy_depth=group.get("hierarchy_depth"),
                index=index,
            )
        )
    return tuple(result)


def harmonized_mapping(values: Iterable[HarmonizedValue]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    identities: set[tuple[str, int]] = set()
    for value in values:
        if not isinstance(value, HarmonizedValue):
            raise TypeError("harmonized values must be HarmonizedValue instances")
        identity = (value.field, value.index)
        if identity in identities:
            raise _model_error(
                f"duplicate harmonized value group hz_{value.field}{value.suffix}"
            )
        identities.add(identity)
        result.update(value.to_mapping())
    return result


def named_harmonized_rows(
    values: Iterable[HarmonizedValue], *, name_key: str = "name"
) -> list[dict[str, Any]]:
    return [
        {name_key: key, "value": item}
        for value in values
        for key, item in value.to_mapping().items()
    ]


def iter_harmonized_values(
    value: Mapping[str, Any] | Sequence[Mapping[str, Any]] | None,
) -> tuple[HarmonizedValue, ...]:
    """Read harmonized groups from an object or a named/tag-value row list."""

    if isinstance(value, Mapping):
        return parse_harmonized_mapping(value)
    if not isinstance(value, (list, tuple)):
        return ()

    flattened: dict[str, Any] = {}
    for row in value:
        if not isinstance(row, Mapping):
            continue
        key = row.get("name", row.get("tag"))
        if not is_harmonized_key(key):
            continue
        if key in flattened:
            raise _model_error(f"duplicate harmonized row {key!r}")
        flattened[str(key)] = row.get("value")
    return parse_harmonized_mapping(flattened)


def next_harmonized_index(
    values: Iterable[HarmonizedValue], *, field: str
) -> int:
    indexes = {value.index for value in values if value.field == field}
    index = 0
    while index in indexes:
        index += 1
    return index


def append_harmonized_value(
    destination: MutableMapping[str, Any] | MutableSequence[MutableMapping[str, Any]],
    value: HarmonizedValue,
    *,
    name_key: str | None = None,
) -> HarmonizedValue:
    """Validate and append one harmonized group without replacing raw evidence.

    Mapping destinations receive flat ``hz_*`` members. Sequence destinations
    receive named rows and therefore require ``name_key`` to be ``"name"`` or
    ``"tag"``. Exact semantic duplicates are idempotent; distinct values for
    the same field receive the next available aligned parenthesized index.
    """

    if not isinstance(value, HarmonizedValue):
        raise TypeError("value must be a HarmonizedValue")
    if isinstance(destination, MutableMapping):
        if name_key is not None:
            raise _model_error("name_key is only valid for named-row destinations")
        existing = iter_harmonized_values(destination)
    elif isinstance(destination, MutableSequence):
        if name_key not in {"name", "tag"}:
            raise _model_error(
                "name_key must be 'name' or 'tag' for named-row destinations"
            )
        existing = iter_harmonized_values(destination)
    else:
        raise TypeError("destination must be a mutable mapping or named-row sequence")

    identity = (
        value.field,
        value.value,
        value.term_source_ref,
        value.term_accession_number,
        value.hierarchy_depth,
    )
    for item in existing:
        if (
            item.field,
            item.value,
            item.term_source_ref,
            item.term_accession_number,
            item.hierarchy_depth,
        ) == identity:
            return item

    appended = HarmonizedValue(
        field=value.field,
        value=value.value,
        term_source_ref=value.term_source_ref,
        term_accession_number=value.term_accession_number,
        hierarchy_depth=value.hierarchy_depth,
        index=next_harmonized_index(existing, field=value.field),
    )
    if isinstance(destination, MutableMapping):
        destination.update(appended.to_mapping())
    else:
        destination.extend(named_harmonized_rows([appended], name_key=name_key))
    return appended


def harmonized_value_mappings(
    value: Mapping[str, Any] | Sequence[Mapping[str, Any]] | None,
) -> tuple[dict[str, Any], ...]:
    return tuple(item.to_annotation_mapping() for item in iter_harmonized_values(value))


__all__ = [
    "HarmonizedValue",
    "append_harmonized_value",
    "harmonized_mapping",
    "harmonized_value_mappings",
    "is_harmonized_key",
    "iter_harmonized_values",
    "named_harmonized_rows",
    "next_harmonized_index",
    "parse_harmonized_key",
    "parse_harmonized_mapping",
]
