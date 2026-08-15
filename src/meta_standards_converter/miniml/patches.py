# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""MSC-owned harmonization patch and retained provenance contracts.

Patch schema 3.1 adds bounded occurrence provenance to the schema 3.0
harmonized-value additions.  Applied package-local fragments use extension
schema 1.0 under ``extensions.msc_harmonization``.
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
import hashlib
import json
import re
from types import MappingProxyType
from typing import Any, Iterable, Mapping
import unicodedata

from .codec import MINiMLCodec
from .harmonization import HarmonizedValue, append_harmonized_value
from .model import MINIML_SCHEMA_VERSION, MINiMLModelError, MINiMLPackage


PATCH_SCHEMA_VERSION = "3.1"
PATCH_EXTENSION_SCHEMA_VERSION = "1.0"
PATCH_EXTENSION_KEY = "msc_harmonization"
_SHA256 = re.compile(r"[0-9a-f]{64}")
_MATCH_KINDS = {"exact_value", "exact_span", "interpreted"}
_SOURCE_EVIDENCE_KEYS = {
    "target_id",
    "source_field",
    "source_label",
    "source_field_path",
    "source_value_path",
    "derivation",
    "match_kind",
}
_SOURCE_LIMITS = {
    "target_id": 160,
    "source_field": 128,
    "source_label": 1000,
    "source_field_path": 2000,
    "source_value_path": 2000,
    "derivation": 128,
}


def _plain(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(key): _plain(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_plain(item) for item in value]
    return deepcopy(value)


def canonical_miniml_document(
    value: MINiMLPackage | Mapping[str, Any] | Iterable[Mapping[str, Any]],
) -> dict[str, Any] | list[dict[str, Any]]:
    """Strictly decode and encode one package or a non-empty package sequence."""

    codec = MINiMLCodec()
    if isinstance(value, MINiMLPackage):
        return codec.encode(codec.decode(value).package)
    if isinstance(value, Mapping):
        return codec.encode(codec.decode(value).package)
    if isinstance(value, (list, tuple)):
        if not value:
            raise MINiMLModelError("MINiML package collection must not be empty")
        return codec.encode_many(codec.decode_many(value).packages)
    raise TypeError("MINiML document must be a package mapping or package sequence")


def _without_retained_patches(value: Any) -> Any:
    result = deepcopy(value)
    packages = result if isinstance(result, list) else [result]
    for package in packages:
        if not isinstance(package, dict):
            continue
        extensions = package.get("extensions")
        if not isinstance(extensions, dict):
            continue
        extensions.pop(PATCH_EXTENSION_KEY, None)
        if not extensions:
            package.pop("extensions", None)
    return result


def miniml_source_fingerprint(
    value: MINiMLPackage | Mapping[str, Any] | Iterable[Mapping[str, Any]],
) -> str:
    """Hash canonical biological content, excluding only retained patches."""

    canonical = _without_retained_patches(canonical_miniml_document(value))
    payload = json.dumps(
        canonical, sort_keys=True, separators=(",", ":"), ensure_ascii=False
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _validated_pointer(value: Any, *, name: str) -> str:
    if not isinstance(value, str) or not value.startswith("/"):
        raise MINiMLModelError(f"{name} must be a JSON pointer")
    if len(value) > 2000:
        raise MINiMLModelError(f"{name} exceeds 2000 characters")
    return value


def _validated_source_evidence(value: Any) -> dict[str, str]:
    if not isinstance(value, Mapping) or set(value) != _SOURCE_EVIDENCE_KEYS:
        raise MINiMLModelError(
            "source_evidence requires target_id, source_field, source_label, "
            "source_field_path, source_value_path, derivation and match_kind"
        )
    result: dict[str, str] = {}
    for key in _SOURCE_EVIDENCE_KEYS:
        item = value.get(key)
        if not isinstance(item, str) or not item.strip():
            raise MINiMLModelError(f"source_evidence {key} must be nonblank")
        item = item.strip()
        limit = _SOURCE_LIMITS.get(key)
        if limit is not None and len(item) > limit:
            raise MINiMLModelError(f"source_evidence {key} exceeds {limit} characters")
        result[key] = item
    if result["match_kind"] not in _MATCH_KINDS:
        raise MINiMLModelError(
            "source_evidence match_kind must be exact_value, exact_span or interpreted"
        )
    _validated_pointer(result["source_field_path"], name="source_field_path")
    _validated_pointer(result["source_value_path"], name="source_value_path")
    return result


def _validated_harmonized_value(value: Any) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise MINiMLModelError("MINiML harmonized_value must be an object")
    allowed = {
        "field",
        "value",
        "term_source_ref",
        "term_accession_number",
        "hierarchy_depth",
    }
    if set(value) - allowed:
        raise MINiMLModelError("MINiML harmonized_value contains unsupported fields")
    typed = HarmonizedValue(
        field=str(value.get("field") or ""),
        value=value.get("value"),
        term_source_ref=value.get("term_source_ref"),
        term_accession_number=value.get("term_accession_number"),
        hierarchy_depth=value.get("hierarchy_depth"),
    )
    result = {"field": typed.field, "value": typed.value}
    for key in ("term_source_ref", "term_accession_number", "hierarchy_depth"):
        item = getattr(typed, key)
        if item is not None:
            result[key] = item
    return result


def _validated_operation(value: Any, *, schema_version: str) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise MINiMLModelError("MINiML harmonization operation must be an object")
    allowed = {"path", "harmonized_value"}
    if schema_version == PATCH_SCHEMA_VERSION:
        allowed.add("source_evidence")
    if set(value) - allowed or not {"path", "harmonized_value"}.issubset(value):
        raise MINiMLModelError(
            "MINiML harmonization operation requires path and harmonized_value"
        )
    result: dict[str, Any] = {
        "path": _validated_pointer(value.get("path"), name="operation path"),
        "harmonized_value": _validated_harmonized_value(
            value.get("harmonized_value")
        ),
    }
    if "source_evidence" in value:
        result["source_evidence"] = _validated_source_evidence(
            value["source_evidence"]
        )
    return result


@dataclass(frozen=True)
class MINiMLHarmonizationPatch:
    """Immutable document-level harmonized additions and source provenance."""

    base_sha256: str
    adds: tuple[Mapping[str, Any], ...]
    schema_version: str = PATCH_SCHEMA_VERSION
    miniml_schema_version: str = MINIML_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if self.schema_version not in {"3.0", PATCH_SCHEMA_VERSION}:
            raise MINiMLModelError(
                "MINiML harmonization patch schema_version must be '3.0' or '3.1'"
            )
        if self.miniml_schema_version != MINIML_SCHEMA_VERSION:
            raise MINiMLModelError("MINiML harmonization patch requires MINiML 3.0")
        fingerprint = str(self.base_sha256 or "").lower()
        if _SHA256.fullmatch(fingerprint) is None:
            raise MINiMLModelError(
                "MINiML harmonization patch base_sha256 must be SHA-256"
            )
        if not isinstance(self.adds, (tuple, list)):
            raise TypeError("MINiML harmonization patch adds must be a sequence")
        normalized: list[Mapping[str, Any]] = []
        for operation in self.adds:
            item = _validated_operation(operation, schema_version=self.schema_version)
            nested = dict(item)
            nested["harmonized_value"] = MappingProxyType(
                nested["harmonized_value"]
            )
            if "source_evidence" in nested:
                nested["source_evidence"] = MappingProxyType(
                    nested["source_evidence"]
                )
            normalized.append(MappingProxyType(nested))
        object.__setattr__(self, "base_sha256", fingerprint)
        object.__setattr__(self, "adds", tuple(normalized))

    @classmethod
    def from_mapping(cls, value: Mapping[str, Any]) -> "MINiMLHarmonizationPatch":
        if not isinstance(value, Mapping):
            raise TypeError("MINiML harmonization patch must be an object")
        allowed = {
            "schema_version",
            "miniml_schema_version",
            "base_sha256",
            "adds",
        }
        if set(value) - allowed:
            raise MINiMLModelError("MINiML harmonization patch contains unsupported fields")
        adds = value.get("adds")
        if not isinstance(adds, list):
            raise MINiMLModelError("MINiML harmonization patch adds must be a list")
        return cls(
            base_sha256=value.get("base_sha256"),
            adds=tuple(adds),
            schema_version=value.get("schema_version"),
            miniml_schema_version=value.get("miniml_schema_version"),
        )

    @property
    def patch_id(self) -> str:
        payload = json.dumps(
            self.to_mapping(), sort_keys=True, separators=(",", ":"), ensure_ascii=False
        ).encode("utf-8")
        return hashlib.sha256(payload).hexdigest()

    def matches(self, document: Any) -> bool:
        return self.base_sha256 == miniml_source_fingerprint(document)

    def to_mapping(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "miniml_schema_version": self.miniml_schema_version,
            "base_sha256": self.base_sha256,
            "adds": [_plain(item) for item in self.adds],
        }


def _unescape(segment: str) -> str:
    return segment.replace("~1", "/").replace("~0", "~")


def _resolve_pointer(document: Any, pointer: str) -> Any:
    current = document
    for raw in pointer.split("/")[1:]:
        segment = _unescape(raw)
        if isinstance(current, Mapping):
            if segment not in current:
                raise MINiMLModelError(f"MINiML harmonization path does not exist: {pointer}")
            current = current[segment]
        elif isinstance(current, list):
            if not segment.isdecimal() or int(segment) >= len(current):
                raise MINiMLModelError(f"MINiML harmonization path does not exist: {pointer}")
            current = current[int(segment)]
        else:
            raise MINiMLModelError(f"MINiML harmonization path does not exist: {pointer}")
    return current


def _resolve_pointer_parent(document: Any, pointer: str) -> tuple[Any, str | int]:
    parts = pointer.split("/")[1:]
    if not parts:
        raise MINiMLModelError("MINiML harmonization path must not target the root")
    parent_path = "/" + "/".join(parts[:-1]) if len(parts) > 1 else ""
    parent = document if not parent_path else _resolve_pointer(document, parent_path)
    final = _unescape(parts[-1])
    return parent, int(final) if isinstance(parent, list) and final.isdecimal() else final


def _normalized_text(value: Any) -> str:
    return " ".join(unicodedata.normalize("NFKC", str(value)).split()).casefold()


def _validate_source_claim(document: Any, operation: Mapping[str, Any]) -> None:
    evidence = operation.get("source_evidence")
    if not isinstance(evidence, Mapping):
        return
    _resolve_pointer(document, str(evidence["source_field_path"]))
    source = _resolve_pointer(document, str(evidence["source_value_path"]))
    kind = evidence["match_kind"]
    if kind == "interpreted":
        return
    source_text = _normalized_text(source)
    label = _normalized_text(evidence["source_label"])
    valid = source_text == label
    if kind == "exact_span":
        valid = re.search(rf"(?<!\w){re.escape(label)}(?!\w)", source_text) is not None
    if not valid:
        raise MINiMLModelError(
            f"source evidence {kind} claim is not present at "
            f"{evidence['source_value_path']}"
        )


def _apply_add(document: Any, operation: Mapping[str, Any]) -> None:
    path = str(operation["path"])
    destination = _resolve_pointer(document, path)
    if not isinstance(destination, dict):
        raise MINiMLModelError(f"MINiML harmonization path is not an object: {path}")
    raw = operation["harmonized_value"]
    value = HarmonizedValue(
        field=str(raw["field"]),
        value=raw["value"],
        term_source_ref=raw.get("term_source_ref"),
        term_accession_number=raw.get("term_accession_number"),
        hierarchy_depth=raw.get("hierarchy_depth"),
    )
    parent, final = _resolve_pointer_parent(document, path)
    name_key = "name" if "name" in destination else "tag" if "tag" in destination else None
    if isinstance(parent, list) and isinstance(final, int) and name_key is not None:
        original_length = len(parent)
        append_harmonized_value(parent, value, name_key=name_key)
        if len(parent) == original_length:
            return
        added = parent[original_length:]
        del parent[original_length:]
        insertion = final + 1
        while insertion < len(parent):
            row = parent[insertion]
            key = row.get(name_key) if isinstance(row, Mapping) else None
            if not (isinstance(key, str) and key.startswith("hz_")):
                break
            insertion += 1
        parent[insertion:insertion] = added
        return
    append_harmonized_value(destination, value)


def _rebase_pointer(pointer: str, package_index: int, *, multi: bool) -> str:
    if not multi:
        return pointer
    prefix = f"/{package_index}"
    if pointer == prefix:
        return ""
    if not pointer.startswith(prefix + "/"):
        raise MINiMLModelError(
            f"operation path {pointer!r} is outside package {package_index}"
        )
    return pointer[len(prefix) :]


def _partition_operations(
    operations: Iterable[Mapping[str, Any]], *, package_count: int
) -> dict[int, list[dict[str, Any]]]:
    multi = package_count > 1
    result: dict[int, list[dict[str, Any]]] = {}
    for operation in operations:
        path = str(operation["path"])
        if multi:
            first = path.split("/", 2)[1] if path.startswith("/") else ""
            if not first.isdecimal() or int(first) >= package_count:
                raise MINiMLModelError(
                    "multi-package harmonization paths require a package index"
                )
            package_index = int(first)
        else:
            package_index = 0
        local = _plain(operation)
        local["path"] = _rebase_pointer(path, package_index, multi=multi)
        evidence = local.get("source_evidence")
        if isinstance(evidence, dict):
            for key in ("source_field_path", "source_value_path"):
                evidence[key] = _rebase_pointer(
                    evidence[key], package_index, multi=multi
                )
        result.setdefault(package_index, []).append(local)
    return result


def _package_patch_ids(package: Mapping[str, Any]) -> set[str]:
    return {
        str(fragment.get("patch_id"))
        for fragment in iter_harmonization_patches(package)
        if fragment.get("patch_id")
    }


def apply_miniml_harmonization_patch(
    document: MINiMLPackage | Mapping[str, Any] | Iterable[Mapping[str, Any]],
    patch: MINiMLHarmonizationPatch | Mapping[str, Any],
) -> dict[str, Any] | list[dict[str, Any]]:
    """Apply and retain one immutable patch without replacing raw metadata."""

    typed = (
        patch
        if isinstance(patch, MINiMLHarmonizationPatch)
        else MINiMLHarmonizationPatch.from_mapping(patch)
    )
    canonical = canonical_miniml_document(document)
    packages = canonical if isinstance(canonical, list) else [canonical]
    partitions = _partition_operations(typed.adds, package_count=len(packages))
    if partitions and all(
        typed.patch_id in _package_patch_ids(packages[index]) for index in partitions
    ):
        return canonical
    if not typed.matches(canonical):
        raise MINiMLModelError(
            "MINiML harmonization patch base document fingerprint differs"
        )
    for operation in typed.adds:
        _validate_source_claim(canonical, operation)

    result = deepcopy(canonical)
    result_packages = result if isinstance(result, list) else [result]
    for package_index, operations in sorted(partitions.items()):
        package = result_packages[package_index]
        source_sha256 = miniml_source_fingerprint(package)
        for operation in operations:
            _apply_add(package, operation)
        result_sha256 = miniml_source_fingerprint(package)
        extensions = package.setdefault("extensions", {})
        retained = extensions.setdefault(
            PATCH_EXTENSION_KEY,
            {"schema_version": PATCH_EXTENSION_SCHEMA_VERSION, "patches": []},
        )
        retained["patches"].append(
            {
                "schema_version": PATCH_EXTENSION_SCHEMA_VERSION,
                "patch_schema_version": typed.schema_version,
                "patch_id": typed.patch_id,
                "package_index": package_index,
                "source_sha256": source_sha256,
                "result_sha256": result_sha256,
                "status": "applied",
                "operations": operations,
            }
        )
    return canonical_miniml_document(result)


def validate_harmonization_extension_mapping(value: Any) -> None:
    if not isinstance(value, Mapping) or set(value) != {"schema_version", "patches"}:
        raise MINiMLModelError(
            "extensions.msc_harmonization requires schema_version and patches"
        )
    if value.get("schema_version") != PATCH_EXTENSION_SCHEMA_VERSION:
        raise MINiMLModelError("msc_harmonization schema_version must be '1.0'")
    patches = value.get("patches")
    if not isinstance(patches, (list, tuple)):
        raise MINiMLModelError("msc_harmonization patches must be a list")
    for fragment in patches:
        required = {
            "schema_version",
            "patch_schema_version",
            "patch_id",
            "package_index",
            "source_sha256",
            "result_sha256",
            "status",
            "operations",
        }
        if not isinstance(fragment, Mapping) or set(fragment) != required:
            raise MINiMLModelError("malformed retained harmonization patch fragment")
        if fragment["schema_version"] != PATCH_EXTENSION_SCHEMA_VERSION:
            raise MINiMLModelError("retained patch schema_version must be '1.0'")
        if fragment["patch_schema_version"] not in {"3.0", PATCH_SCHEMA_VERSION}:
            raise MINiMLModelError("unsupported retained patch schema version")
        for key in ("patch_id", "source_sha256", "result_sha256"):
            if not isinstance(fragment[key], str) or _SHA256.fullmatch(fragment[key]) is None:
                raise MINiMLModelError(f"retained patch {key} must be SHA-256")
        if isinstance(fragment["package_index"], bool) or not isinstance(
            fragment["package_index"], int
        ) or fragment["package_index"] < 0:
            raise MINiMLModelError("retained patch package_index must be non-negative")
        if fragment["status"] != "applied":
            raise MINiMLModelError("retained patch status must be applied")
        if not isinstance(fragment["operations"], (list, tuple)):
            raise MINiMLModelError("retained patch operations must be a list")
        for operation in fragment["operations"]:
            _validated_operation(
                operation, schema_version=str(fragment["patch_schema_version"])
            )


def _package_mapping(package: MINiMLPackage | Mapping[str, Any]) -> dict[str, Any]:
    if isinstance(package, MINiMLPackage):
        return package.to_mapping()
    if not isinstance(package, Mapping):
        raise TypeError("package must be a MINiML package")
    return MINiMLCodec.encode(MINiMLCodec().decode(package).package)


def iter_harmonization_patches(
    package: MINiMLPackage | Mapping[str, Any],
) -> tuple[dict[str, Any], ...]:
    mapping = _package_mapping(package)
    return _harmonization_patches_from_mapping(mapping)


def _harmonization_patches_from_mapping(
    mapping: Mapping[str, Any],
) -> tuple[dict[str, Any], ...]:
    extensions = mapping.get("extensions", {})
    retained = extensions.get(PATCH_EXTENSION_KEY) if isinstance(extensions, Mapping) else None
    if retained is None:
        return ()
    validate_harmonization_extension_mapping(retained)
    return tuple(_plain(item) for item in retained["patches"])


def iter_harmonization_operations(
    package: MINiMLPackage | Mapping[str, Any],
    *,
    sample: str | int | None = None,
    field: str | None = None,
) -> tuple[dict[str, Any], ...]:
    mapping = _package_mapping(package)
    return _harmonization_operations_from_mapping(
        mapping,
        sample=sample,
        field=field,
    )


def _harmonization_operations_from_mapping(
    mapping: Mapping[str, Any],
    *,
    sample: str | int | None = None,
    field: str | None = None,
) -> tuple[dict[str, Any], ...]:
    sample_index: int | None = None
    if isinstance(sample, int):
        sample_index = sample
    elif isinstance(sample, str):
        for index, item in enumerate(mapping.get("sample", [])):
            if isinstance(item, Mapping) and item.get("iid") == sample:
                sample_index = index
                break
    result: list[dict[str, Any]] = []
    for fragment in _harmonization_patches_from_mapping(mapping):
        for operation in fragment["operations"]:
            path = str(operation["path"])
            if sample_index is not None and not path.startswith(f"/sample/{sample_index}/"):
                continue
            if field is not None and operation["harmonized_value"].get("field") != field:
                continue
            item = _plain(operation)
            item["patch_id"] = fragment["patch_id"]
            item["package_index"] = fragment["package_index"]
            result.append(item)
    return tuple(result)


def harmonization_provenance_index(
    package: MINiMLPackage | Mapping[str, Any],
) -> dict[str, dict[str, tuple[dict[str, Any], ...]]]:
    mapping = _package_mapping(package)
    samples = mapping.get("sample", [])
    buckets: dict[str, dict[str, list[dict[str, Any]]]] = {}
    for operation in _harmonization_operations_from_mapping(mapping):
        match = re.match(r"^/sample/([0-9]+)(?:/|$)", str(operation["path"]))
        if match is None:
            continue
        index = int(match.group(1))
        if index >= len(samples) or not isinstance(samples[index], Mapping):
            continue
        sample_id = str(samples[index].get("iid") or index)
        value = operation["harmonized_value"]
        field = str(value["field"])
        entry = {
            "path": operation["path"],
            "patch_id": operation["patch_id"],
            "value": value["value"],
            "id": value.get("term_accession_number"),
            "ontology": value.get("term_source_ref"),
            "hierarchy_depth": value.get("hierarchy_depth"),
        }
        if isinstance(operation.get("source_evidence"), Mapping):
            entry.update(_plain(operation["source_evidence"]))
        buckets.setdefault(sample_id, {}).setdefault(field, []).append(entry)
    return {
        sample_id: {field: tuple(items) for field, items in fields.items()}
        for sample_id, fields in buckets.items()
    }


__all__ = [
    "MINiMLHarmonizationPatch",
    "PATCH_EXTENSION_KEY",
    "PATCH_EXTENSION_SCHEMA_VERSION",
    "PATCH_SCHEMA_VERSION",
    "apply_miniml_harmonization_patch",
    "canonical_miniml_document",
    "harmonization_provenance_index",
    "iter_harmonization_operations",
    "iter_harmonization_patches",
    "miniml_source_fingerprint",
    "validate_harmonization_extension_mapping",
]
