# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""XSD-derived, compatibility-first model for MINiML-compatible JSON."""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, field, fields, replace
from functools import wraps
import json
import os
from pathlib import Path
import re
from typing import Any, Callable, Mapping, TypeVar
from uuid import uuid4

from .harmonization import (
    HarmonizedValue,
    harmonized_mapping,
    is_harmonized_key,
    iter_harmonized_values,
    parse_harmonized_mapping,
)


MINIML_SCHEMA_VERSION = "3.0"


class MINiMLModelError(ValueError):
    """A package cannot be represented by the stable MINiML JSON model."""


class _FrozenJSONMapping(Mapping[str, Any]):
    """Recursively immutable storage for open MINiML extension fields."""

    def __init__(self, value: Mapping[str, Any] | None = None) -> None:
        self._data = {
            str(key): _freeze_json(item) for key, item in (value or {}).items()
        }

    def __getitem__(self, key: str) -> Any:
        return self._data[key]

    def __iter__(self):
        return iter(self._data)

    def __len__(self) -> int:
        return len(self._data)

    def __deepcopy__(self, memo):
        return self


def _freeze_json(value: Any) -> Any:
    if isinstance(value, Mapping):
        return _FrozenJSONMapping(value)
    if isinstance(value, (list, tuple)):
        return tuple(_freeze_json(item) for item in value)
    return deepcopy(value)


def _deep_freeze_constructor(cls):
    """Make direct dataclass construction as immutable as mapping decoding."""
    original = cls.__init__

    @wraps(original)
    def immutable_init(self, *args, **kwargs):
        original(self, *args, **kwargs)
        for model_field in fields(self):
            value = getattr(self, model_field.name)
            if isinstance(value, Mapping):
                frozen = _FrozenJSONMapping(value)
            elif isinstance(value, (list, tuple)):
                frozen = tuple(_freeze_json(item) for item in value)
            else:
                continue
            object.__setattr__(self, model_field.name, frozen)

    cls.__init__ = immutable_init
    return cls


@_deep_freeze_constructor
@dataclass(frozen=True)
class MINiMLValidationIssue:
    path: str
    code: str
    message: str
    severity: str = "warning"


def _reject_unknown(data: Mapping[str, Any], known: set[str], path: str) -> None:
    unknown = sorted(str(key) for key in data if key not in known)
    if unknown:
        label = "MINiML package" if path == "MINiML package" else path
        raise MINiMLModelError(f"unsupported {label} field: {unknown[0]}")


@_deep_freeze_constructor
@dataclass(frozen=True)
class _OccurrenceHarmonizedValue:
    field: str
    value: str
    term_source_ref: str | None = None
    term_accession_number: str | None = None
    hierarchy_depth: int | None = None
    index: int = 0

    @classmethod
    def from_mapping(cls, value: Any) -> "_OccurrenceHarmonizedValue":
        data = _mapping(value, "annotation")
        known = {
            "field", "value", "term_source_ref", "term_accession_number",
            "hierarchy_depth", "index",
        }
        _reject_unknown(data, known, "annotation")
        field_value = str(data.get("field", "")).strip()
        label = str(data.get("value", "")).strip()
        if not field_value or not label:
            raise MINiMLModelError("annotation requires nonblank field and value")
        depth = data.get("hierarchy_depth")
        if depth is not None and (isinstance(depth, bool) or not isinstance(depth, int) or depth < 0):
            raise MINiMLModelError("annotation hierarchy_depth must be a non-negative integer")
        return cls(
            field_value,
            label,
            data.get("term_source_ref"),
            data.get("term_accession_number"),
            depth,
            int(data.get("index", 0)),
        )

    def to_mapping(self) -> dict[str, Any]:
        result = {"field": self.field, "value": self.value}
        for key in ("term_source_ref", "term_accession_number", "hierarchy_depth"):
            _put(result, key, getattr(self, key))
        if self.index:
            result["index"] = self.index
        return result


@_deep_freeze_constructor
@dataclass(frozen=True)
class NamedComment:
    name: str
    value: str

    def __post_init__(self) -> None:
        if not isinstance(self.name, str) or not self.name.strip():
            raise MINiMLModelError("comment requires a nonblank name")
        if not isinstance(self.value, str):
            raise MINiMLModelError("comment value must be a string")

    @classmethod
    def from_mapping(cls, value: Any) -> "NamedComment":
        data = _mapping(value, "comment")
        _reject_unknown(data, {"name", "value"}, "comment")
        name = str(data.get("name", "")).strip()
        if not name:
            raise MINiMLModelError("comment requires a nonblank name")
        return cls(name, str(data.get("value", "")))

    def to_mapping(self) -> dict[str, Any]:
        return {"name": self.name, "value": self.value}


def _annotations(value: Any) -> tuple[_OccurrenceHarmonizedValue, ...]:
    return _objects(value, _OccurrenceHarmonizedValue.from_mapping, "annotations")


def _harmonized_annotations(
    data: Mapping[str, Any], known: set[str], path: str
) -> tuple[_OccurrenceHarmonizedValue, ...]:
    unknown = {str(key) for key in data if key not in known}
    unsupported = sorted(key for key in unknown if not is_harmonized_key(key))
    if unsupported:
        raise MINiMLModelError(f"unsupported {path} field: {unsupported[0]}")
    values = parse_harmonized_mapping(data)
    return tuple(
        _OccurrenceHarmonizedValue(
            field=value.field,
            value=str(value.value),
            term_source_ref=value.term_source_ref,
            term_accession_number=value.term_accession_number,
            hierarchy_depth=value.hierarchy_depth,
            index=value.index,
        )
        for value in values
    )


def _harmonized_wire(
    values: tuple[_OccurrenceHarmonizedValue, ...],
) -> dict[str, Any]:
    return harmonized_mapping(
        HarmonizedValue(
            field=value.field,
            value=value.value,
            term_source_ref=value.term_source_ref,
            term_accession_number=value.term_accession_number,
            hierarchy_depth=value.hierarchy_depth,
            index=value.index,
        )
        for value in values
    )


def _comments(value: Any) -> tuple[NamedComment, ...]:
    return _objects(value, NamedComment.from_mapping, "comments")


@_deep_freeze_constructor
@dataclass(frozen=True)
class OntologyValue:
    value: str
    term_source_ref: str | None = None
    term_accession_number: str | None = None
    annotations: tuple[_OccurrenceHarmonizedValue, ...] = ()

    def __post_init__(self) -> None:
        if not isinstance(self.value, str):
            raise MINiMLModelError("ontology value must be a string")
        if not isinstance(self.annotations, tuple) or not all(
            isinstance(item, _OccurrenceHarmonizedValue) for item in self.annotations
        ):
            raise MINiMLModelError("ontology value annotations must be typed annotations")

    @classmethod
    def from_value(cls, value: Any) -> "OntologyValue":
        if isinstance(value, Mapping):
            known = {"value", "term_source_ref", "term_accession_number"}
            return cls(
                str(value.get("value", "")),
                value.get("term_source_ref"),
                value.get("term_accession_number"),
                _harmonized_annotations(value, known, "ontology value"),
            )
        if isinstance(value, (str, int, float, bool)):
            return cls(str(value))
        raise MINiMLModelError("ontology value must be a scalar or object")

    def to_mapping(self) -> dict[str, Any]:
        result = {"value": self.value}
        for key in ("term_source_ref", "term_accession_number"):
            _put(result, key, getattr(self, key))
        result.update(_harmonized_wire(self.annotations))
        return result


@_deep_freeze_constructor
@dataclass(frozen=True)
class NamedValue:
    name: str
    value: Any
    term_source_ref: str | None = None
    term_accession_number: str | None = None
    unit: OntologyValue | None = None
    annotations: tuple[_OccurrenceHarmonizedValue, ...] = ()
    comments: tuple[NamedComment, ...] = ()
    qualifier: str | None = None
    unit_type: str | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.name, str) or not self.name.strip():
            raise MINiMLModelError("named value requires a nonblank name")
        if is_harmonized_key(self.name):
            from .harmonization import parse_harmonized_key

            _field, role, _index = parse_harmonized_key(self.name)
            if role == "hierarchy_depth":
                if (
                    isinstance(self.value, bool)
                    or not isinstance(self.value, int)
                    or self.value < 0
                ):
                    raise MINiMLModelError(
                        "harmonized hierarchy depth row must contain a non-negative integer"
                    )
            elif not isinstance(self.value, str):
                raise MINiMLModelError("harmonized named values must be strings")
        elif not isinstance(self.value, str):
            raise MINiMLModelError("named value value must be a string")
        if self.unit is not None and not isinstance(self.unit, OntologyValue):
            raise MINiMLModelError("named value unit must be an ontology value")
        if not isinstance(self.annotations, tuple) or not all(
            isinstance(item, _OccurrenceHarmonizedValue) for item in self.annotations
        ):
            raise MINiMLModelError("named value annotations must be typed annotations")
        if not isinstance(self.comments, tuple) or not all(
            isinstance(item, NamedComment) for item in self.comments
        ):
            raise MINiMLModelError("named value comments must be typed comments")

    @classmethod
    def from_mapping(cls, value: Any) -> "NamedValue":
        data = _mapping(value, "named value")
        known = {
            "name", "value", "term_source_ref", "term_accession_number", "unit",
            "comments", "qualifier", "unit_type",
        }
        name = str(data.get("name", "")).strip()
        if not name:
            raise MINiMLModelError("named value requires a nonblank name")
        raw_value = data.get("value", "")
        return cls(
            name,
            raw_value if is_harmonized_key(name) else str(raw_value),
            data.get("term_source_ref"),
            data.get("term_accession_number"),
            None if data.get("unit") is None else OntologyValue.from_value(data["unit"]),
            _harmonized_annotations(data, known, "named value"),
            _comments(data.get("comments")),
            data.get("qualifier"),
            data.get("unit_type"),
        )

    def to_mapping(self) -> dict[str, Any]:
        result = {"name": self.name, "value": self.value}
        for key in (
            "term_source_ref", "term_accession_number", "unit",
            "comments", "qualifier", "unit_type",
        ):
            _put(result, key, getattr(self, key))
        result.update(_harmonized_wire(self.annotations))
        return result


@_deep_freeze_constructor
@dataclass(frozen=True)
class SourceDocument:
    kind: str
    name: str
    uri: str | None = None
    sha256: str | None = None
    media_type: str | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.kind, str) or not self.kind.strip():
            raise MINiMLModelError("source document requires nonblank kind and name")
        if not isinstance(self.name, str) or not self.name.strip():
            raise MINiMLModelError("source document requires nonblank kind and name")
        if self.sha256 is not None:
            if not isinstance(self.sha256, str) or not re.fullmatch(r"[0-9a-fA-F]{64}", self.sha256):
                raise MINiMLModelError("source document sha256 must contain 64 hexadecimal characters")
            object.__setattr__(self, "sha256", self.sha256.lower())
        if self.media_type is not None and (
            not isinstance(self.media_type, str) or not self.media_type.strip()
        ):
            raise MINiMLModelError("source document media_type must be nonblank")

    @classmethod
    def from_mapping(cls, value: Any) -> "SourceDocument":
        data = _mapping(value, "source document")
        _reject_unknown(data, {"kind", "name", "uri", "sha256", "media_type"}, "source document")
        kind = str(data.get("kind", "")).strip()
        name = str(data.get("name", "")).strip()
        if not kind or not name:
            raise MINiMLModelError("source document requires nonblank kind and name")
        digest = data.get("sha256")
        if digest is not None and not re.fullmatch(r"[0-9a-fA-F]{64}", str(digest)):
            raise MINiMLModelError("source document sha256 must contain 64 hexadecimal characters")
        return cls(
            kind,
            name,
            data.get("uri"),
            None if digest is None else str(digest).lower(),
            data.get("media_type"),
        )

    def to_mapping(self) -> dict[str, Any]:
        result = {"kind": self.kind, "name": self.name}
        _put(result, "uri", self.uri)
        _put(result, "sha256", self.sha256)
        _put(result, "media_type", self.media_type)
        return result


@_deep_freeze_constructor
@dataclass(frozen=True)
class SourceInfo:
    format: str
    version: str | None = None
    schema_location: str | None = None
    documents: tuple[SourceDocument, ...] = ()

    def __post_init__(self) -> None:
        if not isinstance(self.format, str) or not self.format.strip():
            raise MINiMLModelError("source requires a nonblank format")
        if not isinstance(self.documents, tuple) or not all(
            isinstance(item, SourceDocument) for item in self.documents
        ):
            raise MINiMLModelError("source documents must be typed source documents")
        names = [item.name for item in self.documents]
        if len(names) != len(set(names)):
            raise MINiMLModelError("source document names must be unique")

    @classmethod
    def from_mapping(cls, value: Any) -> "SourceInfo":
        data = _mapping(value, "source")
        known = {"format", "version", "schema_location", "documents"}
        _reject_unknown(data, known, "source")
        format_value = str(data.get("format", "")).strip()
        if not format_value:
            raise MINiMLModelError("source requires a nonblank format")
        documents = _objects(data.get("documents"), SourceDocument.from_mapping, "source.documents")
        names = [item.name for item in documents]
        if len(names) != len(set(names)):
            raise MINiMLModelError("source document names must be unique")
        return cls(format_value, data.get("version"), data.get("schema_location"), documents)

    def to_mapping(self) -> dict[str, Any]:
        result = {"format": self.format}
        for key in ("version", "schema_location", "documents"):
            _put(result, key, getattr(self, key))
        return result


@_deep_freeze_constructor
@dataclass(frozen=True)
class PubMedPublication:
    pubmed_id: str
    doi: str | None = None
    author_list: str | None = None
    title: str | None = None
    status: str | None = None
    status_term_source_ref: str | None = None
    status_term_accession_number: str | None = None
    extras: Mapping[str, Any] = field(default_factory=dict)

    @classmethod
    def from_mapping(cls, value: Any) -> "PubMedPublication":
        data = _mapping(value, "pubmed_publication")
        known = {"pubmed_id", "doi", "author_list", "title", "status", "status_term_source_ref", "status_term_accession_number"}
        return cls(str(data.get("pubmed_id", "")), data.get("doi"), data.get("author_list"), data.get("title"), data.get("status"), data.get("status_term_source_ref"), data.get("status_term_accession_number"), _extras(data, known))

    def to_mapping(self) -> dict[str, Any]:
        result: dict[str, Any] = {
            "pubmed_id": self.pubmed_id,
            "doi": self.doi,
            "author_list": self.author_list,
            "title": self.title,
            "status": self.status,
            "status_term_source_ref": self.status_term_source_ref,
            "status_term_accession_number": self.status_term_accession_number,
        }
        return _record(result, self.extras)


@_deep_freeze_constructor
@dataclass(frozen=True)
class FASTQFile:
    uri: str | None = None
    filename: str | None = None
    md5: str | None = None
    bytes: str | None = None
    extras: Mapping[str, Any] = field(default_factory=dict)

    @classmethod
    def from_mapping(cls, value: Any) -> "FASTQFile":
        data = _mapping(value, "fastq_file")
        known = {"uri", "filename", "md5", "bytes"}
        return cls(data.get("uri"), data.get("filename"), data.get("md5"), data.get("bytes"), _extras(data, known))

    def to_mapping(self) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key in ("uri", "filename", "md5", "bytes"):
            _put(result, key, getattr(self, key))
        return _record(result, self.extras)


@_deep_freeze_constructor
@dataclass(frozen=True)
class SRARun:
    run: str | None = None
    study: str | None = None
    experiment: str | None = None
    sample: str | None = None
    biosample: str | None = None
    geo_sample: str | None = None
    library_layout: str | None = None
    library_selection: str | None = None
    library_source: str | None = None
    library_strategy: str | None = None
    scan_name: str | None = None
    instrument_model: str | None = None
    fastq_files: tuple[FASTQFile, ...] = ()
    submitted_file_name: str | None = None
    md5: str | None = None
    read_lengths: tuple[Any, ...] = ()
    extras: Mapping[str, Any] = field(default_factory=dict)

    @classmethod
    def from_mapping(cls, value: Any) -> "SRARun":
        data = _mapping(value, "sra_run")
        known = {"run", "study", "experiment", "sample", "biosample", "geo_sample", "library_layout", "library_selection", "library_source", "library_strategy", "scan_name", "instrument_model", "fastq_files", "submitted_file_name", "md5", "read_lengths"}
        return cls(data.get("run"), data.get("study"), data.get("experiment"), data.get("sample"), data.get("biosample"), data.get("geo_sample"), data.get("library_layout"), data.get("library_selection"), data.get("library_source"), data.get("library_strategy"), data.get("scan_name"), data.get("instrument_model"), _objects(data.get("fastq_files"), FASTQFile.from_mapping, "sra_run.fastq_files"), data.get("submitted_file_name"), data.get("md5"), tuple(_items(data.get("read_lengths"))), _extras(data, known))

    def to_mapping(self) -> dict[str, Any]:
        result: dict[str, Any] = {"run": self.run, "study": self.study}
        for key in ("experiment", "sample", "biosample", "geo_sample", "library_layout", "library_selection", "library_source", "library_strategy", "scan_name", "instrument_model", "fastq_files", "submitted_file_name", "md5", "read_lengths"):
            _put(result, key, getattr(self, key))
        return _record(result, self.extras)


def _mapping(value: Any, path: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise MINiMLModelError(f"{path} must be an object")
    return value


def _items(value: Any) -> list[Any]:
    if value is None:
        return []
    return list(value) if isinstance(value, (list, tuple)) else [value]


def _extras(data: Mapping[str, Any], known: set[str]) -> Mapping[str, Any]:
    return _FrozenJSONMapping({key: value for key, value in data.items() if key not in known})


def _put(result: dict[str, Any], key: str, value: Any) -> None:
    if value is None or value == ():
        return
    if isinstance(value, tuple):
        result[key] = [_plain(item) for item in value]
    else:
        result[key] = _plain(value)


def _plain(value: Any) -> Any:
    if hasattr(value, "to_mapping"):
        return value.to_mapping()
    if isinstance(value, Mapping):
        return {str(key): _plain(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_plain(item) for item in value]
    return deepcopy(value)


def _record(result: dict[str, Any], extras: Mapping[str, Any]) -> dict[str, Any]:
    return {**_plain(extras), **result}


T = TypeVar("T")


def _objects(value: Any, cls: Callable[[Any], T], path: str) -> tuple[T, ...]:
    result = []
    for index, item in enumerate(_items(value)):
        if not isinstance(item, Mapping):
            raise MINiMLModelError(f"{path}[{index}] must be an object")
        result.append(cls(item))
    return tuple(result)


@_deep_freeze_constructor
@dataclass(frozen=True)
class Accession:
    value: str
    database: str | None = None
    extras: Mapping[str, Any] = field(default_factory=dict)

    @classmethod
    def from_value(cls, value: Any) -> "Accession":
        if isinstance(value, Mapping):
            return cls(
                value=str(value.get("value", "")),
                database=None if value.get("database") is None else str(value["database"]),
                extras=_extras(value, {"value", "database"}),
            )
        if isinstance(value, (str, int)):
            return cls(str(value))
        raise MINiMLModelError("accession must be a string or object")

    def to_mapping(self) -> dict[str, Any]:
        result = {"value": self.value}
        _put(result, "database", self.database)
        return _record(result, self.extras)


@_deep_freeze_constructor
@dataclass(frozen=True)
class Reference:
    ref: str
    position: str | None = None
    extras: Mapping[str, Any] = field(default_factory=dict)

    @classmethod
    def from_mapping(cls, value: Any) -> "Reference":
        data = _mapping(value, "reference")
        ref = data.get("ref")
        if not isinstance(ref, (str, int)) or not str(ref).strip():
            raise MINiMLModelError("reference ref must be a nonblank string")
        return cls(
            ref=str(ref),
            position=None if data.get("position") is None else str(data["position"]),
            extras=_extras(data, {"ref", "position"}),
        )

    def to_mapping(self) -> dict[str, Any]:
        result = {"ref": self.ref}
        _put(result, "position", self.position)
        return _record(result, self.extras)


@_deep_freeze_constructor
@dataclass(frozen=True)
class Status:
    submission_date: str | None = None
    release_date: str | None = None
    last_update_date: str | None = None
    comments: tuple[Any, ...] = ()
    database: str | None = None
    extras: Mapping[str, Any] = field(default_factory=dict)

    @classmethod
    def from_mapping(cls, value: Any) -> "Status":
        data = _mapping(value, "status")
        return cls(
            submission_date=data.get("submission_date"),
            release_date=data.get("release_date"),
            last_update_date=data.get("last_update_date"),
            comments=tuple(_items(data.get("comment"))),
            database=data.get("database"),
            extras=_extras(data, {"submission_date", "release_date", "last_update_date", "comment", "database"}),
        )

    def to_mapping(self) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key in ("submission_date", "release_date", "last_update_date", "database"):
            _put(result, key, getattr(self, key))
        _put(result, "comment", self.comments)
        return _record(result, self.extras)


@_deep_freeze_constructor
@dataclass(frozen=True)
class SupplementLink:
    value: str
    type: str | None = None
    checksum: str | None = None
    build: str | None = None
    extras: Mapping[str, Any] = field(default_factory=dict)

    @classmethod
    def from_value(cls, value: Any) -> "SupplementLink":
        if isinstance(value, Mapping):
            return cls(
                value=str(value.get("value", "")), type=value.get("type"),
                checksum=value.get("checksum"), build=value.get("build"),
                extras=_extras(value, {"value", "type", "checksum", "build"}),
            )
        return cls(value=str(value))

    def to_mapping(self) -> dict[str, Any]:
        result = {"value": self.value}
        for key in ("type", "checksum", "build"):
            _put(result, key, getattr(self, key))
        return _record(result, self.extras)


@_deep_freeze_constructor
@dataclass(frozen=True)
class Organism:
    value: str
    taxid: str | None = None
    term_source_ref: str | None = None
    term_accession_number: str | None = None
    annotations: tuple[_OccurrenceHarmonizedValue, ...] = ()
    extras: Mapping[str, Any] = field(default_factory=dict)

    @classmethod
    def from_value(cls, value: Any) -> "Organism":
        if isinstance(value, Mapping):
            known = {
                "value", "taxid", "term_source_ref", "term_accession_number",
            }
            return cls(
                str(value.get("value", "")),
                None if value.get("taxid") is None else str(value["taxid"]),
                value.get("term_source_ref"),
                value.get("term_accession_number"),
                _harmonized_annotations(value, known, "organism"),
                _extras(value, {
                    *known, *[key for key in value if is_harmonized_key(key)],
                }),
            )
        return cls(str(value))

    def to_mapping(self) -> dict[str, Any]:
        result = {"value": self.value}
        for key in ("taxid", "term_source_ref", "term_accession_number"):
            _put(result, key, getattr(self, key))
        result.update(_harmonized_wire(self.annotations))
        return _record(result, self.extras)


@_deep_freeze_constructor
@dataclass(frozen=True)
class Relation:
    type: str
    target: str
    comment: str | None = None
    extras: Mapping[str, Any] = field(default_factory=dict)

    @classmethod
    def from_mapping(cls, value: Any) -> "Relation":
        data = _mapping(value, "relation")
        return cls(str(data.get("type", "")), str(data.get("target", "")), data.get("comment"), _extras(data, {"type", "target", "comment"}))

    def to_mapping(self) -> dict[str, Any]:
        result = {"type": self.type, "target": self.target}
        _put(result, "comment", self.comment)
        return _record(result, self.extras)


@_deep_freeze_constructor
@dataclass(frozen=True)
class Address:
    lines: tuple[Any, ...] = ()
    city: str | None = None
    state: str | None = None
    province: str | None = None
    zip_code: str | None = None
    postal_code: str | None = None
    country: str | None = None
    extras: Mapping[str, Any] = field(default_factory=dict)

    @classmethod
    def from_mapping(cls, value: Any) -> "Address":
        data = _mapping(value, "address")
        return cls(tuple(_items(data.get("line"))), data.get("city"), data.get("state"), data.get("province"), data.get("zip_code"), data.get("postal_code"), data.get("country"), _extras(data, {"line", "city", "state", "province", "zip_code", "postal_code", "country"}))

    def to_mapping(self) -> dict[str, Any]:
        result: dict[str, Any] = {}
        _put(result, "line", self.lines)
        for key in ("city", "state", "province", "zip_code", "postal_code", "country"):
            _put(result, key, getattr(self, key))
        return _record(result, self.extras)


@_deep_freeze_constructor
@dataclass(frozen=True)
class Person:
    first: str | None = None
    middle: str | None = None
    last: str | None = None
    extras: Mapping[str, Any] = field(default_factory=dict)

    @classmethod
    def from_mapping(cls, value: Any) -> "Person":
        data = _mapping(value, "person")
        return cls(data.get("first"), data.get("middle"), data.get("last"), _extras(data, {"first", "middle", "last"}))

    def to_mapping(self) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key in ("first", "middle", "last"):
            _put(result, key, getattr(self, key))
        return _record(result, self.extras)


@_deep_freeze_constructor
@dataclass(frozen=True)
class Characteristics:
    name: str
    value: Any
    term_source_ref: str | None = None
    term_accession_number: str | None = None
    unit: OntologyValue | None = None
    annotations: tuple[_OccurrenceHarmonizedValue, ...] = ()
    comments: tuple[NamedComment, ...] = ()
    qualifier: str | None = None
    unit_type: str | None = None

    @classmethod
    def from_value(cls, value: Any) -> "Characteristics":
        named = NamedValue.from_mapping(value)
        return cls(
            named.name,
            named.value,
            named.term_source_ref,
            named.term_accession_number,
            named.unit,
            named.annotations,
            named.comments,
            named.qualifier,
            named.unit_type,
        )

    def to_mapping(self) -> dict[str, Any]:
        return NamedValue(
            self.name,
            self.value,
            self.term_source_ref,
            self.term_accession_number,
            self.unit,
            self.annotations,
            self.comments,
            self.qualifier,
            self.unit_type,
        ).to_mapping()


@_deep_freeze_constructor
@dataclass(frozen=True)
class InstrumentModel:
    predefined: str | None = None
    other: str | None = None
    value: str | None = None
    extras: Mapping[str, Any] = field(default_factory=dict)

    @classmethod
    def from_value(cls, value: Any) -> "InstrumentModel":
        if isinstance(value, Mapping):
            return cls(value.get("predefined"), value.get("other"), value.get("value"), _extras(value, {"predefined", "other", "value"}))
        return cls(value=str(value))

    def to_mapping(self) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key in ("predefined", "other", "value"):
            _put(result, key, getattr(self, key))
        return _record(result, self.extras)


@_deep_freeze_constructor
@dataclass(frozen=True)
class DataColumn:
    name: str | None = None
    type: str | None = None
    unit: str | None = None
    description: str | None = None
    link_prefix: str | None = None
    link_suffix: str | None = None
    link_delimiter: str | None = None
    position: str | None = None
    extras: Mapping[str, Any] = field(default_factory=dict)

    @classmethod
    def from_mapping(cls, value: Any) -> "DataColumn":
        data = _mapping(value, "data column")
        known = {"name", "type", "unit", "description", "link_prefix", "link_suffix", "link_delimiter", "position"}
        return cls(*(data.get(key) for key in ("name", "type", "unit", "description", "link_prefix", "link_suffix", "link_delimiter", "position")), extras=_extras(data, known))

    def to_mapping(self) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key in ("name", "type", "unit", "description", "link_prefix", "link_suffix", "link_delimiter", "position"):
            _put(result, key, getattr(self, key))
        return _record(result, self.extras)


@_deep_freeze_constructor
@dataclass(frozen=True)
class TableData:
    value: str
    rows: str | None = None
    extras: Mapping[str, Any] = field(default_factory=dict)

    @classmethod
    def from_value(cls, value: Any) -> "TableData":
        if isinstance(value, Mapping):
            return cls(str(value.get("value", "")), None if value.get("rows") is None else str(value["rows"]), _extras(value, {"value", "rows"}))
        return cls(str(value))

    def to_mapping(self) -> dict[str, Any]:
        result = {"value": self.value}
        _put(result, "rows", self.rows)
        return _record(result, self.extras)


@_deep_freeze_constructor
@dataclass(frozen=True)
class DataTable:
    external_file: SupplementLink | None = None
    title: str | None = None
    columns: tuple[DataColumn, ...] = ()
    internal_data: TableData | None = None
    external_data: TableData | None = None
    extras: Mapping[str, Any] = field(default_factory=dict)

    @classmethod
    def from_mapping(cls, value: Any) -> "DataTable":
        data = _mapping(value, "data table")
        return cls(
            None if data.get("external_file") is None else SupplementLink.from_value(data["external_file"]),
            data.get("title"),
            _objects(data.get("column"), DataColumn.from_mapping, "data_table.column"),
            None if data.get("internal_data") is None else TableData.from_value(data["internal_data"]),
            None if data.get("external_data") is None else TableData.from_value(data["external_data"]),
            _extras(data, {"external_file", "title", "column", "internal_data", "external_data"}),
        )

    def to_mapping(self) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key in ("external_file", "title", "columns", "internal_data", "external_data"):
            _put(result, "column" if key == "columns" else key, getattr(self, key))
        return _record(result, self.extras)


@_deep_freeze_constructor
@dataclass(frozen=True)
class Channel:
    source: OntologyValue | None = None
    organisms: tuple[Organism, ...] = ()
    characteristics: tuple[Characteristics, ...] = ()
    biomaterial_providers: tuple[Any, ...] = ()
    treatment_protocol: str | None = None
    growth_protocol: str | None = None
    molecule: OntologyValue | None = None
    extract_protocol: str | None = None
    label: str | None = None
    label_protocol: str | None = None
    annotations: tuple[_OccurrenceHarmonizedValue, ...] = ()
    extras: Mapping[str, Any] = field(default_factory=dict)

    @classmethod
    def from_mapping(cls, value: Any) -> "Channel":
        data = _mapping(value, "channel")
        known = {"source", "organism", "characteristics", "biomaterial_provider", "treatment_protocol", "growth_protocol", "molecule", "extract_protocol", "label", "label_protocol", "extensions"}
        characteristic_values = _items(data.get("characteristics"))
        # Named hz rows form one occurrence-local evidence stream. Validate the
        # stream as a whole so no ID/ontology/depth companion can be orphaned.
        iter_harmonized_values(
            [item for item in characteristic_values if isinstance(item, Mapping)]
        )
        return cls(
            None if data.get("source") is None else OntologyValue.from_value(data["source"]), tuple(Organism.from_value(item) for item in _items(data.get("organism"))),
            tuple(Characteristics.from_value(item) for item in characteristic_values),
            tuple(_items(data.get("biomaterial_provider"))), data.get("treatment_protocol"),
            data.get("growth_protocol"), None if data.get("molecule") is None else OntologyValue.from_value(data["molecule"]), data.get("extract_protocol"),
            data.get("label"), data.get("label_protocol"), _harmonized_annotations(data, known, "channel"), _FrozenJSONMapping(_mapping(data.get("extensions", {}), "channel.extensions")),
        )

    def to_mapping(self) -> dict[str, Any]:
        result: dict[str, Any] = {}
        mapping = {"organisms": "organism", "biomaterial_providers": "biomaterial_provider"}
        for key in ("source", "organisms", "characteristics", "biomaterial_providers", "treatment_protocol", "growth_protocol", "molecule", "extract_protocol", "label", "label_protocol"):
            _put(result, mapping.get(key, key), getattr(self, key))
        result.update(_harmonized_wire(self.annotations))
        if self.extras:
            result["extensions"] = _plain(self.extras)
        return result


@_deep_freeze_constructor
@dataclass(frozen=True)
class Variable:
    factor: str | None = None
    type: OntologyValue | None = None
    description: str | None = None
    sample_ref: tuple[Reference, ...] = ()
    position: str | None = None
    extras: Mapping[str, Any] = field(default_factory=dict)

    @classmethod
    def from_mapping(cls, value: Any) -> "Variable":
        data = _mapping(value, "variable")
        return cls(data.get("factor"), None if data.get("type") is None else OntologyValue.from_value(data["type"]), data.get("description"), _objects(data.get("sample_ref"), Reference.from_mapping, "variable.sample_ref"), data.get("position"), _extras(data, {"factor", "type", "description", "sample_ref", "position"}))

    def to_mapping(self) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key in ("factor", "type", "description", "sample_ref", "position"):
            _put(result, key, getattr(self, key))
        return _record(result, self.extras)


@_deep_freeze_constructor
@dataclass(frozen=True)
class Repeat:
    factor: str | None = None
    sample_ref: tuple[Reference, ...] = ()
    position: str | None = None
    extras: Mapping[str, Any] = field(default_factory=dict)

    @classmethod
    def from_mapping(cls, value: Any) -> "Repeat":
        data = _mapping(value, "repeat")
        return cls(data.get("factor"), _objects(data.get("sample_ref"), Reference.from_mapping, "repeat.sample_ref"), data.get("position"), _extras(data, {"factor", "sample_ref", "position"}))

    def to_mapping(self) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key in ("factor", "sample_ref", "position"):
            _put(result, key, getattr(self, key))
        return _record(result, self.extras)


def _accessions(value: Any) -> tuple[Accession, ...]:
    return tuple(Accession.from_value(item) for item in _items(value))


def _statuses(value: Any) -> tuple[Status, ...]:
    return _objects(value, Status.from_mapping, "status")


def _refs(value: Any, path: str) -> tuple[Reference, ...]:
    return _objects(value, Reference.from_mapping, path)


def _relations(value: Any) -> tuple[Relation, ...]:
    return _objects(value, Relation.from_mapping, "relation")


def _links(value: Any) -> tuple[SupplementLink, ...]:
    return tuple(SupplementLink.from_value(item) for item in _items(value))


@_deep_freeze_constructor
@dataclass(frozen=True)
class Database:
    iid: str | None = None
    name: str | None = None
    public_id: str | None = None
    organization_ref: Reference | None = None
    organization: str | None = None
    web_link: str | None = None
    email: str | None = None
    extras: Mapping[str, Any] = field(default_factory=dict)

    @classmethod
    def from_mapping(cls, value: Any) -> "Database":
        data = _mapping(value, "database")
        return cls(data.get("iid"), data.get("name"), data.get("public_id"), None if data.get("organization_ref") is None else Reference.from_mapping(data["organization_ref"]), data.get("organization"), data.get("web_link"), data.get("email"), _extras(data, {"iid", "name", "public_id", "organization_ref", "organization", "web_link", "email"}))

    def to_mapping(self) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key in ("iid", "name", "public_id", "organization_ref", "organization", "web_link", "email"):
            _put(result, key, getattr(self, key))
        return _record(result, self.extras)


@_deep_freeze_constructor
@dataclass(frozen=True)
class Organization:
    iid: str | None = None
    name: str | None = None
    address: Address | str | None = None
    extras: Mapping[str, Any] = field(default_factory=dict)

    @classmethod
    def from_mapping(cls, value: Any) -> "Organization":
        data = _mapping(value, "organization")
        address = data.get("address")
        if isinstance(address, Mapping):
            address = Address.from_mapping(address)
        return cls(data.get("iid"), data.get("name"), address, _extras(data, {"iid", "name", "address"}))

    def to_mapping(self) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key in ("iid", "name", "address"):
            _put(result, key, getattr(self, key))
        return _record(result, self.extras)


@_deep_freeze_constructor
@dataclass(frozen=True)
class Contributor:
    iid: str | None = None
    person: Person | None = None
    organization: str | None = None
    company: str | None = None
    email: str | None = None
    phone: str | None = None
    fax: str | None = None
    laboratory: str | None = None
    department: str | None = None
    address: Address | str | None = None
    organization_ref: Reference | None = None
    web_link: str | None = None
    roles: tuple[OntologyValue, ...] = ()
    extras: Mapping[str, Any] = field(default_factory=dict)

    @classmethod
    def from_mapping(cls, value: Any) -> "Contributor":
        data = _mapping(value, "contributor")
        known = {"iid", "person", "organization", "company", "email", "phone", "fax", "laboratory", "department", "address", "organization_ref", "web_link", "roles", "extensions"}
        _reject_unknown(data, known, "contributor")
        address = data.get("address")
        if isinstance(address, Mapping):
            address = Address.from_mapping(address)
        return cls(data.get("iid"), None if data.get("person") is None else Person.from_mapping(data["person"]), data.get("organization"), data.get("company"), data.get("email"), data.get("phone"), data.get("fax"), data.get("laboratory"), data.get("department"), address, None if data.get("organization_ref") is None else Reference.from_mapping(data["organization_ref"]), data.get("web_link"), tuple(OntologyValue.from_value(item) for item in _items(data.get("roles"))), _FrozenJSONMapping(_mapping(data.get("extensions", {}), "contributor.extensions")))

    def to_mapping(self) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key in ("iid", "person", "organization", "company", "email", "phone", "fax", "laboratory", "department", "address", "organization_ref", "web_link", "roles"):
            _put(result, key, getattr(self, key))
        if self.extras:
            result["extensions"] = _plain(self.extras)
        return result


@_deep_freeze_constructor
@dataclass(frozen=True)
class Platform:
    iid: str | None = None
    accessions: tuple[Accession, ...] = ()
    statuses: tuple[Status, ...] = ()
    title: str | None = None
    technology: str | None = None
    distribution: str | None = None
    organisms: tuple[Organism, ...] = ()
    manufacturer: str | None = None
    manufacture_protocol: str | None = None
    catalog_number: str | None = None
    support: str | None = None
    coating: str | None = None
    description: str | None = None
    web_links: tuple[Any, ...] = ()
    pubmed_ids: tuple[Any, ...] = ()
    citations: tuple[Any, ...] = ()
    contributor_ref: tuple[Reference, ...] = ()
    contributors: tuple[Contributor, ...] = ()
    contact_ref: tuple[Reference, ...] = ()
    contacts: tuple[Contributor, ...] = ()
    supplementary_data: tuple[SupplementLink, ...] = ()
    relations: tuple[Relation, ...] = ()
    data_table: DataTable | None = None
    extras: Mapping[str, Any] = field(default_factory=dict)

    @classmethod
    def from_mapping(cls, value: Any) -> "Platform":
        data = _mapping(value, "platform")
        known = {"iid", "accession", "status", "title", "technology", "distribution", "organism", "manufacturer", "manufacture_protocol", "catalog_number", "support", "coating", "description", "web_link", "pubmed_id", "citation", "contributor_ref", "contributor", "contact_ref", "contact", "supplementary_data", "relation", "data_table"}
        return cls(data.get("iid"), _accessions(data.get("accession")), _statuses(data.get("status")), data.get("title"), data.get("technology"), data.get("distribution"), tuple(Organism.from_value(x) for x in _items(data.get("organism"))), data.get("manufacturer"), data.get("manufacture_protocol"), data.get("catalog_number"), data.get("support"), data.get("coating"), data.get("description"), tuple(_items(data.get("web_link"))), tuple(_items(data.get("pubmed_id"))), tuple(_items(data.get("citation"))), _refs(data.get("contributor_ref"), "platform.contributor_ref"), _objects(data.get("contributor"), Contributor.from_mapping, "platform.contributor"), _refs(data.get("contact_ref"), "platform.contact_ref"), _objects(data.get("contact"), Contributor.from_mapping, "platform.contact"), _links(data.get("supplementary_data")), _relations(data.get("relation")), None if data.get("data_table") is None else DataTable.from_mapping(data["data_table"]), _extras(data, known))

    def to_mapping(self) -> dict[str, Any]:
        result: dict[str, Any] = {}
        keys = {"accessions": "accession", "statuses": "status", "organisms": "organism", "web_links": "web_link", "pubmed_ids": "pubmed_id", "citations": "citation", "contributors": "contributor", "contacts": "contact", "relations": "relation"}
        for key in ("iid", "accessions", "statuses", "title", "technology", "distribution", "organisms", "manufacturer", "manufacture_protocol", "catalog_number", "support", "coating", "description", "web_links", "pubmed_ids", "citations", "contributor_ref", "contributors", "contact_ref", "contacts", "supplementary_data", "relations", "data_table"):
            _put(result, keys.get(key, key), getattr(self, key))
        return _record(result, self.extras)


@_deep_freeze_constructor
@dataclass(frozen=True)
class Sample:
    iid: str | None = None
    accessions: tuple[Accession, ...] = ()
    statuses: tuple[Status, ...] = ()
    title: str | None = None
    type: str | None = None
    anchor: str | None = None
    tag_length: str | None = None
    tag_count: str | None = None
    channel_count: str | None = None
    channels: tuple[Channel, ...] = ()
    hybridization_protocol: str | None = None
    scan_protocol: str | None = None
    description: str | None = None
    data_processing: str | None = None
    platform_ref: Reference | None = None
    library_strategy: str | None = None
    library_source: str | None = None
    library_selection: str | None = None
    instrument_model: InstrumentModel | None = None
    barcode: str | None = None
    contact_ref: tuple[Reference, ...] = ()
    contacts: tuple[Contributor, ...] = ()
    supplementary_data: tuple[SupplementLink, ...] = ()
    raw_data: tuple[SupplementLink, ...] = ()
    relations: tuple[Relation, ...] = ()
    data_table: DataTable | None = None
    sra_accessions: tuple[str, ...] = ()
    ena_accessions: tuple[str, ...] = ()
    sra_runs: tuple[SRARun, ...] = ()
    sra_runs_present: bool = False
    extras: Mapping[str, Any] = field(default_factory=dict)

    @classmethod
    def from_mapping(cls, value: Any) -> "Sample":
        data = _mapping(value, "sample")
        known = {"iid", "accession", "status", "title", "type", "anchor", "tag_length", "tag_count", "channel_count", "channel", "hybridization_protocol", "scan_protocol", "description", "data_processing", "platform_ref", "library_strategy", "library_source", "library_selection", "instrument_model", "barcode", "contact_ref", "contact", "supplementary_data", "raw_data", "relation", "data_table", "sra_accession", "ena_accession", "sra_run"}
        return cls(data.get("iid"), _accessions(data.get("accession")), _statuses(data.get("status")), data.get("title"), data.get("type"), data.get("anchor"), data.get("tag_length"), data.get("tag_count"), data.get("channel_count"), _objects(data.get("channel"), Channel.from_mapping, "sample.channel"), data.get("hybridization_protocol"), data.get("scan_protocol"), data.get("description"), data.get("data_processing"), None if data.get("platform_ref") is None else Reference.from_mapping(data["platform_ref"]), data.get("library_strategy"), data.get("library_source"), data.get("library_selection"), None if data.get("instrument_model") is None else InstrumentModel.from_value(data["instrument_model"]), data.get("barcode"), _refs(data.get("contact_ref"), "sample.contact_ref"), _objects(data.get("contact"), Contributor.from_mapping, "sample.contact"), _links(data.get("supplementary_data")), _links(data.get("raw_data")), _relations(data.get("relation")), None if data.get("data_table") is None else DataTable.from_mapping(data["data_table"]), tuple(str(item) for item in _items(data.get("sra_accession"))), tuple(str(item) for item in _items(data.get("ena_accession"))), _objects(data.get("sra_run"), SRARun.from_mapping, "sample.sra_run"), "sra_run" in data, _extras(data, known))

    def to_mapping(self) -> dict[str, Any]:
        result: dict[str, Any] = {}
        keys = {"accessions": "accession", "statuses": "status", "channels": "channel", "contacts": "contact", "relations": "relation"}
        keys.update({"sra_accessions": "sra_accession", "ena_accessions": "ena_accession", "sra_runs": "sra_run"})
        for key in ("iid", "accessions", "statuses", "title", "type", "anchor", "tag_length", "tag_count", "channel_count", "channels", "hybridization_protocol", "scan_protocol", "description", "data_processing", "platform_ref", "library_strategy", "library_source", "library_selection", "instrument_model", "barcode", "contact_ref", "contacts", "supplementary_data", "raw_data", "relations", "data_table", "sra_accessions", "ena_accessions", "sra_runs"):
            _put(result, keys.get(key, key), getattr(self, key))
        if self.sra_runs_present and "sra_run" not in result:
            result["sra_run"] = []
        return _record(result, self.extras)


ASSAY_NODE_KINDS = {
    "source", "sample", "extract", "labeled_extract", "hybridization", "assay",
    "scan", "normalization", "array_data_file", "derived_array_data_file",
    "array_data_matrix_file", "derived_array_data_matrix_file", "image_file",
}


@_deep_freeze_constructor
@dataclass(frozen=True)
class Protocol:
    name: str
    type: OntologyValue | None = None
    description: str | None = None
    parameters: tuple[str, ...] = ()
    hardware: tuple[str, ...] = ()
    software: tuple[str, ...] = ()
    contacts: tuple[str, ...] = ()
    performers: tuple[str, ...] = ()
    comments: tuple[NamedComment, ...] = ()

    def __post_init__(self) -> None:
        if not isinstance(self.name, str) or not self.name.strip():
            raise MINiMLModelError("protocol requires a nonblank name")

    @classmethod
    def from_mapping(cls, value: Any) -> "Protocol":
        data = _mapping(value, "protocol")
        known = {
            "name", "type", "description", "parameters", "hardware", "software",
            "contacts", "performers", "comments",
        }
        _reject_unknown(data, known, "protocol")
        name = str(data.get("name", "")).strip()
        if not name:
            raise MINiMLModelError("protocol requires a nonblank name")
        return cls(
            name,
            None if data.get("type") is None else OntologyValue.from_value(data["type"]),
            data.get("description"),
            tuple(str(item) for item in _items(data.get("parameters"))),
            tuple(str(item) for item in _items(data.get("hardware"))),
            tuple(str(item) for item in _items(data.get("software"))),
            tuple(str(item) for item in _items(data.get("contacts"))),
            tuple(str(item) for item in _items(data.get("performers"))),
            _comments(data.get("comments")),
        )

    def to_mapping(self) -> dict[str, Any]:
        result = {"name": self.name}
        for key in (
            "type", "description", "parameters", "hardware", "software", "contacts",
            "performers", "comments",
        ):
            _put(result, key, getattr(self, key))
        return result


@_deep_freeze_constructor
@dataclass(frozen=True)
class ProtocolApplication:
    protocol_ref: str
    parameter_values: tuple[NamedValue, ...] = ()
    performer: str | None = None
    date: str | None = None
    comments: tuple[NamedComment, ...] = ()
    kind: str = field(default="protocol_application", init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.protocol_ref, str) or not self.protocol_ref.strip():
            raise MINiMLModelError("protocol application requires protocol_ref")

    @classmethod
    def from_mapping(cls, value: Any) -> "ProtocolApplication":
        data = _mapping(value, "protocol application")
        known = {"kind", "protocol_ref", "parameter_values", "performer", "date", "comments"}
        _reject_unknown(data, known, "protocol application")
        reference = str(data.get("protocol_ref", "")).strip()
        if not reference:
            raise MINiMLModelError("protocol application requires protocol_ref")
        return cls(
            reference,
            _objects(data.get("parameter_values"), NamedValue.from_mapping, "parameter_values"),
            data.get("performer"),
            data.get("date"),
            _comments(data.get("comments")),
        )

    def to_mapping(self) -> dict[str, Any]:
        result = {"kind": self.kind, "protocol_ref": self.protocol_ref}
        for key in ("parameter_values", "performer", "date", "comments"):
            _put(result, key, getattr(self, key))
        return result


@_deep_freeze_constructor
@dataclass(frozen=True)
class AssayNode:
    kind: str
    name: str
    sample_ref: str | None = None
    characteristics: tuple[NamedValue, ...] = ()
    factor_values: tuple[NamedValue, ...] = ()
    provider: str | None = None
    material_type: OntologyValue | None = None
    description: str | None = None
    label: OntologyValue | None = None
    technology_type: OntologyValue | None = None
    array_design_ref: Reference | None = None
    link: SupplementLink | None = None
    comments: tuple[NamedComment, ...] = ()

    def __post_init__(self) -> None:
        if self.kind not in ASSAY_NODE_KINDS:
            raise MINiMLModelError(f"unsupported assay node kind: {self.kind}")
        if not isinstance(self.name, str) or not self.name.strip():
            raise MINiMLModelError("assay node requires a nonblank name")

    @classmethod
    def from_mapping(cls, value: Any) -> "AssayNode":
        data = _mapping(value, "assay node")
        known = {
            "kind", "name", "sample_ref", "characteristics", "factor_values",
            "provider", "material_type", "description", "label", "technology_type",
            "array_design_ref", "link", "comments",
        }
        _reject_unknown(data, known, "assay node")
        kind = str(data.get("kind", "")).strip()
        name = str(data.get("name", "")).strip()
        if kind not in ASSAY_NODE_KINDS:
            raise MINiMLModelError(f"unsupported assay node kind: {kind}")
        if not name:
            raise MINiMLModelError("assay node requires a nonblank name")
        return cls(
            kind,
            name,
            data.get("sample_ref"),
            _objects(data.get("characteristics"), NamedValue.from_mapping, "assay node characteristics"),
            _objects(data.get("factor_values"), NamedValue.from_mapping, "assay node factor_values"),
            data.get("provider"),
            None if data.get("material_type") is None else OntologyValue.from_value(data["material_type"]),
            data.get("description"),
            None if data.get("label") is None else OntologyValue.from_value(data["label"]),
            None if data.get("technology_type") is None else OntologyValue.from_value(data["technology_type"]),
            None if data.get("array_design_ref") is None else Reference.from_mapping(data["array_design_ref"]),
            None if data.get("link") is None else SupplementLink.from_value(data["link"]),
            _comments(data.get("comments")),
        )

    def to_mapping(self) -> dict[str, Any]:
        result = {"kind": self.kind, "name": self.name}
        for key in (
            "sample_ref", "characteristics", "factor_values", "provider", "material_type",
            "description", "label", "technology_type", "array_design_ref", "link", "comments",
        ):
            _put(result, key, getattr(self, key))
        return result


AssayStep = AssayNode | ProtocolApplication


def _assay_step(value: Any) -> AssayStep:
    data = _mapping(value, "assay step")
    return (
        ProtocolApplication.from_mapping(data)
        if data.get("kind") == "protocol_application"
        else AssayNode.from_mapping(data)
    )


@_deep_freeze_constructor
@dataclass(frozen=True)
class AssayPath:
    steps: tuple[AssayStep, ...]
    document: str | None = None
    comments: tuple[NamedComment, ...] = ()

    def __post_init__(self) -> None:
        if not self.steps:
            raise MINiMLModelError("assay path requires at least one step")

    @classmethod
    def from_mapping(cls, value: Any) -> "AssayPath":
        data = _mapping(value, "assay path")
        _reject_unknown(data, {"document", "steps", "comments"}, "assay path")
        steps = tuple(_assay_step(item) for item in _items(data.get("steps")))
        if not steps:
            raise MINiMLModelError("assay path requires at least one step")
        return cls(steps, data.get("document"), _comments(data.get("comments")))

    def to_mapping(self) -> dict[str, Any]:
        result: dict[str, Any] = {"steps": [_plain(item) for item in self.steps]}
        _put(result, "document", self.document)
        _put(result, "comments", self.comments)
        return result


@_deep_freeze_constructor
@dataclass(frozen=True)
class Series:
    iid: str | None = None
    accessions: tuple[Accession, ...] = ()
    statuses: tuple[Status, ...] = ()
    title: str | None = None
    pubmed_ids: tuple[Any, ...] = ()
    citations: tuple[Any, ...] = ()
    web_links: tuple[Any, ...] = ()
    summary: str | None = None
    overall_design: str | None = None
    types: tuple[Any, ...] = ()
    contributor_ref: tuple[Reference, ...] = ()
    contributors: tuple[Contributor, ...] = ()
    contact_ref: tuple[Reference, ...] = ()
    contacts: tuple[Contributor, ...] = ()
    sample_ref: tuple[Reference, ...] = ()
    variables: tuple[Variable, ...] = ()
    repeats: tuple[Repeat, ...] = ()
    supplementary_data: tuple[SupplementLink, ...] = ()
    relations: tuple[Relation, ...] = ()
    data_tables: tuple[DataTable, ...] = ()
    pubmed_publications: tuple[PubMedPublication, ...] = ()
    experiment_date: str | None = None
    protocols: tuple[Protocol, ...] = ()
    assay_paths: tuple[AssayPath, ...] = ()
    quality_controls: tuple[OntologyValue, ...] = ()
    replicate_types: tuple[OntologyValue, ...] = ()
    normalization_types: tuple[OntologyValue, ...] = ()
    comments: tuple[NamedComment, ...] = ()
    extras: Mapping[str, Any] = field(default_factory=dict)

    @classmethod
    def from_mapping(cls, value: Any) -> "Series":
        data = _mapping(value, "series")
        known = {
            "iid", "accession", "status", "title", "pubmed_id", "citation", "web_link",
            "summary", "overall_design", "type", "contributor_ref", "contributor",
            "contact_ref", "contact", "sample_ref", "variable", "repeats",
            "supplementary_data", "relation", "data_table", "pubmed_publication",
            "experiment_date", "protocols", "assay_paths", "quality_controls",
            "replicate_types", "normalization_types", "comments", "extensions",
        }
        _reject_unknown(data, known, "series")
        return cls(
            iid=data.get("iid"),
            accessions=_accessions(data.get("accession")),
            statuses=_statuses(data.get("status")),
            title=data.get("title"),
            pubmed_ids=tuple(_items(data.get("pubmed_id"))),
            citations=tuple(_items(data.get("citation"))),
            web_links=tuple(_items(data.get("web_link"))),
            summary=data.get("summary"),
            overall_design=data.get("overall_design"),
            types=tuple(OntologyValue.from_value(item) for item in _items(data.get("type"))),
            contributor_ref=_refs(data.get("contributor_ref"), "series.contributor_ref"),
            contributors=_objects(data.get("contributor"), Contributor.from_mapping, "series.contributor"),
            contact_ref=_refs(data.get("contact_ref"), "series.contact_ref"),
            contacts=_objects(data.get("contact"), Contributor.from_mapping, "series.contact"),
            sample_ref=_refs(data.get("sample_ref"), "series.sample_ref"),
            variables=_objects(data.get("variable"), Variable.from_mapping, "series.variable"),
            repeats=_objects(data.get("repeats"), Repeat.from_mapping, "series.repeats"),
            supplementary_data=_links(data.get("supplementary_data")),
            relations=_relations(data.get("relation")),
            data_tables=_objects(data.get("data_table"), DataTable.from_mapping, "series.data_table"),
            pubmed_publications=_objects(data.get("pubmed_publication"), PubMedPublication.from_mapping, "series.pubmed_publication"),
            experiment_date=data.get("experiment_date"),
            protocols=_objects(data.get("protocols"), Protocol.from_mapping, "series.protocols"),
            assay_paths=_objects(data.get("assay_paths"), AssayPath.from_mapping, "series.assay_paths"),
            quality_controls=tuple(OntologyValue.from_value(item) for item in _items(data.get("quality_controls"))),
            replicate_types=tuple(OntologyValue.from_value(item) for item in _items(data.get("replicate_types"))),
            normalization_types=tuple(OntologyValue.from_value(item) for item in _items(data.get("normalization_types"))),
            comments=_comments(data.get("comments")),
            extras=_FrozenJSONMapping(_mapping(data.get("extensions", {}), "series.extensions")),
        )

    def to_mapping(self) -> dict[str, Any]:
        result: dict[str, Any] = {}
        keys = {"accessions": "accession", "statuses": "status", "pubmed_ids": "pubmed_id", "citations": "citation", "web_links": "web_link", "types": "type", "contributors": "contributor", "contacts": "contact", "variables": "variable", "relations": "relation", "data_tables": "data_table"}
        keys["pubmed_publications"] = "pubmed_publication"
        for key in ("iid", "accessions", "statuses", "title", "pubmed_ids", "citations", "web_links", "summary", "overall_design", "types", "contributor_ref", "contributors", "contact_ref", "contacts", "sample_ref", "variables", "repeats", "supplementary_data", "relations", "data_tables", "pubmed_publications", "experiment_date", "protocols", "assay_paths", "quality_controls", "replicate_types", "normalization_types", "comments"):
            _put(result, keys.get(key, key), getattr(self, key))
        if self.extras:
            result["extensions"] = _plain(self.extras)
        return result


TECHNOLOGIES = {"high-throughput sequencing", "in situ oligonucleotide", "spotted oligonucleotide", "mixed spotted oligonucleotide", "spotted DNA/cDNA", "spotted peptide or protein", "antibody", "tissue", "oligonucleotide beads", "MS", "SAGE NlaIII", "SAGE Sau3A", "SAGE RsaI", "SARST", "MPSS", "RT-PCR", "other"}
SAMPLE_TYPES = {"RNA", "genomic", "protein", "SAGE", "MPSS", "SARST", "mixed", "other", "SRA"}
MOLECULES = {"genomic DNA", "polyA RNA", "total RNA", "cytoplasmic RNA", "nuclear RNA", "protein", "other"}
VARIABLE_FACTORS = {"dose", "time", "tissue", "strain", "gender", "cell line", "development stage", "age", "agent", "cell type", "infection", "isolate", "metabolism", "shock", "stress", "temperature", "speciman", "disease state", "protocol", "growth protocol", "other", "genotype/variation", "species", "individual"}


@_deep_freeze_constructor
@dataclass(frozen=True)
class MINiMLPackage(Mapping[str, Any]):
    series: Series
    source: SourceInfo
    databases: tuple[Database, ...] = ()
    organizations: tuple[Organization, ...] = ()
    contributors: tuple[Contributor, ...] = ()
    platforms: tuple[Platform, ...] = ()
    samples: tuple[Sample, ...] = ()
    extensions: Mapping[str, Any] = field(default_factory=dict)
    miniml_schema_version: str = MINIML_SCHEMA_VERSION

    def __post_init__(self) -> None:
        extensions = _plain(self.extensions)
        for key, item in self.series.extras.items():
            if key == "msc_harmonization":
                raise MINiMLModelError(
                    "msc_harmonization is reserved for package extensions"
                )
            plain_item = _plain(item)
            if key in extensions and extensions[key] != plain_item:
                raise MINiMLModelError(
                    f"conflicting package and series extension {key!r}"
                )
            extensions[key] = plain_item
        object.__setattr__(self, "extensions", _FrozenJSONMapping(extensions))
        if self.series.extras:
            object.__setattr__(self, "series", replace(self.series, extras={}))

    def __getitem__(self, key: str) -> Any:
        return self.to_mapping()[key]

    def __iter__(self):
        return iter(self.to_mapping())

    def __len__(self) -> int:
        return len(self.to_mapping())

    @classmethod
    def from_mapping(cls, value: Mapping[str, Any]) -> "MINiMLPackage":
        data = _mapping(value, "MINiML package")
        version = data.get("miniml_schema_version")
        if version == "2.0":
            from .migration import MINiMLV2Migrator

            return MINiMLV2Migrator().migrate(data).package
        if version != MINIML_SCHEMA_VERSION:
            raise MINiMLModelError(
                "runtime decoding requires MSC MINiML schema version '3.0'"
            )
        known = {
            "miniml_schema_version", "source", "database", "organization",
            "contributor", "platform", "sample", "series", "extensions",
        }
        _reject_unknown(data, known, "MINiML package")
        series = Series.from_mapping(data.get("series"))
        if not (isinstance(series.iid, str) and series.iid.strip()) and not any(item.value.strip() for item in series.accessions):
            raise MINiMLModelError("series requires iid or accession")
        root_extensions = _plain(
            _mapping(data.get("extensions", {}), "extensions")
        )
        for key, item in series.extras.items():
            if key == "msc_harmonization":
                raise MINiMLModelError(
                    "msc_harmonization is reserved for package extensions"
                )
            plain_item = _plain(item)
            if key in root_extensions and root_extensions[key] != plain_item:
                raise MINiMLModelError(
                    f"conflicting package and series extension {key!r}"
                )
            root_extensions[key] = plain_item
        retained = root_extensions.get("msc_harmonization")
        if retained is not None:
            from .patches import validate_harmonization_extension_mapping

            validate_harmonization_extension_mapping(retained)
        package = cls(
            series=series,
            source=SourceInfo.from_mapping(data.get("source")),
            databases=_objects(data.get("database"), Database.from_mapping, "database"),
            organizations=_objects(data.get("organization"), Organization.from_mapping, "organization"),
            contributors=_objects(data.get("contributor"), Contributor.from_mapping, "contributor"),
            platforms=_objects(data.get("platform"), Platform.from_mapping, "platform"),
            samples=_objects(data.get("sample"), Sample.from_mapping, "sample"),
            extensions=_FrozenJSONMapping(root_extensions),
        )
        package._raise_structural_errors()
        return package

    @classmethod
    def load(cls, path: str | Path) -> "MINiMLPackage":
        with Path(path).open(encoding="utf-8") as handle:
            return cls.from_mapping(json.load(handle))

    def to_mapping(self) -> dict[str, Any]:
        extensions = _plain(self.extensions)
        for key, item in self.series.extras.items():
            if key == "msc_harmonization":
                raise MINiMLModelError(
                    "msc_harmonization is reserved for package extensions"
                )
            plain_item = _plain(item)
            if key in extensions and extensions[key] != plain_item:
                raise MINiMLModelError(
                    f"conflicting package and series extension {key!r}"
                )
            extensions[key] = plain_item
        series = self.series.to_mapping()
        series.pop("extensions", None)
        result: dict[str, Any] = {
            "miniml_schema_version": self.miniml_schema_version,
            "source": self.source.to_mapping(),
            "database": [_plain(item) for item in self.databases],
            "organization": [_plain(item) for item in self.organizations],
            "contributor": [_plain(item) for item in self.contributors],
            "platform": [_plain(item) for item in self.platforms],
            "sample": [_plain(item) for item in self.samples],
            "series": series,
        }
        if extensions:
            result["extensions"] = extensions
        return result

    def dump(self, path: str | Path) -> None:
        destination = Path(path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        temporary = destination.with_name(f".{destination.name}.{uuid4().hex}.tmp")
        try:
            with temporary.open("x", encoding="utf-8") as handle:
                json.dump(self.to_mapping(), handle, indent=2, sort_keys=True, ensure_ascii=False)
                handle.write("\n")
                handle.flush()
                os.fsync(handle.fileno())
            os.replace(temporary, destination)
        finally:
            temporary.unlink(missing_ok=True)

    def _raise_structural_errors(self) -> None:
        for name, entities in (("database", self.databases), ("organization", self.organizations), ("contributor", self.contributors), ("platform", self.platforms), ("sample", self.samples)):
            seen: set[str] = set()
            for entity in entities:
                iid = getattr(entity, "iid", None)
                if not isinstance(iid, str) or not iid.strip():
                    continue
                if iid in seen:
                    raise MINiMLModelError(f"duplicate {name} iid: {iid}")
                seen.add(iid)

    def validate(self) -> tuple[MINiMLValidationIssue, ...]:
        issues: list[MINiMLValidationIssue] = []
        def warn(path: str, code: str, message: str) -> None:
            issues.append(MINiMLValidationIssue(path, code, message))

        sample_ids = {item.iid for item in self.samples if item.iid}
        platform_ids = {item.iid for item in self.platforms if item.iid}
        contributor_ids = {item.iid for item in self.contributors if item.iid}
        database_ids = {item.iid for item in self.databases if item.iid}
        organization_ids = {item.iid for item in self.organizations if item.iid}

        def check_refs(refs: tuple[Reference, ...], known: set[str], path: str) -> None:
            for index, ref in enumerate(refs):
                if ref.ref not in known:
                    warn(f"{path}/{index}/ref", "unresolved_reference", f"reference {ref.ref!r} is not present in this package")

        check_refs(self.series.sample_ref, sample_ids, "/series/sample_ref")
        check_refs(self.series.contributor_ref + self.series.contact_ref, contributor_ids, "/series/contributors")
        for index, variable in enumerate((*self.series.variables, *self.series.repeats)):
            check_refs(variable.sample_ref, {ref.ref for ref in self.series.sample_ref}, f"/series/subsets/{index}/sample_ref")
        for index, sample in enumerate(self.samples):
            if sample.platform_ref and sample.platform_ref.ref not in platform_ids:
                warn(f"/sample/{index}/platform_ref/ref", "unresolved_reference", f"platform {sample.platform_ref.ref!r} is not present in this package")
            check_refs(sample.contact_ref, contributor_ids, f"/sample/{index}/contact_ref")
            if sample.channel_count is not None and str(sample.channel_count).isdigit() and int(sample.channel_count) != len(sample.channels):
                warn(f"/sample/{index}/channel_count", "channel_count_mismatch", "channel_count does not match channel array length")
            if sample.type and sample.type not in SAMPLE_TYPES:
                warn(f"/sample/{index}/type", "xsd_enumeration", f"{sample.type!r} is outside the MINiML 0.5.4 vocabulary")
            for channel_index, channel in enumerate(sample.channels):
                molecule = channel.molecule.value if channel.molecule else None
                if molecule and molecule not in MOLECULES:
                    warn(f"/sample/{index}/channel/{channel_index}/molecule", "xsd_enumeration", f"{molecule!r} is outside the MINiML 0.5.4 vocabulary")
            self._validate_links(sample.supplementary_data + sample.raw_data, f"/sample/{index}", warn)
        for index, platform in enumerate(self.platforms):
            if platform.technology and platform.technology not in TECHNOLOGIES:
                warn(f"/platform/{index}/technology", "xsd_enumeration", f"{platform.technology!r} is outside the MINiML 0.5.4 vocabulary")
            check_refs(platform.contributor_ref + platform.contact_ref, contributor_ids, f"/platform/{index}/contributors")
            self._validate_links(platform.supplementary_data, f"/platform/{index}", warn)
        for index, database in enumerate(self.databases):
            if database.organization_ref and database.organization_ref.ref not in organization_ids:
                warn(f"/database/{index}/organization_ref/ref", "unresolved_reference", f"organization {database.organization_ref.ref!r} is not present in this package")
        for entity_name, entities in (("series", (self.series,)), ("sample", self.samples), ("platform", self.platforms)):
            titles: set[str] = set()
            for index, entity in enumerate(entities):
                title = getattr(entity, "title", None)
                if title and title in titles:
                    warn(f"/{entity_name}/{index}/title", "xsd_uniqueness", f"duplicate {entity_name} title {title!r}")
                if title:
                    titles.add(title)
                for accession_index, accession in enumerate(getattr(entity, "accessions", ())):
                    if accession.database and accession.database not in database_ids:
                        warn(f"/{entity_name}/{index}/accession/{accession_index}/database", "unresolved_reference", f"database {accession.database!r} is not present in this package")
        for index, variable in enumerate(self.series.variables):
            if variable.factor and variable.factor not in VARIABLE_FACTORS:
                warn(f"/series/variable/{index}/factor", "xsd_enumeration", f"{variable.factor!r} is outside the MINiML 0.5.4 vocabulary")
        protocol_names = [item.name for item in self.series.protocols]
        if len(protocol_names) != len(set(protocol_names)):
            duplicate = next(name for name in protocol_names if protocol_names.count(name) > 1)
            raise MINiMLModelError(f"duplicate protocol name: {duplicate}")
        known_protocols = set(protocol_names)
        known_documents = {item.name for item in self.source.documents if item.kind.casefold() == "sdrf"}
        for path in self.series.assay_paths:
            if path.document and known_documents and path.document not in known_documents:
                raise MINiMLModelError(f"unknown assay path document: {path.document}")
            for step in path.steps:
                if isinstance(step, ProtocolApplication):
                    if step.protocol_ref not in known_protocols:
                        warn("/series/assay_paths/protocol_ref", "external_protocol_reference", f"protocol {step.protocol_ref!r} is not declared in this package")
                    continue
                if step.sample_ref and step.sample_ref not in sample_ids:
                    warn("/series/assay_paths/sample_ref", "external_sample_reference", f"sample {step.sample_ref!r} is not declared in this package")
        self._validate_links(self.series.supplementary_data, "/series", warn)
        return tuple(issues)

    @staticmethod
    def _validate_links(links: tuple[SupplementLink, ...], path: str, warn: Callable[[str, str, str], None]) -> None:
        for index, link in enumerate(links):
            if link.checksum and not re.fullmatch(r"[0-9a-fA-F]{32}", link.checksum):
                warn(f"{path}/supplementary_data/{index}/checksum", "xsd_checksum", "checksum is not a 32-character MD5 value")
