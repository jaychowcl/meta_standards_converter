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
from dataclasses import dataclass, field
from importlib.resources import files
import json
import os
from pathlib import Path
import re
from typing import Any, Callable, Mapping, TypeVar
from uuid import uuid4


MINIML_SCHEMA_VERSION = "1.0"


class MINiMLModelError(ValueError):
    """A package cannot be represented by the stable MINiML JSON model."""


class FrozenJSONMapping(Mapping[str, Any]):
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
        return FrozenJSONMapping(value)
    if isinstance(value, (list, tuple)):
        return tuple(_freeze_json(item) for item in value)
    return deepcopy(value)


@dataclass(frozen=True)
class MINiMLValidationIssue:
    path: str
    code: str
    message: str
    severity: str = "warning"


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


def miniml_schema_path() -> Path:
    return Path(str(files("meta_standards_converter.miniml").joinpath("miniml-package-v1.schema.json")))


def _mapping(value: Any, path: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise MINiMLModelError(f"{path} must be an object")
    return value


def _items(value: Any) -> list[Any]:
    if value is None:
        return []
    return list(value) if isinstance(value, (list, tuple)) else [value]


def _extras(data: Mapping[str, Any], known: set[str]) -> Mapping[str, Any]:
    return FrozenJSONMapping({key: value for key, value in data.items() if key not in known})


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


@dataclass(frozen=True)
class Organism:
    value: str
    taxid: str | None = None
    extras: Mapping[str, Any] = field(default_factory=dict)

    @classmethod
    def from_value(cls, value: Any) -> "Organism":
        if isinstance(value, Mapping):
            return cls(str(value.get("value", "")), None if value.get("taxid") is None else str(value["taxid"]), _extras(value, {"value", "taxid"}))
        return cls(str(value))

    def to_mapping(self) -> dict[str, Any]:
        result = {"value": self.value}
        _put(result, "taxid", self.taxid)
        return _record(result, self.extras)


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


@dataclass(frozen=True)
class Characteristics:
    value: str
    tag: str | None = None
    extras: Mapping[str, Any] = field(default_factory=dict)

    @classmethod
    def from_value(cls, value: Any) -> "Characteristics":
        if isinstance(value, Mapping):
            return cls(str(value.get("value", "")), value.get("tag"), _extras(value, {"value", "tag"}))
        return cls(str(value))

    def to_mapping(self) -> dict[str, Any]:
        result = {"value": self.value}
        _put(result, "tag", self.tag)
        return _record(result, self.extras)


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


@dataclass(frozen=True)
class Channel:
    source: str | None = None
    organisms: tuple[Organism, ...] = ()
    characteristics: tuple[Characteristics, ...] = ()
    biomaterial_providers: tuple[Any, ...] = ()
    treatment_protocol: str | None = None
    growth_protocol: str | None = None
    molecule: str | None = None
    extract_protocol: str | None = None
    label: str | None = None
    label_protocol: str | None = None
    position: str | None = None
    extras: Mapping[str, Any] = field(default_factory=dict)

    @classmethod
    def from_mapping(cls, value: Any) -> "Channel":
        data = _mapping(value, "channel")
        known = {"source", "organism", "characteristics", "biomaterial_provider", "treatment_protocol", "growth_protocol", "molecule", "extract_protocol", "label", "label_protocol", "position"}
        return cls(
            data.get("source"), tuple(Organism.from_value(item) for item in _items(data.get("organism"))),
            tuple(Characteristics.from_value(item) for item in _items(data.get("characteristics"))),
            tuple(_items(data.get("biomaterial_provider"))), data.get("treatment_protocol"),
            data.get("growth_protocol"), data.get("molecule"), data.get("extract_protocol"),
            data.get("label"), data.get("label_protocol"), data.get("position"), _extras(data, known),
        )

    def to_mapping(self) -> dict[str, Any]:
        result: dict[str, Any] = {}
        mapping = {"organisms": "organism", "biomaterial_providers": "biomaterial_provider"}
        for key in ("source", "organisms", "characteristics", "biomaterial_providers", "treatment_protocol", "growth_protocol", "molecule", "extract_protocol", "label", "label_protocol", "position"):
            _put(result, mapping.get(key, key), getattr(self, key))
        return _record(result, self.extras)


@dataclass(frozen=True)
class Variable:
    factor: str | None = None
    description: str | None = None
    sample_ref: tuple[Reference, ...] = ()
    position: str | None = None
    extras: Mapping[str, Any] = field(default_factory=dict)

    @classmethod
    def from_mapping(cls, value: Any) -> "Variable":
        data = _mapping(value, "variable")
        return cls(data.get("factor"), data.get("description"), _objects(data.get("sample_ref"), Reference.from_mapping, "variable.sample_ref"), data.get("position"), _extras(data, {"factor", "description", "sample_ref", "position"}))

    def to_mapping(self) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key in ("factor", "description", "sample_ref", "position"):
            _put(result, key, getattr(self, key))
        return _record(result, self.extras)


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
    position: str | None = None
    extras: Mapping[str, Any] = field(default_factory=dict)

    @classmethod
    def from_mapping(cls, value: Any) -> "Contributor":
        data = _mapping(value, "contributor")
        known = {"iid", "person", "organization", "company", "email", "phone", "fax", "laboratory", "department", "address", "organization_ref", "web_link", "position"}
        address = data.get("address")
        if isinstance(address, Mapping):
            address = Address.from_mapping(address)
        return cls(data.get("iid"), None if data.get("person") is None else Person.from_mapping(data["person"]), data.get("organization"), data.get("company"), data.get("email"), data.get("phone"), data.get("fax"), data.get("laboratory"), data.get("department"), address, None if data.get("organization_ref") is None else Reference.from_mapping(data["organization_ref"]), data.get("web_link"), data.get("position"), _extras(data, known))

    def to_mapping(self) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key in ("iid", "person", "organization", "company", "email", "phone", "fax", "laboratory", "department", "address", "organization_ref", "web_link", "position"):
            _put(result, key, getattr(self, key))
        return _record(result, self.extras)


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
    extras: Mapping[str, Any] = field(default_factory=dict)

    @classmethod
    def from_mapping(cls, value: Any) -> "Series":
        data = _mapping(value, "series")
        known = {"iid", "accession", "status", "title", "pubmed_id", "citation", "web_link", "summary", "overall_design", "type", "contributor_ref", "contributor", "contact_ref", "contact", "sample_ref", "variable", "repeats", "supplementary_data", "relation", "data_table", "pubmed_publication"}
        return cls(data.get("iid"), _accessions(data.get("accession")), _statuses(data.get("status")), data.get("title"), tuple(_items(data.get("pubmed_id"))), tuple(_items(data.get("citation"))), tuple(_items(data.get("web_link"))), data.get("summary"), data.get("overall_design"), tuple(_items(data.get("type"))), _refs(data.get("contributor_ref"), "series.contributor_ref"), _objects(data.get("contributor"), Contributor.from_mapping, "series.contributor"), _refs(data.get("contact_ref"), "series.contact_ref"), _objects(data.get("contact"), Contributor.from_mapping, "series.contact"), _refs(data.get("sample_ref"), "series.sample_ref"), _objects(data.get("variable"), Variable.from_mapping, "series.variable"), _objects(data.get("repeats"), Repeat.from_mapping, "series.repeats"), _links(data.get("supplementary_data")), _relations(data.get("relation")), _objects(data.get("data_table"), DataTable.from_mapping, "series.data_table"), _objects(data.get("pubmed_publication"), PubMedPublication.from_mapping, "series.pubmed_publication"), _extras(data, known))

    def to_mapping(self) -> dict[str, Any]:
        result: dict[str, Any] = {}
        keys = {"accessions": "accession", "statuses": "status", "pubmed_ids": "pubmed_id", "citations": "citation", "web_links": "web_link", "types": "type", "contributors": "contributor", "contacts": "contact", "variables": "variable", "relations": "relation", "data_tables": "data_table"}
        keys["pubmed_publications"] = "pubmed_publication"
        for key in ("iid", "accessions", "statuses", "title", "pubmed_ids", "citations", "web_links", "summary", "overall_design", "types", "contributor_ref", "contributors", "contact_ref", "contacts", "sample_ref", "variables", "repeats", "supplementary_data", "relations", "data_tables", "pubmed_publications"):
            _put(result, keys.get(key, key), getattr(self, key))
        return _record(result, self.extras)


TECHNOLOGIES = {"high-throughput sequencing", "in situ oligonucleotide", "spotted oligonucleotide", "mixed spotted oligonucleotide", "spotted DNA/cDNA", "spotted peptide or protein", "antibody", "tissue", "oligonucleotide beads", "MS", "SAGE NlaIII", "SAGE Sau3A", "SAGE RsaI", "SARST", "MPSS", "RT-PCR", "other"}
SAMPLE_TYPES = {"RNA", "genomic", "protein", "SAGE", "MPSS", "SARST", "mixed", "other", "SRA"}
MOLECULES = {"genomic DNA", "polyA RNA", "total RNA", "cytoplasmic RNA", "nuclear RNA", "protein", "other"}
VARIABLE_FACTORS = {"dose", "time", "tissue", "strain", "gender", "cell line", "development stage", "age", "agent", "cell type", "infection", "isolate", "metabolism", "shock", "stress", "temperature", "speciman", "disease state", "protocol", "growth protocol", "other", "genotype/variation", "species", "individual"}


@dataclass(frozen=True)
class MINiMLPackage(Mapping[str, Any]):
    series: Series
    databases: tuple[Database, ...] = ()
    organizations: tuple[Organization, ...] = ()
    contributors: tuple[Contributor, ...] = ()
    platforms: tuple[Platform, ...] = ()
    samples: tuple[Sample, ...] = ()
    version: str | None = None
    schema_location: str | None = None
    mage_tab: Mapping[str, Any] | None = None
    extras: Mapping[str, Any] = field(default_factory=dict)
    miniml_schema_version: str = MINIML_SCHEMA_VERSION

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
        if version is not None and version != MINIML_SCHEMA_VERSION:
            raise MINiMLModelError(f"unsupported MINiML schema version: {version!r}")
        series = Series.from_mapping(data.get("series"))
        if not (isinstance(series.iid, str) and series.iid.strip()) and not any(item.value.strip() for item in series.accessions):
            raise MINiMLModelError("series requires iid or accession")
        package = cls(
            series=series,
            databases=_objects(data.get("database"), Database.from_mapping, "database"),
            organizations=_objects(data.get("organization"), Organization.from_mapping, "organization"),
            contributors=_objects(data.get("contributor"), Contributor.from_mapping, "contributor"),
            platforms=_objects(data.get("platform"), Platform.from_mapping, "platform"),
            samples=_objects(data.get("sample"), Sample.from_mapping, "sample"),
            version=None if data.get("version") is None else str(data["version"]),
            schema_location=data.get("schema_location"),
            mage_tab=None if data.get("mage_tab") is None else FrozenJSONMapping(_mapping(data["mage_tab"], "mage_tab")),
            extras=_extras(data, {"miniml_schema_version", "version", "schema_location", "database", "organization", "contributor", "platform", "sample", "series", "mage_tab"}),
        )
        package._raise_structural_errors()
        return package

    @classmethod
    def load(cls, path: str | Path) -> "MINiMLPackage":
        with Path(path).open(encoding="utf-8") as handle:
            return cls.from_mapping(json.load(handle))

    def to_mapping(self) -> dict[str, Any]:
        result: dict[str, Any] = {
            "miniml_schema_version": self.miniml_schema_version,
            "database": [_plain(item) for item in self.databases],
            "organization": [_plain(item) for item in self.organizations],
            "contributor": [_plain(item) for item in self.contributors],
            "platform": [_plain(item) for item in self.platforms],
            "sample": [_plain(item) for item in self.samples],
            "series": self.series.to_mapping(),
        }
        _put(result, "version", self.version)
        _put(result, "schema_location", self.schema_location)
        _put(result, "mage_tab", self.mage_tab)
        return _record(result, self.extras)

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
                if channel.molecule and channel.molecule not in MOLECULES:
                    warn(f"/sample/{index}/channel/{channel_index}/molecule", "xsd_enumeration", f"{channel.molecule!r} is outside the MINiML 0.5.4 vocabulary")
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
        self._validate_links(self.series.supplementary_data, "/series", warn)
        return tuple(issues)

    @staticmethod
    def _validate_links(links: tuple[SupplementLink, ...], path: str, warn: Callable[[str, str, str], None]) -> None:
        for index, link in enumerate(links):
            if link.checksum and not re.fullmatch(r"[0-9a-fA-F]{32}", link.checksum):
                warn(f"{path}/supplementary_data/{index}/checksum", "xsd_checksum", "checksum is not a 32-character MD5 value")
