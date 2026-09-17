# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Bounded recognition and deterministic logical-input discovery."""

from __future__ import annotations

import json
import os
import re
from collections.abc import Mapping
from dataclasses import replace
from pathlib import Path
from urllib.parse import urlsplit

from meta_standards_converter.miniml import MINiMLPackage
from meta_standards_converter.sources.archive_support import accession_kind
from .contracts import InputSpec, InputError

ALIASES = {
    "sra": "sra_accession",
    "ena": "ena_accession",
    "geo": "geo_accession",
    "msc_miniml": "miniml",
    "atlas_v1": "atlas",
}
KINDS = {
    "sra_records",
    "ena_records",
    "geo_accession",
    "ae_accession",
    "sra_accession",
    "ena_accession",
    "miniml",
    "json",
    "atlas",
    "curator",
    "geo_xml",
    "geo_archive",
    "magetab",
    "sra_xml",
    "ena_xml",
    "h5ad",
    "anndata",
    "matrix",
    "fastq",
}
CACHE_NAMES = {
    ".cache",
    ".processed",
    ".git",
    ".dev",
    ".artifact-bundles",
    "__pycache__",
    "work",
}


def kind_name(value):
    if value is None:
        return None
    value = ALIASES.get(value, value)
    if value not in KINDS:
        raise ValueError(f"Unknown input type: {value}")
    return value


def is_url(value):
    return isinstance(value, str) and bool(
        re.match(r"^[a-zA-Z][a-zA-Z0-9+.-]*://", value)
    )


def package_value(value):
    return (
        isinstance(value, MINiMLPackage)
        or isinstance(value, Mapping)
        and "miniml_schema_version" in value
    )


def metadata_kind(value):
    if isinstance(value, MINiMLPackage):
        return "miniml"
    if isinstance(value, (tuple, list)):
        return "miniml" if value and all(package_value(p) for p in value) else None
    if isinstance(value, Mapping):
        matches = []
        if "miniml_schema_version" in value:
            matches.append("miniml")
        if "miniml_json" in value:
            matches.append("curator")
        if "atlas" in value or "datasets" in value:
            matches.append("atlas")
        if len(matches) > 1:
            raise InputError(
                "ambiguous_input", "Conflicting metadata envelope signatures"
            )
        return matches[0] if matches else None
    return None


def is_anndata(value):
    # Inspect the type without touching .X: backed arrays must stay lazy during probing.
    return any(
        cls.__name__ == "AnnData" and cls.__module__.split(".")[0] == "anndata"
        for cls in type(value).__mro__
    )


def read_text(path, limit):
    with Path(path).open("rb") as stream:
        data = stream.read(limit + 1)
    if len(data) > limit:
        raise InputError("input_too_large", "Metadata exceeds the input byte limit")
    return data.decode("utf-8-sig")


def xml_kind(text):
    # Full bounded secure parsing is performed by validate/load, never an XML regex parser.
    from meta_standards_converter.xml_safety import parse_xml

    root = parse_xml(text, max_bytes=max(1, len(text.encode("utf-8"))))
    name = root.tag.rsplit("}", 1)[-1].upper()
    if name in {"MINIML", "MINIML_SET"}:
        return "geo_xml"
    if name in {"EXPERIMENT_PACKAGE_SET", "EXPERIMENT_PACKAGE"}:
        return "sra_xml"
    if name in {
        "STUDY_SET",
        "PROJECT_SET",
        "SAMPLE_SET",
        "EXPERIMENT_SET",
        "RUN_SET",
        "STUDY",
        "PROJECT",
        "SAMPLE",
        "EXPERIMENT",
        "RUN",
    }:
        return "ena_xml"
    raise InputError(
        "unsupported_xml", "XML root does not identify a supported source format"
    )


def tenx_member(path):
    from meta_standards_converter.expression.readers import _tenx_member

    name = Path(urlsplit(str(path)).path).name.lower()
    if name.endswith(".gz"):
        name = name[:-3]
    roles = {
        "matrix.mtx": "matrix",
        "barcodes.tsv": "barcodes",
        "features.tsv": "features",
        "genes.tsv": "features",
    }
    return ("", roles[name]) if name in roles else _tenx_member(str(path))


def tenx_directory(path):
    path = Path(path)
    if not path.is_dir():
        return False
    members = [tenx_member(member) for member in path.iterdir()]
    return (
        len(members) == 3
        and all(member and member[0] == "" for member in members)
        and {member[1] for member in members} == {"matrix", "features", "barcodes"}
    )


def detect(value, runtime, limit):
    from meta_standards_converter.sources.archive_support import StudyRecords

    if isinstance(value, StudyRecords):
        return (
            "sra_records"
            if any(
                root.tag in {"EXPERIMENT_PACKAGE_SET", "EXPERIMENT_PACKAGE"}
                for root in value.xml
            )
            else "ena_records"
        )
    kind = metadata_kind(value)
    if kind:
        return kind
    if is_anndata(value):
        return "anndata"
    if isinstance(value, (list, tuple)):
        kinds = {detect(v, runtime, limit) for v in value}
        if len(kinds) == 1 and kinds <= {"sra_xml", "ena_xml", "fastq"}:
            return kinds.pop()
        raise InputError(
            "ambiguous_bundle", "Specify the format and roles of compound input sources"
        )
    if not isinstance(value, (str, os.PathLike)):
        raise InputError("unsupported_input", "Unsupported input object")
    token = str(value)
    if token.lstrip().startswith(("{", "[", "<")):
        raise InputError(
            "explicit_type_required",
            "Inline JSON/XML requires force_in_type or InputSpec.in_type",
        )
    if is_url(token):
        name = urlsplit(token).path.lower()
        if ".idf." in name or ".sdrf." in name:
            return "magetab"
        if name.endswith(".json"):
            return "json"
        if name.endswith(".h5ad"):
            return "h5ad"
        if name.endswith((".tgz", ".tar.gz")):
            return "geo_archive"
        if name.endswith(".xml"):
            return "remote_xml"
        return file_kind(name)
    path = Path(token)
    if path.exists():
        if path.is_dir():
            if path.is_symlink():
                raise InputError(
                    "directory_symlink", "Directory symlinks are not traversed"
                )
            if tenx_directory(path):
                return "matrix"
            raise InputError(
                "ambiguous_directory", "Directory must be expanded into logical inputs"
            )
        if not path.is_file():
            raise InputError("unsupported_input", "Input must be a regular file")
        name = path.name.lower()
        if name.endswith(
            (".h5ad", ".h5", ".fastq", ".fq", ".fastq.gz", ".fq.gz", ".tgz", ".tar.gz")
        ):
            return file_kind(name)
        with path.open("rb") as stream:
            head = stream.read(min(4096, limit))
        stripped = head.lstrip(b"\xef\xbb\xbf \r\n\t")
        if stripped.startswith((b"{", b"[")):
            obj = json.loads(read_text(path, limit))
            return metadata_kind(obj) or "json"
        if stripped.startswith(b"<"):
            return xml_kind(read_text(path, limit))
        if ".idf." in name or ".sdrf." in name or b"SDRF File\t" in head:
            return "magetab"
        return file_kind(name)
    if isinstance(value, os.PathLike):
        raise InputError("missing_file", "Input path does not exist")
    normalized = token.strip().upper()
    if re.fullmatch(r"GSE\d+", normalized):
        return "geo_accession"
    if re.fullmatch(r"(?:E-[A-Z]+-\d+|S-[A-Z]+\d+)", normalized):
        return "ae_accession"
    try:
        accession, _ = accession_kind(token)
    except ValueError:
        if (
            isinstance(value, os.PathLike)
            or "/" in token
            or "\\" in token
            or Path(token).suffix
        ):
            raise InputError("missing_file", "Input path does not exist")
        raise InputError("unsupported_input", "Unrecognized accession or input")
    if accession.startswith(("SR", "PRJNA", "SAMN")):
        return "sra_accession"
    if accession.startswith(("ER", "PRJEB", "SAMEA")):
        return "ena_accession"
    return runtime["insdc_default"] + "_accession"


def file_kind(name):
    lower = name.lower()
    if lower.endswith((".tgz", ".tar.gz")):
        return "geo_archive"
    for suffix in (".gz", ".bz2", ".xz", ".zip"):
        if lower.endswith(suffix):
            lower = lower[: -len(suffix)]
            break
    if lower.endswith(".h5ad"):
        return "h5ad"
    if lower.endswith((".fastq", ".fq")):
        return "fastq"
    if lower.endswith((".h5", ".mtx", ".csv", ".tsv", ".txt")):
        return "matrix"
    raise InputError(
        "unsupported_input", "File content does not identify a supported format"
    )


def idf_references(path, limit):
    import csv

    rows = csv.reader(read_text(path, limit).splitlines(), delimiter="\t")
    return [
        value
        for row in rows
        if row and row[0].strip().casefold() == "sdrf file"
        for value in row[1:]
        if value
    ]


def idf_candidates(sdrf, limit):
    path = Path(sdrf).resolve()
    return [
        idf
        for idf in sorted(path.parent.glob("*.idf.txt"))
        if path
        in {
            (idf.parent / ref).resolve()
            for ref in idf_references(idf, limit)
            if not is_url(ref)
        }
    ]


def manifest_specs(value, limit):
    if isinstance(value, Mapping):
        payload, base = value, Path.cwd()
    else:
        path = Path(value)
        payload, base = json.loads(read_text(path, limit)), path.resolve().parent
    if (
        not isinstance(payload, Mapping)
        or set(payload) != {"schema_version", "inputs"}
        or payload["schema_version"] != "1.0"
    ):
        raise ValueError("Input manifest requires schema_version 1.0 and inputs")
    if not isinstance(payload["inputs"], list) or not payload["inputs"]:
        raise ValueError("Manifest inputs must be nonempty")
    result, seen = [], set()

    def resolve(v):
        if isinstance(v, list):
            return [resolve(x) for x in v]
        if isinstance(v, Mapping):
            return v
        if isinstance(v, str) and not is_url(v):
            if v.lstrip().startswith(("{", "[", "<")):
                return v
            if re.fullmatch(
                r"(?:GSE\d+|[SED]R[PRXS]\d+|PRJ(?:NA|EB|DB|DA)\d+|SAM(?:N|EA|D)\d+|E-[A-Z]+-\d+|S-[A-Z]+\d+)",
                v,
                re.I,
            ):
                return v
            return str(base / v)
        return v

    for entry in payload["inputs"]:
        allowed = {
            "id",
            "sources",
            "in_type",
            "companions",
            "metadata",
            "input_options",
            "output_options",
        }
        if (
            not isinstance(entry, Mapping)
            or set(entry) - allowed
            or not {"id", "sources"} <= set(entry)
        ):
            raise ValueError("Invalid manifest input entry")
        if not isinstance(entry["id"], str) or not entry["id"] or entry["id"] in seen:
            raise ValueError("Manifest input IDs must be unique nonempty strings")
        if not isinstance(entry["sources"], list) or not entry["sources"]:
            raise ValueError("Manifest sources must be a nonempty list")
        seen.add(entry["id"])
        fields = dict(entry)
        fields["sources"] = resolve(fields["sources"])
        if len(fields["sources"]) == 1:
            fields["sources"] = fields["sources"][0]
        if not isinstance(fields.get("companions", {}), Mapping):
            raise ValueError("Manifest companions must be an object")
        fields["companions"] = {
            k: resolve(v) for k, v in fields.get("companions", {}).items()
        }
        if isinstance(fields.get("metadata"), str):
            fields["metadata"] = resolve(fields["metadata"])
        for group in ("input_options", "output_options"):
            fields[group] = dict(fields.get(group, {}))
            for key in (
                "asset_manifest",
                "fasta",
                "gtf",
                "gff",
                "params_file",
                "nextflow_config",
                "work_dir",
                "processed_checkpoint_dir",
                "evidence_dir",
                "sdrf_sources",
            ):
                if key in fields[group]:
                    fields[group][key] = resolve(fields[group][key])
        result.append(InputSpec(**fields))
    return result


def expand(value, runtime, excluded=()):
    if isinstance(value, InputSpec):
        specs = [value]
    elif isinstance(value, (list, tuple)) and not metadata_kind(value):
        specs = [v if isinstance(v, InputSpec) else InputSpec(v) for v in value]
    elif value is None:
        specs = []
    else:
        specs = [InputSpec(value)]
    result = []
    excluded = [Path(p).resolve() for p in excluded if p is not None]

    def children(spec, directory):
        paths = sorted(directory.iterdir(), key=lambda p: p.name)
        consumed = set()
        grouped = {}

        matrix_groups = {}
        xml_groups = {}
        for path in paths:
            if not path.is_file() or any(
                path.resolve() == p or path.resolve().is_relative_to(p)
                for p in excluded
            ):
                continue
            member = tenx_member(path)
            if member:
                prefix, role = member
                matrix_groups.setdefault(prefix, {}).setdefault(role, []).append(path)
            if path.suffix.lower() == ".xml":
                try:
                    kind = xml_kind(read_text(path, 32 * 1024 * 1024))
                    if kind in {"sra_xml", "ena_xml"}:
                        xml_groups.setdefault(kind, []).append(path)
                except (OSError, ValueError, UnicodeError):
                    pass
        for members in matrix_groups.values():
            if all(
                len(members.get(role, [])) == 1
                for role in ("matrix", "features", "barcodes")
            ):
                primary = members["matrix"][0]
                grouped[primary] = replace(
                    spec,
                    sources=primary,
                    id=None,
                    companions={
                        **spec.companions,
                        "features": str(members["features"][0]),
                        "barcodes": str(members["barcodes"][0]),
                    },
                )
                consumed.update(
                    members[role][0].resolve() for role in ("features", "barcodes")
                )
        for kind, members in xml_groups.items():
            grouped[members[0]] = replace(
                spec, sources=members, id=None, in_type=spec.in_type
            )
            consumed.update(p.resolve() for p in members[1:])
        # IDF declarations bind their SDRFs before standalone discovery.
        for path in paths:
            if path.is_file() and path.name.lower().endswith(".idf.txt"):
                try:
                    for line in read_text(path, 8 * 1024 * 1024).splitlines():
                        parts = line.split("\t")
                        if parts[0].strip().casefold() == "sdrf file":
                            consumed.update(
                                (path.parent / p.strip('"')).resolve()
                                for p in parts[1:]
                                if p
                            )
                except (OSError, ValueError, UnicodeError):
                    pass
        for path in paths:
            if path.resolve() in consumed or path.name in CACHE_NAMES:
                continue
            if any(
                path.resolve() == p or path.resolve().is_relative_to(p)
                for p in excluded
            ):
                continue
            if path.is_symlink() and path.is_dir():
                continue
            if path.is_dir() and not path.is_symlink() and not tenx_directory(path):
                if runtime["recursive"]:
                    children(replace(spec, id=None), path)
                continue
            result.append(grouped.get(path, replace(spec, sources=path, id=None)))

    for spec in specs:
        value = spec.sources
        if (
            isinstance(value, (str, os.PathLike))
            and not is_url(str(value))
            and not str(value).lstrip().startswith(("{", "[", "<"))
        ):
            path = Path(value)
            if path.is_dir() and not path.is_symlink() and not tenx_directory(path):
                before = len(result)
                children(spec, path)
                if len(result) == before and (not excluded or not any(path.iterdir())):
                    result.append(spec)
                continue
        result.append(spec)
    return result
