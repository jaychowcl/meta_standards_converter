# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Focused metadata input loaders; parsers retain ownership of format semantics."""
import hashlib
import json
import re
import tarfile
from pathlib import Path
from meta_standards_converter.sources.json import JSONPackageSource
from .contracts import InputError, LoadedInput, Diagnostic
from .discovery import is_url, metadata_kind, read_text, xml_kind

def values(source):
    return list(source) if isinstance(source, (tuple, list)) else [source]


def safe_name(value):
    if (
        not isinstance(value, str)
        or value in {"", ".", ".."}
        or not re.fullmatch(r"[A-Za-z0-9_.-]+", value)
    ):
        raise InputError(
            "unsafe_identity",
            "Dataset/output identity must be a safe filename component",
        )
    return value


def canonical_identity(loaded):
    payload = [p.to_mapping() for g in loaded.groups for p in g.packages]
    return hashlib.sha256(
        json.dumps(
            payload, sort_keys=True, ensure_ascii=False, separators=(",", ":")
        ).encode()
    ).hexdigest()



def metadata_input(packages, provider=None, diagnostics=(), *, enriched=False):
    metadata = JSONPackageSource().decode(packages, newly_parsed=True)
    return LoadedInput(
        metadata=metadata,
        provider=provider,
        enrichment_applied=enriched,
        content_id=canonical_identity(metadata),
        diagnostics=list(diagnostics),
    )

def load_geo_xml(text, context):
    if xml_kind(text) != "geo_xml":
        raise InputError("type_mismatch", "Expected GEO MINiML XML")
    from meta_standards_converter.miniml.geo_parser import GEOParser

    parser = context.service(
        "geo_parser", lambda: GEOParser(resource_profile=context.profile)
    )
    return metadata_input(
        parser.parse(
            text, remove_empty=context.input_options.get("remove_empty", True)
        ),
        "geo",
    )


def load_json(kind, spec, context):
    source, options = spec.sources, context.input_options
    payload = source
    path = None
    if isinstance(source, (str, Path)):
        if isinstance(source, str) and source.lstrip().startswith(("{", "[")):
            if len(source.encode()) > context.profile.max_xml_bytes:
                raise InputError(
                    "input_too_large", "JSON exceeds input byte limit"
                )
            payload = json.loads(source)
        else:
            path = context.localize(source)
            with Path(path).open("rb") as stream:
                raw = stream.read(context.profile.max_xml_bytes + 1)
            if len(raw) > context.profile.max_xml_bytes:
                raise InputError("input_too_large", "JSON exceeds input byte limit")
            payload = json.loads(raw)
    found = metadata_kind(payload)
    if kind != "json" and found != kind:
        raise InputError(
            "type_mismatch", f"Input does not satisfy the {kind} format"
        )
    service = context.service("package_source", JSONPackageSource)
    loaded = (service.decode_file(path, raw) if path and hasattr(service, "decode_file") else
              service.decode(payload, fallback=Path(path).stem if path else "input"))
    return LoadedInput(
        metadata=loaded,
        source_path=path,
        content_id=canonical_identity(loaded),
        diagnostics=(
            [
                Diagnostic("source_partial", warning, "load", "warning")
                for warning in loaded.warnings
            ]
            if found == "atlas"
            else []
        ),
    )


def load_magetab(kind, spec, context):
    source, options = spec.sources, context.input_options
    from meta_standards_converter.sources.metadata import MAGETabMetadataSource

    sources = values(source)
    if any(
        isinstance(v, (str, Path)) and str(v).lstrip().startswith("<")
        for v in sources
    ):
        raise InputError(
            "type_mismatch", "MAGE-TAB requires an IDF and associated SDRFs"
        )
    explicit = spec.companions.get("sdrf", options.get("sdrf_sources"))
    if len(sources) > 1:
        idfs = [v for v in sources if ".idf." in str(v).lower()]
        if len(idfs) != 1:
            raise InputError(
                "ambiguous_bundle", "MAGE-TAB bundle requires exactly one IDF"
            )
        primary = idfs[0]
        explicit = [str(v) for v in sources if v != primary]
    else:
        primary = sources[0]
    from .discovery import idf_candidates, idf_references

    if ".sdrf." in str(primary).lower():
        explicit = [str(primary)]
        if spec.companions.get("idf"):
            primary = spec.companions["idf"]
        elif not is_url(str(primary)):
            candidates = idf_candidates(primary, context.profile.max_xml_bytes)
            if len(candidates) != 1:
                raise InputError(
                    "ambiguous_companion",
                    "SDRF requires one explicitly bound or uniquely associated IDF",
                )
            primary = candidates[0]
        else:
            raise InputError(
                "missing_companion",
                "Remote SDRF requires an explicit IDF companion",
            )
    elif explicit is None and not is_url(str(primary)):
        for ref in idf_references(primary, context.profile.max_xml_bytes):
            if (
                not is_url(ref)
                and len(
                    idf_candidates(
                        Path(primary).parent / ref,
                        context.profile.max_xml_bytes,
                    )
                )
                > 1
            ):
                raise InputError(
                    "ambiguous_companion",
                    "Multiple IDFs reference the same SDRF; bind the intended inputs explicitly",
                )
    if isinstance(explicit, (str, Path)):
        explicit = [str(explicit)]
    converter = context.service(
        "ae_source",
        lambda: MAGETabMetadataSource(
            resource_profile=context.profile,
            allowed_hosts=context.runtime["allowed_hosts"],
        ),
    )
    packages = context.session.load(converter, str(primary), sdrf_sources=explicit)
    return metadata_input(packages, provider="biostudies")


def load_geo_archive(kind, spec, context):
    source, options = spec.sources, context.input_options
    sources = values(source)
    if len(sources) != 1:
        raise InputError(
            "ambiguous_bundle", "A GEO archive input requires one archive"
        )
    path = context.localize(sources[0])
    profile = context.profile
    if Path(path).stat().st_size > profile.max_compressed_archive_bytes:
        raise InputError("input_too_large", "GEO archive is too large")
    from meta_standards_converter.sources.geo import (
        _normalise_archive_member_name,
    )

    texts, names, expanded = [], set(), 0
    with tarfile.open(path, mode="r|gz") as archive:
        for index, member in enumerate(archive):
            if index >= 10000:
                raise InputError("invalid_archive", "Too many archive members")
            name = _normalise_archive_member_name(member.name)
            if name in names or not (member.isfile() or member.isdir()):
                raise InputError(
                    "invalid_archive", "Unsafe or duplicate archive member"
                )
            names.add(name)
            expanded += member.size
            if expanded > profile.max_expanded_archive_bytes:
                raise InputError(
                    "input_too_large", "Expanded archive is too large"
                )
            if member.isfile() and name.lower().endswith(".xml"):
                if member.size > profile.max_xml_bytes:
                    raise InputError(
                        "input_too_large", "Archive XML is too large"
                    )
                texts.append(
                    archive.extractfile(member)
                    .read(profile.max_xml_bytes + 1)
                    .decode("utf-8-sig")
                )
    if len(texts) != 1:
        raise InputError(
            "ambiguous_archive",
            "GEO family archive requires exactly one XML document",
        )
    return load_geo_xml(texts[0], context)
