# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Registered, defensive input adapters. Domain parsers remain the owners of meaning."""

from __future__ import annotations
import hashlib
import json
import re
import tarfile
from collections.abc import Mapping
from pathlib import Path
from xml.etree import ElementTree as ET

from meta_standards_converter.miniml import MINiMLCodec, MINiMLPackage
from meta_standards_converter.sources.json import JSONPackageSource
from .contracts import Diagnostic, InputError, LoadedInput
from .discovery import (
    detect,
    metadata_kind,
    read_text,
    xml_kind,
    is_url,
    is_anndata,
    KINDS,
)


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


class InputHandler:
    """One registered format with recognition, validation and loading boundaries."""

    def __init__(self, kind):
        self.kind = kind

    def probe(self, source, context):
        return (
            detect(source, context.runtime, context.profile.max_xml_bytes) == self.kind
        )

    def validate(self, spec, context):
        kind, source = self.kind, spec.sources
        if kind.endswith("_accession") and isinstance(source, Path):
            if not source.exists():
                raise InputError("missing_file", "Input path does not exist")
            raise InputError(
                "type_mismatch",
                "Accession handlers require an accession string, not a Path object",
            )
        for value in values(source):
            if (
                isinstance(value, (str, Path))
                and not is_url(str(value))
                and not str(value).lstrip().startswith(("{", "[", "<"))
            ):
                path = Path(value)
                if path.is_symlink() and path.is_dir():
                    raise InputError(
                        "directory_symlink", "Directory symlinks are not traversed"
                    )
        if kind.endswith("_accession"):
            token = str(source).strip().upper()
            if kind == "geo_accession" and not re.fullmatch(r"GSE\d+", token):
                raise InputError("invalid_accession", "Expected a GEO Series accession")
            if kind == "ae_accession" and not re.fullmatch(
                r"(?:E-[A-Z]+-\d+|S-[A-Z]+\d+)", token
            ):
                raise InputError(
                    "invalid_accession", "Expected an ArrayExpress/BioStudies accession"
                )
            if kind in {"sra_accession", "ena_accession"}:
                from meta_standards_converter.sources.archive_support import (
                    accession_kind,
                )

                try:
                    accession_kind(token)
                except ValueError:
                    raise InputError(
                        "invalid_accession", "Expected a supported INSDC accession"
                    )
        if kind == "anndata" and not is_anndata(source):
            raise InputError("type_mismatch", "Expected an AnnData object")

    def load(self, spec, context):
        self.validate(spec, context)
        kind, source, options = self.kind, spec.sources, context.input_options
        if kind in {"sra_records", "ena_records"}:
            from copy import deepcopy
            from meta_standards_converter.sources.archive_support import (
                StudyRecords,
                accession_kind,
                identifier,
            )
            from meta_standards_converter.miniml.sra_parser import SRAParser
            from meta_standards_converter.miniml.ena_parser import ENAParser

            if not isinstance(source, StudyRecords):
                raise InputError(
                    "type_mismatch", "Expected a native StudyRecords object"
                )
            accession_kind(source.seed.study)
            records = deepcopy(source)
            expected = {records.seed.study, records.seed.primary}
            references = {
                identifier(node)
                for root in records.xml
                for node in root.findall(".//EXPERIMENT/STUDY_REF")
            }
            references.update(
                row.get("secondary_study_accession") or row.get("study_accession")
                for name, rows in records.indexed.items()
                if name in {"read_run", "read_experiment"}
                for row in rows
            )
            references.discard(None)
            if not references or not references <= expected:
                raise InputError(
                    "missing_study_binding",
                    "Native records require matching explicit study references",
                )
            parser = SRAParser() if kind == "sra_records" else ENAParser()
            package = parser.parse(records)
            diagnostics = [
                Diagnostic("source_partial", issue, "load", "warning")
                for issue in records.issues
            ]
            return self._metadata(package, kind.split("_")[0], diagnostics)
        if kind in {"miniml", "json", "atlas", "curator"}:
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
                    payload = json.loads(read_text(path, context.profile.max_xml_bytes))
            found = metadata_kind(payload)
            if kind != "json" and found != kind:
                raise InputError(
                    "type_mismatch", f"Input does not satisfy the {kind} format"
                )
            loaded = context.service("package_source", JSONPackageSource).decode(
                payload, fallback=Path(path).stem if path else "input"
            )
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
        if kind.endswith("_accession"):
            return self._accession(spec, context)
        if kind == "magetab":
            from meta_standards_converter.converters.ae2json import AE2JSONConverter

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
                "ae2json",
                lambda: AE2JSONConverter(
                    resource_profile=context.profile,
                    source_hosts=context.runtime["allowed_hosts"],
                ),
            )
            packages = converter.convert(str(primary), sdrf_sources=explicit)
            return self._metadata(packages, provider="biostudies")
        if kind == "geo_archive":
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
            return self._geo_xml(texts[0], context)
        if kind == "geo_xml":
            return self._geo_xml(context.text(source), context)
        if kind in {"sra_xml", "ena_xml"}:
            return self._archive_xml(spec, context)
        if kind in {"h5ad", "anndata", "matrix", "fastq"}:
            from .expression import load_expression

            return load_expression(kind, spec, context)
        raise InputError("unsupported_input", "No loader registered for input")

    def _metadata(self, packages, provider=None, diagnostics=(), *, enriched=False):
        metadata = JSONPackageSource().decode(packages)
        return LoadedInput(
            metadata=metadata,
            provider=provider,
            enrichment_applied=enriched,
            content_id=canonical_identity(metadata),
            diagnostics=list(diagnostics),
        )

    def _geo_xml(self, text, context):
        if xml_kind(text) != "geo_xml":
            raise InputError("type_mismatch", "Expected GEO MINiML XML")
        from meta_standards_converter.miniml.geo_parser import GEOParser

        parser = context.service(
            "geo_parser", lambda: GEOParser(resource_profile=context.profile)
        )
        return self._metadata(
            parser.parse(
                text, remove_empty=context.input_options.get("remove_empty", True)
            ),
            "geo",
        )

    def _accession(self, spec, context):
        token = str(spec.sources).strip().upper()
        options = context.input_options
        if self.kind == "geo_accession":
            from meta_standards_converter.converters.geo2json import GEO2JSONConverter

            converter = context.service(
                "geo2json", lambda: GEO2JSONConverter(resource_profile=context.profile)
            )
            return self._metadata(
                converter.convert(token, **{**options, "enrich": False, "related_series": False}),
                "geo",
                enriched=False,
            )
        if self.kind == "ae_accession":
            from meta_standards_converter.converters.ae2json import AE2JSONConverter

            converter = context.service(
                "ae2json", lambda: AE2JSONConverter(resource_profile=context.profile)
            )
            return self._metadata(converter.convert(token, **options), "biostudies")
        from meta_standards_converter.converters.sra2json import SRA2JSONConverter
        from meta_standards_converter.converters.ena2json import ENA2JSONConverter

        provider = self.kind.split("_")[0]
        cls = SRA2JSONConverter if provider == "sra" else ENA2JSONConverter
        converter = context.service(
            provider + "2json", lambda: cls(resource_profile=context.profile)
        )
        result = converter.convert(token, **options)
        diagnostics = [
            Diagnostic("source_" + s.status, f"{s.study}: {issue}", "load", "warning")
            for s in result.studies
            for issue in s.issues
        ]
        if not result.packages:
            raise InputError("source_failed", "No native study package was retrieved")
        if not result.ok:
            diagnostics.append(
                Diagnostic(
                    "source_partial", "Native import was incomplete", "load", "warning"
                )
            )
        return self._metadata(result.packages, provider, diagnostics, enriched=True)

    def _archive_xml(self, spec, context):
        from meta_standards_converter.xml_safety import parse_xml
        from meta_standards_converter.sources.archive_support import (
            StudyRecords,
            StudySeed,
            identifier,
        )
        from meta_standards_converter.miniml.sra_parser import SRAParser
        from meta_standards_converter.miniml.ena_parser import ENAParser

        roots = []
        for value in values(spec.sources):
            text = context.text(value)
            found = xml_kind(text)
            if (
                self.kind == "ena_xml"
                and found != "ena_xml"
                or self.kind == "sra_xml"
                and found not in {"sra_xml", "ena_xml"}
            ):
                raise InputError(
                    "type_mismatch", "XML does not match selected archive handler"
                )
            root = parse_xml(text, max_bytes=context.profile.max_xml_bytes)
            if root.tag in {
                "EXPERIMENT_PACKAGE",
                "STUDY",
                "PROJECT",
                "SAMPLE",
                "EXPERIMENT",
                "RUN",
            }:
                parent = ET.Element(root.tag + "_SET")
                parent.append(root)
                root = parent
            roots.append(root)
        studies = sorted(
            {
                identifier(n)
                for r in roots
                for n in r.findall(".//EXPERIMENT/STUDY_REF")
                if identifier(n)
            }
        )
        if not studies:
            raise InputError(
                "missing_study_binding",
                "Archive XML requires experiments with explicit study references",
            )
        # Reject conflicting identities before any parser can keep an arbitrary first record.
        records_by_id = {}
        for root in roots:
            for node in root.iter():
                if node.tag not in {
                    "STUDY",
                    "PROJECT",
                    "SAMPLE",
                    "EXPERIMENT",
                    "RUN",
                } or not identifier(node):
                    continue
                key = (node.tag, identifier(node))
                fingerprint = ET.canonicalize(
                    ET.tostring(node, encoding="unicode"), strip_text=True
                )
                if key in records_by_id and records_by_id[key][0] != fingerprint:
                    raise InputError(
                        "conflicting_records",
                        "Archive bundle contains conflicting records for one identity",
                    )
                records_by_id[key] = (fingerprint, node)
        packages, diagnostics = [], []
        parser = SRAParser() if self.kind == "sra_xml" else ENAParser()
        for study in studies:
            primary = study
            selected_roots = roots
            if self.kind == "ena_xml":
                nodes = [node for _fingerprint, node in records_by_id.values()]
                study_node = next(
                    (
                        n
                        for n in nodes
                        if n.tag in {"STUDY", "PROJECT"}
                        and study
                        in {
                            identifier(n),
                            *(x.text for x in n.findall("IDENTIFIERS/*")),
                        }
                    ),
                    None,
                )
                if study_node is not None:
                    projects = [
                        n.text
                        for n in study_node.findall("IDENTIFIERS/*")
                        if n.text and n.text.startswith("PRJ")
                    ]
                    if projects:
                        primary = projects[0]
                experiments = [
                    n
                    for n in nodes
                    if n.tag == "EXPERIMENT"
                    and identifier(n.find("STUDY_REF")) in {study, primary}
                ]
                experiment_ids = {identifier(n) for n in experiments}
                sample_ids = {
                    identifier(n.find("DESIGN/SAMPLE_DESCRIPTOR")) for n in experiments
                }
                sample_ids.update(
                    identifier(m)
                    for n in experiments
                    for m in n.findall("DESIGN/SAMPLE_DESCRIPTOR/POOL/*")
                )
                selected = ET.Element("ROOT")
                if study_node is not None:
                    selected.append(study_node)
                for node in nodes:
                    aliases = {
                        identifier(node),
                        *(x.text for x in node.findall("IDENTIFIERS/*")),
                    }
                    if (
                        node in experiments
                        or node.tag == "SAMPLE"
                        and aliases & sample_ids
                        or node.tag == "RUN"
                        and identifier(node.find("EXPERIMENT_REF")) in experiment_ids
                    ):
                        selected.append(node)
                selected_roots = [selected]
            records = StudyRecords(StudySeed(study, primary), xml=selected_roots)
            package = parser.parse(records)
            packages.append(package)
            diagnostics.extend(
                Diagnostic("source_partial", issue, "load", "warning")
                for issue in records.issues
            )
        return self._metadata(packages, self.kind.split("_")[0], diagnostics)


def default_handlers():
    return {kind: InputHandler(kind) for kind in sorted(KINDS)}
