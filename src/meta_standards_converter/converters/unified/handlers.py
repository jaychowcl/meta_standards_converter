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
from meta_standards_converter.conversion_logging import measured
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


from .loaders import (values, safe_name, canonical_identity, metadata_input,
                      load_json, load_magetab, load_geo_archive, load_geo_xml)


class InputHandler:
    """One registered format with recognition, validation and loading boundaries."""

    def __init__(self, kind):
        self.kind = kind

    _metadata = staticmethod(metadata_input)

    def probe(self, source, context):
        return (
            detect(source, context.runtime, context.profile.max_xml_bytes) == self.kind
        )

    @measured("input_validation")
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

    @measured("load")
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
            return metadata_input(package, kind.split("_")[0], diagnostics)
        if kind in {"miniml", "json", "atlas", "curator"}:
            return load_json(kind, spec, context)
        if kind.endswith("_accession"):
            return self._accession(spec, context)
        if kind == "magetab":
            return load_magetab(kind, spec, context)
        if kind == "geo_archive":
            return load_geo_archive(kind, spec, context)
        if kind == "geo_xml":
            return load_geo_xml(context.text(source), context)
        if kind in {"sra_xml", "ena_xml"}:
            return self._archive_xml(spec, context)
        if kind in {"h5ad", "anndata", "matrix", "fastq"}:
            from .expression import load_expression

            return load_expression(kind, spec, context)
        raise InputError("unsupported_input", "No loader registered for input")

    def _accession(self, spec, context):
        token = str(spec.sources).strip().upper()
        if self.kind == "geo_accession":
            packages = context.session.load(context.geo_source(), token, remove_empty=context.input_options.get("remove_empty", True))
            return metadata_input(packages, "geo")
        if self.kind == "ae_accession":
            from meta_standards_converter.sources.errors import MetadataSourceError
            try:
                packages = context.session.load(context.ae_source(), token, **context.input_options)
                return metadata_input(packages, "biostudies")
            except MetadataSourceError as error:
                if error.code not in error.FALLBACK_CODES:
                    raise
                return self._ae_fallback(token, context, error)
        provider = self.kind.split("_")[0]
        return self._native(token, provider, context, **context.input_options)

    @staticmethod
    def _native_source(provider, context):
        from meta_standards_converter.sources.metadata import NativeMetadataSource
        return context.service(provider + "_source", lambda: NativeMetadataSource(provider, resource_profile=context.profile))

    def _ae_fallback(self, token, context, original):
        from meta_standards_converter.metadata.source_fallback import annotate
        from .preparation import study_ids
        from dataclasses import replace
        original_diagnostic = Diagnostic(original.code, str(original), 'load', 'warning')
        try:
            source = self._native_source('ena', context)
            mapping = context.session.retrieve(('ae_fallback', id(source), token),
                lambda: source.resolve_arrayexpress(token))
            loaded = self._native(mapping.study, 'ena', context)
            groups = []
            for group in loaded.metadata.groups:
                packages = []
                for package in group.packages:
                    if not set(study_ids(package)) & {mapping.study, mapping.project}:
                        raise InputError('fallback_invalid', 'ENA metadata does not match the verified study; supply the ArrayExpress IDF and SDRF files.')
                    packages.append(annotate(package, token, mapping, original.code))
                groups.append(replace(group, packages=tuple(packages)))
            loaded.metadata = replace(loaded.metadata, groups=tuple(groups))
            loaded.content_id = canonical_identity(loaded.metadata)
            loaded.diagnostics.insert(0, original_diagnostic)
            return loaded
        except Exception as error:
            # Keep the original condition and the independent fallback failure.
            # The final exception retains its own retry eligibility.
            error._msc_source_diagnostics = [original_diagnostic, Diagnostic(
                'fallback_failed', 'Verified ENA fallback could not provide usable metadata. Supply the ArrayExpress IDF and SDRF files.', 'load')]
            raise

    def _native(self, token, provider, context, **options):
        source = self._native_source(provider, context)
        result = context.session.load(source, token, publications=context.preparation_policy.allows("publications"), **options)
        diagnostics = [
            Diagnostic("source_" + s.status, f"{s.study}: {issue}", "load", "warning")
            for s in result.studies
            for issue in s.issues
        ]
        from ..archive_results import native_failure
        error = native_failure(getattr(s, 'error', None) for s in result.studies)
        if not result.packages:
            if error is not None:
                raise error
            from meta_standards_converter.metadata.retrieval_budget import _active
            budget = _active.get()
            detail = f'; metadata budget exhausted: {budget.events}' if budget and budget.exhausted else ''
            raise InputError("source_failed", "No native study package was retrieved" + detail)
        if not result.ok:
            diagnostics.append(
                Diagnostic(
                    "source_partial", "Native import was incomplete", "load", "warning"
                )
            )
        loaded = metadata_input(result.packages, provider, diagnostics, enriched=True)
        loaded._native_failure = error
        return loaded

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
        return metadata_input(packages, self.kind.split("_")[0], diagnostics)


def default_handlers():
    return {kind: InputHandler(kind) for kind in sorted(KINDS)}
