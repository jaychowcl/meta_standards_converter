# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Unified request orchestration; handlers and domain services perform conversion."""

from __future__ import annotations
import os
import tempfile
from collections.abc import Mapping
from contextlib import ExitStack
from pathlib import Path
from urllib.parse import urlsplit

from meta_standards_converter.runtime_contracts import get_resource_profile
from meta_standards_converter.retrieval import RetrievalPolicy, RetrievalService
from meta_standards_converter.sources.json import SourceLoadResult
from .contracts import (
    ConversionBatchResult,
    ConversionItemResult,
    Diagnostic,
    InputError,
)
from .discovery import (
    detect,
    expand,
    kind_name,
    manifest_specs,
    is_url,
    read_text,
    xml_kind,
)
from .handlers import default_handlers
from .options import TARGETS, INPUT, OUTPUT, RUNTIME_DEFAULTS, settings
from .publication import publish_jsons, publish_tables, publish_magetab
from .routes import select_route
from .requests import (AUTO, STANDARD, normalize, input_settings, output_settings,
                       preparation_options, merge_input_options, apply_preparation_overrides)


def _missing_scientific_dependency(exc):
    scientific = {"anndata", "numpy", "pandas", "scipy", "h5py", "scanpy"}
    return any(
        isinstance(error, ImportError)
        and (error.name or "").split(".")[0] in scientific
        for error in (exc, exc.__cause__)
    )


def source_label(source):
    if isinstance(source, (str, os.PathLike)):
        value = str(source)
        if value.lstrip().startswith(("{", "[", "<")):
            return "inline document"
        if is_url(value):
            url = urlsplit(value)
            return f'{url.scheme}://{url.hostname or ""}{url.path}'
        return value
    if isinstance(source, (tuple, list)):
        return "input collection"
    return type(source).__name__


class ExecutionContext:
    def __init__(self, converter, runtime, target, outfile, outdir, stack):
        self.converter = converter
        self.runtime = runtime
        self.target = target
        self.outfile = Path(outfile) if outfile is not None else None
        self.outdir = Path(outdir) if outdir is not None else None
        self.destination = self.outfile or self.outdir
        self.stack = stack
        self.reserved = {Path(p).resolve() for p in runtime["reserved_paths"]}
        self.input_options = {}
        self.output_options = {}
        self.profile = get_resource_profile(
            runtime["resource_profile"], overrides=runtime["resource_overrides"]
        )
        self.policy = RetrievalPolicy(
            resource_profile=self.profile,
            allowed_hosts=frozenset(runtime["allowed_hosts"]),
        )
        self._retrieval = None
        self._services = {}
        from .preparation import MetadataPreparation
        self.preparer = MetadataPreparation(self)
        self.stage = "input"

    def service(self, key, factory):
        # Default instances are operation-local. Injected instances retain caller ownership.
        if key in self.converter.services:
            return self.converter.services[key]
        if key not in self._services:
            direct = self.converter.services.get("geo2ae")
            if direct is not None and key == "geo2json":
                from meta_standards_converter.converters.geo2json import GEO2JSONConverter
                self._services[key] = GEO2JSONConverter(
                    geo_fetcher=direct.geo_fetcher, parser=direct.parser,
                    enricher=direct.enricher, resource_profile=self.profile)
            elif direct is not None and key == "json2ae":
                from meta_standards_converter.converters.json2ae import JSON2AEConverter
                self._services[key] = JSON2AEConverter(
                    enricher=direct.enricher, ae_constructor=direct.ae_constructor)
            else:
                self._services[key] = factory()
        return self._services[key]

    def localize(self, source):
        source = str(source)
        if is_url(source):
            if self._retrieval is None:
                cache = self.stack.enter_context(
                    tempfile.TemporaryDirectory(prefix="msc-input-")
                )
                self._retrieval = self.service(
                    "retrieval", lambda: RetrievalService(cache, policy=self.policy)
                )
            return self._retrieval.localize(source)
        path = Path(source)
        if not path.exists():
            raise InputError("missing_file", "Input path does not exist")
        return str(path)

    def text(self, source):
        if isinstance(source, str) and source.lstrip().startswith("<"):
            if len(source.encode()) > self.profile.max_xml_bytes:
                raise InputError("input_too_large", "XML exceeds input byte limit")
            return source
        return read_text(self.localize(source), self.profile.max_xml_bytes)


class Converter:
    """Route validated input objects to supported destinations.

    ``services`` injects existing converters or sources by their documented names.
    ``handlers`` adds/replaces handlers implementing kind/probe/validate/load.
    """

    def __init__(self, *, services=None, handlers=None):
        self.services = dict(services or {})
        self.handlers = default_handlers()
        self.custom_handlers = tuple(handlers or ())
        if handlers:
            self.handlers.update({h.kind: h for h in self.custom_handlers})

    def _detect(self, source, context):
        matches = [
            handler.kind
            for handler in self.custom_handlers
            if handler.probe(source, context)
        ]
        if len(matches) > 1:
            raise InputError(
                "ambiguous_input", "Multiple registered handlers recognize the input"
            )
        return (
            matches[0]
            if matches
            else detect(source, context.runtime, context.profile.max_xml_bytes)
        )

    def convert(
        self,
        input=None,
        outdir=None,
        *,
        out_type,
        in_type=AUTO,
        enrichment=STANDARD,
        options=None,
        force_in_type=None,
        outfile=None,
        input_manifest=None,
        input_options=None,
        output_options=None,
        runtime_options=None,
    ):
        if out_type not in TARGETS:
            raise ValueError(f"Unsupported output type: {out_type}")
        forced, outfile, input_manifest, input_options, output_options, runtime_options = normalize(
            target=out_type, in_type=in_type, enrichment=enrichment, options=options,
            force_in_type=force_in_type, outfile=outfile, input_manifest=input_manifest,
            input_options=input_options, output_options=output_options,
            runtime_options=runtime_options)
        if outfile is not None and outdir is not None:
            raise ValueError("outfile and outdir are mutually exclusive")
        if outfile is not None and out_type == "magetab":
            raise ValueError("MAGE-TAB requires outdir")
        runtime = RUNTIME_DEFAULTS | settings(
            runtime_options, RUNTIME_DEFAULTS, "runtime options"
        )
        base_input = input_settings(input_options)
        base_output = settings(output_options, OUTPUT[out_type], "output options")
        result = ConversionBatchResult()
        with ExitStack() as stack:
            from meta_standards_converter.metadata.preparation_scope import preparation_session
            stack.enter_context(preparation_session())
            context = ExecutionContext(self, runtime, out_type, outfile, outdir, stack)
            declared = (
                manifest_specs(input_manifest, context.profile.max_xml_bytes)
                if input_manifest is not None
                else []
            )
            companion_paths = [
                v
                for spec in declared
                for value in spec.companions.values()
                for v in (value if isinstance(value, (list, tuple)) else [value])
                if isinstance(v, (str, Path)) and not is_url(str(v))
            ]
            declared_paths = [
                v
                for spec in declared
                for v in (
                    spec.sources
                    if isinstance(spec.sources, (list, tuple))
                    else [spec.sources]
                )
                if isinstance(v, (str, Path))
                and not is_url(str(v))
                and not str(v).lstrip().startswith(("{", "[", "<"))
            ]
            specs = declared + expand(
                input,
                runtime,
                output_directory=outdir,
                excluded=(
                    outfile,
                    input_manifest if isinstance(input_manifest, (str, Path)) else None,
                    *companion_paths,
                    *declared_paths,
                    *runtime["reserved_paths"],
                ),
            )
            # Manifest-owned paths must not be discovered again from the accompanying directory.
            owned = {
                str(Path(v).resolve())
                for s in declared
                for v in (
                    s.sources if isinstance(s.sources, (list, tuple)) else [s.sources]
                )
                if isinstance(v, (str, Path)) and not is_url(str(v))
            }
            specs = declared + [
                s
                for s in specs[len(declared) :]
                if not (
                    isinstance(s.sources, (str, Path))
                    and not is_url(str(s.sources))
                    and str(s.sources)[:1] not in "{[<"
                    and str(Path(s.sources).resolve()) in owned
                )
            ]
            if not specs:
                raise ValueError("At least one input or manifest entry is required")
            if outfile is not None and len(specs) > 1:
                raise ValueError("Batch inputs require outdir")
            planned = []
            ids = set()
            for index, spec in enumerate(specs, 1):
                if not isinstance(spec.companions, Mapping):
                    raise TypeError("companions must be a mapping")
                allowed_companions = {"idf", "sdrf", "features", "barcodes"}
                if set(spec.companions) - allowed_companions:
                    raise ValueError("Unknown companion role")
                chosen = kind_name(spec.in_type) or forced
                inp = merge_input_options(base_input, spec.input_options)
                out = settings(
                    base_output | output_settings(spec.output_options, out_type, partial=True),
                    OUTPUT[out_type],
                    "output options",
                )
                inp, out = apply_preparation_overrides(inp, out, spec.input_options, spec.output_options)
                inp, out, preparation = preparation_options(inp, out)
                error = None
                kind = chosen
                if kind is None:
                    try:
                        kind = self._detect(spec.sources, context)
                    except Exception as exc:
                        error = exc
                if kind in INPUT:
                    settings(inp, INPUT[kind], f"{kind} input options")
                    companion_roles = {
                        "magetab": {"idf", "sdrf"},
                        "matrix": {"features", "barcodes"},
                    }.get(kind, set())
                    if set(spec.companions) - companion_roles:
                        raise ValueError(f"Irrelevant companions for {kind}")
                    if spec.metadata is not None and kind not in {
                        "matrix",
                        "fastq",
                        "h5ad",
                        "anndata",
                    }:
                        raise ValueError(
                            f"Explicit metadata is not applicable to {kind}"
                        )
                if kind in {"h5ad", "anndata"} and out_type in {"h5ad", "obs"}:
                    settings(
                        out,
                        {"include_var", "include_uns"} if out_type == "obs" else set(),
                        "standalone expression options",
                    )
                if (
                    kind == "matrix"
                    and spec.metadata is None
                    and out_type in {"h5ad", "obs"}
                ):
                    settings(
                        out,
                        (
                            {"matrix_orientation", "include_var", "include_uns"}
                            if out_type == "obs"
                            else {"matrix_orientation"}
                        ),
                        "standalone matrix options",
                    )
                    if (
                        "orientation" in inp
                        and "matrix_orientation" in out
                        and inp["orientation"] != out["matrix_orientation"]
                    ):
                        raise ValueError("Conflicting matrix orientations")
                    if "matrix_orientation" in out:
                        inp["orientation"] = out["matrix_orientation"]
                item_id = spec.id or f"input_{index}"
                if not isinstance(item_id, str) or not item_id or item_id in ids:
                    raise ValueError("Input IDs must be unique nonempty strings")
                ids.add(item_id)
                planned.append((spec, kind, chosen, inp, out, error, item_id, preparation))
            seen = {}
            stopped = False
            for spec, kind, chosen, inp, out, error, item_id, preparation in planned:
                item = ConversionItemResult(
                    item_id,
                    source_label(spec.sources),
                    in_type=kind,
                    forced=bool(chosen),
                )
                result.items.append(item)
                sources = (
                    spec.sources
                    if isinstance(spec.sources, (list, tuple))
                    else [spec.sources]
                )
                item.origins = tuple(source_label(source) for source in sources)
                if stopped:
                    item.status = "skipped"
                    item.execution = "skipped"
                    item.diagnostics.append(
                        Diagnostic("fail_fast", "Not attempted after previous failure")
                    )
                    continue
                context.preparation_policy = preparation
                context.input_options = inp
                context.output_options = out
                context.stage = "input"
                try:
                    if error is not None:
                        raise error
                    if isinstance(spec.sources, (str, os.PathLike)) and not str(
                        spec.sources
                    ).lstrip().startswith(("{", "[", "<")):
                        key = (
                            (
                                str(Path(spec.sources).resolve())
                                if not is_url(str(spec.sources))
                                else str(spec.sources)
                            ),
                            kind,
                        )
                        configuration = repr((inp, out, preparation, spec.companions, spec.metadata))
                        if key in seen:
                            if seen[key] != configuration:
                                raise InputError(
                                    "conflicting_input",
                                    "The same source has conflicting settings or bindings",
                                )
                            item.status = "skipped"
                            item.execution = "skipped"
                            item.diagnostics.append(
                                Diagnostic(
                                    "duplicate_input",
                                    "Duplicate input was not converted twice",
                                    "discovery",
                                    "warning",
                                )
                            )
                            continue
                        seen[key] = configuration
                    if kind == "remote_xml":
                        kind = xml_kind(context.text(spec.sources))
                        item.in_type = kind
                    if kind not in self.handlers:
                        raise InputError(
                            "unsupported_input", "No registered input handler"
                        )
                    context.stage = "load"
                    self.handlers[kind].validate(spec, context)
                    from meta_standards_converter.metadata.preparation_scope import loading_source
                    with loading_source(preparation):
                        loaded = self.handlers[kind].load(spec, context)
                    sources = (
                        spec.sources
                        if isinstance(spec.sources, (list, tuple))
                        else [spec.sources]
                    )
                    loaded.origins = tuple(source_label(source) for source in sources)
                    loaded.companions = dict(spec.companions)
                    item.provider = loaded.provider
                    item.content_id = loaded.content_id
                    item.origins = loaded.origins
                    item.companions = {
                        key: (
                            [source_label(v) for v in value]
                            if isinstance(value, (list, tuple))
                            else source_label(value)
                        )
                        for key, value in spec.companions.items()
                    }
                    route, item.route = select_route(loaded, kind, out_type)
                    item.diagnostics.extend(loaded.diagnostics)
                    if loaded.metadata is not None:
                        context.stage = "validation"
                        issues = loaded.metadata.diagnostics
                        item.diagnostics.extend(
                            Diagnostic(
                                (
                                    "metadata_warning"
                                    if d.severity == "warning"
                                    else "invalid_metadata"
                                ),
                                f"{d.path}: {d.message}",
                                "validation",
                                d.severity,
                            )
                            for d in issues
                        )
                        if any(d.severity == "error" for d in issues):
                            raise InputError(
                                "invalid_metadata",
                                "Metadata failed structural or reference validation",
                            )
                        from .validation import validate_metadata
                        validate_metadata(loaded.metadata)
                        item.dataset_ids = tuple(
                            g.dataset_id for g in loaded.metadata.groups
                        )
                        item.diagnostics.extend(
                            Diagnostic("source_warning", w, "load", "warning")
                            for w in loaded.metadata.warnings
                            if not any(d.message == w for d in loaded.diagnostics)
                        )
                    context.stage = "preparation"
                    before = len(loaded.diagnostics)
                    loaded = context.preparer.prepare(loaded)
                    item.provider = loaded.provider
                    item.preparation = loaded.preparation
                    item.diagnostics.extend(loaded.diagnostics[before:])
                    if loaded.metadata is not None:
                        item.dataset_ids = tuple(g.dataset_id for g in loaded.metadata.groups)
                        validate_metadata(loaded.metadata)
                    context.stage = "conversion"
                    from meta_standards_converter.metadata.preparation_scope import exporting_prepared
                    with exporting_prepared():
                        self._export(loaded, context, item, route)
                    item.execution = "succeeded"
                    incomplete = any(
                        d.code
                        in {
                            "source_partial",
                            "source_failed",
                            "conversion_partial",
                            "sample_skipped",
                        }
                        for d in item.diagnostics
                    )
                    invalid = any(
                        d.code == "projection_invalid" for d in item.diagnostics
                    )
                    item.completeness = "partial" if incomplete else "complete"
                    item.validation = "invalid" if invalid else "valid"
                    item.status = "partial" if incomplete or invalid else "complete"
                except Exception as exc:
                    item.status = "failed"
                    item.execution = "failed"
                    item.completeness = "partial" if item.artifacts else "unknown"
                    if isinstance(exc, InputError):
                        code, message = exc.code, str(exc)
                    elif _missing_scientific_dependency(exc):
                        code, message = (
                            "missing_optional_dependency",
                            "Expression conversion requires optional dependencies; install "
                            "meta-standards-converter[h5ad] or use the full Docker image.",
                        )
                    elif type(exc).__name__ in {
                        "MINiMLModelError",
                        "MINiMLCompatibilityError",
                        "AtlasV1Error",
                    }:
                        code, message = "invalid_metadata", str(exc)
                    elif type(exc).__name__ in {
                        "TabularProjectionError",
                        "AnnDataProjectionError",
                        "MAGETabValidationError",
                    }:
                        code, message = "projection_invalid", str(exc)
                    else:
                        code, message = (
                            "conversion_failed",
                            f"Conversion failed ({type(exc).__name__})",
                        )
                    item.diagnostics.append(Diagnostic(code, message, context.stage))
                    if code in {
                        "type_mismatch",
                        "invalid_accession",
                        "invalid_metadata",
                        "projection_invalid",
                    }:
                        item.validation = "invalid"
                    stopped = runtime["fail_fast"]
        return result

    def _export(self, loaded, context, item, route):
        target = context.target
        options = context.output_options
        if route.service == "direct_magetab":
            publish_magetab(loaded.direct_magetab, context, item)
            return
        if route.service == "expression":
            from .expression import export_expression

            export_expression(loaded, context, item)
            return
        if loaded.metadata is None:
            raise InputError(
                "metadata_required",
                "Destination requires supplied or valid embedded MSC metadata",
            )
        groups = loaded.metadata.groups
        if route.service == "json":
            publish_jsons(groups, context, item)
            return
        if route.service == "magetab":
            from meta_standards_converter.converters.json2ae import JSON2AEConverter

            converter = context.service("json2ae", JSON2AEConverter)
            # Accession chains have already applied their native import enrichment.
            effective = {**options, "enrich": False}
            tables = converter.convert_loaded(loaded.metadata, **effective)
            publish_magetab(tables, context, item)
            return
        from meta_standards_converter.converters.json2tsv import JSON2TSVConverter

        converter = context.service("json2tsv", JSON2TSVConverter)
        aggregate = options.get("aggregate", False)
        selection = (
            [(item.id, loaded.metadata)]
            if aggregate
            else [
                (
                    g.dataset_id,
                    SourceLoadResult(
                        (g,), loaded.metadata.warnings, loaded.metadata.diagnostics
                    ),
                )
                for g in groups
            ]
        )
        tables = {}
        for name, source in selection:
            table = converter.project_loaded(
                source, **{k: v for k, v in options.items() if k != "aggregate"}
            )
            item.diagnostics.extend(
                Diagnostic("projection_warning", w, "projection", "warning")
                for w in table.warnings
            )
            item.diagnostics.extend(
                Diagnostic("projection_invalid", e, "projection") for e in table.errors
            )
            tables[name] = table
        publish_tables(tables, context, item)
