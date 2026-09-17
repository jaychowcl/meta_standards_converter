# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Expression input adapters with lazy scientific imports and explicit processing."""

from __future__ import annotations
import copy
import gzip
import json
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path

from .contracts import Diagnostic, InputError, LoadedInput
from .handlers import canonical_identity, safe_name, values
from .publication import reserve


@dataclass
class ExpressionInput:
    value: object
    kind: str
    name: str
    asset: object = None
    orientation: str = "auto"


def _metadata(value, context):
    from meta_standards_converter.sources.json import JSONPackageSource

    source = context.service("package_source", JSONPackageSource)
    if isinstance(value, (str, Path)):
        if isinstance(value, str) and value.lstrip().startswith(("{", "[")):
            return source.decode(json.loads(value))
        from .discovery import read_text

        return source.decode(
            json.loads(
                read_text(context.localize(value), context.profile.max_xml_bytes)
            )
        )
    return source.decode(value)


def load_expression(kind, spec, context):
    from meta_standards_converter.expression.assets import Asset
    from meta_standards_converter.sources.json import JSONPackageSource

    metadata = _metadata(spec.metadata, context) if spec.metadata is not None else None
    sources = values(spec.sources)
    if kind == "fastq":
        if (
            context.target in {"h5ad", "obs"}
            and not context.runtime["allow_processing"]
        ):
            raise InputError(
                "processing_required",
                "Raw reads require allow_processing=True and processing configuration",
            )
        if metadata is None:
            raise InputError(
                "metadata_required",
                "Raw reads require explicit sample/library metadata",
            )
        for read in sources:
            path = context.localize(read)
            opener = gzip.open if str(path).lower().endswith(".gz") else open
            with opener(path, "rt", encoding="ascii") as stream:
                record = [stream.readline(1024 * 1024).rstrip("\r\n") for _ in range(4)]
            header, sequence, separator, quality = record
            if (
                not header.startswith("@")
                or not separator.startswith("+")
                or not sequence
                or len(sequence) != len(quality)
                or any(len(line) >= 1024 * 1024 - 2 for line in record)
            ):
                raise InputError(
                    "invalid_fastq",
                    "Raw read input does not have a valid bounded FASTQ record",
                )
    if kind != "fastq" and len(sources) != 1:
        raise InputError(
            "ambiguous_bundle",
            "Expression input requires one primary asset and named companions",
        )
    source = sources[0]
    name = spec.id or (
        Path(str(source)).stem if isinstance(source, (str, Path)) else "input"
    )
    asset = None
    if kind in {"matrix", "fastq"}:
        if kind == "matrix":
            from .discovery import is_url, file_kind, tenx_member
            from meta_standards_converter.expression.readers import underlying_suffix

            if not is_url(str(source)):
                context.localize(source)
            if not Path(str(source)).is_dir() and file_kind(str(source)) != "matrix":
                raise InputError("type_mismatch", "Expected a supported matrix input")
            if underlying_suffix(str(source)) == ".mtx" and not all(
                spec.companions.get(k) for k in ("features", "barcodes")
            ):
                if is_url(str(source)):
                    raise InputError(
                        "missing_companion",
                        "Remote MEX requires features and barcodes companions",
                    )
                member = tenx_member(source)
                prefix = member[0] if member else ""
                matches = {"features": [], "barcodes": []}
                for candidate in Path(source).parent.iterdir():
                    found = tenx_member(candidate)
                    if found and found[0] == prefix and found[1] in matches:
                        matches[found[1]].append(candidate)
                if any(len(paths) > 1 for paths in matches.values()):
                    raise InputError(
                        "ambiguous_companion",
                        "Multiple feature or barcode files match the matrix",
                    )
                if any(not paths for paths in matches.values()):
                    raise InputError(
                        "missing_companion",
                        "MEX requires features and barcodes companions",
                    )
        sample_ids = (
            list(
                dict.fromkeys(
                    s.iid
                    for g in metadata.groups
                    for p in g.packages
                    for s in p.samples
                )
            )
            if metadata
            else []
        )
        if metadata and len(sample_ids) != 1:
            raise InputError(
                "sample_binding_required",
                "A standalone matrix/read bundle requires one explicitly bound metadata sample",
            )
        scope = sample_ids[0] if sample_ids else name
        source = str(source)
        asset = Asset(
            scope,
            source,
            "raw" if kind == "fastq" else "matrix",
            source="cli",
            features_path=spec.companions.get("features"),
            barcodes_path=spec.companions.get("barcodes"),
            orientation=context.input_options.get("orientation", "auto"),
            members=tuple({"uri": str(v)} for v in sources) if kind == "fastq" else (),
        )
        if kind == "fastq" or metadata is not None:
            return LoadedInput(
                metadata=metadata,
                assets=(asset,),
                content_id=canonical_identity(metadata),
            )
    if kind == "h5ad":
        source = context.localize(source)
        if str(source).lower().endswith(".gz"):
            if (
                Path(source).stat().st_size
                > context.profile.max_compressed_archive_bytes
            ):
                raise InputError(
                    "input_too_large", "Compressed H5AD exceeds input byte limit"
                )
            directory = context.stack.enter_context(
                tempfile.TemporaryDirectory(prefix="msc-h5ad-input-")
            )
            expanded = Path(directory) / "input.h5ad"
            total = 0
            with gzip.open(source, "rb") as stream, expanded.open("xb") as output:
                while block := stream.read(1024 * 1024):
                    total += len(block)
                    if total > context.profile.max_expanded_archive_bytes:
                        raise InputError(
                            "input_too_large", "Expanded H5AD exceeds input byte limit"
                        )
                    output.write(block)
            source = str(expanded)
        import anndata

        obj = anndata.read_h5ad(source, backed="r")
        try:
            if not obj.n_obs or not obj.n_vars:
                raise InputError("empty_matrix", "H5AD has an empty axis")
            embedded = obj.uns.get("msc_miniml")
            if metadata is None and embedded is not None:
                metadata = _embedded(embedded)
        finally:
            obj.file.close()
    elif kind == "anndata":
        if not source.n_obs or not source.n_vars:
            raise InputError("empty_matrix", "AnnData has an empty axis")
        embedded = source.uns.get("msc_miniml")
        if metadata is None and embedded is not None:
            metadata = _embedded(embedded)
    return LoadedInput(
        metadata=metadata,
        expression=ExpressionInput(
            source, kind, name, asset, context.input_options.get("orientation", "auto")
        ),
        content_id=canonical_identity(metadata) if metadata else None,
    )


def _embedded(value):
    from meta_standards_converter.sources.json import JSONPackageSource

    if (
        not isinstance(value, dict)
        or str(value.get("schema_version")) != "1.0"
        or not isinstance(value.get("packages_json"), str)
    ):
        raise InputError(
            "invalid_embedded_metadata", "Unsupported MSC MINiML transport in AnnData"
        )
    return JSONPackageSource().decode(json.loads(value["packages_json"]))


def _read_expression(expression, context, *, metadata_only=False):
    from meta_standards_converter.expression.memory import (
        _available_memory_bytes,
        _estimate_asset_memory_bytes,
    )
    from meta_standards_converter.expression.assets import Asset

    if expression.kind != "anndata" and not (
        metadata_only and expression.kind == "h5ad"
    ):
        path = context.localize(expression.value)
        estimator = context.service(
            "memory_estimator", lambda: _estimate_asset_memory_bytes
        )
        available = context.service("available_memory", lambda: _available_memory_bytes)
        estimate = estimator(path, expression.asset or Asset("input", path, "h5ad"))
        limit = min(
            context.profile.max_in_memory_matrix_bytes,
            int(available() * context.profile.available_memory_fraction),
        )
        if estimate <= 0 or limit <= 0 or estimate > limit:
            raise InputError(
                "memory_limit", "Expression input exceeds memory admission limits"
            )
    if expression.kind == "anndata":
        obj = expression.value
        if metadata_only:
            return obj, False
        return obj.to_memory() if getattr(obj, "isbacked", False) else obj.copy(), False
    if expression.kind == "h5ad":
        import anndata

        obj = anndata.read_h5ad(expression.value, backed="r" if metadata_only else None)
        return obj, metadata_only
    from meta_standards_converter.expression.readers import ProcessedAssetReader

    reader = context.service("reader", ProcessedAssetReader)
    from meta_standards_converter.expression.readers import underlying_suffix

    if expression.orientation == "auto" and underlying_suffix(
        expression.asset.path
    ) in {".csv", ".tsv", ".txt"}:
        raise InputError(
            "orientation_required",
            "Delimited matrices require an explicit matrix orientation",
        )
    obj = reader.read(
        expression.asset,
        orientation=expression.orientation,
        localize=lambda path, md5=None: context.localize(path),
    )
    return obj, False


def _publish_expression(obj, path, context):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".msc-h5ad-", dir=path.parent) as temporary:
        staged = Path(temporary) / path.name
        obj.write_h5ad(staged)
        with staged.open("rb") as stream:
            os.fsync(stream.fileno())
        if context.runtime["overwrite"]:
            os.replace(staged, path)
        else:
            os.link(staged, path)


def _standalone(loaded, context, item):
    expression = loaded.expression
    include_var = context.output_options.get("include_var", False)
    include_uns = context.output_options.get("include_uns", False)
    name = safe_name(expression.name) if context.outdir else expression.name
    obj, close = _read_expression(
        expression, context, metadata_only=context.target == "obs"
    )
    try:
        if context.target == "h5ad":
            if loaded.metadata is not None:
                obj.uns["msc_miniml"] = {
                    "schema_version": "1.0",
                    "packages_json": json.dumps(
                        [
                            p.to_mapping()
                            for g in loaded.metadata.groups
                            for p in g.packages
                        ],
                        sort_keys=True,
                    ),
                }
            item.payload = obj
            if context.destination:
                path = context.outfile or context.outdir / (name + ".h5ad")
                reserve([path], context)
                _publish_expression(obj, path, context)
                item.artifacts[path.name] = str(path)
            return
        payload = {"obs": obj.obs.copy()}
        if include_var:
            payload["var"] = obj.var.copy()
        if include_uns:
            payload["uns"] = copy.deepcopy(dict(obj.uns))
        item.payload = payload
        if not context.destination:
            return
        if context.outfile and (include_var or include_uns):
            raise InputError(
                "output_cardinality", "Observation sidecars require outdir"
            )
        if context.outfile:
            from .publication import atomic_write

            reserve([context.outfile], context)
            atomic_write(
                context.outfile,
                lambda stream: payload["obs"].rename_axis("cell_id").to_csv(stream),
                context.runtime["overwrite"],
            )
            item.artifacts[context.outfile.name] = str(context.outfile)
            return
        from meta_standards_converter.expression.components import (
            _publish_bundle,
            _json_value,
        )

        destinations = {
            "obs": context.outdir / (name + ".obs.csv"),
            "manifest": context.outdir / (name + ".json2obs.json"),
        }
        if include_var:
            destinations["var"] = context.outdir / (name + ".var.csv")
        if include_uns:
            destinations["uns"] = context.outdir / (name + ".uns.json")
        reserve(destinations.values(), context)
        with tempfile.TemporaryDirectory(prefix="msc-components-") as temporary:
            staged = {
                key: Path(temporary) / path.name for key, path in destinations.items()
            }
            payload["obs"].rename_axis("cell_id").to_csv(
                staged["obs"], lineterminator="\n"
            )
            if include_var:
                payload["var"].rename_axis("feature_id").to_csv(
                    staged["var"], lineterminator="\n"
                )
            if include_uns:
                staged["uns"].write_text(
                    json.dumps(
                        {
                            "schema_version": "1.0",
                            "values": _json_value(payload["uns"]),
                        },
                        ensure_ascii=False,
                    )
                    + "\n"
                )
            staged["manifest"].write_text(
                json.dumps(
                    {
                        "operation": "anndata_metadata",
                        "expression_integration": "none",
                        "artifacts": {k: str(v) for k, v in destinations.items()},
                        "rows": len(payload["obs"]),
                    },
                    sort_keys=True,
                )
                + "\n"
            )
            published = _publish_bundle(
                staged, destinations, overwrite=context.runtime["overwrite"]
            )
        item.artifacts.update({k: str(v) for k, v in destinations.items()})
        item.artifacts["bundle_pointer"] = str(published.pointer_path)
    finally:
        if close:
            obj.file.close()


def export_expression(loaded, context, item):
    if loaded.expression is not None:
        _standalone(loaded, context, item)
        return
    if loaded.metadata is None:
        raise InputError("metadata_required", "Expression processing requires metadata")
    if context.outdir is None:
        raise InputError(
            "destination_required", "Metadata-to-expression catalogues require outdir"
        )
    if len(loaded.metadata.groups) > 1:
        _export_groups(loaded, context, item)
        return
    from meta_standards_converter.converters.json2h5ad import JSON2H5ADConverter
    from meta_standards_converter.expression.assets import AssetManifest
    from meta_standards_converter.expression.components import AnnDataComponentExporter
    from meta_standards_converter.expression.catalogue import BatchConversionResult

    converter = context.service(
        "json2h5ad",
        lambda: JSON2H5ADConverter(
            retrieval_policy=context.policy, resource_profile=context.profile
        ),
    )
    options = {
        k: v
        for k, v in context.output_options.items()
        if k not in {"include_var", "include_uns"}
    }
    caller_assets = list(options.pop("explicit_assets", ()) or ())
    explicit = caller_assets + list(loaded.assets)
    manifest = AssetManifest()
    explicit.extend(
        manifest.parse_spec(v) for v in options.get("asset_specs", ()) or ()
    )
    if options.get("asset_manifest"):
        explicit.extend(manifest.load(options["asset_manifest"]))
    # Plan once for preflight; the underlying converter retains its own authoritative planning.
    destinations = []
    groups = loaded.metadata.groups
    for group in groups:
        safe_name(group.dataset_id)
        planned = converter.planner.plan(
            list(group.packages),
            explicit_assets=explicit,
            force_reprocess=options.get("force_reprocess", False),
        )
        if (
            any(a.kind == "raw" for a in planned.values())
            and not context.runtime["allow_processing"]
        ):
            raise InputError(
                "processing_required",
                "Selected raw assets require allow_processing=True",
            )
        if any(a.kind == "raw" for a in planned.values()):
            from meta_standards_converter.expression.references import ReferenceResolver

            resolver = (
                getattr(converter.pipeline_runner, "reference_resolver", None)
                or ReferenceResolver()
            )
            try:
                resolver.resolve(
                    [p.to_mapping() for p in group.packages],
                    **{k: options.get(k) for k in ("genome", "fasta", "gtf", "gff")},
                    accept_inferred=options.get("accept_inferred_reference", False),
                )
            except ValueError as exc:
                raise InputError("reference_required", str(exc)) from exc
        root = context.outdir / group.dataset_id if len(groups) > 1 else context.outdir
        if context.target == "h5ad":
            destinations.extend(
                root / (safe_name(sample) + ".h5ad") for sample in planned
            )
            destinations.append(root / (group.dataset_id + ".json2h5ad.json"))
        else:
            if context.output_options.get("include_var") and len(planned) != 1:
                raise InputError(
                    "output_cardinality", "A var export requires one sample per dataset"
                )
            destinations.extend(
                [
                    root / (group.dataset_id + ".obs.csv"),
                    root / (group.dataset_id + ".json2obs.json"),
                ]
            )
            if context.output_options.get("include_var"):
                destinations.append(root / (group.dataset_id + ".var.csv"))
            if context.output_options.get("include_uns"):
                destinations.append(root / (group.dataset_id + ".uns.json"))
    reserve(destinations, context)
    options["explicit_assets"] = caller_assets + list(loaded.assets)
    if context.target == "h5ad":
        assembly = context.outdir
    else:
        assembly = (
            Path(
                context.stack.enter_context(
                    tempfile.TemporaryDirectory(prefix="msc-obs-")
                )
            )
            / "assembly"
        )
    converted = converter.convert_loaded(
        loaded.metadata,
        out=str(assembly),
        source_json=loaded.source_path,
        overwrite=context.runtime["overwrite"],
        **options,
    )
    conversions = (
        converted.conversions
        if isinstance(converted, BatchConversionResult)
        else {converted.study_accession: converted}
    )
    for failure in getattr(converted, "failures", ()):
        item.diagnostics.append(
            Diagnostic("conversion_partial", str(failure), "expression", "warning")
        )
    outputs = {}
    for dataset, conversion in conversions.items():
        for failure in conversion.failures:
            item.diagnostics.append(
                Diagnostic("sample_skipped", str(failure), "expression", "warning")
            )
        for error in conversion.errors:
            item.diagnostics.append(
                Diagnostic("projection_invalid", str(error), "expression")
            )
        for warning in conversion.warnings:
            item.diagnostics.append(
                Diagnostic("projection_warning", str(warning), "expression", "warning")
            )
        if context.target == "h5ad":
            outputs[dataset] = conversion
            item.artifacts.update(
                {
                    f"{dataset}/{key}.h5ad": path
                    for key, path in conversion.sample_h5ads.items()
                }
            )
            if conversion.manifest_path:
                item.artifacts[dataset + "/manifest"] = conversion.manifest_path
            item.artifacts.update(
                {
                    dataset + "/retained/" + str(n): p
                    for n, p in enumerate(conversion.retained_h5ads)
                }
            )
        else:
            try:
                target = context.outdir / dataset if len(groups) > 1 else context.outdir
                exported = context.service(
                    "components", AnnDataComponentExporter
                ).export(
                    conversion,
                    target,
                    include_var=context.output_options.get("include_var", False),
                    include_uns=context.output_options.get("include_uns", False),
                    overwrite=context.runtime["overwrite"],
                )
                outputs[dataset] = exported
                item.artifacts.update(
                    {
                        dataset + "/" + k: str(v)
                        for k, v in exported.to_dict()["artifacts"].items()
                    }
                )
            except Exception as exc:
                item.diagnostics.append(
                    Diagnostic(
                        "conversion_partial",
                        f"Component export failed ({type(exc).__name__})",
                        "publication",
                        "warning",
                    )
                )
    if not outputs:
        raise InputError("conversion_failed", "No expression outputs were produced")
    item.payload = outputs


def _export_groups(loaded, context, item):
    """Keep study requirements, pipelines and publication independent within a source."""
    from dataclasses import replace
    from .contracts import ConversionItemResult

    outputs = {}
    for group in loaded.metadata.groups:
        child = copy.copy(context)
        outcome = ConversionItemResult(item.id, item.source)
        try:
            child.outdir = context.outdir / safe_name(group.dataset_id)
            child.destination = child.outdir
            source = replace(loaded, metadata=replace(loaded.metadata, groups=(group,)))
            export_expression(source, child, outcome)
            outputs.update(outcome.payload or {})
        except Exception as exc:
            message = (
                str(exc)
                if isinstance(exc, InputError)
                else f"Conversion failed ({type(exc).__name__})"
            )
            item.diagnostics.append(
                Diagnostic(
                    "conversion_partial",
                    f"{group.dataset_id}: {message}",
                    child.stage,
                    "warning",
                )
            )
            if isinstance(exc, InputError):
                item.diagnostics.append(
                    Diagnostic(exc.code, message, child.stage, "warning")
                )
            if context.runtime["fail_fast"]:
                break
        finally:
            item.diagnostics.extend(outcome.diagnostics)
            item.artifacts.update(outcome.artifacts)
    item.payload = outputs
    if not outputs:
        raise InputError("conversion_failed", "No dataset produced expression outputs")
