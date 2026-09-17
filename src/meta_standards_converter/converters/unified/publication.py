# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Publication adapters with exact paths and explicit atomic boundaries."""

import csv
import json
import os
import tempfile
from pathlib import Path
from uuid import uuid4
from .contracts import InputError
from .handlers import safe_name


def reserve(paths, context):
    context.stage = "publication"
    resolved = [Path(p).resolve() for p in paths]
    if len(set(resolved)) != len(resolved) or any(
        p in context.reserved for p in resolved
    ):
        raise InputError(
            "output_collision", "Multiple inputs would publish the same destination"
        )
    if not context.runtime["overwrite"] and any(p.exists() for p in resolved):
        raise InputError(
            "output_exists", "Destination already exists; overwrite is disabled"
        )
    context.reserved.update(resolved)


def atomic_write(path, writer, overwrite):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name("." + path.name + "." + uuid4().hex + ".tmp")
    try:
        with temporary.open("x", encoding="utf-8", newline="") as stream:
            writer(stream)
            stream.flush()
            os.fsync(stream.fileno())
        if overwrite:
            os.replace(temporary, path)
        else:
            os.link(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def publish_jsons(groups, context, item):
    from meta_standards_converter.miniml import MINiMLCodec

    aggregate = context.output_options.get("aggregate", False)
    packages = [p for g in groups for p in g.packages]
    item.payload = packages
    if not context.destination:
        return
    if context.outfile and len(groups) > 1 and not aggregate:
        raise InputError(
            "output_cardinality",
            "Multiple datasets need outdir or explicit aggregation",
        )
    if aggregate:
        paths = [
            (
                context.outfile or context.outdir / (safe_name(item.id) + ".json"),
                packages,
            )
        ]
    else:
        paths = [
            (
                context.outfile or context.outdir / (safe_name(g.dataset_id) + ".json"),
                g.packages,
            )
            for g in groups
        ]
    reserve([p for p, _ in paths], context)
    for path, values in paths:
        from meta_standards_converter.converters.archive_results import publish_json

        publish_json(
            path, MINiMLCodec().encode_many(values), context.runtime["overwrite"]
        )
        item.artifacts[str(path.name)] = str(path)


def publish_tables(tables, context, item):
    item.payload = tables
    if not context.destination:
        return
    if context.outfile and len(tables) > 1:
        raise InputError(
            "output_cardinality", "Multiple tables require outdir or aggregation"
        )
    paths = [
        (
            context.outfile or context.outdir / (safe_name(key) + "." + context.target),
            table,
        )
        for key, table in tables.items()
    ]
    reserve([p for p, _ in paths], context)
    for path, table in paths:

        def write(stream):
            writer = csv.DictWriter(
                stream,
                fieldnames=table.columns,
                delimiter="\t" if context.target == "tsv" else ",",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(table.rows)

        atomic_write(path, write, context.runtime["overwrite"])
        item.artifacts[path.name] = str(path)


def publish_magetab(tables, context, item):
    from meta_standards_converter.magetab.writer import MAGETabWriter
    from meta_standards_converter.expression.components import _publish_bundle

    writer = MAGETabWriter()
    names = [safe_name(writer._magetab_accession(t)) for t in tables]
    item.dataset_ids = tuple(names)
    item.payload = tables
    if not context.destination:
        return
    if context.outfile:
        raise InputError(
            "output_cardinality", "MAGE-TAB requires outdir for its IDF/SDRF pair"
        )
    destinations = [
        (context.outdir / (name + ".idf.txt"), context.outdir / (name + ".sdrf.txt"))
        for name in names
    ]
    reserve([p for pair in destinations for p in pair], context)
    for table, (idf, sdrf) in zip(tables, destinations):
        with tempfile.TemporaryDirectory(prefix="msc-magetab-") as temporary:
            writer.write(table, temporary)
            published = _publish_bundle(
                {
                    "idf": Path(temporary) / idf.name,
                    "sdrf": Path(temporary) / sdrf.name,
                },
                {"idf": idf, "sdrf": sdrf},
                overwrite=context.runtime["overwrite"],
            )
        item.artifacts[idf.name] = str(idf)
        item.artifacts[sdrf.name] = str(sdrf)
        item.artifacts[idf.stem + ".bundle_pointer"] = str(published.pointer_path)
