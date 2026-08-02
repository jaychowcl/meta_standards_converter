# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""
Command line interface for parsed MINiML or canonical Atlas v1 JSON to H5AD.
"""

import argparse
import json
import logging
import sys

from meta_standards_converter.cli.common import add_logging_arguments, configure_logging
from meta_standards_converter.converters import JSONDataOutputOrchestrator


logger = logging.getLogger(__name__)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Convert parsed MINiML or canonical Atlas v1 JSON files to H5AD."
        )
    )
    parser.add_argument(
        "json_path",
        nargs="+",
        help=(
            "Parsed MINiML or canonical Atlas v1 JSON file(s), "
            "for example GSE234602.json."
        ),
    )
    parser.add_argument(
        "--out",
        "--outdir",
        dest="outdir",
        default=".",
        help="Directory for generated H5AD files. Defaults to the current directory.",
    )
    parser.add_argument("--asset-manifest", help="CSV/TSV mapping GEO accessions to local or remote assets.")
    parser.add_argument(
        "--asset",
        action="append",
        default=[],
        metavar="ACCESSION=PATH",
        help="Explicit H5AD, matrix, or FASTQ asset. Repeat for multiple assets.",
    )
    parser.add_argument(
        "--force-reprocess",
        action="store_true",
        help="Ignore processed assets and rebuild every eligible sample from raw FASTQs.",
    )
    parser.add_argument(
        "--pipeline",
        choices=("auto", "scrnaseq", "rnaseq"),
        default="auto",
        help="nf-core pipeline for raw inputs. Defaults to metadata-based selection.",
    )
    reference = parser.add_argument_group("reference")
    reference.add_argument("--genome", help="nf-core genome key, for example GRCh38.")
    reference.add_argument("--fasta", help="Custom reference genome FASTA path.")
    reference.add_argument("--gtf", help="Custom reference annotation GTF path.")
    reference.add_argument("--gff", help="Custom reference annotation GFF3 path.")
    reference.add_argument(
        "--accept-inferred-reference",
        action="store_true",
        help="Allow a supported reference inferred from GEO organism metadata.",
    )
    workflow = parser.add_argument_group("nf-core execution")
    workflow.add_argument("--profile", default="docker", help="Nextflow profile. Defaults to docker.")
    workflow.add_argument("--revision", help="Override the pinned nf-core pipeline revision.")
    workflow.add_argument("--params-file", help="Additional nf-core JSON parameters.")
    workflow.add_argument("--nextflow-config", help="Nextflow resource/infrastructure config path.")
    workflow.add_argument("--work-dir", help="Nextflow work directory.")
    workflow.add_argument("--resume", action="store_true", help="Resume from the Nextflow cache.")
    parser.add_argument("--overwrite", action="store_true", help="Replace existing normalized outputs.")
    parser.add_argument(
        "--allow-invalid",
        action="store_true",
        help="Write outputs carrying projector-reported validation errors.",
    )
    parser.add_argument(
        "--matrix-orientation",
        choices=("auto", "genes-by-observations", "observations-by-genes"),
        default="auto",
        help="Orientation for generic delimited matrices.",
    )
    parser.add_argument(
        "--use-harmonization-overrides",
        action="store_true",
        help="Apply an Agentic Curator harmonization override profile.",
    )
    add_logging_arguments(parser)
    return parser


def main(argv=None) -> int:
    args = _parser().parse_args(argv)
    configure_logging(args, stream=sys.stderr)
    orchestrator = JSONDataOutputOrchestrator()
    failed = False
    summaries = []
    logger.debug(
        "Starting json2h5ad CLI with %d JSON file(s), out=%s",
        len(args.json_path),
        args.outdir,
    )

    for json_path in args.json_path:
        logger.info("%s: H5AD conversion started", json_path)
        try:
            convert_options = {}
            for name, value, default in (
                ("asset_manifest", args.asset_manifest, None),
                ("asset_specs", args.asset, []),
                ("force_reprocess", args.force_reprocess, False),
                ("pipeline", args.pipeline, "auto"),
                ("genome", args.genome, None),
                ("fasta", args.fasta, None),
                ("gtf", args.gtf, None),
                ("gff", args.gff, None),
                ("accept_inferred_reference", args.accept_inferred_reference, False),
                ("profile", args.profile, "docker"),
                ("revision", args.revision, None),
                ("params_file", args.params_file, None),
                ("nextflow_config", args.nextflow_config, None),
                ("work_dir", args.work_dir, None),
                ("resume", args.resume, False),
                ("overwrite", args.overwrite, False),
                ("allow_invalid", args.allow_invalid, False),
                ("matrix_orientation", args.matrix_orientation, "auto"),
            ):
                if value != default:
                    convert_options[name] = value
            if args.use_harmonization_overrides:
                convert_options["use_harmonization_overrides"] = True
            conversion = orchestrator.export_h5ad(
                json_path, outdir=args.outdir, **convert_options
            )
        except Exception as error:
            failed = True
            logger.exception("%s: H5AD conversion failed", json_path)
            summaries.append(
                {"source": json_path, "status": "failed", "error": str(error)}
            )
            continue

        if getattr(conversion, "partial", False):
            failed = True
            logger.error("%s: conversion completed partially: %s", json_path, conversion)
        else:
            logger.info("%s: converted to %s", json_path, conversion)

        summaries.append(conversion.to_dict())

    print(
        json.dumps(
            {
                "operation": "h5ad",
                "status": "partial" if failed else "complete",
                "datasets": summaries,
            },
            ensure_ascii=False,
            sort_keys=True,
        )
    )

    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
