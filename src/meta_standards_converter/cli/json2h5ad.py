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

from meta_standards_converter.cli.common import (
    add_logging_arguments,
    configure_logging,
    record_safe_cli_error,
)
from meta_standards_converter.converters import (
    JSON2H5ADConverter,
    JSONDataOutputOrchestrator,
)
from meta_standards_converter.retrieval import RetrievalPolicy
from meta_standards_converter.runtime_contracts import get_resource_profile


logger = logging.getLogger(__name__)


def _resource_override(value: str) -> tuple[str, int | float]:
    name, separator, raw_value = value.partition("=")
    if not separator or not name.strip() or not raw_value.strip():
        raise argparse.ArgumentTypeError(
            "resource override must use FIELD=VALUE"
        )
    try:
        parsed: int | float
        parsed = (
            float(raw_value)
            if name.strip() == "disk_headroom_fraction"
            else int(raw_value)
        )
    except ValueError as error:
        raise argparse.ArgumentTypeError(
            "resource override VALUE must be numeric"
        ) from error
    return name.strip(), parsed


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
    workflow.add_argument(
        "--processed-checkpoint-dir",
        help="Persistent directory for resumable processed-sample checkpoints.",
    )
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
    resources = parser.add_argument_group("resource and retrieval policy")
    resources.add_argument(
        "--resource-profile",
        choices=("standard", "large"),
        default="standard",
        help="Typed resource envelope. Defaults to standard.",
    )
    resources.add_argument(
        "--resource-override",
        action="append",
        default=[],
        type=_resource_override,
        metavar="FIELD=VALUE",
        help="Override one typed profile field; repeat for multiple fields.",
    )
    resources.add_argument(
        "--asset-host",
        action="append",
        default=[],
        help="Explicitly allow one additional exact remote asset hostname.",
    )
    add_logging_arguments(parser)
    return parser


def main(argv=None) -> int:
    args = _parser().parse_args(argv)
    configure_logging(args, stream=sys.stderr)
    if (
        args.resource_profile != "standard"
        or args.resource_override
        or args.asset_host
    ):
        resource_profile = get_resource_profile(
            args.resource_profile,
            overrides=dict(args.resource_override),
        )
        retrieval_policy = RetrievalPolicy(
            resource_profile=resource_profile,
            allowed_hosts=frozenset(args.asset_host),
        )
        orchestrator = JSONDataOutputOrchestrator(
            h5ad_converter=JSON2H5ADConverter(
                retrieval_policy=retrieval_policy,
            )
        )
    else:
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
                ("processed_checkpoint_dir", args.processed_checkpoint_dir, None),
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
            safe_error = record_safe_cli_error(
                logger,
                error,
                location=json_path,
                stage="h5ad_conversion",
            )
            summaries.append(
                {
                    "source": safe_error.location,
                    "status": "failed",
                    "error": safe_error.to_dict(),
                }
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
