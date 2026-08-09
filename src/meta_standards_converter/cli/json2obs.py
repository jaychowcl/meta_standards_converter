# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Command line interface for JSON plus expression assets to AnnData metadata."""

from __future__ import annotations

import argparse
import json
import logging
import sys

from meta_standards_converter.cli.common import (
    add_logging_arguments,
    configure_logging,
    record_safe_cli_error,
)
from meta_standards_converter.converters import JSONDataOutputOrchestrator

logger = logging.getLogger(__name__)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Aggregate AnnData observation metadata without combining expression matrices."
    )
    parser.add_argument("json_path", nargs="+")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--include-var", action="store_true")
    parser.add_argument("--include-uns", action="store_true")
    parser.add_argument("--asset-manifest")
    parser.add_argument("--asset", action="append", default=[])
    parser.add_argument("--force-reprocess", action="store_true")
    parser.add_argument("--pipeline", choices=("auto", "scrnaseq", "rnaseq"), default="auto")
    parser.add_argument("--genome")
    parser.add_argument("--fasta")
    parser.add_argument("--gtf")
    parser.add_argument("--gff")
    parser.add_argument("--accept-inferred-reference", action="store_true")
    parser.add_argument("--profile", default="docker")
    parser.add_argument("--revision")
    parser.add_argument("--params-file")
    parser.add_argument("--nextflow-config")
    parser.add_argument("--work-dir")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--processed-checkpoint-dir")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--allow-invalid", action="store_true")
    parser.add_argument(
        "--use-harmonization-overrides",
        action="store_true",
        help="Apply an Agentic Curator harmonization override profile.",
    )
    parser.add_argument(
        "--matrix-orientation",
        choices=("auto", "genes-by-observations", "observations-by-genes"),
        default="auto",
    )
    add_logging_arguments(parser)
    return parser


def main(argv=None) -> int:
    args = _parser().parse_args(argv)
    configure_logging(args, stream=sys.stderr)
    orchestrator = JSONDataOutputOrchestrator()
    summaries = []
    failed = False
    for source in args.json_path:
        try:
            convert_options = dict(
                outdir=args.outdir,
                include_var=args.include_var,
                include_uns=args.include_uns,
                asset_manifest=args.asset_manifest,
                asset_specs=args.asset,
                force_reprocess=args.force_reprocess,
                pipeline=args.pipeline,
                genome=args.genome,
                fasta=args.fasta,
                gtf=args.gtf,
                gff=args.gff,
                accept_inferred_reference=args.accept_inferred_reference,
                profile=args.profile,
                revision=args.revision,
                params_file=args.params_file,
                nextflow_config=args.nextflow_config,
                work_dir=args.work_dir,
                resume=args.resume,
                processed_checkpoint_dir=args.processed_checkpoint_dir,
                overwrite=args.overwrite,
                allow_invalid=args.allow_invalid,
                matrix_orientation=args.matrix_orientation,
            )
            if args.use_harmonization_overrides:
                convert_options["use_harmonization_overrides"] = True
            result = orchestrator.export_anndata_metadata(source, **convert_options)
        except Exception as error:
            failed = True
            safe_error = record_safe_cli_error(
                logger,
                error,
                location=source,
                stage="observation_export",
            )
            summaries.append(
                {
                    "source": safe_error.location,
                    "status": "failed",
                    "error": safe_error.to_dict(),
                }
            )
            continue
        failed = failed or result.partial
        summaries.append(result.to_dict())
    print(
        json.dumps(
            {
                "operation": "anndata_metadata",
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
