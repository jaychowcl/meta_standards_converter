# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Command line interface for parsed MINiML or Atlas v1 JSON to TSV."""

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
from meta_standards_converter.converters import JSON2TSVConverter

logger = logging.getLogger(__name__)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Convert parsed MINiML or canonical Atlas v1 JSON files to a sample manifest."
    )
    parser.add_argument("json_path", nargs="+")
    parser.add_argument("--out", "--outdir", dest="outdir", default=".")
    parser.add_argument("--format", choices=("tsv", "csv"), default="tsv")
    parser.add_argument("--allow-invalid", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
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
    orchestrator = JSON2TSVConverter()
    failed = False
    summaries = []
    for value in args.json_path:
        try:
            export_options = dict(
                outdir=args.outdir,
                output_format=args.format,
                allow_invalid=args.allow_invalid,
                overwrite=args.overwrite,
            )
            if args.use_harmonization_overrides:
                export_options["use_harmonization_overrides"] = True
            result = orchestrator.export_manifest(value, **export_options)
        except Exception as error:
            failed = True
            safe_error = record_safe_cli_error(
                logger,
                error,
                location=value,
                stage="manifest_export",
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
                "operation": "manifest",
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
