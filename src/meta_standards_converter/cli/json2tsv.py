# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Command line interface for parsed MINiML or Atlas v2 JSON to TSV."""

from __future__ import annotations

import argparse
import json
import logging

from meta_standards_converter.cli.common import (
    add_logging_arguments,
    configure_logging,
)
from meta_standards_converter.converters import JSONDataOutputOrchestrator

logger = logging.getLogger(__name__)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Convert parsed MINiML or canonical Atlas v2 JSON files to a sample manifest."
    )
    parser.add_argument("json_path", nargs="+")
    parser.add_argument("--out", "--outdir", dest="outdir", default=".")
    parser.add_argument("--format", choices=("tsv", "csv"), default="tsv")
    parser.add_argument("--allow-invalid", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    add_logging_arguments(parser)
    return parser


def main(argv=None) -> int:
    args = _parser().parse_args(argv)
    configure_logging(args)
    orchestrator = JSONDataOutputOrchestrator()
    failed = False
    summaries = []
    for value in args.json_path:
        try:
            result = orchestrator.export_manifest(
                value,
                outdir=args.outdir,
                output_format=args.format,
                allow_invalid=args.allow_invalid,
                overwrite=args.overwrite,
            )
        except Exception:
            failed = True
            logger.exception("%s: TSV conversion failed", value)
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
