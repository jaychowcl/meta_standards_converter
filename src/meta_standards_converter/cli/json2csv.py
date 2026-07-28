# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Command line interface for parsed or Atlas JSON to CSV conversion."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from meta_standards_converter.cli.common import (
    add_logging_arguments,
    configure_logging,
)
from meta_standards_converter.converters.json2tabular import JSON2CSVConverter

logger = logging.getLogger(__name__)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Convert parsed MINiML or Atlas JSON files to CSV."
    )
    parser.add_argument("json_path", nargs="+")
    parser.add_argument("--out", default=".")
    parser.add_argument("--allow-invalid", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    add_logging_arguments(parser)
    return parser


def main(argv=None) -> int:
    args = _parser().parse_args(argv)
    configure_logging(args)
    converter = JSON2CSVConverter()
    failed = False
    output_dir = Path(args.out)
    for value in args.json_path:
        source = Path(value)
        try:
            result = converter.convert_source(
                source,
                output_dir / f"{source.stem}.csv",
                allow_invalid=args.allow_invalid,
                overwrite=args.overwrite,
            )
        except Exception:
            failed = True
            logger.exception("%s: CSV conversion failed", value)
            continue
        failed = failed or result.partial
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
