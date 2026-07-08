# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""
Command line interface for parsed MINiML JSON to H5AD conversion.
"""

import argparse
import logging

from meta_standards_converter.cli.common import add_logging_arguments, configure_logging
from meta_standards_converter.converters.json2h5ad import json2h5ad


logger = logging.getLogger(__name__)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Convert one or more parsed MINiML JSON files to H5AD files."
    )
    parser.add_argument(
        "json_path",
        nargs="+",
        help="Parsed MINiML JSON file(s), for example GSE234602.json.",
    )
    parser.add_argument(
        "--out",
        default=".",
        help="Directory for generated H5AD files. Defaults to the current directory.",
    )
    add_logging_arguments(parser)
    return parser


def main(argv=None) -> int:
    args = _parser().parse_args(argv)
    configure_logging(args)
    converter = json2h5ad()
    failed = False
    logger.debug(
        "Starting json2h5ad CLI with %d JSON file(s), out=%s",
        len(args.json_path),
        args.out,
    )

    for json_path in args.json_path:
        logger.info("%s: H5AD conversion started", json_path)
        try:
            h5ad_path = converter.convert(json_path=json_path, out=args.out)
        except Exception:
            failed = True
            logger.exception("%s: H5AD conversion failed", json_path)
            continue

        logger.info("%s: converted to %s", json_path, h5ad_path)

    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
