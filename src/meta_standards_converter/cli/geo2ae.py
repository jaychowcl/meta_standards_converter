# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""
Command line interface for GEO to ArrayExpress MAGE-TAB conversion.
"""

import argparse
import logging

from meta_standards_converter.cli.common import (
    add_logging_arguments,
    add_platform_handler_arguments,
    configure_logging,
    print_platform_handlers,
)
from meta_standards_converter.converters.geo2ae import geo2ae

logger = logging.getLogger(__name__)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Convert one or more GEO Series accessions to ArrayExpress MAGE-TAB files."
    )
    parser.add_argument(
        "gse",
        nargs="*",
        help="GEO Series accession(s), for example GSE234602.",
    )
    parser.add_argument(
        "--related",
        "--related-series",
        "--get-related-series",
        action="store_true",
        dest="related_series",
        help="Include related GEO super/subseries where available.",
    )
    cleanup_group = parser.add_mutually_exclusive_group()
    cleanup_group.add_argument(
        "--remove-empty",
        action="store_true",
        dest="remove_empty",
        default=True,
        help="Remove empty parsed MINiML fields before conversion. This is the default.",
    )
    cleanup_group.add_argument(
        "--keep-empty",
        action="store_false",
        dest="remove_empty",
        help="Preserve empty parsed MINiML fields before conversion.",
    )
    parser.add_argument(
        "--out",
        default=".",
        help="Directory for generated IDF and SDRF files. Defaults to the current directory.",
    )
    add_platform_handler_arguments(parser)
    add_logging_arguments(parser)
    return parser


def main(argv=None) -> int:
    parser = _parser()
    args = parser.parse_args(argv)
    if args.list_platform_handlers:
        print_platform_handlers()
        return 0
    if not args.gse:
        parser.error("the following arguments are required: gse")
    configure_logging(args)
    converter = geo2ae()
    failed = False
    logger.debug(
        "Starting geo2ae CLI with %d accession(s), related_series=%s, remove_empty=%s, out=%s",
        len(args.gse),
        args.related_series,
        args.remove_empty,
        args.out,
    )

    for gse in args.gse:
        logger.info("%s: conversion started", gse)
        try:
            convert_options = dict(
                gse=gse,
                related_series=args.related_series,
                remove_empty=args.remove_empty,
                out=args.out,
            )
            if args.platform_handler:
                convert_options["platform_handler"] = args.platform_handler
            magetabs = converter.convert(**convert_options)
        except Exception:
            failed = True
            logger.exception("%s: conversion failed", gse)
            continue

        logger.info("%s: converted %d MAGE-TAB output(s) to %s", gse, len(magetabs), args.out)

    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
