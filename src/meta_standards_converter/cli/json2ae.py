# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Command line interface for parsed MINiML JSON to MAGE-TAB conversion."""

import argparse
import logging

from meta_standards_converter.cli.common import (
    add_logging_arguments,
    add_platform_handler_arguments,
    configure_logging,
    print_platform_handlers,
)
from meta_standards_converter.converters.json2ae import json2ae


logger = logging.getLogger(__name__)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Convert one or more parsed MINiML JSON files to ArrayExpress MAGE-TAB files."
    )
    parser.add_argument(
        "json_path",
        nargs="*",
        help="Parsed MINiML JSON file(s), for example GSE234602.json.",
    )
    parser.add_argument(
        "--no-enrich",
        action="store_false",
        dest="enrich",
        default=True,
        help="Skip PubMed/SRA enrichment and convert the JSON exactly as supplied.",
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
    if not args.json_path:
        parser.error("the following arguments are required: json_path")
    configure_logging(args)
    converter = json2ae()
    failed = False
    logger.debug(
        "Starting json2ae CLI with %d JSON file(s), enrich=%s, out=%s",
        len(args.json_path),
        args.enrich,
        args.out,
    )

    for json_path in args.json_path:
        logger.info("%s: conversion started", json_path)
        try:
            convert_options = dict(
                json_path=json_path,
                enrich=args.enrich,
                out=args.out,
            )
            if args.platform_handler:
                convert_options["platform_handler"] = args.platform_handler
            magetabs = converter.convert(**convert_options)
        except Exception:
            failed = True
            logger.exception("%s: conversion failed", json_path)
            continue

        logger.info(
            "%s: converted %d MAGE-TAB output(s) to %s",
            json_path,
            len(magetabs),
            args.out,
        )

    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
