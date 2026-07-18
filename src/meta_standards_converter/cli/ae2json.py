# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Command line interface for MAGE-TAB to parsed JSON conversion."""

import argparse
import logging

from meta_standards_converter.cli.common import add_logging_arguments, configure_logging
from meta_standards_converter.converters.ae2json import ae2json


logger = logging.getLogger(__name__)


def _parser():
    parser = argparse.ArgumentParser(
        description="Convert local, remote, or BioStudies MAGE-TAB metadata to parsed JSON."
    )
    parser.add_argument(
        "source",
        nargs="+",
        help="IDF path, HTTP(S) IDF URL, or ArrayExpress/BioStudies accession.",
    )
    parser.add_argument(
        "--sdrf",
        action="append",
        default=None,
        help="Explicit SDRF path or HTTP(S) URL. Repeat for multiple SDRFs; requires one source.",
    )
    parser.add_argument("--out", default=".", help="Directory for generated JSON files.")
    add_logging_arguments(parser)
    return parser


def main(argv=None):
    parser = _parser()
    args = parser.parse_args(argv)
    if args.sdrf and len(args.source) != 1:
        parser.error("--sdrf overrides require exactly one source")
    configure_logging(args)
    converter = ae2json()
    failed = False
    logger.debug("Starting ae2json CLI with %d source(s), out=%s", len(args.source), args.out)
    for source in args.source:
        logger.info("%s: conversion started", source)
        try:
            packages = converter.convert(
                source=source,
                out=args.out,
                sdrf_sources=args.sdrf,
            )
        except Exception:
            failed = True
            logger.exception("%s: conversion failed", source)
            continue
        logger.info("%s: converted %d package(s) to %s", source, len(packages), args.out)
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
