# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""
Shared command-line helpers.
"""

import argparse
import logging
import sys


def add_logging_arguments(parser: argparse.ArgumentParser) -> None:
    verbosity_group = parser.add_mutually_exclusive_group()
    verbosity_group.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Increase logging verbosity. Use -v for INFO and -vv for DEBUG.",
    )
    verbosity_group.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Only emit ERROR logs.",
    )
    parser.add_argument(
        "--log-file",
        help="Optional file path to write logs.",
    )


def log_level(args) -> int:
    if args.quiet:
        return logging.ERROR
    if args.verbose >= 2:
        return logging.DEBUG
    if args.verbose == 1:
        return logging.INFO
    return logging.WARNING


def configure_logging(args) -> None:
    level = log_level(args)
    package_logger = logging.getLogger("meta_standards_converter")
    for handler in package_logger.handlers[:]:
        package_logger.removeHandler(handler)
        handler.close()
    package_logger.setLevel(logging.DEBUG)
    package_logger.propagate = False

    stdout_handler = logging.StreamHandler(sys.stdout)
    stdout_handler.setLevel(level)
    stdout_handler.setFormatter(logging.Formatter("%(levelname)s %(name)s: %(message)s"))
    package_logger.addHandler(stdout_handler)

    if args.log_file:
        file_handler = logging.FileHandler(args.log_file, mode="w")
        file_handler.setLevel(level)
        file_handler.setFormatter(
            logging.Formatter("%(asctime)s %(levelname)s %(name)s: %(message)s")
        )
        package_logger.addHandler(file_handler)
