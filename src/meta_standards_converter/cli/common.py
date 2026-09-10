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

from meta_standards_converter.magetab.constructor import PLATFORM_HANDLER_KEYS
from meta_standards_converter.runtime_contracts import (
    ResourceProfile,
    SafeErrorEnvelope,
    get_resource_profile,
)


def parse_resource_override(value: str) -> tuple[str, int | float]:
    name, separator, raw_value = value.partition("=")
    if not separator or not name.strip() or not raw_value.strip():
        raise argparse.ArgumentTypeError(
            "resource override must use FIELD=VALUE"
        )
    try:
        parsed: int | float = (
            float(raw_value)
            if name.strip() == "disk_headroom_fraction"
            else int(raw_value)
        )
    except ValueError as error:
        raise argparse.ArgumentTypeError(
            "resource override VALUE must be numeric"
        ) from error
    return name.strip(), parsed


def add_resource_profile_arguments(parser) -> None:
    parser.add_argument(
        "--resource-profile",
        choices=("standard", "large"),
        default="standard",
        help="Typed resource envelope. Defaults to standard.",
    )
    parser.add_argument(
        "--resource-override",
        action="append",
        default=[],
        type=parse_resource_override,
        metavar="FIELD=VALUE",
        help="Override one typed profile field; repeat for multiple fields.",
    )


def configured_resource_profile(args, parser) -> ResourceProfile:
    try:
        return get_resource_profile(
            args.resource_profile,
            overrides=dict(args.resource_override),
        )
    except ValueError as error:
        parser.error(str(error))


def add_platform_handler_arguments(parser: argparse.ArgumentParser) -> None:
    handler_group = parser.add_mutually_exclusive_group()
    handler_group.add_argument(
        "--platform-handler",
        choices=PLATFORM_HANDLER_KEYS,
        help="Force IDF and SDRF generation through the selected platform handler.",
    )
    handler_group.add_argument(
        "--list-platform-handlers",
        action="store_true",
        help="List available platform handler keys and exit.",
    )


def print_platform_handlers() -> None:
    for name in PLATFORM_HANDLER_KEYS:
        print(name)


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


def configure_logging(args, *, stream=None) -> None:
    level = log_level(args)
    package_logger = logging.getLogger("meta_standards_converter")
    for handler in package_logger.handlers[:]:
        package_logger.removeHandler(handler)
        handler.close()
    package_logger.setLevel(logging.DEBUG)
    package_logger.propagate = False

    handler = logging.StreamHandler(sys.stdout if stream is None else stream)
    handler.setLevel(level)
    handler.setFormatter(logging.Formatter("%(levelname)s %(name)s: %(message)s"))
    package_logger.addHandler(handler)

    if args.log_file:
        file_handler = logging.FileHandler(args.log_file, mode="w")
        file_handler.setLevel(level)
        file_handler.setFormatter(
            logging.Formatter("%(asctime)s %(levelname)s %(name)s: %(message)s")
        )
        package_logger.addHandler(file_handler)


def record_safe_cli_error(
    logger: logging.Logger,
    error: BaseException,
    *,
    location: str,
    stage: str,
    provider: str | None = None,
) -> SafeErrorEnvelope:
    """Log and return a durable error without serializing its raw message."""

    envelope = SafeErrorEnvelope.from_exception(
        error,
        provider=provider,
        location=location,
        stage=stage,
    )
    logger.error(
        "%s: %s failed error_type=%s correlation_id=%s",
        envelope.location or "input",
        stage,
        envelope.error_type,
        envelope.correlation_id,
    )
    return envelope
