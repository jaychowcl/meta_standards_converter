# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Command-line interface for study-scoped ENA conversion."""

from __future__ import annotations

import argparse
import logging

from meta_standards_converter.cli.common import (
    add_logging_arguments,
    add_resource_profile_arguments,
    configured_resource_profile,
    configure_logging,
    record_safe_cli_error,
)
from meta_standards_converter.converters.ena2json import ena2json


logger = logging.getLogger(__name__)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Convert INSDC accessions through ENA to study-scoped MINiML JSON."
    )
    parser.add_argument("accession", nargs="+", help="INSDC project, study, sample, experiment, or run accession(s).")
    enrichment = parser.add_mutually_exclusive_group()
    enrichment.add_argument("--enrich-geo", action="store_true", help="Explicitly merge a linked GEO study.")
    enrichment.add_argument("--enrich-ae", action="store_true", help="Explicitly merge a linked ArrayExpress study.")
    parser.add_argument("--out", default=".", help="Directory for generated JSON files.")
    add_resource_profile_arguments(parser.add_argument_group("resource policy"))
    add_logging_arguments(parser)
    return parser


def main(argv=None) -> int:
    parser = _parser()
    args = parser.parse_args(argv)
    configure_logging(args)
    converter = ena2json(resource_profile=configured_resource_profile(args, parser))
    failed = False
    for accession in args.accession:
        try:
            converter.convert(
                accession,
                enrich_geo=args.enrich_geo,
                enrich_ae=args.enrich_ae,
                out=args.out,
            )
        except Exception as error:
            failed = True
            record_safe_cli_error(
                logger,
                error,
                location=accession,
                stage="ena_json_conversion",
                provider="ena",
            )
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
