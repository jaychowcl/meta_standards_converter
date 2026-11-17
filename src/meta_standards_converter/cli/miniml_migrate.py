# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Explicit legacy MSC MINiML to canonical 3.0 migration command."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from meta_standards_converter.miniml import MINiMLCodec
from meta_standards_converter.conversion_logging import operation, stage
from .common import add_logging_arguments, configure_logging
from .logging_output import cli_entrypoint


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Import legacy unversioned/1.0 source JSON directly as MSC MINiML 3.0. MINiML 2.0 is unsupported."
    )
    parser.add_argument("source", help="Legacy MINiML JSON file.")
    parser.add_argument("destination", help="Destination for MSC MINiML 3.0 JSON.")
    add_logging_arguments(parser)
    return parser


@cli_entrypoint
def main(argv=None) -> int:
    args = _parser().parse_args(argv)
    configure_logging(args)
    with operation('MINiMLCodec.migrate_v1') as run:
        with stage('load'):
            payload = json.loads(Path(args.source).read_text(encoding='utf-8'))
        with stage('validation'):
            values = payload if isinstance(payload, list) else [payload]
            results = [MINiMLCodec.migrate_v1(value) for value in values]
            encoded = [result.package.to_mapping() for result in results]
            output = encoded if isinstance(payload, list) else encoded[0]
        with stage('publication'):
            from meta_standards_converter.converters.archive_results import publish_json
            publish_json(args.destination, output, overwrite=True)
        run.statistics.update(packages_migrated=len(results), artifact_paths={'json': args.destination})
    print(json.dumps({
        "packages_migrated": len(results),
        "diagnostics": [
            {
                "path": item.path,
                "code": item.code,
                "message": item.message,
                "severity": item.severity,
            }
            for result in results
            for item in result.diagnostics
        ],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
