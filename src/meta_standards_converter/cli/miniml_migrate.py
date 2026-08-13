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


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Migrate MINiML 1.x or 2.0 JSON to MSC MINiML 3.0."
    )
    parser.add_argument("source", help="Legacy MINiML JSON file.")
    parser.add_argument("destination", help="Destination for MSC MINiML 3.0 JSON.")
    return parser


def main(argv=None) -> int:
    args = _parser().parse_args(argv)
    payload = json.loads(Path(args.source).read_text(encoding="utf-8"))
    values = payload if isinstance(payload, list) else [payload]
    results = [
        (
            MINiMLCodec.migrate_v2(value)
            if value.get("miniml_schema_version") == "2.0"
            else MINiMLCodec.migrate_v1(value)
        )
        for value in values
    ]
    encoded = [result.package.to_mapping() for result in results]
    output = encoded if isinstance(payload, list) else encoded[0]
    Path(args.destination).write_text(
        json.dumps(output, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
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
