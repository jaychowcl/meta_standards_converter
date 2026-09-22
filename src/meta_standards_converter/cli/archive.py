# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Argument/report handling shared by the independent native commands."""
import argparse
import logging
from pathlib import Path
from .common import add_logging_arguments, add_resource_profile_arguments, configured_resource_profile, configure_logging
from ..converters.archive_results import publish_json


def parser_for(provider):
    parser = argparse.ArgumentParser(description=f'Convert {provider} accessions to study-scoped MSC MINiML JSON.')
    parser.add_argument('accession', nargs='+')
    parser.add_argument('--out', default='.')
    parser.add_argument('--enrich-from-geo-ae', action='store_true')
    parser.add_argument('--include-peer', action='store_true')
    parser.add_argument('--report', metavar='PATH')
    parser.add_argument('--evidence-dir', metavar='DIR')
    parser.add_argument('--overwrite', action='store_true')
    add_logging_arguments(parser)
    add_resource_profile_arguments(parser.add_argument_group('resource policy'))
    return parser


def run_cli(parser, converter_type, argv):
    args = parser.parse_args(argv)
    configure_logging(args)
    converter = converter_type(resource_profile=configured_resource_profile(args, parser))
    if args.evidence_dir and hasattr(converter.source, "http"):
        converter.source.http.evidence_dir = Path(args.evidence_dir)
    from ..conversion_logging import (
        operation, stage, labels, emit, count, result_diagnostics, outcome, statistics,
    )
    seen, reports, failed = set(), [], False
    from ..sources.archive_support import Resolution
    with operation(type(converter).__name__) as run:
        emit("batch_planned", inputs=len(args.accession))
        resolutions = []
        with stage("discovery", inputs=len(args.accession)):
            for accession in args.accession:
                with labels(item_id=accession, source=accession):
                    try:
                        resolutions.append(converter.source.resolve(accession))
                    except Exception as error:
                        resolutions.append(Resolution(issues=[type(error).__name__]))
        seeds = list(dict.fromkeys(seed for resolution in resolutions for seed in resolution.studies))
        totals = dict.fromkeys(("complete", "partial", "failed", "skipped"), 0)
        for index, (accession, resolution) in enumerate(zip(args.accession, resolutions), 1):
            with labels(item_id=str(index), source=accession):
                emit("item_started", status="running", index=index, inputs=len(args.accession))
                result = converter.convert(accession, _resolution=resolution, _filename_seeds=seeds, out=args.out, enrich_from_geo_ae=args.enrich_from_geo_ae,
                    include_peer=args.include_peer, evidence_dir=args.evidence_dir, overwrite=args.overwrite,
                    seen_studies=seen)
                reports.append(result.to_mapping())
                failed |= not result.ok
                status = "skipped" if not result.studies and resolution.studies else outcome(result)
                result_diagnostics(result)
                totals[status] += 1
                details = statistics(result)
                for key in ("studies", "packages", "samples", "runs"):
                    if key in details:
                        count(key, details[key])
                emit("item_finished", status=status, index=index, inputs=len(args.accession),
                     **{key: value for key, value in details.items() if key not in totals})
        if args.report:
            try:
                publish_json(args.report, {'imports': reports}, args.overwrite)
            except OSError as error:
                logging.getLogger(__name__).error('Report publication failed: %s', type(error).__name__)
                failed = True
        run.statistics.update(totals | {"inputs": len(args.accession)})
        run.status = "failed" if totals["failed"] == len(args.accession) else "partial" if failed else "complete"
    return int(failed)
