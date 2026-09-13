"""Argument/report handling shared by the independent native commands."""
import argparse
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
    seen, reports, failed = set(), [], False
    for accession in args.accession:
        result = converter.convert(accession, out=args.out, enrich_from_geo_ae=args.enrich_from_geo_ae,
            include_peer=args.include_peer, evidence_dir=args.evidence_dir, overwrite=args.overwrite,
            seen_studies=seen)
        reports.append(result.to_mapping())
        failed |= not result.ok
    if args.report:
        publish_json(args.report, {'imports': reports}, args.overwrite)
    return int(failed)
