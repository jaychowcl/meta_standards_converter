# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Terminal interface to the unified Converter; domain behavior stays in the API."""

import argparse
from dataclasses import asdict, is_dataclass
import json
import logging
import os
from pathlib import Path
import sys

from meta_standards_converter.converters import Converter
from meta_standards_converter.converters.unified.discovery import (
    ALIASES, KINDS, expand, idf_references, is_url, manifest_specs,
)
from meta_standards_converter.converters.unified.options import TARGETS, RUNTIME_DEFAULTS
from meta_standards_converter.converters.unified.publication import atomic_write
from meta_standards_converter.runtime_contracts import get_resource_profile
from .common import (
    add_logging_arguments, add_platform_handler_arguments,
    add_replacement_profile_arguments, configure_logging, parse_resource_override,
    print_platform_handlers, record_safe_cli_error, replacement_profile_from_args,
)

logger = logging.getLogger(__name__)


def _boolean(group, name, help, *, aliases=(), negative_aliases=()):
    choice = group.add_mutually_exclusive_group()
    dest = name.replace('-', '_')
    choice.add_argument('--' + name, *aliases, dest=dest, action='store_true',
                        default=argparse.SUPPRESS, help=help)
    choice.add_argument('--no-' + name, *negative_aliases, dest=dest,
                        action='store_false', default=argparse.SUPPRESS,
                        help='Disable ' + name.replace('-', ' ') + '.')


def _parser():
    parser = argparse.ArgumentParser(
        prog='msc-convert',
        description='Convert accessions, files, URLs or directories through the unified Converter.',
        epilog='Defaults: automatic input detection, standard enrichment, study expansion on, '
               'output in the current directory, overwrite off. Per-input manifest settings '
               'override command-level defaults. Raw processing requires --allow-processing.',
    )
    parser.add_argument('inputs', nargs='*', metavar='INPUT')
    parser.add_argument('--out-type', choices=sorted(TARGETS))
    parser.add_argument('--in-type', '--force-in-type', choices=['auto', *sorted(KINDS | ALIASES.keys())], default=argparse.SUPPRESS)
    parser.add_argument('--input-manifest', default=argparse.SUPPRESS, help='Version 1.0 JSON input manifest; may accompany positional inputs.')
    destination = parser.add_mutually_exclusive_group()
    destination.add_argument('--out', '--outdir', '-o', dest='outdir', help='Output directory (default: current directory).')
    destination.add_argument('--outfile', help='Exact single-file output; incompatible with MAGE-TAB and batch inputs.')
    listing = parser.add_mutually_exclusive_group()
    listing.add_argument('--list-input-types', action='store_true')
    listing.add_argument('--list-output-types', action='store_true')
    add_platform_handler_arguments(parser)
    preparation = parser.add_argument_group('input preparation')
    preparation.add_argument('--enrichment', choices=('standard', 'curators', 'off'), default=argparse.SUPPRESS)
    _boolean(preparation, 'expand-studies', 'Expand verified related studies (default: enabled).', aliases=('--related', '--related-series'), negative_aliases=('--no-related',))
    for name, help in (
        ('enrich', 'Compatibility enrichment switch; prefer --enrichment.'),
        ('include-peer', 'Allow peer archive enrichment.'),
        ('enrich-from-geo-ae', 'Allow linked GEO/ArrayExpress metadata enrichment.'),
        ('remove-empty', 'Remove empty parsed GEO fields.'),
    ):
        _boolean(preparation, name, help)
    preparation.add_argument('--evidence-dir', default=argparse.SUPPRESS)
    preparation.add_argument('--sdrf-source', dest='sdrf_sources', action='append', default=argparse.SUPPRESS)
    preparation.add_argument('--orientation', choices=('auto', 'genes-by-observations', 'observations-by-genes'), default=argparse.SUPPRESS, help='Orientation of an input matrix.')
    output = parser.add_argument_group('output construction')
    for name in ('aggregate', 'allow-invalid', 'include-var', 'include-uns'):
        _boolean(output, name, {'aggregate': 'Aggregate datasets within one logical metadata input.', 'allow-invalid': 'Allow outputs with projector-reported validation errors.', 'include-var': 'Include feature metadata in OBS output.', 'include-uns': 'Include unstructured metadata in OBS output.'}[name])
    add_replacement_profile_arguments(output)
    assets = parser.add_argument_group('expression assets and processing')
    assets.add_argument('--asset', dest='asset_specs', action='append', metavar='ACCESSION=PATH_OR_URL', default=argparse.SUPPRESS)
    assets.add_argument('--asset-manifest', default=argparse.SUPPRESS, help='Existing CSV/TSV asset manifest.')
    explicit = assets.add_mutually_exclusive_group()
    explicit.add_argument('--explicit-assets', help='JSON array of Asset field mappings.')
    explicit.add_argument('--explicit-assets-file', help='JSON file containing an array of Asset field mappings.')
    assets.add_argument('--matrix-orientation', choices=('auto', 'genes-by-observations', 'observations-by-genes'), default=argparse.SUPPRESS)
    assets.add_argument('--pipeline', choices=('auto', 'rnaseq', 'scrnaseq'), default=argparse.SUPPRESS)
    for name in ('allow-processing', 'force-reprocess', 'accept-inferred-reference', 'resume', 'force-memory'):
        _boolean(assets, name, {'allow-processing': 'Permit raw-read pipeline execution; references are still required.', 'force-reprocess': 'Rebuild from raw reads instead of using processed assets.', 'accept-inferred-reference': 'Allow a supported reference inferred from organism metadata.', 'resume': 'Resume from workflow and processing checkpoints.', 'force-memory': 'On resume, bypass the profile ceiling but retain the available-memory safety ceiling.'}[name])
    _boolean(assets, 'allow-unverified-combination', 'Deprecated compatibility option; never combines matrices.')
    for name in ('genome', 'fasta', 'gtf', 'gff', 'revision', 'params-file', 'nextflow-config', 'work-dir', 'processed-checkpoint-dir'):
        assets.add_argument('--' + name, default=argparse.SUPPRESS)
    assets.add_argument('--execution-profile', '--profile', dest='profile', default=argparse.SUPPRESS, help='Nextflow execution profile.')
    runtime = parser.add_argument_group('runtime and retrieval policy')
    for name in ('overwrite', 'fail-fast', 'recursive'):
        _boolean(runtime, name, {'overwrite': 'Replace existing conversion outputs and report/log files.', 'fail-fast': 'Skip later inputs after the first unsuccessful input.', 'recursive': 'Discover inputs recursively, excluding caches and nested outputs.'}[name])
    runtime.add_argument('--resource-profile', choices=('standard', 'large'), default=argparse.SUPPRESS)
    runtime.add_argument('--resource-override', action='append', type=parse_resource_override, default=argparse.SUPPRESS, metavar='FIELD=VALUE')
    runtime.add_argument('--insdc-default', choices=('ena', 'sra'), default=argparse.SUPPRESS)
    runtime.add_argument('--allowed-host', '--asset-host', dest='allowed_hosts', action='append', default=argparse.SUPPRESS)
    runtime.add_argument('--reserved-path', dest='reserved_paths', action='append', default=argparse.SUPPRESS, help='Reserve an additional path against discovery/publication.')
    reports = parser.add_argument_group('results and logging')
    reports.add_argument('--report-json', metavar='PATH_OR_MINUS', help='Write a structured result summary; use - for JSON-only stdout.')
    add_logging_arguments(reports)
    return parser


def _asset_objects(args):
    raw = args.explicit_assets
    if args.explicit_assets_file is not None:
        raw = Path(args.explicit_assets_file).read_text(encoding='utf-8')
    if raw is None:
        return None
    from meta_standards_converter.expression.assets import Asset
    values = json.loads(raw)
    if not isinstance(values, list) or not all(isinstance(v, dict) for v in values):
        raise ValueError('Explicit assets must be a JSON array of Asset field mappings')
    return [Asset(**v) for v in values]


def _options(args, parser):
    cli_only = {'inputs', 'out_type', 'in_type', 'outdir', 'outfile', 'enrichment',
                'list_input_types', 'list_output_types', 'list_platform_handlers',
                'replacement_profile', 'replacement_profile_file', 'explicit_assets',
                'explicit_assets_file', 'report_json', 'verbose', 'quiet', 'log_file'}
    options = {k: v for k, v in vars(args).items() if k not in cli_only and v is not None}
    overrides = options.pop('resource_override', None)
    if overrides is not None:
        options['resource_overrides'] = {}
        for key, value in overrides:
            if key in options['resource_overrides'] and options['resource_overrides'][key] != value:
                parser.error('Conflicting resource override: ' + key)
            options['resource_overrides'][key] = value
    profile = replacement_profile_from_args(args, parser)
    if profile is not None:
        options['replacement_profile'] = profile
    assets = _asset_objects(args)
    if assets is not None:
        options['explicit_assets'] = assets
    return options


def _paths(value):
    """Conservative local-path inventory for protecting auxiliary output files."""
    if is_dataclass(value):
        yield from _paths(asdict(value))
    elif isinstance(value, dict):
        for child in value.values():
            yield from _paths(child)
    elif isinstance(value, (list, tuple)):
        for child in value:
            yield from _paths(child)
    elif isinstance(value, (str, Path)) and not is_url(str(value)):
        yield Path(value).resolve()


def _auxiliary_paths(args, options, parser):
    auxiliary = [Path(p).resolve() for p in (args.report_json if args.report_json != '-' else None, args.log_file) if p]
    if not auxiliary:
        return
    if len(auxiliary) != len(set(auxiliary)):
        parser.error('Report and log must have different paths')
    protected = set(_paths(args.inputs))
    # Manifests own sources, companions, metadata and path-bearing options.
    profile = get_resource_profile(options.get('resource_profile', 'standard'), overrides=options.get('resource_overrides'))
    specs = expand(args.inputs, RUNTIME_DEFAULTS | {k: v for k, v in options.items() if k in RUNTIME_DEFAULTS},
                   output_directory=args.outdir or (None if args.outfile else '.'))
    if options.get('input_manifest'):
        specs += manifest_specs(options['input_manifest'], profile.max_xml_bytes)
    for spec in specs:
        protected.update(_paths(vars(spec)))
        for path in _paths(spec.sources):
            if path.is_file() and path.name.lower().endswith('.idf.txt'):
                protected.update((path.parent / ref).resolve() for ref in idf_references(path, profile.max_xml_bytes))
    if options.get('asset_specs') or options.get('asset_manifest'):
        from meta_standards_converter.expression.assets import AssetManifest
        manifest = AssetManifest()
        assets = [manifest.parse_spec(value) for value in options.get('asset_specs', ())]
        if options.get('asset_manifest'):
            assets.extend(manifest.load(options['asset_manifest']))
        protected.update(_paths(assets))
    protected.update(_paths({k: v for k, v in options.items() if k != 'reserved_paths'}))
    protected.update(_paths([args.replacement_profile_file, args.explicit_assets_file]))
    if args.outfile:
        protected.add(Path(args.outfile).resolve())
    for path in auxiliary:
        output_root = Path(args.outdir or '.').resolve() if not args.outfile else Path(args.outfile).resolve()
        if output_root.is_relative_to(path) or '.artifact-bundles' in path.parts:
            parser.error('Report/log path conflicts with an output directory or internal artifact bundle')
        if path in protected:
            parser.error('Report/log path conflicts with an input, configuration or exact output')
        if path.exists() and (not path.is_file() or not options.get('overwrite', False)):
            parser.error('Report/log destination exists; choose another path or --overwrite')
        parent = path.parent
        while not parent.exists():
            parent = parent.parent
        if not parent.is_dir() or not os.access(parent, os.W_OK):
            parser.error('Report/log destination parent is not writable')
    options['reserved_paths'] = list(options.get('reserved_paths', ())) + [str(p) for p in auxiliary]


def _human_summary(result):
    print('Conversion: ' + result.status)
    for item in result.items:
        print(f'{item.id}: {item.status} ({item.source})')
        for name, path in item.artifacts.items():
            print(f'  {name}: {path}')
        for diagnostic in item.diagnostics:
            print(f'  {diagnostic.severity} [{diagnostic.code}]: {diagnostic.message}')


def main(argv=None):
    parser = _parser()
    args = parser.parse_args(argv)
    if args.list_input_types:
        for kind in sorted(KINDS):
            aliases = ', '.join(k for k, v in ALIASES.items() if v == kind)
            print(kind + (f' (aliases: {aliases})' if aliases else ''))
        print('In-memory-only objects require the Python API; use files or manifests in the CLI.')
        return 0
    if args.list_output_types:
        print('\n'.join(sorted(TARGETS)))
        return 0
    if args.list_platform_handlers:
        print_platform_handlers()
        return 0
    if not args.out_type:
        parser.error('--out-type is required')
    if not args.inputs and not hasattr(args, 'input_manifest'):
        parser.error('At least one INPUT or --input-manifest is required')
    try:
        options = _options(args, parser)
        _auxiliary_paths(args, options, parser)
    except (ValueError, TypeError, OSError) as error:
        parser.error(str(error))
    kwargs = {'out_type': args.out_type, 'outdir': args.outdir or (None if args.outfile else '.'), 'options': options}
    for name in ('in_type', 'enrichment'):
        if hasattr(args, name):
            kwargs[name] = getattr(args, name)
    if args.outfile:
        options['outfile'] = args.outfile
    try:
        if args.log_file:
            Path(args.log_file).parent.mkdir(parents=True, exist_ok=True)
        configure_logging(args, stream=sys.stderr)
        result = Converter().convert(args.inputs or None, **kwargs)
    except (ValueError, TypeError) as error:
        parser.error(str(error))
    except Exception as error:
        record_safe_cli_error(logger, error, location='msc-convert', stage='conversion')
        return 1
    if args.report_json != '-':
        _human_summary(result)
    if args.report_json:
        try:
            if args.report_json == '-':
                print(json.dumps(result.to_dict(), ensure_ascii=False))
            else:
                report_path = Path(args.report_json).resolve()
                if any(report_path == Path(path).resolve() for item in result.items for path in item.artifacts.values()):
                    raise ValueError('Report destination is a published conversion artifact')
                atomic_write(args.report_json, lambda stream: json.dump(result.to_dict(), stream, ensure_ascii=False, indent=2), options.get('overwrite', False))
        except Exception as error:
            record_safe_cli_error(logger, error, location='msc-convert', stage='report')
            return 1
    return 0 if result.status == 'complete' else 1


if __name__ == '__main__':
    raise SystemExit(main())
