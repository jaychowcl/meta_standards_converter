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
import os
from pathlib import Path

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

    parser.add_argument("--log-format", choices=("text", "jsonl"), default="text", help="Log file format (default: text).")
    parser.add_argument("--log-file-level", choices=("DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"), default="DEBUG", help="File threshold, independent of console verbosity (default: DEBUG).")
    parser.add_argument("--progress", choices=("auto", "plain", "off"), default="auto", help="Terminal progress; quiet disables it (default: auto).")


def log_level(args) -> int:
    if args.quiet:
        return logging.ERROR
    if args.verbose >= 2:
        return logging.DEBUG
    if args.verbose == 1:
        return logging.INFO
    return logging.WARNING


def configure_logging(args, *, stream=None) -> None:
    from .logging_output import ConsoleHandler, EventFormatter, SafeFileHandler, register_configuration
    preflight_reporting(args)
    package_logger = logging.getLogger("meta_standards_converter")
    previous = (package_logger.level, package_logger.propagate)
    for handler in package_logger.handlers[:]:
        if getattr(handler, "_msc_owned", False):
            package_logger.removeHandler(handler)
            handler.close()
    package_logger.setLevel(logging.DEBUG)
    package_logger.propagate = False
    total = None
    for name in ("inputs", "gse", "json_path", "source", "accession"):
        value = getattr(args, name, None)
        if isinstance(value, (list, tuple)):
            total = len(value)
            break
    handler = ConsoleHandler(sys.stderr if stream is None else stream, threshold=log_level(args),
        mode=getattr(args, "progress", "auto"), quiet=args.quiet, total=total)
    handlers = [handler]
    package_logger.addHandler(handler)
    register_configuration(args, package_logger, handlers, previous)
    if args.log_file:
        Path(args.log_file).parent.mkdir(parents=True, exist_ok=True)
        file_handler = SafeFileHandler(args.log_file, mode="w", encoding="utf-8")
        file_handler._msc_owned = True
        file_handler.setLevel(getattr(args, "log_file_level", "DEBUG"))
        file_handler.setFormatter(EventFormatter(jsonl=getattr(args, "log_format", "text") == "jsonl", timestamps=True))
        package_logger.addHandler(file_handler)
        handlers.append(file_handler)


def record_safe_cli_error(
    logger: logging.Logger,
    error: BaseException,
    *,
    location: str,
    stage: str,
    provider: str | None = None,
) -> SafeErrorEnvelope:
    """Log and return a durable error without serializing its raw message."""

    from dataclasses import replace
    from meta_standards_converter.diagnostics import diagnose
    from meta_standards_converter.converters.unified.outcomes import failure
    if stage in {'report', 'logging'}:
        from .logging_output import _session
        session = _session.get()
        if session is not None:
            session['reporting_failed'] = True
    envelope, _ = failure(error, stage=stage, provider=provider, item_id=None)
    safe_location = SafeErrorEnvelope.from_exception(error, location=location).location
    envelope = replace(envelope, location=envelope.location or safe_location)
    diagnostic = getattr(error, '_msc_diagnostic', None) or diagnose(
        error, stage=envelope.stage or stage, provider=provider, location=location)
    error._msc_diagnostic = diagnostic
    if not getattr(error, '_msc_logged_envelope', None):
        from .logging_output import _session
        message = (f'{diagnostic.stage} [{diagnostic.code}]: {diagnostic.message} '
                   f'Action: {diagnostic.action} [correlation_id={diagnostic.correlation_id}]')
        session = _session.get()
        if session is not None and not session.get('handlers'):
            # Preflight failed before installing any command-owned handler.
            sys.stderr.write('ERROR ' + message + '\n')
        else:
            try:
                logger.error('%s', message)
            except Exception:
                from meta_standards_converter.operational_events import _report_sink_failure
                _report_sink_failure()
                sys.stderr.write('ERROR ' + message + '\n')
        error._msc_logged_envelope = envelope
    return envelope


def preflight_reporting(args):
    """Protect all command inputs/configuration before an auxiliary file opens."""
    from meta_standards_converter.diagnostics import RequestError
    from .convert import _paths
    auxiliary = [Path(p).resolve() for p in (
        getattr(args, 'log_file', None), getattr(args, 'report', None),
        getattr(args, 'report_json', None)) if p and p != '-']
    if not auxiliary:
        return
    if len(set(auxiliary)) != len(auxiliary):
        raise RequestError('Report and log must have different paths.')
    names = ('inputs', 'json_path', 'source', 'input_manifest', 'asset_manifest',
             'replacement_profile_file', 'explicit_assets_file', 'fasta', 'gtf', 'gff',
             'params_file', 'nextflow_config', 'destination', 'outfile')
    protected = set(_paths([getattr(args, name, None) for name in names]))
    assets = getattr(args, 'asset', None) or []
    if assets or getattr(args, 'asset_manifest', None):
        from meta_standards_converter.expression.assets import AssetManifest
        manifest = AssetManifest()
        try:
            entries = [manifest.parse_spec(value) for value in assets]
            if getattr(args, 'asset_manifest', None):
                entries.extend(manifest.load(args.asset_manifest))
            protected.update(_paths(entries))
        except (ValueError, OSError) as error:
            raise RequestError('Cannot validate asset paths for log/report protection; check the asset manifest.') from error
    directories = [Path(getattr(args, name)).resolve() for name in ('outdir', 'out')
                   if getattr(args, name, None)]
    for path in auxiliary:
        if path in protected or '.artifact-bundles' in path.parts or any(root.is_relative_to(path) for root in directories):
            raise RequestError('Report/log path conflicts with an input, configuration or conversion destination.')
        if path.exists() and (not path.is_file() or not getattr(args, 'overwrite', False)):
            raise RequestError('Report/log destination exists; choose another path or explicitly enable overwrite.')
        parent = path.parent
        while not parent.exists():
            parent = parent.parent
        if not parent.is_dir() or not os.access(parent, os.W_OK):
            raise RequestError('Report/log destination parent is not writable; choose another destination.')
    # A permission hint alone cannot establish that report publication can be
    # initialized (ACLs, mounts and capacity can disagree with os.access).
    # Probe a detached file now; never truncate the requested report itself.
    import tempfile
    from meta_standards_converter.diagnostics import OperationalError
    for value in (getattr(args, 'report', None), getattr(args, 'report_json', None)):
        if value and value != '-':
            path = Path(value).resolve()
            try:
                path.parent.mkdir(parents=True, exist_ok=True)
                with tempfile.TemporaryFile(dir=path.parent):
                    pass
            except OSError as error:
                raise OperationalError('reporting_failed', location=path.name) from error
    from meta_standards_converter.operational_events import _reserved_outputs
    from .logging_output import _session
    if _session.get() is not None:
        _reserved_outputs.set(frozenset(auxiliary))


def add_replacement_profile_arguments(parser):
    """Add explicit, mutually exclusive export replacement policy inputs."""
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--replacement-profile", help="Replacement profile as inline JSON; supplying it activates replacements.")
    group.add_argument("--replacement-profile-file", help="Path to a replacement profile JSON object.")


def replacement_profile_from_args(args, parser):
    """Parse once before conversion; malformed policy input must not publish files."""
    import json
    from pathlib import Path
    try:
        value = args.replacement_profile
        if args.replacement_profile_file is not None:
            value = Path(args.replacement_profile_file).read_text(encoding="utf-8")
        if value is None:
            return None
        profile = json.loads(value)
        if not isinstance(profile, dict):
            raise ValueError("replacement profile must be a JSON object")
        return profile
    except (OSError, UnicodeError, ValueError) as error:
        parser.error(str(error))
