# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Controlled operational diagnostics, independent of console presentation.

Domain validators retain their own scientific codes. This catalogue describes
operational consequences and recovery; it never decides retry or completeness.
Optional scientific dependencies are not imported to classify their failures.
"""
import errno
import json
import re
import xml.etree.ElementTree as ET
from dataclasses import replace


CATALOGUE = {
    'missing_file': ('An input or required companion file was not found.', 'Check the input path and supply the required files.'),
    'permission_denied': ('A required file or directory is not accessible.', 'Check permissions or choose an accessible destination.'),
    'disk_full': ('Storage is full; the requested write could not finish.', 'Free disk space or choose another destination, then retry.'),
    'disk_budget_exceeded': ('The destination has insufficient disk headroom.', 'Free disk space or choose a destination with sufficient capacity.'),
    'output_exists': ('The output already exists.', 'Choose another destination or explicitly enable overwrite.'),
    'invalid_json': ('The input is not valid JSON.', 'Correct the JSON syntax at the reported position and retry.'),
    'invalid_xml': ('The input is not valid XML.', 'Correct the XML syntax or obtain a fresh source file.'),
    'invalid_encoding': ('The input cannot be decoded using its required encoding.', 'Supply a correctly encoded source file.'),
    'invalid_metadata': ('Metadata failed validation.', 'Correct the reported metadata fields or regenerate from the original source.'),
    'projection_invalid': ('The requested output failed validation.', 'Correct the source annotations or projector diagnostics before exporting.'),
    'missing_optional_dependency': ('Expression conversion requires optional dependencies.', 'Install meta-standards-converter[h5ad] or use the full Docker image.'),
    'missing_runtime': ('A required processing tool or container runtime is unavailable.', 'Install and configure the required processing runtime, then retry.'),
    'processing_failed': ('The expression processing pipeline failed.', 'Inspect the processing log and correct the pipeline inputs or runtime before retrying.'),
    'processing_output_missing': ('Processing did not produce the required expression output.', 'Inspect the processing log and check the selected pipeline and inputs.'),
    'memory_limit': ('The operation exceeds the available memory budget.', 'Use a smaller input or a host with more memory; review the resource profile.'),
    'retrieval_budget_exhausted': ('Metadata retrieval exhausted the configured budget.', 'Increase the retrieval budget or supply local metadata files.'),
    'retrieval_size_limit': ('A download, archive or cache exceeds its resource limit.', 'Review the resource profile or supply a smaller local input.'),
    'unsafe_input': ('The input violates the retrieval or XML safety policy.', 'Use a trusted input satisfying the documented retrieval and XML restrictions.'),
    'cache_integrity': ('Cached data failed its integrity check.', 'Replace the affected cache entry from its original source and retry.'),
    'provider_transport_failure': ('The provider request could not complete.', 'Retry after the reported cooldown; usable partial output can be retained.'),
    'provider_http_failure': ('The provider returned an unsuccessful HTTP response.', 'Check provider availability and the reported status and retry eligibility.'),
    'provider_response_invalid': ('The provider returned malformed or incomplete metadata.', 'Retry later or supply verified local metadata.'),
    'publication_failed': ('An output could not be published.', 'Check the destination and retain any artifacts already reported as written.'),
    'reporting_failed': ('A requested log or report could not be written.', 'Check the reporting destination; conversion artifacts already written remain available.'),
    'conversion_failed': ('An unexpected internal error stopped this operation.', 'Report this diagnostic and correlation ID with a DEBUG log.'),
    'source_partial': ('Optional source metadata is incomplete.', 'Retain usable output; inspect the source diagnostic before retrying.'),
    'source_fallback': ('An alternative verified source supplied partial metadata.', 'Supply the original source files to restore missing source-specific metadata.'),
}


class RequestError(ValueError):
    """Explicit request validation, distinguishable from an internal ValueError."""


class OperationalError(RuntimeError):
    """Controlled operational failure retaining a catalogue code and safe location."""
    def __init__(self, code, *, location=None):
        self.code = code
        self.location = location
        super().__init__(CATALOGUE[code][0])


def _causes(error):
    seen = set()
    while error is not None and id(error) not in seen:
        seen.add(id(error))
        yield error
        error = error.__cause__ or (None if error.__suppress_context__ else error.__context__)


def _safe_label(value):
    return re.sub(r'[^A-Za-z0-9_./:< >-]', '_', str(value))[:160]


def _frames(error):
    """Bounded locations only: never read locals, source lines or exception text."""
    causes = []
    for cause in _causes(error):
        frames = []
        tb = cause.__traceback__
        while tb is not None:
            frames.append({'module': _safe_label(tb.tb_frame.f_globals.get('__name__', 'unknown')),
                           'function': _safe_label(tb.tb_frame.f_code.co_name), 'line': tb.tb_lineno})
            frames = frames[-20:]
            tb = tb.tb_next
        causes.append({'error_type': _safe_label(type(cause).__name__), 'frames': frames})
        if len(causes) == 5:
            break
    return causes


def _classify(error, stage):
    from .converters.unified.contracts import InputError
    from .sources.errors import MetadataSourceError
    from .runtime_contracts import DiskBudgetError
    from .metadata.retrieval_budget import _BudgetExhausted
    from .miniml import MINiMLModelError
    from .miniml.codec import MINiMLCompatibilityError
    from .atlas_v1.reader import AtlasV1Error
    from .magetab.validation import MAGETabValidationError
    from .magetab.semantics import MAGETabModelError
    from .metadata.projection.tabular import TabularProjectionError
    from .metadata.projection.anndata import AnnDataProjectionError
    from .retrieval import RetrievalSecurityError, RetrievalSizeError, CacheIntegrityError
    from .xml_safety import UnsafeXMLDocumentError, XMLSizeLimitError
    import requests
    if isinstance(error, (InputError, MetadataSourceError, OperationalError)):
        return error.code, str(error)
    if isinstance(error, _BudgetExhausted): return 'retrieval_budget_exhausted', None
    if isinstance(error, (MemoryError,)): return 'memory_limit', None
    if isinstance(error, DiskBudgetError): return 'disk_budget_exceeded', None
    if isinstance(error, (RetrievalSecurityError, UnsafeXMLDocumentError)): return 'unsafe_input', None
    if isinstance(error, (RetrievalSizeError, XMLSizeLimitError)): return 'retrieval_size_limit', None
    if isinstance(error, CacheIntegrityError): return 'cache_integrity', None
    scientific = {'anndata', 'numpy', 'pandas', 'scipy', 'h5py', 'scanpy'}
    if any(isinstance(c, ImportError) and (c.name or '').split('.')[0] in scientific for c in _causes(error)):
        return 'missing_optional_dependency', None
    if isinstance(error, (MINiMLModelError, MINiMLCompatibilityError, AtlasV1Error)):
        return 'invalid_metadata', None
    if isinstance(error, (TabularProjectionError, AnnDataProjectionError, MAGETabValidationError, MAGETabModelError)):
        return 'projection_invalid', None
    if isinstance(error, json.JSONDecodeError): return 'invalid_json', None
    if isinstance(error, ET.ParseError): return 'invalid_xml', None
    if isinstance(error, UnicodeError): return 'invalid_encoding', None
    if isinstance(error, FileNotFoundError): return 'missing_file', None
    if isinstance(error, PermissionError): return 'permission_denied', None
    if isinstance(error, FileExistsError): return 'output_exists', None
    if isinstance(error, OSError) and error.errno == errno.ENOSPC: return 'disk_full', None
    if isinstance(error, requests.HTTPError): return 'provider_http_failure', None
    if isinstance(error, (requests.RequestException, TimeoutError, ConnectionError)): return 'provider_transport_failure', None
    if isinstance(error, OSError):
        return ('reporting_failed' if stage in {'report', 'logging'} else 'publication_failed'), None
    return 'conversion_failed', None


def diagnose(error, *, stage, provider=None, item_id=None, location=None, severity='error'):
    """Describe an exception safely, sharing its existing recovery envelope."""
    from .converters.unified.contracts import Diagnostic
    from .converters.unified.outcomes import failure
    from .runtime_contracts import SafeErrorEnvelope
    from .conversion_logging import emit
    envelope, _ = failure(error, stage=stage, provider=provider, item_id=item_id)
    code, message = _classify(error, stage)
    issue = next((getattr(c, '_msc_provider_issue', None) for c in _causes(error)
                  if getattr(c, '_msc_provider_issue', None) is not None), None)
    if issue is not None and code in {'invalid_json', 'invalid_xml', 'invalid_encoding', 'conversion_failed'}:
        code, message = issue.reason, None
    default, action = CATALOGUE.get(code, ('The requested operation could not complete.',
                                         'Review the reported input condition and correct it before retrying.'))
    message = message or default
    if isinstance(error, json.JSONDecodeError):
        message += f' Line {error.lineno}, column {error.colno}.'
    elif isinstance(error, ET.ParseError) and getattr(error, 'position', None):
        message += f' Line {error.position[0]}, column {error.position[1]}.'
    if code == 'invalid_metadata':
        from .miniml.codec import MINiMLCompatibilityError
        if isinstance(error, MINiMLCompatibilityError):
            details = [f'{_safe_label(d.code)} at {_safe_label(d.path)}' for d in error.diagnostics[:5]]
            message += ' ' + '; '.join(details)
    if envelope.http_status is not None and f'HTTP {envelope.http_status}' not in message:
        message += f' HTTP {envelope.http_status}.'
    if code == 'conversion_failed' and not getattr(error, '_msc_debug_recorded', False):
        emit('internal_failure', level='DEBUG', correlation_id=envelope.correlation_id, causes=_frames(error))
        error._msc_debug_recorded = True
    safe_location = SafeErrorEnvelope.from_exception(error, location=location or getattr(error, 'location', None) or getattr(error, 'filename', None)).location
    return Diagnostic(code, message, stage, severity, action=action,
                      provider=envelope.provider, operation=issue.operation if issue else None,
                      location=safe_location, correlation_id=envelope.correlation_id)


def enrich(diagnostic):
    """Supply catalogue guidance/context without changing domain severity or codes."""
    issue = diagnostic.provider_issue
    return replace(diagnostic,
        action=diagnostic.action or CATALOGUE.get(diagnostic.code, (None, None))[1],
        provider=diagnostic.provider or (issue.error.provider if issue else None),
        operation=diagnostic.operation or (issue.operation if issue else None),
        correlation_id=diagnostic.correlation_id or (issue.error.correlation_id if issue else None))
