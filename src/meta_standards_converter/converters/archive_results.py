# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Native import outcomes and independent atomic file publication."""
from meta_standards_converter.conversion_logging import measured, count, progress
from meta_standards_converter.diagnostics import diagnostic_mapping
from dataclasses import dataclass, field
from pathlib import Path
import json
import os
import re
from uuid import uuid4
from meta_standards_converter.miniml import MINiMLPackage


@dataclass
class StudyImportOutcome:
    study: str
    primary: str
    status: str
    package: MINiMLPackage | None = None
    output: str | None = None
    issues: list[str] = field(default_factory=list)
    error: Exception | None = field(default=None, repr=False, compare=False)
    diagnostics: list = field(default_factory=list)

    def to_mapping(self):
        return {**diagnostic_mapping(self), 'study': self.study, 'primary': self.primary, 'status': self.status,
                'output': self.output, 'issues': list(self.issues)}


@dataclass
class ArchiveImportResult:
    accession: str
    studies: list[StudyImportOutcome] = field(default_factory=list)

    @property
    def packages(self):
        return [s.package for s in self.studies if s.package is not None]

    @property
    def ok(self):
        return all(s.status == 'complete' for s in self.studies)

    def to_mapping(self):
        return {'accession': self.accession, 'studies': [s.to_mapping() for s in self.studies]}


class NativeProviderFailure(RuntimeError):
    """Copy-safe, redacted provider failure retained by operation-local memoization."""


def native_failure(errors):
    """Preserve a transient provider failure and the latest applicable cooldown."""
    from .unified.outcomes import failure
    from ..runtime_contracts import RetryCategory
    from ..helpers.metadata_stream import BatchPayloadFailure
    candidates = []
    for error in errors:
        if error is None:
            continue
        if isinstance(error, BatchPayloadFailure) and error.__cause__ is not None:
            error = error.__cause__
        envelope, retry_at = failure(error, stage='load', provider=None, item_id=None)
        candidates.append((envelope.retry_category == RetryCategory.RETRYABLE, retry_at or 0, error))
    if not candidates:
        return None
    transient, retry_at, original = max(candidates, key=lambda item: item[:2])
    if isinstance(original, NativeProviderFailure):
        return original
    error = NativeProviderFailure(f'Native provider retrieval failed ({type(original).__name__})')
    # Keep the redacted provider identity across the copy-safe native wrapper.
    # Retaining the original exception would also retain response bodies/URLs.
    cause, seen = original, set()
    while cause is not None and id(cause) not in seen:
        seen.add(id(cause))
        issue = getattr(cause, '_msc_provider_issue', None)
        if issue is not None:
            error._msc_provider_issue = issue
            break
        cause = cause.__cause__
    from ..diagnostics import diagnose
    error._msc_preserved_diagnostic = diagnose(original, stage='load')
    if hasattr(original, '_msc_failure'):
        error._msc_failure = original._msc_failure
    error.retry_at = retry_at or None
    error.request_trace = {'errors': [{'retry_category': 'retryable' if transient else 'terminal'}]}
    response = getattr(original, 'response', None)
    status = getattr(response, 'status_code', None)
    if isinstance(status, int):
        from types import SimpleNamespace
        error.response = SimpleNamespace(status_code=status)
    return error


def _load_native(source, parser, seed, *, _outcome=None, publications=None):
    """Retain a strictly valid core before optional provider discovery can stop."""
    from inspect import signature
    from ..metadata.retrieval_budget import _active, _phase, _budget_stage, _BudgetExhausted, _coverage
    from ..miniml import MINiMLCodec
    from ..miniml.codec import MINiMLCompatibilityError
    from ..miniml.native_checkpoint import _checkpoint, _complete, _valid_native
    budget, candidate, captured, checkpoint = _active.get(), None, None, None
    def core(records):
        nonlocal candidate, captured, checkpoint
        captured = records
        if not records.xml and not any(records.indexed.values()): return
        try:
            checkpoint = _checkpoint(parser,records)
            candidate = checkpoint.package if checkpoint is not None else MINiMLCodec().decode(parser.parse(records),strict=True).package
        except MINiMLCompatibilityError: pass
        _phase.set('native_supplementary')  # Restored by the enclosing stage scope.
    try:
        with _budget_stage('native_core'):
            options = {'_on_core': core} if budget and '_on_core' in signature(source.fetch).parameters else {}
            if publications is not None:
                options["publications"] = publications
            records = source.fetch(seed, **options)
            captured = records
        if not records.xml and not any(records.indexed.values()):
            error = native_failure(getattr(records, 'failures', []))
            if error is not None:
                raise error
            raise ValueError('No native metadata records retrieved')
        if checkpoint is not None:
            from ..miniml.model import MINiMLModelError
            try:
                return _complete(checkpoint,parser,records), records.issues
            except (MINiMLCompatibilityError, MINiMLModelError) as error:
                from ..miniml.native_checkpoint import _projection_failed
                records.issues.append(f'native supplementary projection: {type(error).__name__}; valid native candidate retained')
                records.failures.append(error)
                return _projection_failed(checkpoint,error), records.issues
        return parser.parse(records), records.issues
    except _BudgetExhausted as exc:
        if candidate is None: raise
        try:
            candidate = (_valid_native(_complete(checkpoint,parser,captured)) if checkpoint is not None
                         else MINiMLCodec().decode(parser.parse(captured), strict=True).package)
        except Exception:
            pass  # Supplementary records cannot invalidate the accepted core.
        return _coverage(candidate, budget), [*captured.issues, str(exc) + '; valid native candidate retained']
    finally:
        if _outcome is not None and captured is not None:
            _outcome.error = native_failure(getattr(captured, 'failures', []))


@measured("publication")
def publish_json(path, value, overwrite=False, *, sort_keys=False):
    from ..operational_events import _guard_output
    _guard_output([path])
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name('.' + path.name + '.' + uuid4().hex + '.tmp')
    try:
        with temporary.open('x', encoding='utf-8') as handle:
            json.dump(value, handle, indent=2, ensure_ascii=False, sort_keys=sort_keys)
            handle.write('\n')
            handle.flush()
            os.fsync(handle.fileno())
        if overwrite:
            os.replace(temporary, path)
        else:
            os.link(temporary, path)  # atomic no-replace, including competing imports
        count("artifacts_committed")
        progress(1, 1, "artifacts")
    finally:
        temporary.unlink(missing_ok=True)


def output_name(seed, seeds):
    value = seed.primary
    if sum(other.primary == seed.primary for other in seeds) > 1:
        value += '__' + seed.study
    if not re.fullmatch(r'[A-Za-z0-9_.-]+', value):
        raise ValueError('Unsafe study filename')
    return value + '.json'
