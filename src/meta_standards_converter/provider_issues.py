# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Structured recovery facts at remote I/O boundaries, scoped to a conversion.

No warning-text parsing or global exception tracing. Local parsing and semantic
disagreements are outside these boundaries. Existing exceptions still propagate.
"""
from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass
from functools import wraps
from json import JSONDecodeError
import xml.etree.ElementTree as ET
from urllib.parse import parse_qs, urlencode, urlsplit

from .runtime_contracts import RetryCategory, SafeErrorEnvelope


@dataclass(frozen=True)
class ProviderIssue:
    reason: str
    operation: str
    error: SafeErrorEnvelope
    retry_at: float | None = None

    def to_dict(self):
        value = {'reason': self.reason, 'operation': self.operation, 'error': self.error.to_dict()}
        if self.retry_at is not None:
            value['retry_at'] = self.retry_at
        return value


_issues = ContextVar('msc_provider_issues', default=None)
_operation = ContextVar('msc_provider_operation', default='metadata_retrieval')
_provider = ContextVar('msc_issue_provider', default=None)


def snapshot():
    return dict(_issues.get() or {})


def restore(issues):
    current = _issues.get()
    if current is not None:
        current.update(issues)


def recover(error):
    """Call only when an alternative actually recovered the requested evidence."""
    issue = getattr(error, '_msc_provider_issue', None)
    if issue is not None:
        restore({issue.error.correlation_id: (issue, True)})


def record(error, provider=None):
    from dataclasses import replace
    import requests
    from .converters.unified.outcomes import failure
    from .sources.errors import MetadataSourceError
    issue = getattr(error, '_msc_provider_issue', None)
    cause = error.__cause__
    seen = set()
    while issue is None and cause is not None and id(cause) not in seen:
        seen.add(id(cause))
        issue = getattr(cause, '_msc_provider_issue', None)
        cause = cause.__cause__
    if issue is None:
        envelope, retry_at = failure(error, stage='retrieval', provider=provider, item_id=None)
        if isinstance(error, MetadataSourceError):
            reason = error.code
        elif isinstance(error, (JSONDecodeError, ET.ParseError, UnicodeDecodeError)):
            reason = 'provider_response_invalid'
            envelope = replace(envelope, retry_category=RetryCategory.RETRYABLE)
        elif envelope.http_status is not None:
            reason = 'provider_http_failure'
        elif envelope.retry_category is RetryCategory.RETRYABLE or isinstance(error, requests.RequestException):
            reason = 'provider_transport_failure'
        else:
            return None
        issue = ProviderIssue(reason, _operation.get(), envelope, retry_at)
        error._msc_provider_issue = issue
    current = _issues.get()
    if current is not None and issue.error.correlation_id not in current:
        current[issue.error.correlation_id] = (issue, False)
        from .conversion_logging import emit
        emit('provider_issue', level='WARNING', **issue.to_dict())
    return issue


def _remote_provider(args, kwargs):
    url = kwargs.get('url') or (args[1] if len(args) > 1 and isinstance(args[1], str) else '')
    host = urlsplit(url).hostname or ''
    if '/ena/' in url:
        return 'ena'
    if 'biostudies' in url:
        return 'biostudies'
    if '/biosamples/' in url:
        return 'biosamples'
    if host == 'pmc.ncbi.nlm.nih.gov':
        return 'pubmed'
    if 'eutils.ncbi' in host:
        params = kwargs.get('params', args[2] if len(args) > 2 else None)
        # requests accepts mappings, ordered pairs and already encoded queries.
        # Provider attribution must not narrow that existing transport contract.
        try:
            query = (params.decode('ascii') if isinstance(params, bytes) else
                     params if isinstance(params, str) else urlencode(params or {}, doseq=True))
            databases = parse_qs(query).get('db') or parse_qs(urlsplit(url).query).get('db')
            return databases[0] if databases else 'ncbi'
        except (TypeError, ValueError, UnicodeError):
            return 'ncbi'
    if 'geo' in url:
        return 'geo'
    return {'ena_portal': 'ena', 'ncbi_eutils': 'ncbi'}.get(
        getattr(args[0], 'service', None) if args else None, 'provider')


def provider_boundary(*, provider=None, operation=None):
    """Record exhausted transport/parse failures, never intermediate attempts."""
    def decorate(call):
        @wraps(call)
        def wrapped(*args, **kwargs):
            provider_token = _provider.set(provider(*args, **kwargs) if callable(provider) else provider or _provider.get() or _remote_provider(args, kwargs))
            operation_token = _operation.set(operation(*args, **kwargs) if callable(operation) else operation or _operation.get())
            try:
                return call(*args, **kwargs)
            except Exception as error:
                record(error, _provider.get())
                raise
            finally:
                _operation.reset(operation_token)
                _provider.reset(provider_token)
        return wrapped
    return decorate


@contextmanager
def collect_provider_issues(item):
    token = _issues.set({})
    try:
        yield
    finally:
        from .converters.unified.contracts import Diagnostic
        found = snapshot()
        _issues.reset(token)
        unresolved = []
        for issue, recovered in found.values():
            item.diagnostics.append(Diagnostic(
                'provider_recovered' if recovered else 'source_partial',
                f'{issue.error.provider}: {issue.operation}: {issue.reason}',
                'retrieval', 'info' if recovered else 'warning', provider_issue=issue))
            if not recovered:
                unresolved.append(issue)
        if unresolved:
            if item.execution == 'succeeded':
                item.completeness = item.status = 'partial'
            if item.error is not None and item.error.error_type == 'InputError' and any(
                d.code == 'source_failed' for d in item.diagnostics
            ):
                cause = next((i for i in unresolved if i.error.retry_category is RetryCategory.RETRYABLE), unresolved[0])
                item.error, item.retry_at = cause.error, cause.retry_at
