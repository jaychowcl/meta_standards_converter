# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Optional metadata limits shared by nested services within one facade input."""
from contextlib import contextmanager
from contextvars import ContextVar
from copy import deepcopy
import time
from threading import RLock

_active = ContextVar('msc_metadata_budget', default=None)
_phase = ContextVar('msc_metadata_retrieval_phase', default='discovery')
_cancelled = ContextVar('msc_retrieval_cancelled', default=None)


class _RetrievalCancelled(BaseException):
    """Internal worker shutdown; never swallow as an optional provider failure."""


def _check_cancelled():
    event = _cancelled.get()
    if event is not None and event.is_set():
        raise _RetrievalCancelled()


class _BudgetExhausted(RuntimeError):
    pass


class _MetadataBudget:
    def __init__(self, timeout_seconds, max_requests, clock):
        self.clock = clock
        self.deadline = None if timeout_seconds is None else clock() + timeout_seconds
        self.max_requests = max_requests
        self.requests = 0
        self.events = []
        self.lock = RLock()

    @property
    def exhausted(self):
        return bool(self.events)

    def remaining(self):
        self.check()
        return None if self.deadline is None else max(0, self.deadline - self.clock())

    def stop(self, reason):
        from ..conversion_logging import count, emit
        event = {'stage': _phase.get(), 'reason': reason, 'requests': self.requests}
        if event not in self.events:
            self.events.append(event)
            count('metadata_budget_exhaustions')
            emit('metadata_budget_exhausted', level='WARNING', **event)
        raise _BudgetExhausted(f'metadata budget exhausted during {_phase.get()}: {reason}; {self.requests} HTTP attempts')

    def check(self, *, request=False):
        with self.lock:
            if self.exhausted: self.stop(self.events[0]['reason'])
            if self.deadline is not None and self.clock() >= self.deadline: self.stop('deadline reached')
            if request and self.max_requests is not None and self.requests >= self.max_requests: self.stop('request allowance reached')

    def attempt(self, timeout):
        with self.lock:
            self.check(request=True)
            remaining = self.remaining()
            self.requests += 1
        from ..conversion_logging import count
        count('metadata_http_attempts')
        if remaining is None: return timeout
        from urllib3.util import Timeout
        if isinstance(timeout, Timeout):
            connect, read = timeout.connect_timeout, timeout.read_timeout
        elif isinstance(timeout, (list, tuple)): connect, read = timeout
        else: connect = read = timeout
        return Timeout(total=remaining, connect=min(connect or remaining, remaining), read=min(read or remaining, remaining))


@contextmanager
def _budget_scope(timeout_seconds=None, max_requests=None, *, clock=time.monotonic):
    current = _active.get()
    if current is not None:
        yield current
        return
    budget = None if timeout_seconds is None and max_requests is None else _MetadataBudget(timeout_seconds, max_requests, clock)
    token = _active.set(budget)
    try: yield budget
    finally: _active.reset(token)


@contextmanager
def _budget_stage(name):
    token = _phase.set(name)
    try: yield
    finally: _phase.reset(token)


@contextmanager
def _without_budget():
    token = _active.set(None)
    try: yield
    finally: _active.reset(token)


def _coverage(package, budget):
    if budget is None or not budget.events: return package
    from ..miniml import MINiMLCodec
    data = package.to_mapping()
    records = data.setdefault('extensions', {}).setdefault('insdc', {'version': '1.0', 'records': []})['records']
    record = {'provider': data.get('source', {}).get('format', 'MSC'), 'kind': 'retrieval_coverage',
              'accession': data['series']['iid'], 'metadata': {'complete': False, 'unfinished': deepcopy(budget.events)}}
    if record not in records: records.append(record)
    return MINiMLCodec().decode(data, strict=True).package


def _coverage_issues(package):
    for record in package.to_mapping().get('extensions', {}).get('insdc', {}).get('records', []):
        if record.get('kind') == 'retrieval_coverage' and isinstance(record.get('metadata'), dict) and record['metadata'].get('complete') is False:
            yield f"recorded metadata budget coverage incomplete for {record.get('accession')}: {record['metadata'].get('unfinished', [])}"
