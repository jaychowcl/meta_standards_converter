# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Operation-scoped policy for shared services invoked by the unified facade.

Outside these scopes legacy converter behavior is unchanged. ContextVar tokens
are always restored, including when a provider or an injected service fails.
"""
from contextlib import contextmanager
from contextvars import ContextVar
from ..conversion_logging import measured

_loading = ContextVar("msc_preparation_loading", default=None)
_exporting = ContextVar("msc_preparation_exporting", default=False)


@contextmanager
def loading_source(policy):
    token = _loading.set(policy)
    try:
        yield
    finally:
        _loading.reset(token)


def converter_enrichment_enabled():
    return _loading.get() is None


def source_publications_enabled():
    policy = _loading.get()
    return policy is None or policy.allows("publications")


@contextmanager
def exporting_prepared():
    token = _exporting.set(True)
    try:
        yield
    finally:
        _exporting.reset(token)


def exporter_retrieval_enabled():
    return not _exporting.get()


_memo = ContextVar("msc_preparation_memo", default=None)


@contextmanager
def preparation_session():
    token = _memo.set({})
    try:
        yield
    finally:
        _memo.reset(token)


@measured('copying')
def _candidate_copy(value):
    """Copy mutable shells/evidence, sharing only construction-certified branches."""
    from copy import deepcopy
    from dataclasses import is_dataclass
    from ..miniml.model import _certified_immutable
    memo, visited = {}, set()
    def remember(node):
        if id(node) in visited:
            return
        visited.add(id(node))
        if _certified_immutable(node):
            memo[id(node)] = node
        elif type(node) in (list, tuple):
            for child in node:
                remember(child)
        elif type(node) is dict:
            for child in node.values():
                remember(child)
        elif is_dataclass(node) and not isinstance(node, type):
            for child in vars(node).values():
                remember(child)
    remember(value)
    return deepcopy(value, memo)


def retrieve_once(key, call):
    """Reuse evidence only within one facade invocation; legacy calls are unchanged."""
    from meta_standards_converter.provider_issues import snapshot, restore
    memo = _memo.get()
    if memo is None:
        return call()
    from meta_standards_converter.conversion_logging import count, emit
    if key in memo:
        count("preparation_cache_hits")
        emit("preparation_cache_hit", level="DEBUG")
    if key not in memo:
        before = snapshot()
        try:
            value = call()
            from .retrieval_budget import _active
            if _active.get() and _active.get().exhausted: return value
            memo[key] = (
                True,
                _candidate_copy(value),
                {k: v for k, v in snapshot().items() if before.get(k) != v},
            )
        except Exception as exc:
            from .retrieval_budget import _BudgetExhausted
            if isinstance(exc, _BudgetExhausted): raise
            memo[key] = (False, exc, {k: v for k, v in snapshot().items() if before.get(k) != v})
    ok, value, issues = memo[key]
    restore(issues)
    if not ok:
        raise value
    return _candidate_copy(value)


def convert_source(converter, accession, **options):
    key = ("converter", id(converter), accession, _loading.get(),
           tuple(sorted((k, repr(v)) for k, v in options.items())))
    return retrieve_once(key, lambda: converter.convert(accession, **options))
