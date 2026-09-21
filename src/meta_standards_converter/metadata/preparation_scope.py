"""Operation-scoped policy for shared services invoked by the unified facade.

Outside these scopes legacy converter behavior is unchanged. ContextVar tokens
are always restored, including when a provider or an injected service fails.
"""
from contextlib import contextmanager
from contextvars import ContextVar

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


def retrieve_once(key, call):
    """Reuse evidence only within one facade invocation; legacy calls are unchanged."""
    from copy import deepcopy
    memo = _memo.get()
    if memo is None:
        return call()
    if key not in memo:
        try:
            memo[key] = (True, deepcopy(call()))
        except Exception as exc:
            memo[key] = (False, exc)
    ok, value = memo[key]
    if not ok:
        raise value
    return deepcopy(value)


def convert_source(converter, accession, **options):
    key = ("converter", id(converter), accession, _loading.get(),
           tuple(sorted((k, repr(v)) for k, v in options.items())))
    return retrieve_once(key, lambda: converter.convert(accession, **options))
