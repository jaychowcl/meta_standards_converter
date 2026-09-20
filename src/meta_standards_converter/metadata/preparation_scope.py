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
