# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Explicit per-conversion evidence cache, including failures and provider issues."""
from ..conversion_logging import measured

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



class PreparationSession:
    def __init__(self):
        self.memo = {}

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.memo.clear()

    def retrieve(self, key, call):
        from meta_standards_converter.provider_issues import snapshot, restore
        memo = self.memo
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

    def load(self, source, accession, **options):
        key = ("source", id(source), accession,
               tuple(sorted((k, repr(v)) for k, v in options.items())))
        return self.retrieve(key, lambda: source.load(accession, **options))


def retrieve_once(key, call, *, session=None):
    return call() if session is None else session.retrieve(key, call)
