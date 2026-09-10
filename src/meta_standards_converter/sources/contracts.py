# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Explicit retrieval and telemetry contracts for conversion consumers."""
from dataclasses import dataclass
from typing import Protocol, Any

@dataclass(frozen=True)
class RequestMetrics:
    provider_attempts: int = 0
    retry_count: int = 0
    rate_wait_seconds: float = 0.0

class MetricsProvider(Protocol):
    def metrics(self) -> RequestMetrics: ...

class INSDCClient(Protocol):
    def fetch_sra_xml(self, nrx: str) -> Any: ...
    def fetch_ena_file_report(self, accession: str) -> list: ...

class PubMedClient(Protocol):
    def pubmed_summary(self, pubmed_id: str) -> tuple: ...

def request_metrics(*requesters) -> RequestMetrics:
    """Sum explicit, distinct requesters; never inspect nested collaborators."""
    unique = {id(r): r for r in requesters if r is not None}
    totals = [0, 0, 0.0]
    for r in unique.values():
        m = r.metrics() if callable(getattr(r, 'metrics', None)) else r
        for i, key in enumerate(('provider_attempts', 'retry_count', 'rate_wait_seconds')):
            value = getattr(m, key, 0)
            if isinstance(value, (int, float)):
                totals[i] += value
    return RequestMetrics(int(totals[0]), int(totals[1]), float(totals[2]))


class GEOXMLParser(Protocol):
    """Parse supplied XML without performing retrieval."""
    def parse(self, miniml: str, remove_empty: bool = False) -> Any: ...


class MAGETabSourceResolver(Protocol):
    def resolve(self, source: str, sdrf_sources: list[str] | None = None) -> Any: ...


class PackageLoader(Protocol):
    def load(self, json_path: str) -> Any: ...
