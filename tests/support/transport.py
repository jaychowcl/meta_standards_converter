# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Shared transport builders and deterministic collaborators."""

from __future__ import annotations

from unittest.mock import Mock
from meta_standards_converter.helpers.request_helper import RateLimitedRequester, RequestSettings
from meta_standards_converter.retrieval import RetrievalPolicy
from meta_standards_converter.runtime_contracts import get_resource_profile
from meta_standards_converter.converters.sra2json import SRA2JSONConverter
from meta_standards_converter.converters.ena2json import ENA2JSONConverter
from meta_standards_converter.sources.archive_support import Resolution, StudySeed
from tests.support.archive import fixture_records


def response(chunks):
    def stream(**kwargs):
        for item in chunks:
            if isinstance(item, BaseException):
                raise item
            yield item
    return Mock(status_code=200, headers={}, iter_content=stream)


def client(get, retries=3):
    return RateLimitedRequester('ena_portal', RequestSettings(request_delay=0, max_retries=retries),
        get=get, host_gate=Mock(wait=Mock(return_value=0)), random_value=lambda: 0)


class _Response:
    def __init__(
        self,
        *,
        status_code: int = 200,
        headers: dict[str, str] | None = None,
        chunks: tuple[bytes, ...] = (b"abc",),
    ) -> None:
        self.status_code = status_code
        self.headers = headers or {"Content-Length": str(sum(map(len, chunks)))}
        self._chunks = chunks
        self.closed = False

    def iter_content(self, chunk_size: int):
        yield from self._chunks

    def raise_for_status(self) -> None:
        if self.status_code >= 400:
            raise RuntimeError(f"HTTP {self.status_code}")

    def close(self) -> None:
        self.closed = True


class _Session:
    def __init__(self, responses: list[_Response]) -> None:
        self.responses = list(responses)
        self.calls: list[tuple[str, dict]] = []

    def get(self, url: str, **kwargs):
        self.calls.append((url, kwargs))
        return self.responses.pop(0)


def _resolver(host: str, port: int, *args, **kwargs):
    return [(2, 1, 6, "", ("93.184.216.34", port))]


def _policy(**overrides) -> RetrievalPolicy:
    values = {
        "resource_profile": get_resource_profile(
            "standard",
            overrides={
                "max_matrix_bytes": 8,
                "max_aggregate_download_bytes": 12,
                "max_cache_bytes": 32,
            },
        ),
        "allowed_hosts": frozenset({"data.example.org"}),
        "resolver": _resolver,
        "disk_preflight": lambda path, required_bytes, headroom_fraction: None,
    }
    values.update(overrides)
    return RetrievalPolicy(**values)


class HTTP:
    def __init__(self, handler):
        self.handler = handler
        self.calls = []
    def get(self, url, params=None, fmt='xml'):
        self.calls.append((url, params or {}, fmt))
        return self.handler(url, params or {}, fmt)


class Source:
    def __init__(self, provider='sra', partial=False):
        self.provider, self.partial = provider, partial
        self.calls = []
    def resolve(self, accession):
        primary = 'SRP250911' if self.provider == 'sra' else 'PRJNA609050'
        return Resolution([StudySeed('SRP250911', primary), StudySeed('SRP2', 'SRP2')])
    def fetch(self, seed):
        self.calls.append(seed.study)
        if seed.study == 'SRP2':
            raise OSError('failed study')
        records = fixture_records(self.provider)
        if self.partial:
            records.issues.append('sample retrieval incomplete')
        return records
