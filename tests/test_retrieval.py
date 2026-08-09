# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from meta_standards_converter.retrieval import (
    CacheIntegrityError,
    RetrievalPolicy,
    RetrievalSecurityError,
    RetrievalService,
    RetrievalSizeError,
)
from meta_standards_converter.runtime_contracts import get_resource_profile


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


def test_retrieval_rejects_file_urls_and_unapproved_hosts(tmp_path) -> None:
    service = RetrievalService(tmp_path, policy=_policy(), session=_Session([]))

    with pytest.raises(RetrievalSecurityError, match="file URLs are disabled"):
        service.localize("file:///etc/passwd")
    with pytest.raises(RetrievalSecurityError, match="host is not allowed"):
        service.localize("https://attacker.example/asset.h5ad")


def test_retrieval_rejects_private_or_mixed_dns_answers(tmp_path) -> None:
    def private_resolver(host: str, port: int, *args, **kwargs):
        return [(2, 1, 6, "", ("127.0.0.1", port))]

    service = RetrievalService(
        tmp_path,
        policy=_policy(resolver=private_resolver),
        session=_Session([]),
    )

    with pytest.raises(RetrievalSecurityError, match="non-public address"):
        service.localize("https://data.example.org/asset.h5ad")


def test_retrieval_revalidates_redirect_targets(tmp_path) -> None:
    session = _Session(
        [
            _Response(
                status_code=302,
                headers={"Location": "https://127.0.0.1/private"},
                chunks=(),
            )
        ]
    )
    service = RetrievalService(tmp_path, policy=_policy(), session=session)

    with pytest.raises(RetrievalSecurityError, match="host is not allowed"):
        service.localize("https://data.example.org/asset.h5ad")

    assert len(session.calls) == 1
    assert session.calls[0][1]["allow_redirects"] is False


def test_retrieval_enforces_declared_and_streamed_object_limits(tmp_path) -> None:
    declared = _Session([_Response(headers={"Content-Length": "9"})])
    with pytest.raises(RetrievalSizeError, match="exceeds the 8 byte object limit"):
        RetrievalService(tmp_path, policy=_policy(), session=declared).localize(
            "https://data.example.org/declared.h5ad"
        )

    streamed = _Session(
        [_Response(headers={}, chunks=(b"1234", b"5678", b"9"))]
    )
    with pytest.raises(RetrievalSizeError, match="exceeds the 8 byte object limit"):
        RetrievalService(tmp_path, policy=_policy(), session=streamed).localize(
            "https://data.example.org/chunked.h5ad"
        )

    assert not list(tmp_path.glob("*.stage"))


def test_retrieval_uses_connect_read_timeouts_and_integrity_sidecar(tmp_path) -> None:
    session = _Session([_Response(chunks=(b"abc", b"123"))])
    service = RetrievalService(tmp_path, policy=_policy(), session=session)
    signed_url = "https://data.example.org/data.h5ad?token=canary-secret"

    localized = Path(
        service.localize(
            signed_url,
            md5="e99a18c428cb38d5f260853678922e03",
        )
    )

    assert localized.read_bytes() == b"abc123"
    assert session.calls[0][1]["timeout"] == (10, 60)
    sidecar = json.loads(
        localized.with_suffix(localized.suffix + ".metadata.json").read_text(
            encoding="utf-8"
        )
    )
    assert sidecar["origin"] == "https://data.example.org/data.h5ad"
    assert "canary-secret" not in json.dumps(sidecar)
    assert sidecar["sha256"] == hashlib.sha256(b"abc123").hexdigest()
    assert sidecar["byte_count"] == 6


def test_retrieval_validates_cached_bytes_before_reuse(tmp_path) -> None:
    session = _Session([_Response(chunks=(b"abc",))])
    service = RetrievalService(tmp_path, policy=_policy(), session=session)
    value = "https://data.example.org/data.h5ad"
    localized = Path(service.localize(value))
    localized.write_bytes(b"tampered")

    with pytest.raises(CacheIntegrityError, match="SHA-256"):
        service.localize(value)

    assert len(session.calls) == 1


def test_retrieval_enforces_aggregate_run_limit(tmp_path) -> None:
    session = _Session(
        [
            _Response(chunks=(b"12345678",)),
            _Response(chunks=(b"12345",)),
        ]
    )
    service = RetrievalService(tmp_path, policy=_policy(), session=session)

    service.localize("https://data.example.org/one.h5ad")
    with pytest.raises(RetrievalSizeError, match="aggregate download limit"):
        service.localize("https://data.example.org/two.h5ad")


def test_retrieval_uses_bounded_ncbi_range_fallback_under_same_policy(tmp_path) -> None:
    session = _Session(
        [
            _Response(status_code=403, headers={}, chunks=()),
            _Response(
                status_code=206,
                headers={"Content-Range": "bytes 0-2/6"},
                chunks=(b"abc",),
            ),
            _Response(
                status_code=206,
                headers={"Content-Range": "bytes 3-5/6"},
                chunks=(b"123",),
            ),
        ]
    )
    policy = _policy(
        allowed_hosts=frozenset({"ftp.ncbi.nlm.nih.gov"}),
    )
    service = RetrievalService(tmp_path, policy=policy, session=session)

    localized = Path(
        service.localize(
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE1/matrix.h5ad"
        )
    )

    assert localized.read_bytes() == b"abc123"
    assert [call[1].get("headers", {}).get("Range") for call in session.calls] == [
        None,
        "bytes=0-7",
        "bytes=3-5",
    ]
    assert all(call[1]["allow_redirects"] is False for call in session.calls)


def test_ncbi_range_fallback_closes_first_response_when_size_is_rejected(
    tmp_path,
) -> None:
    rejected = _Response(
        status_code=206,
        headers={"Content-Range": "bytes 0-7/9"},
        chunks=(b"12345678",),
    )
    service = RetrievalService(
        tmp_path,
        policy=_policy(allowed_hosts=frozenset({"ftp.ncbi.nlm.nih.gov"})),
        session=_Session(
            [_Response(status_code=403, headers={}, chunks=()), rejected]
        ),
    )

    with pytest.raises(RetrievalSizeError, match="object limit"):
        service.localize(
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE1/matrix.h5ad"
        )

    assert rejected.closed is True
