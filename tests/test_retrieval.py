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
from datetime import datetime, timezone
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


def test_unknown_length_retrieval_reserves_object_budget_before_streaming(
    tmp_path,
) -> None:
    class ObservedResponse(_Response):
        stream_started = False

        def iter_content(self, chunk_size: int):
            self.stream_started = True
            yield b"12345678"

    (tmp_path / "existing.bin").write_bytes(b"x" * 25)
    response = ObservedResponse(headers={"X-Unknown-Length": "true"})
    service = RetrievalService(
        tmp_path,
        policy=_policy(),
        session=_Session([response]),
    )

    with pytest.raises(RetrievalSizeError, match="cache byte limit"):
        service.localize(
            "https://data.example.org/chunked.h5ad",
            max_bytes=8,
        )

    assert response.stream_started is False


def test_streaming_retrieval_scans_cache_and_preflights_disk_once(tmp_path) -> None:
    disk_preflights: list[int] = []

    class CountingService(RetrievalService):
        cache_scans = 0

        def _cache_usage_bytes(self) -> int:
            self.cache_scans += 1
            return super()._cache_usage_bytes()

    policy = _policy(
        disk_preflight=lambda path, required_bytes, headroom_fraction: (
            disk_preflights.append(required_bytes)
        )
    )
    service = CountingService(
        tmp_path,
        policy=policy,
        session=_Session(
            [
                _Response(
                    headers={"X-Unknown-Length": "true"},
                    chunks=(b"a",) * 8,
                )
            ]
        ),
    )

    localized = Path(
        service.localize(
            "https://data.example.org/chunked.h5ad",
            max_bytes=8,
        )
    )

    assert localized.read_bytes() == b"a" * 8
    assert service.cache_scans == 1
    assert disk_preflights == [8]


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


def test_cached_reuse_updates_last_use_without_exposing_origin_query(tmp_path) -> None:
    session = _Session([_Response(chunks=(b"abc",))])
    service = RetrievalService(tmp_path, policy=_policy(), session=session)
    value = "https://data.example.org/data.h5ad?token=canary-secret"
    localized = Path(service.localize(value))
    sidecar_path = localized.with_suffix(localized.suffix + ".metadata.json")
    initial = json.loads(sidecar_path.read_text(encoding="utf-8"))

    service.localize(value)

    reused = json.loads(sidecar_path.read_text(encoding="utf-8"))
    assert reused["last_used_at"] >= initial["last_used_at"]
    assert "canary-secret" not in json.dumps(reused)


def test_cache_retention_is_dry_run_verified_active_safe_and_quarantined(
    tmp_path,
) -> None:
    now = datetime(2026, 8, 10, tzinfo=timezone.utc)

    def asset(name: str, content: bytes, last_used_at: str) -> Path:
        path = tmp_path / name
        path.write_bytes(content)
        sidecar = path.with_suffix(path.suffix + ".metadata.json")
        sidecar.write_text(
            json.dumps(
                {
                    "contract_version": "1.0",
                    "origin": f"https://data.example.org/{name}",
                    "byte_count": len(content),
                    "sha256": hashlib.sha256(content).hexdigest(),
                    "md5": None,
                    "fetched_at": last_used_at,
                    "last_used_at": last_used_at,
                    "resource_profile": "standard",
                }
            ),
            encoding="utf-8",
        )
        return path

    candidate = asset("old.h5ad", b"old", "2026-01-01T00:00:00Z")
    active = asset("active.h5ad", b"active", "2026-01-01T00:00:00Z")
    retained = asset("new.h5ad", b"new", "2026-08-09T00:00:00Z")
    corrupt = asset("corrupt.h5ad", b"corrupt", "2026-01-01T00:00:00Z")
    corrupt.write_bytes(b"tampered")
    service = RetrievalService(tmp_path, policy=_policy(), session=_Session([]))

    dry_run = service.retention_report(
        max_age_seconds=30 * 24 * 60 * 60,
        min_retained_assets=1,
        active_paths=[active],
        now=now,
    )

    by_name = {item["name"]: item for item in dry_run["items"]}
    assert dry_run["dry_run"] is True
    assert by_name["old.h5ad"]["action"] == "candidate"
    assert by_name["active.h5ad"]["reason"] == "active_reference"
    assert by_name["new.h5ad"]["reason"] == "minimum_retained"
    assert by_name["corrupt.h5ad"]["reason"] == "integrity_failed"
    assert all(path.exists() for path in (candidate, active, retained, corrupt))

    applied = service.retention_report(
        max_age_seconds=30 * 24 * 60 * 60,
        min_retained_assets=1,
        active_paths=[active],
        now=now,
        apply=True,
    )
    applied_by_name = {item["name"]: item for item in applied["items"]}

    assert applied["dry_run"] is False
    assert applied_by_name["old.h5ad"]["action"] == "quarantined"
    assert not candidate.exists()
    assert Path(applied_by_name["old.h5ad"]["quarantine_path"]).is_file()
    assert active.exists() and retained.exists() and corrupt.exists()


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
