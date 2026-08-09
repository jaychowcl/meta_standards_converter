# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================

"""Policy-enforced, integrity-checked retrieval for large curation assets."""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass, field
from datetime import datetime, timezone
import hashlib
import ipaddress
import json
import os
from pathlib import Path
import re
import socket
import tempfile
import threading
from typing import Any, Callable
import urllib.request
from urllib.parse import urljoin, urlparse, urlsplit, urlunsplit

import requests

try:  # POSIX cache coordination; non-POSIX retains process-local safety.
    import fcntl
except ImportError:  # pragma: no cover - exercised on non-POSIX platforms.
    fcntl = None

from meta_standards_converter.runtime_contracts import (
    ResourceProfile,
    get_resource_profile,
    require_disk_headroom,
)


DEFAULT_PROVIDER_HOST_SUFFIXES = frozenset(
    {
        "ncbi.nlm.nih.gov",
        "ebi.ac.uk",
        "biostudies.org",
    }
)
_REDIRECT_STATUSES = frozenset({301, 302, 303, 307, 308})
_CONTENT_RANGE = re.compile(r"^bytes (\d+)-(\d+)/(\d+)$")
_RANGE_CHUNK_BYTES = 16 * 1024 * 1024
_FALLBACK_CACHE_LOCK = threading.RLock()


class RetrievalError(RuntimeError):
    """Base class for policy-enforced retrieval failures."""


class RetrievalSecurityError(RetrievalError):
    """A URL violated scheme, host, address, credential, or redirect policy."""


class RetrievalSizeError(RetrievalError):
    """An object, aggregate run, cache, or disk ceiling would be exceeded."""


class CacheIntegrityError(RetrievalError):
    """A cached file is missing or fails its immutable integrity metadata."""


@dataclass(frozen=True)
class RetrievalPolicy:
    resource_profile: ResourceProfile = field(default_factory=get_resource_profile)
    allowed_hosts: frozenset[str] = frozenset()
    allowed_host_suffixes: frozenset[str] = DEFAULT_PROVIDER_HOST_SUFFIXES
    allowed_schemes: frozenset[str] = frozenset({"https", "ftp"})
    allow_file_urls: bool = False
    ranged_fallback_hosts: frozenset[str] = frozenset(
        {"ftp.ncbi.nlm.nih.gov"}
    )
    resolver: Callable[..., list[Any]] = socket.getaddrinfo
    disk_preflight: Callable[..., Any] = require_disk_headroom

    def validate_url(self, value: str) -> None:
        parsed = urlsplit(value)
        scheme = parsed.scheme.casefold()
        if scheme == "file" and not self.allow_file_urls:
            raise RetrievalSecurityError("file URLs are disabled by retrieval policy.")
        if scheme not in self.allowed_schemes:
            raise RetrievalSecurityError(
                f"URL scheme {scheme!r} is not allowed by retrieval policy."
            )
        if parsed.username is not None or parsed.password is not None:
            raise RetrievalSecurityError("URL userinfo credentials are not allowed.")
        host = (parsed.hostname or "").casefold().rstrip(".")
        if not host:
            raise RetrievalSecurityError("Remote asset URL requires a hostname.")
        if not self._host_allowed(host):
            raise RetrievalSecurityError(f"Remote asset host is not allowed: {host}")
        port = parsed.port or (443 if scheme == "https" else 21)
        try:
            answers = self.resolver(host, port, type=socket.SOCK_STREAM)
        except OSError as error:
            raise RetrievalSecurityError(
                f"Remote asset host could not be resolved: {host}"
            ) from error
        addresses = {
            str(answer[4][0]).split("%", 1)[0]
            for answer in answers
            if len(answer) >= 5 and answer[4]
        }
        if not addresses:
            raise RetrievalSecurityError(
                f"Remote asset host has no resolved addresses: {host}"
            )
        for address in addresses:
            try:
                parsed_address = ipaddress.ip_address(address)
            except ValueError as error:
                raise RetrievalSecurityError(
                    f"Remote asset host returned an invalid address: {host}"
                ) from error
            if not parsed_address.is_global:
                raise RetrievalSecurityError(
                    f"Remote asset host resolves to a non-public address: {address}"
                )

    def _host_allowed(self, host: str) -> bool:
        if host in {value.casefold().rstrip(".") for value in self.allowed_hosts}:
            return True
        return any(
            host == suffix or host.endswith(f".{suffix}")
            for suffix in (
                value.casefold().lstrip(".").rstrip(".")
                for value in self.allowed_host_suffixes
            )
        )


class RetrievalService:
    """Retrieve remote assets without crossing configured trust/resource bounds."""

    def __init__(
        self,
        cache_dir: str | Path,
        *,
        policy: RetrievalPolicy | None = None,
        session: Any | None = None,
        urlopen: Callable[..., Any] | None = None,
    ) -> None:
        self.cache_dir = Path(cache_dir)
        self.policy = policy or RetrievalPolicy()
        self.session = session or requests.Session()
        self.urlopen = urlopen or urllib.request.urlopen
        self._downloaded_bytes = 0

    def localize(
        self,
        value: str,
        *,
        md5: str | None = None,
        max_bytes: int | None = None,
    ) -> str:
        parsed = urlparse(value)
        if parsed.scheme == "":
            return value
        if parsed.scheme == "file":
            self.policy.validate_url(value)
            return parsed.path
        self.policy.validate_url(value)
        object_limit = max_bytes or self.policy.resource_profile.max_matrix_bytes
        if object_limit <= 0:
            raise ValueError("max_bytes must be positive")
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        destination = self._destination(value)
        with self._cache_lock():
            if destination.exists():
                self._verify_cached(destination, origin=value, expected_md5=md5)
                return str(destination)
            if parsed.scheme in {"http", "https"}:
                self._download_http(
                    value,
                    destination=destination,
                    object_limit=object_limit,
                    expected_md5=md5,
                )
            else:
                self._download_ftp(
                    value,
                    destination=destination,
                    object_limit=object_limit,
                    expected_md5=md5,
                )
        return str(destination)

    def _download_http(
        self,
        value: str,
        *,
        destination: Path,
        object_limit: int,
        expected_md5: str | None,
    ) -> None:
        current = value
        response = None
        for redirect_count in range(self.policy.resource_profile.max_redirects + 1):
            self.policy.validate_url(current)
            response = self.session.get(
                current,
                stream=True,
                allow_redirects=False,
                timeout=(
                    self.policy.resource_profile.connect_timeout_seconds,
                    self.policy.resource_profile.read_timeout_seconds,
                ),
            )
            if response.status_code not in _REDIRECT_STATUSES:
                break
            location = response.headers.get("Location")
            response.close()
            if not location:
                raise RetrievalSecurityError("Redirect response omitted Location.")
            if redirect_count >= self.policy.resource_profile.max_redirects:
                raise RetrievalSecurityError("Remote asset exceeded redirect limit.")
            target = urljoin(current, location)
            if urlsplit(target).scheme.casefold() != urlsplit(current).scheme.casefold():
                raise RetrievalSecurityError("Cross-scheme redirects are not allowed.")
            self.policy.validate_url(target)
            current = target
        if response is None:
            raise RetrievalError("Remote asset request did not produce a response.")
        parsed_current = urlsplit(current)
        if (
            response.status_code == 403
            and parsed_current.scheme.casefold() == "https"
            and (parsed_current.hostname or "").casefold().rstrip(".")
            in {
                host.casefold().rstrip(".")
                for host in self.policy.ranged_fallback_hosts
            }
        ):
            response.close()
            self._download_http_ranges(
                current,
                destination=destination,
                object_limit=object_limit,
                expected_md5=expected_md5,
                origin=value,
            )
            return
        try:
            response.raise_for_status()
            declared = _content_length(response.headers)
            self._preflight_declared(declared, object_limit=object_limit)
            self._write_stream(
                response.iter_content(chunk_size=1024 * 1024),
                destination=destination,
                origin=value,
                object_limit=object_limit,
                declared_bytes=declared,
                expected_md5=expected_md5,
            )
        finally:
            response.close()

    def _download_http_ranges(
        self,
        value: str,
        *,
        destination: Path,
        object_limit: int,
        expected_md5: str | None,
        origin: str,
    ) -> None:
        first = self._range_response(
            value,
            start=0,
            end=min(_RANGE_CHUNK_BYTES, object_limit) - 1,
        )
        try:
            first_start, _first_end, total = _content_range(first.headers)
            if first_start != 0:
                raise RetrievalSizeError(
                    "NCBI ranged asset response is not a valid initial range."
                )
            self._preflight_declared(total, object_limit=object_limit)
        except BaseException:
            first.close()
            raise

        def chunks():
            response = first
            expected_start = 0
            expected_total = total
            while expected_start < expected_total:
                returned_start, returned_end, returned_total = _content_range(
                    response.headers
                )
                if (
                    returned_start != expected_start
                    or returned_end < returned_start
                    or returned_end >= returned_total
                    or returned_total != expected_total
                ):
                    response.close()
                    raise RetrievalSizeError(
                        "NCBI ranged asset response is not contiguous."
                    )
                received = 0
                try:
                    for chunk in response.iter_content(chunk_size=1024 * 1024):
                        if chunk:
                            received += len(chunk)
                            yield chunk
                finally:
                    response.close()
                expected_bytes = returned_end - returned_start + 1
                if received != expected_bytes:
                    raise RetrievalSizeError(
                        "NCBI ranged asset body does not match Content-Range."
                    )
                expected_start = returned_end + 1
                if expected_start < expected_total:
                    response = self._range_response(
                        value,
                        start=expected_start,
                        end=min(
                            expected_total - 1,
                            expected_start + _RANGE_CHUNK_BYTES - 1,
                        ),
                    )

        self._write_stream(
            chunks(),
            destination=destination,
            origin=origin,
            object_limit=object_limit,
            declared_bytes=total,
            expected_md5=expected_md5,
        )

    def _range_response(self, value: str, *, start: int, end: int) -> Any:
        self.policy.validate_url(value)
        response = self.session.get(
            value,
            headers={"Range": f"bytes={start}-{end}"},
            stream=True,
            allow_redirects=False,
            timeout=(
                self.policy.resource_profile.connect_timeout_seconds,
                self.policy.resource_profile.read_timeout_seconds,
            ),
        )
        if response.status_code != 206:
            response.close()
            raise RetrievalError(
                "NCBI ranged asset request did not return partial content."
            )
        return response

    def _download_ftp(
        self,
        value: str,
        *,
        destination: Path,
        object_limit: int,
        expected_md5: str | None,
    ) -> None:
        self.policy.validate_url(value)
        with self.urlopen(
            value,
            timeout=self.policy.resource_profile.read_timeout_seconds,
        ) as response:
            declared = _content_length(getattr(response, "headers", {}))
            self._preflight_declared(declared, object_limit=object_limit)

            def chunks():
                while True:
                    chunk = response.read(1024 * 1024)
                    if not chunk:
                        break
                    yield chunk

            self._write_stream(
                chunks(),
                destination=destination,
                origin=value,
                object_limit=object_limit,
                declared_bytes=declared,
                expected_md5=expected_md5,
            )

    def _preflight_declared(self, declared: int | None, *, object_limit: int) -> None:
        if declared is None:
            return
        if declared > object_limit:
            raise RetrievalSizeError(
                f"Remote asset declares {declared} bytes and exceeds the "
                f"{object_limit} byte object limit."
            )
        self._check_aggregate(declared)

    def _write_stream(
        self,
        chunks,
        *,
        destination: Path,
        origin: str,
        object_limit: int,
        declared_bytes: int | None,
        expected_md5: str | None,
    ) -> None:
        sha256 = hashlib.sha256()
        md5_digest = hashlib.md5(usedforsecurity=False)
        byte_count = 0
        reserved_bytes = (
            declared_bytes if declared_bytes is not None else object_limit
        )
        self._check_aggregate(reserved_bytes)
        self._check_cache_quota(reserved_bytes)
        self._disk_preflight(reserved_bytes)
        with tempfile.NamedTemporaryFile(
            dir=self.cache_dir,
            prefix=f".{destination.name}.",
            suffix=".stage",
            delete=False,
        ) as handle:
            temporary = Path(handle.name)
            try:
                for chunk in chunks:
                    if not chunk:
                        continue
                    byte_count += len(chunk)
                    self._downloaded_bytes += len(chunk)
                    if byte_count > object_limit:
                        raise RetrievalSizeError(
                            f"Remote asset exceeds the {object_limit} byte object limit."
                        )
                    if (
                        self._downloaded_bytes
                        > self.policy.resource_profile.max_aggregate_download_bytes
                    ):
                        raise RetrievalSizeError(
                            "Remote assets exceed the aggregate download limit."
                        )
                    handle.write(chunk)
                    sha256.update(chunk)
                    md5_digest.update(chunk)
                if declared_bytes is not None and byte_count != declared_bytes:
                    raise RetrievalSizeError(
                        f"Remote asset body has {byte_count} bytes; "
                        f"Content-Length declared {declared_bytes}."
                    )
                if expected_md5 and md5_digest.hexdigest().casefold() != expected_md5.casefold():
                    raise CacheIntegrityError("MD5 checksum mismatch for remote asset.")
                handle.flush()
                os.fsync(handle.fileno())
                os.replace(temporary, destination)
                self._write_sidecar(
                    destination,
                    origin=origin,
                    byte_count=byte_count,
                    sha256=sha256.hexdigest(),
                    md5=md5_digest.hexdigest() if expected_md5 else None,
                )
            finally:
                temporary.unlink(missing_ok=True)

    def _destination(self, value: str) -> Path:
        parsed = urlsplit(value)
        basename = Path(parsed.path).name or "asset"
        safe_basename = "".join(
            character if character.isalnum() or character in "._-" else "_"
            for character in basename
        )
        prefix = hashlib.sha256(value.encode("utf-8")).hexdigest()[:12]
        return self.cache_dir / f"{prefix}-{safe_basename}"

    def _verify_cached(
        self,
        destination: Path,
        *,
        origin: str,
        expected_md5: str | None,
    ) -> None:
        sidecar_path = self._sidecar_path(destination)
        if not sidecar_path.is_file():
            raise CacheIntegrityError(
                f"Cached asset {destination.name} has no integrity sidecar."
            )
        try:
            metadata = json.loads(sidecar_path.read_text(encoding="utf-8"))
        except (OSError, ValueError) as error:
            raise CacheIntegrityError(
                f"Cached asset {destination.name} has an invalid integrity sidecar."
            ) from error
        if metadata.get("origin") != _sanitize_origin(origin):
            raise CacheIntegrityError("Cached asset origin does not match its request.")
        digest = _file_digest(destination, "sha256")
        if digest != metadata.get("sha256"):
            raise CacheIntegrityError("Cached asset SHA-256 integrity check failed.")
        if destination.stat().st_size != metadata.get("byte_count"):
            raise CacheIntegrityError("Cached asset byte count integrity check failed.")
        if expected_md5 and _file_digest(destination, "md5") != expected_md5.casefold():
            raise CacheIntegrityError("Cached asset MD5 integrity check failed.")
        metadata["last_used_at"] = self._utc_now()
        self._write_sidecar_payload(sidecar_path, metadata)

    def _write_sidecar(
        self,
        destination: Path,
        *,
        origin: str,
        byte_count: int,
        sha256: str,
        md5: str | None,
    ) -> None:
        sidecar = self._sidecar_path(destination)
        timestamp = self._utc_now()
        payload = {
            "contract_version": "1.0",
            "origin": _sanitize_origin(origin),
            "byte_count": byte_count,
            "sha256": sha256,
            "md5": md5,
            "fetched_at": timestamp,
            "last_used_at": timestamp,
            "resource_profile": self.policy.resource_profile.name,
        }
        self._write_sidecar_payload(sidecar, payload)

    @staticmethod
    def _write_sidecar_payload(sidecar: Path, payload: dict) -> None:
        temporary = sidecar.with_name(f".{sidecar.name}.{os.getpid()}.stage")
        try:
            temporary.write_text(
                json.dumps(payload, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
            os.replace(temporary, sidecar)
        finally:
            temporary.unlink(missing_ok=True)

    @staticmethod
    def _utc_now() -> str:
        return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")

    def retention_report(
        self,
        *,
        max_age_seconds: float,
        min_retained_assets: int = 1,
        active_paths=(),
        now: datetime | None = None,
        apply: bool = False,
    ) -> dict:
        """Plan or quarantine old verified assets; never delete cache data."""

        if max_age_seconds < 0:
            raise ValueError("max_age_seconds must be non-negative")
        if min_retained_assets < 0:
            raise ValueError("min_retained_assets must be non-negative")
        current = now or datetime.now(timezone.utc)
        if current.tzinfo is None:
            raise ValueError("now must be timezone-aware")
        active = {
            Path(path).resolve(strict=False)
            for path in active_paths
        }
        with self._cache_lock():
            entries = [
                self._retention_entry(path)
                for path in sorted(
                    (
                        path
                        for path in self.cache_dir.iterdir()
                        if path.is_file()
                        and not path.name.startswith(".")
                        and not path.name.endswith(".metadata.json")
                    ),
                    key=lambda path: path.name,
                )
            ]
            valid_by_recency = sorted(
                (entry for entry in entries if entry["integrity"] == "verified"),
                key=lambda entry: entry["last_used_timestamp"],
                reverse=True,
            )
            minimum_retained = {
                entry["path"] for entry in valid_by_recency[:min_retained_assets]
            }
            quarantine_dir = None
            if apply:
                stamp = current.astimezone(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
                quarantine_dir = self.cache_dir / ".quarantine" / stamp

            items = []
            for entry in entries:
                path = Path(entry.pop("path"))
                sidecar = self._sidecar_path(path)
                item = {
                    **entry,
                    "name": path.name,
                    "path": str(path.resolve(strict=False)),
                }
                item.pop("last_used_timestamp", None)
                if entry["integrity"] != "verified":
                    item.update(action="retained", reason="integrity_failed")
                elif path.resolve(strict=False) in active:
                    item.update(action="retained", reason="active_reference")
                elif str(path) in minimum_retained:
                    item.update(action="retained", reason="minimum_retained")
                else:
                    age_seconds = max(
                        0.0,
                        current.timestamp() - entry["last_used_timestamp"],
                    )
                    if age_seconds < max_age_seconds:
                        item.update(action="retained", reason="recently_used")
                    elif not apply:
                        item.update(action="candidate", reason="retention_age")
                    else:
                        quarantine_dir.mkdir(parents=True, exist_ok=True)
                        quarantined = quarantine_dir / path.name
                        quarantined_sidecar = self._sidecar_path(quarantined)
                        os.replace(path, quarantined)
                        try:
                            os.replace(sidecar, quarantined_sidecar)
                        except BaseException:
                            os.replace(quarantined, path)
                            raise
                        item.update(
                            action="quarantined",
                            reason="retention_age",
                            quarantine_path=str(quarantined.resolve()),
                            recovery_path=str(path.resolve(strict=False)),
                        )
                items.append(item)

        return {
            "contract_version": "1.0",
            "dry_run": not apply,
            "max_age_seconds": max_age_seconds,
            "min_retained_assets": min_retained_assets,
            "active_paths": sorted(str(path) for path in active),
            "quarantine_dir": (
                str(quarantine_dir.resolve()) if quarantine_dir is not None else None
            ),
            "items": items,
        }

    def _retention_entry(self, path: Path) -> dict:
        sidecar = self._sidecar_path(path)
        metadata = None
        integrity = "verified"
        try:
            metadata = json.loads(sidecar.read_text(encoding="utf-8"))
            if metadata.get("contract_version") != "1.0":
                raise ValueError("unsupported sidecar contract")
            if path.stat().st_size != metadata.get("byte_count"):
                raise ValueError("byte count mismatch")
            if _file_digest(path, "sha256") != metadata.get("sha256"):
                raise ValueError("sha256 mismatch")
            last_used = self._parse_cache_timestamp(
                metadata.get("last_used_at") or metadata.get("fetched_at")
            )
        except (OSError, TypeError, ValueError, json.JSONDecodeError):
            integrity = "failed"
            last_used = datetime.fromtimestamp(
                path.stat().st_mtime,
                tz=timezone.utc,
            )
        return {
            "path": str(path),
            "origin": metadata.get("origin") if isinstance(metadata, dict) else None,
            "byte_count": path.stat().st_size,
            "sha256": metadata.get("sha256") if isinstance(metadata, dict) else None,
            "last_used_at": last_used.isoformat().replace("+00:00", "Z"),
            "last_used_timestamp": last_used.timestamp(),
            "integrity": integrity,
        }

    @staticmethod
    def _parse_cache_timestamp(value) -> datetime:
        if not isinstance(value, str) or not value:
            raise ValueError("cache timestamp is missing")
        parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
        if parsed.tzinfo is None:
            raise ValueError("cache timestamp must include timezone")
        return parsed.astimezone(timezone.utc)

    def _check_aggregate(self, additional_bytes: int) -> None:
        if (
            self._downloaded_bytes + additional_bytes
            > self.policy.resource_profile.max_aggregate_download_bytes
        ):
            raise RetrievalSizeError("Remote assets exceed the aggregate download limit.")

    def _check_cache_quota(self, additional_bytes: int) -> None:
        cache_bytes = self._cache_usage_bytes()
        if cache_bytes + additional_bytes > self.policy.resource_profile.max_cache_bytes:
            raise RetrievalSizeError("Remote asset exceeds the cache byte limit.")

    def _cache_usage_bytes(self) -> int:
        """Snapshot committed asset bytes once for a locked reservation."""

        return sum(
            path.stat().st_size
            for path in self.cache_dir.iterdir()
            if path.is_file()
            and not path.name.startswith(".")
            and not path.name.endswith(".metadata.json")
        )

    @contextmanager
    def _cache_lock(self):
        """Serialize cache verification, reservation, and publication."""

        lock_path = self.cache_dir / ".retrieval-cache.lock"
        with lock_path.open("a+b") as handle:
            if fcntl is None:
                with _FALLBACK_CACHE_LOCK:
                    yield
                return
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
            try:
                yield
            finally:
                fcntl.flock(handle.fileno(), fcntl.LOCK_UN)

    def _disk_preflight(self, required_bytes: int) -> None:
        self.policy.disk_preflight(
            self.cache_dir,
            required_bytes=required_bytes,
            headroom_fraction=self.policy.resource_profile.disk_headroom_fraction,
        )

    @staticmethod
    def _sidecar_path(destination: Path) -> Path:
        return destination.with_suffix(destination.suffix + ".metadata.json")


class AssetDownloader(RetrievalService):
    """Compatibility facade for the former JSON-to-H5AD downloader."""


def _content_length(headers: Any) -> int | None:
    value = headers.get("Content-Length") if hasattr(headers, "get") else None
    if value in (None, ""):
        return None
    try:
        length = int(value)
    except (TypeError, ValueError) as error:
        raise RetrievalSizeError("Remote asset Content-Length is invalid.") from error
    if length < 0:
        raise RetrievalSizeError("Remote asset Content-Length is invalid.")
    return length


def _content_range(headers: Any) -> tuple[int, int, int]:
    value = headers.get("Content-Range") if hasattr(headers, "get") else None
    match = _CONTENT_RANGE.fullmatch(str(value or ""))
    if match is None:
        raise RetrievalSizeError(
            "NCBI ranged asset response has invalid Content-Range."
        )
    start, end, total = map(int, match.groups())
    if total <= 0 or start < 0 or end < start or end >= total:
        raise RetrievalSizeError(
            "NCBI ranged asset response has invalid Content-Range."
        )
    return start, end, total


def _sanitize_origin(value: str) -> str:
    parsed = urlsplit(value)
    hostname = parsed.hostname or ""
    if ":" in hostname and not hostname.startswith("["):
        hostname = f"[{hostname}]"
    netloc = hostname
    if parsed.port is not None:
        netloc = f"{hostname}:{parsed.port}"
    return urlunsplit((parsed.scheme.casefold(), netloc, parsed.path, "", ""))


def _file_digest(path: Path, algorithm: str) -> str:
    if algorithm == "md5":
        digest = hashlib.md5(usedforsecurity=False)
    else:
        digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest().casefold()


__all__ = [
    "AssetDownloader",
    "CacheIntegrityError",
    "DEFAULT_PROVIDER_HOST_SUFFIXES",
    "RetrievalError",
    "RetrievalPolicy",
    "RetrievalSecurityError",
    "RetrievalService",
    "RetrievalSizeError",
]
