# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""
Shared HTTP request helpers for external metadata services.
"""

from __future__ import annotations

from meta_standards_converter.provider_issues import provider_boundary
from dataclasses import dataclass, field
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime
import fcntl
import hashlib
import json
from meta_standards_converter.conversion_logging import count, emit, current
import logging
import os
from pathlib import Path
import random
import re
import tempfile
import threading
import time
from typing import Any, Callable
from urllib.parse import parse_qsl, urlsplit

import requests

from meta_standards_converter.runtime_contracts import ResourceProfile


logger = logging.getLogger(__name__)

DEFAULT_NCBI_TOOL = "meta_standards_converter"
DEFAULT_NCBI_EMAIL = "jaychowcl@gmail.com"
_NCBI_TOOL_PATTERN = re.compile(r"^[A-Za-z0-9_.-]{1,64}$")
_EMAIL_PATTERN = re.compile(r"^[^@\s]+@[^@\s]+\.[^@\s]+$")


@dataclass(frozen=True)
class NCBIApplicationIdentity:
    """Contactable application identity required by NCBI E-utilities."""

    tool: str = DEFAULT_NCBI_TOOL
    email: str | None = DEFAULT_NCBI_EMAIL
    api_key: str | None = field(
        default_factory=lambda: os.environ.get("NCBI_API_KEY") or None,
        repr=False,
    )
    _warned_missing_email: bool = field(
        default=False,
        init=False,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        if not _NCBI_TOOL_PATTERN.fullmatch(self.tool):
            raise ValueError(
                "NCBI application tool must contain 1-64 letters, digits, "
                "underscores, dots, or hyphens"
            )
        if self.email is not None and not _EMAIL_PATTERN.fullmatch(self.email):
            raise ValueError("NCBI application email must be a valid contact address")
        if self.api_key is not None and (
            not self.api_key.strip() or any(character.isspace() for character in self.api_key)
        ):
            raise ValueError("NCBI API key must be a nonblank token without whitespace")

    def params(self) -> dict[str, str]:
        params = {"tool": self.tool}
        if self.email is not None:
            params["email"] = self.email
        elif not self._warned_missing_email:
            logger.warning(
                "NCBI E-utilities contact email is not configured; "
                "conservative unauthenticated pacing remains active"
            )
            object.__setattr__(self, "_warned_missing_email", True)
        if self.api_key is not None:
            params["api_key"] = self.api_key
        return params


class HostRequestCooldownDeferred(TimeoutError):
    """The shared provider cooldown exceeds the caller's inline wait budget."""

    def __init__(self, key: str, retry_at: float) -> None:
        super().__init__(f"provider cooldown remains active for {key}")
        self.key = key
        self.retry_at = float(retry_at)


class HostRequestGate:
    """Coordinate provider request starts across local processes for one user."""

    CONTRACT_VERSION = 1
    ENVIRONMENT_VARIABLE = "SCIENTIFIC_PROVIDER_GATE_DIR"
    _default = None
    _default_lock = threading.Lock()

    def __init__(self, directory: str | Path | None = None) -> None:
        if directory is not None:
            self.directory = self._validated_directory(Path(directory))
            return

        candidate = self.default_directory()
        try:
            self.directory = self._validated_directory(candidate)
        except OSError:
            fallback = self.fallback_directory()
            if candidate == fallback:
                raise
            self.directory = self._validated_directory(fallback)

    @classmethod
    def default(cls) -> "HostRequestGate":
        with cls._default_lock:
            if cls._default is None:
                cls._default = cls()
            return cls._default

    @classmethod
    def reset_default(cls) -> None:
        """Forget the process singleton without deleting shared pacing state."""

        with cls._default_lock:
            cls._default = None

    @classmethod
    def default_directory(cls) -> Path:
        configured = os.environ.get(cls.ENVIRONMENT_VARIABLE)
        if configured:
            return Path(configured)
        runtime_root = os.environ.get("XDG_RUNTIME_DIR")
        if runtime_root:
            return Path(runtime_root) / "scientific-provider-rate-gates-v1"
        return cls.fallback_directory()

    @staticmethod
    def fallback_directory() -> Path:
        return (
            Path(tempfile.gettempdir())
            / f"scientific-provider-rate-gates-{os.getuid()}"
            / "v1"
        )

    def wait(
        self,
        key: str,
        *,
        min_interval_seconds: int | float,
        max_wait_seconds: int | float | None = None,
        sleep: Callable[[float], None] = time.sleep,
        clock: Callable[[], float] = time.time,
    ) -> float:
        normalized_key = self._normalized_key(key)
        interval = self._nonnegative_seconds(
            min_interval_seconds, name="minimum request interval"
        )
        maximum_wait = (
            None
            if max_wait_seconds is None
            else self._nonnegative_seconds(max_wait_seconds, name="maximum gate wait")
        )
        with self._locked_state(normalized_key) as (handle, state):
            now = float(clock())
            next_allowed = max(
                float(state.get("next_start_at") or 0),
                float(state.get("cooldown_until") or 0),
            )
            delay = max(0.0, next_allowed - now)
            if maximum_wait is not None and delay > maximum_wait:
                raise HostRequestCooldownDeferred(normalized_key, next_allowed)
            if delay:
                sleep(delay)
            granted_at = float(clock())
            state["next_start_at"] = granted_at + interval
            self._write_state(handle, state)
            return delay

    def defer(
        self,
        key: str,
        *,
        delay_seconds: int | float,
        clock: Callable[[], float] = time.time,
    ) -> float:
        normalized_key = self._normalized_key(key)
        delay = self._nonnegative_seconds(delay_seconds, name="provider cooldown")
        with self._locked_state(normalized_key) as (handle, state):
            retry_at = float(clock()) + delay
            state["cooldown_until"] = max(
                float(state.get("cooldown_until") or 0), retry_at
            )
            self._write_state(handle, state)
            return float(state["cooldown_until"])

    def slot(
        self,
        key: str,
        *,
        min_interval_seconds: int | float,
        max_wait_seconds: int | float | None = None,
        sleep: Callable[[float], None] = time.sleep,
        monotonic_clock: Callable[[], float] = time.monotonic,
        wall_clock: Callable[[], float] = time.time,
        cancelled: Callable[[], bool] | None = None,
    ):
        """Return a one-at-a-time lease that also applies start pacing."""

        return _HostRequestSlot(
            self,
            self._normalized_key(key),
            min_interval_seconds=self._nonnegative_seconds(
                min_interval_seconds, name="minimum request interval"
            ),
            max_wait_seconds=(
                None
                if max_wait_seconds is None
                else self._nonnegative_seconds(
                    max_wait_seconds, name="maximum gate wait"
                )
            ),
            sleep=sleep,
            monotonic_clock=monotonic_clock,
            wall_clock=wall_clock,
            cancelled=cancelled,
        )

    @staticmethod
    def retry_after_seconds(
        value: Any,
        *,
        clock: Callable[[], float] = time.time,
    ) -> float | None:
        if value in {None, ""}:
            return None
        try:
            return max(float(value), 0.0)
        except (TypeError, ValueError):
            try:
                retry_at = parsedate_to_datetime(str(value))
            except (TypeError, ValueError, OverflowError):
                return None
            if retry_at.tzinfo is None:
                retry_at = retry_at.replace(tzinfo=timezone.utc)
            now = datetime.fromtimestamp(float(clock()), timezone.utc)
            return max((retry_at - now).total_seconds(), 0.0)

    @staticmethod
    def _normalized_key(key: str) -> str:
        if not isinstance(key, str) or not key.strip():
            raise ValueError("provider gate key must be a non-empty string")
        return key.strip().casefold().rstrip(".")

    @staticmethod
    def _nonnegative_seconds(value: Any, *, name: str) -> float:
        if isinstance(value, bool) or not isinstance(value, (int, float)) or value < 0:
            raise ValueError(f"{name} must be a non-negative number")
        return float(value)

    @staticmethod
    def _validated_directory(directory: Path) -> Path:
        directory.mkdir(mode=0o700, parents=True, exist_ok=True)
        if directory.is_symlink():
            raise RuntimeError("provider gate directory must not be a symbolic link")
        stat = directory.stat()
        if stat.st_uid != os.getuid():
            raise RuntimeError("provider gate directory must be owned by the current user")
        directory.chmod(0o700)
        return directory

    def _locked_state(self, key: str):
        return _LockedGateState(self, key)

    def _state_path(self, key: str) -> Path:
        digest = hashlib.sha256(key.encode("utf-8")).hexdigest()
        return self.directory / f"host-{digest}.json"

    def _slot_path(self, key: str) -> Path:
        digest = hashlib.sha256(key.encode("utf-8")).hexdigest()
        return self.directory / f"slot-{digest}.lock"

    @staticmethod
    def _read_state(handle) -> dict[str, Any]:
        handle.seek(0)
        raw = handle.read()
        if not raw:
            return {
                "contract_version": HostRequestGate.CONTRACT_VERSION,
                "next_start_at": 0.0,
                "cooldown_until": 0.0,
            }
        try:
            value = json.loads(raw.decode("utf-8"))
        except (UnicodeDecodeError, json.JSONDecodeError) as error:
            raise RuntimeError("provider gate state is invalid") from error
        if not isinstance(value, dict) or value.get("contract_version") != 1:
            raise RuntimeError("provider gate state contract is incompatible")
        return value

    @staticmethod
    def _write_state(handle, state: dict[str, Any]) -> None:
        encoded = json.dumps(
            {
                "contract_version": HostRequestGate.CONTRACT_VERSION,
                "next_start_at": float(state.get("next_start_at") or 0),
                "cooldown_until": float(state.get("cooldown_until") or 0),
            },
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")
        handle.seek(0)
        handle.truncate()
        handle.write(encoded)
        handle.flush()
        os.fsync(handle.fileno())


class _LockedGateState:
    def __init__(self, gate: HostRequestGate, key: str) -> None:
        self.gate = gate
        self.key = key
        self.handle = None

    def __enter__(self):
        flags = os.O_CREAT | os.O_RDWR
        if hasattr(os, "O_NOFOLLOW"):
            flags |= os.O_NOFOLLOW
        descriptor = os.open(self.gate._state_path(self.key), flags, 0o600)
        self.handle = os.fdopen(descriptor, "r+b")
        fcntl.flock(self.handle.fileno(), fcntl.LOCK_EX)
        return self.handle, self.gate._read_state(self.handle)

    def __exit__(self, exc_type, exc, traceback) -> None:
        if self.handle is not None:
            fcntl.flock(self.handle.fileno(), fcntl.LOCK_UN)
            self.handle.close()


class _HostRequestSlot:
    POLL_SECONDS = 0.05

    def __init__(
        self,
        gate: HostRequestGate,
        key: str,
        *,
        min_interval_seconds: float,
        max_wait_seconds: float | None,
        sleep: Callable[[float], None],
        monotonic_clock: Callable[[], float],
        wall_clock: Callable[[], float],
        cancelled: Callable[[], bool] | None,
    ) -> None:
        self.gate = gate
        self.key = key
        self.min_interval_seconds = min_interval_seconds
        self.max_wait_seconds = max_wait_seconds
        self.sleep = sleep
        self.monotonic_clock = monotonic_clock
        self.wall_clock = wall_clock
        self.cancelled = cancelled
        self.handle = None

    def __enter__(self) -> float:
        started = float(self.monotonic_clock())
        flags = os.O_CREAT | os.O_RDWR
        if hasattr(os, "O_NOFOLLOW"):
            flags |= os.O_NOFOLLOW
        descriptor = os.open(self.gate._slot_path(self.key), flags, 0o600)
        self.handle = os.fdopen(descriptor, "r+b")
        try:
            while True:
                if self.cancelled is not None and self.cancelled():
                    raise InterruptedError("provider request cancelled before start")
                try:
                    fcntl.flock(
                        self.handle.fileno(),
                        fcntl.LOCK_EX | fcntl.LOCK_NB,
                    )
                    break
                except BlockingIOError:
                    elapsed = float(self.monotonic_clock()) - started
                    if (
                        self.max_wait_seconds is not None
                        and elapsed >= self.max_wait_seconds
                    ):
                        raise HostRequestCooldownDeferred(
                            self.key,
                            float(self.wall_clock()),
                        )
                    remaining = (
                        self.POLL_SECONDS
                        if self.max_wait_seconds is None
                        else min(
                            self.POLL_SECONDS,
                            max(self.max_wait_seconds - elapsed, 0.0),
                        )
                    )
                    self.sleep(remaining)

            lease_wait = float(self.monotonic_clock()) - started
            remaining_wait = (
                None
                if self.max_wait_seconds is None
                else max(self.max_wait_seconds - lease_wait, 0.0)
            )
            rate_wait = self.gate.wait(
                self.key,
                min_interval_seconds=self.min_interval_seconds,
                max_wait_seconds=remaining_wait,
                sleep=self.sleep,
                clock=self.wall_clock,
            )
            return lease_wait + rate_wait
        except BaseException:
            self._release()
            raise

    def __exit__(self, exc_type, exc, traceback) -> None:
        self._release()

    def _release(self) -> None:
        if self.handle is not None:
            fcntl.flock(self.handle.fileno(), fcntl.LOCK_UN)
            self.handle.close()
            self.handle = None


@dataclass(frozen=True)
class RequestSettings:
    timeout: float | tuple[float, float] = 30
    # None selects the policy for the actual URL. Numeric overrides retain the
    # existing caller-controlled spacing contract (including offline clients).
    request_delay: float | None = None
    max_in_flight: int = 2
    max_retries: int = 3
    retry_statuses: frozenset[int] = frozenset(
        {403, 408, 425, 429, 500, 502, 503, 504}
    )
    backoff_base: float = 0.5
    backoff_max: float = 8.0
    max_inline_wait: float = 30.0

    def __post_init__(self) -> None:
        timeout_values = (
            self.timeout
            if isinstance(self.timeout, tuple)
            else (self.timeout,)
        )
        if (
            len(timeout_values) not in {1, 2}
            or any(value <= 0 for value in timeout_values)
        ):
            raise ValueError("timeout must be positive")
        if self.request_delay is not None and self.request_delay < 0:
            raise ValueError("request_delay must be non-negative")
        if self.max_in_flight < 1:
            raise ValueError("max_in_flight must be a positive integer")
        if self.max_retries < 0:
            raise ValueError("max_retries must be non-negative")
        if self.max_inline_wait < 0:
            raise ValueError("max_inline_wait must be non-negative")

    @classmethod
    def from_resource_profile(
        cls,
        profile: ResourceProfile,
        **overrides: Any,
    ) -> "RequestSettings":
        return cls(
            timeout=(
                profile.connect_timeout_seconds,
                profile.read_timeout_seconds,
            ),
            max_in_flight=profile.network_workers,
            **overrides,
        )


DEFAULT_REQUEST_SETTINGS = {
    service: RequestSettings(max_in_flight=2)
    for service in ("ncbi_eutils", "geo_ftp", "ena_portal", "biosamples", "biostudies")
}


class RateLimitedRequester:
    """Apply shared service pacing and process in-flight limits per HTTP host."""

    _host_state = {}
    _state_lock = threading.Lock()

    def __init__(
        self,
        service: str,
        settings: RequestSettings | None = None,
        get: Callable | None = None,
        sleep: Callable[[float], None] | None = None,
        clock: Callable[[], float] | None = None,
        event_emitter=None,
        host_gate: HostRequestGate | None = None,
        random_value: Callable[[], float] = random.random,
    ):
        self.service = service
        self.settings = settings or DEFAULT_REQUEST_SETTINGS.get(service, RequestSettings())
        from .http_transport import get as pooled_get
        self._get = get or pooled_get
        self._sleep = sleep or time.sleep
        self._clock = clock or time.monotonic
        self._wall_clock = clock or time.time
        self._event_emitter = event_emitter
        self._host_gate = host_gate or HostRequestGate.default()
        self._random_value = random_value
        self.provider_attempts = 0
        self.retry_count = 0
        self.rate_wait_seconds = 0.0
        self._metrics_lock = threading.Lock()

    @provider_boundary()
    def get(self, url: str, **kwargs: Any) -> requests.Response:
        """GET one URL, accounting for each retry and redirect under its policy."""
        if kwargs.get('allow_redirects') is False:
            return self._request(url, **kwargs)
        # Requests otherwise follows redirects inside one opaque call. Expose
        # every hop to the same host gate, retry policy and operation allowance.
        from urllib.parse import urljoin
        kwargs['allow_redirects'] = False
        for redirect in range(31):
            response = self._request(url, **kwargs)
            if response.status_code not in {301, 302, 303, 307, 308} or not response.headers.get('Location'):
                return response
            target = urljoin(url, response.headers['Location'])
            response.close()
            if redirect == 30: raise requests.TooManyRedirects('Metadata redirect limit exceeded')
            kwargs.pop('params', None)
            if urlsplit(target).hostname != urlsplit(url).hostname:
                kwargs.pop('auth', None)
                kwargs['headers'] = {k:v for k,v in kwargs.get('headers', {}).items() if k.lower() not in {'authorization', 'cookie'}}
                kwargs.pop('cookies', None)
            url = target

    def consume(self, url, consumer, **kwargs):
        """Consume a streamed body within the existing physical retry loop."""
        response = self.get(url, _consume_response=consumer, **kwargs)
        try:
            response.raise_for_status()
            return response._msc_consumed
        finally:
            # Successful consumers already close their response, including when
            # decoding fails. Error/redirect terminal responses still need closing.
            if not getattr(response, '_msc_body_consumed', False):
                response.close()

    def _request_policy(self, url: str, params: Any) -> tuple[str, float]:
        """Classify by URL, not the source that happened to create this client."""
        parsed = urlsplit(url)
        host = (parsed.hostname or "").lower().rstrip(".")
        key, interval = host, 1.0
        if host == "www.ebi.ac.uk":
            for path, delay in (("ena", 0.1), ("biosamples", 1.0), ("biostudies", 1.0)):
                if parsed.path == "/" + path or parsed.path.startswith("/" + path + "/"):
                    key, interval = host + "/" + path, delay
                    break
        elif host in {"eutils.ncbi.nlm.nih.gov", "www.ncbi.nlm.nih.gov"} and parsed.path.startswith("/entrez/eutils/"):
            key, interval = "eutils.ncbi.nlm.nih.gov", 0.5
            # Match Requests' actual query encoding, including URL parameters,
            # tuple sequences and encoded strings. Environment presence alone
            # does not grant the keyed rate; the key must be sent to E-utilities.
            prepared = requests.PreparedRequest()
            prepared.prepare_url(url, params)
            if any(name == "api_key" and value.strip()
                   for name, value in parse_qsl(urlsplit(prepared.url).query)):
                interval = 0.125
        if self.settings.request_delay is not None:
            interval = self.settings.request_delay
        return key, interval

    def _request(self, url, **kwargs):
        from ..metadata.retrieval_budget import _active, _check_cancelled
        budget = _active.get()
        consumer = kwargs.pop('_consume_response', None)
        response_check = kwargs.pop('_response_check', None)
        if "timeout" not in kwargs:
            kwargs["timeout"] = self.settings.timeout

        response = None
        for attempt in range(self.settings.max_retries + 1):
            _check_cancelled()
            if budget: budget.check(request=True)
            host = (urlsplit(url).hostname or "").lower()
            gate_key, interval = self._request_policy(url, kwargs.get("params"))
            state = self._acquire_host_slot(host)
            started = self._clock()
            logger.debug(
                "HTTP request service=%s host=%s attempt=%s timeout=%s",
                self.service,
                host,
                attempt + 1,
                kwargs.get("timeout"),
            )
            try:
                try:
                    remaining = budget.remaining() if budget else None
                    maximum = self.settings.max_inline_wait
                    if remaining is not None: maximum = min(maximum, remaining) if maximum is not None else remaining
                    try:
                        waited = self._host_gate.wait(
                            gate_key, min_interval_seconds=interval,
                            max_wait_seconds=maximum, sleep=self._cancellable_sleep, clock=self._wall_clock)
                    except HostRequestCooldownDeferred:
                        if remaining is not None and maximum == remaining:
                            budget.stop('deadline cannot accommodate provider rate-limit wait')
                        raise
                    with self._metrics_lock:
                        self.rate_wait_seconds += float(waited or 0)
                    count("rate_wait_seconds", float(waited or 0))
                    effective = dict(kwargs)
                    _check_cancelled()
                    from .metadata_stream import charge_attempt
                    charge_attempt(url, effective.get('params'))
                    if budget: effective['timeout'] = budget.attempt(effective['timeout'])
                    count("provider_attempts")
                    with self._metrics_lock:
                        self.provider_attempts += 1
                    request_started = time.monotonic()
                    response = self._get(url, **effective)
                    response._msc_request_started_at = request_started
                    if response_check is not None:
                        try:
                            response_check(response)
                        except BaseException:
                            response.close()
                            self._emit_request_event(host, 'failed', attempt + 1,
                                                     self._clock() - started, response.status_code)
                            raise
                    if consumer is not None and response.status_code < 400 and not (
                            response.status_code in {301, 302, 303, 307, 308}
                            and response.headers.get('Location')):
                        response._msc_consumed = consumer(response)
                        response._msc_body_consumed = True
                finally:
                    self._release_host_slot(state)
            except (
                requests.exceptions.ConnectionError,
                requests.exceptions.Timeout,
                requests.exceptions.ChunkedEncodingError,
                requests.exceptions.ContentDecodingError,
            ) as error:
                from .metadata_stream import split_payload, _BatchHeaderTimeout
                if split_payload.get() and isinstance(error, requests.exceptions.ReadTimeout):
                    # An availability check in ENASource distinguishes a slow
                    # batch from an unavailable service before splitting it.
                    self._host_gate.defer(gate_key,
                        delay_seconds=self._retry_delay(response=None, attempt=attempt),
                        clock=self._wall_clock)
                    self._emit_request_event(host, 'failed', attempt+1,
                                             self._clock()-started, None)
                    raise _BatchHeaderTimeout('Large XML request headers stalled') from error
                if attempt >= self.settings.max_retries:
                    self._emit_request_event(
                        host,
                        "failed",
                        attempt + 1,
                        self._clock() - started,
                        None,
                    )
                    raise
                delay = self._retry_delay(response=None, attempt=attempt)
                self._host_gate.defer(
                    gate_key,
                    delay_seconds=delay,
                    clock=self._wall_clock,
                )
                count("retries")
                with self._metrics_lock:
                    self.retry_count += 1
                logger.debug(
                    "HTTP retry service=%s host=%s exception=%s "
                    "next_attempt=%s delay_seconds=%.3f",
                    self.service,
                    host,
                    type(error).__name__,
                    attempt + 2,
                    delay,
                )
                continue
            elapsed = self._clock() - started
            logger.debug(
                "HTTP response service=%s host=%s attempt=%s status=%s elapsed_seconds=%.3f",
                self.service,
                host,
                attempt + 1,
                response.status_code,
                elapsed,
            )

            if response.status_code not in self.settings.retry_statuses:
                self._emit_request_event(
                    host, "success" if response.status_code < 400 else "failed",
                    attempt + 1, elapsed, response.status_code,
                )
                return response

            delay = self._retry_delay(response=response, attempt=attempt)
            self._host_gate.defer(
                gate_key,
                delay_seconds=delay,
                clock=self._wall_clock,
            )
            # Even when this caller has exhausted retries, other callers must
            # observe the provider's cooldown before starting their next request.
            if attempt >= self.settings.max_retries:
                self._emit_request_event(host, "failed", attempt + 1, elapsed, response.status_code)
                try:
                    response.raise_for_status()
                except BaseException:
                    response.close()
                    raise
                return response
            count("retries")
            with self._metrics_lock:
                self.retry_count += 1
            logger.debug(
                "HTTP retry service=%s host=%s status=%s next_attempt=%s delay_seconds=%.3f",
                self.service,
                host,
                response.status_code,
                attempt + 2,
                delay,
            )
            close = getattr(response, "close", None)
            if callable(close):
                close()

        return response

    def _emit_request_event(
        self,
        host: str,
        status: str,
        attempts: int,
        elapsed: float,
        status_code: int | None,
    ) -> None:
        emit("provider_request", level="DEBUG", status=status, provider=self.service,
             host=host, attempts=attempts, status_code=status_code, duration_seconds=elapsed)
        if self._event_emitter is not None and current() is None:
            self._event_emitter.emit(
                "provider_request",
                component=self.service,
                level="INFO" if status == "success" else "ERROR",
                status=status,
                duration_seconds=elapsed,
                attributes={"host": host, "attempts": attempts, "status_code": status_code},
            )

    def _cancellable_sleep(self, delay):
        from ..metadata.retrieval_budget import _cancelled, _check_cancelled
        if _cancelled.get() is None:
            return self._sleep(delay)
        while delay > 0:
            _check_cancelled()
            step = min(delay, .1)
            self._sleep(step)
            delay -= step
        _check_cancelled()

    def _acquire_host_slot(self, host: str) -> dict[str, Any]:
        state = self._state_for_host(host)
        with state["condition"]:
            state["max_in_flight"] = min(
                state["max_in_flight"], self.settings.max_in_flight
            )
            while state["in_flight"] >= state["max_in_flight"]:
                from ..metadata.retrieval_budget import _active, _cancelled, _check_cancelled
                _check_cancelled()
                budget = _active.get()
                timeout = budget.remaining() if budget else None
                if _cancelled.get() is not None:
                    timeout = min(timeout, .1) if timeout is not None else .1
                state["condition"].wait(timeout=timeout)
            state["in_flight"] += 1
        return state

    @staticmethod
    def _release_host_slot(state: dict[str, Any]) -> None:
        with state["condition"]:
            state["in_flight"] -= 1
            state["condition"].notify()

    def _retry_delay(
        self, response: requests.Response | None, attempt: int
    ) -> float:
        retry_after = response.headers.get("Retry-After") if response is not None else None
        parsed_retry_after = HostRequestGate.retry_after_seconds(
            retry_after,
            clock=self._wall_clock,
        )
        if parsed_retry_after is not None:
            return parsed_retry_after
        ceiling = min(
            self.settings.backoff_base * (2 ** attempt),
            self.settings.backoff_max,
        )
        jitter = float(self._random_value())
        if not 0 <= jitter <= 1:
            raise ValueError("random_value must return a number between zero and one")
        return ceiling * jitter

    def _state_for_host(self, host: str) -> dict[str, Any]:
        with self._state_lock:
            state = self._host_state.get(host)
            if state is None:
                state = {
                    "condition": threading.Condition(),
                    "last_request_at": None,
                    "in_flight": 0,
                    "max_in_flight": self.settings.max_in_flight,
                }
                self._host_state[host] = state
            return state

    @classmethod
    def reset_service_state(cls) -> None:
        """Reset process host state; retained name preserves the v1 test API."""

        with cls._state_lock:
            cls._host_state = {}
