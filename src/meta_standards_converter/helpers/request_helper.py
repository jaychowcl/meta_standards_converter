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

from dataclasses import dataclass
import logging
import threading
import time
from typing import Any, Callable
from urllib.parse import urlsplit

import requests

from meta_standards_converter.runtime_contracts import ResourceProfile


logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class RequestSettings:
    timeout: float | tuple[float, float] = 30
    request_delay: float = 1.0
    max_in_flight: int = 2
    max_retries: int = 3
    retry_statuses: frozenset[int] = frozenset({429, 500, 502, 503, 504})
    backoff_base: float = 0.5
    backoff_max: float = 8.0

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
        if self.request_delay < 0:
            raise ValueError("request_delay must be non-negative")
        if self.max_in_flight < 1:
            raise ValueError("max_in_flight must be a positive integer")
        if self.max_retries < 0:
            raise ValueError("max_retries must be non-negative")

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
    "ncbi_eutils": RequestSettings(request_delay=0.5, max_in_flight=2),
    "geo_ftp": RequestSettings(request_delay=1.0, max_in_flight=2),
    "ena_portal": RequestSettings(request_delay=1.0, max_in_flight=2),
    "biostudies": RequestSettings(request_delay=1.0, max_in_flight=2),
}


class RateLimitedRequester:
    """Apply process-wide request-start and in-flight limits per HTTP host."""

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
    ):
        self.service = service
        self.settings = settings or DEFAULT_REQUEST_SETTINGS.get(service, RequestSettings())
        self._get = get or requests.get
        self._sleep = sleep or time.sleep
        self._clock = clock or time.monotonic
        self._event_emitter = event_emitter

    def get(self, url: str, **kwargs: Any) -> requests.Response:
        """GET one URL under the process-wide host policy and bounded retries."""
        if "timeout" not in kwargs:
            kwargs["timeout"] = self.settings.timeout

        response = None
        for attempt in range(self.settings.max_retries + 1):
            host = (urlsplit(url).hostname or "").lower()
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
                    response = self._get(url, **kwargs)
                finally:
                    self._release_host_slot(state)
            except (
                requests.exceptions.ConnectionError,
                requests.exceptions.Timeout,
                requests.exceptions.ChunkedEncodingError,
            ) as error:
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
                logger.info(
                    "HTTP retry service=%s host=%s exception=%s "
                    "next_attempt=%s delay_seconds=%.3f",
                    self.service,
                    host,
                    type(error).__name__,
                    attempt + 2,
                    delay,
                )
                self._sleep(delay)
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

            if attempt >= self.settings.max_retries:
                self._emit_request_event(host, "failed", attempt + 1, elapsed, response.status_code)
                response.raise_for_status()
                return response

            delay = self._retry_delay(response=response, attempt=attempt)
            logger.info(
                "HTTP retry service=%s host=%s status=%s next_attempt=%s delay_seconds=%.3f",
                self.service,
                host,
                response.status_code,
                attempt + 2,
                delay,
            )
            self._sleep(delay)

        return response

    def _emit_request_event(
        self,
        host: str,
        status: str,
        attempts: int,
        elapsed: float,
        status_code: int | None,
    ) -> None:
        if self._event_emitter is not None:
            self._event_emitter.emit(
                "provider_request",
                component=self.service,
                level="INFO" if status == "success" else "ERROR",
                status=status,
                duration_seconds=elapsed,
                attributes={"host": host, "attempts": attempts, "status_code": status_code},
            )

    def _acquire_host_slot(self, host: str) -> dict[str, Any]:
        state = self._state_for_host(host)
        with state["condition"]:
            state["request_delay"] = max(
                state["request_delay"], self.settings.request_delay
            )
            state["max_in_flight"] = min(
                state["max_in_flight"], self.settings.max_in_flight
            )
            while state["in_flight"] >= state["max_in_flight"]:
                state["condition"].wait()
            last_request_at = state["last_request_at"]
            now = self._clock()
            if last_request_at is not None:
                wait = state["request_delay"] - (now - last_request_at)
                if wait > 0:
                    self._sleep(wait)
                    now = self._clock()
            state["last_request_at"] = now
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
        if retry_after:
            try:
                return max(float(retry_after), 0)
            except ValueError:
                pass
        return min(self.settings.backoff_base * (2 ** attempt), self.settings.backoff_max)

    def _state_for_host(self, host: str) -> dict[str, Any]:
        with self._state_lock:
            state = self._host_state.get(host)
            if state is None:
                state = {
                    "condition": threading.Condition(),
                    "last_request_at": None,
                    "in_flight": 0,
                    "max_in_flight": self.settings.max_in_flight,
                    "request_delay": self.settings.request_delay,
                }
                self._host_state[host] = state
            return state

    @classmethod
    def reset_service_state(cls) -> None:
        """Reset process host state; retained name preserves the v1 test API."""

        with cls._state_lock:
            cls._host_state = {}
