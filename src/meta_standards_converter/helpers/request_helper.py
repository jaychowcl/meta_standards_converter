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
from typing import Callable
from urllib.parse import urlsplit

import requests


logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class RequestSettings:
    timeout: float = 30
    request_delay: float = 1.0
    max_retries: int = 3
    retry_statuses: frozenset[int] = frozenset({429, 500, 502, 503, 504})
    backoff_base: float = 0.5
    backoff_max: float = 8.0


DEFAULT_REQUEST_SETTINGS = {
    "ncbi_eutils": RequestSettings(request_delay=0.5),
    "geo_ftp": RequestSettings(request_delay=1.0),
    "ena_portal": RequestSettings(request_delay=1.0),
    "biostudies": RequestSettings(request_delay=1.0),
}


class RateLimitedRequester:
    _service_state = {}
    _state_lock = threading.Lock()

    def __init__(
        self,
        service: str,
        settings: RequestSettings | None = None,
        get: Callable | None = None,
        sleep: Callable[[float], None] | None = None,
        clock: Callable[[], float] | None = None,
    ):
        self.service = service
        self.settings = settings or DEFAULT_REQUEST_SETTINGS.get(service, RequestSettings())
        self._get = get or requests.get
        self._sleep = sleep or time.sleep
        self._clock = clock or time.monotonic

    def get(self, url: str, **kwargs):
        if "timeout" not in kwargs:
            kwargs["timeout"] = self.settings.timeout

        response = None
        for attempt in range(self.settings.max_retries + 1):
            self._wait_for_service_slot()
            started = self._clock()
            host = urlsplit(url).hostname or ""
            logger.debug(
                "HTTP request service=%s host=%s attempt=%s timeout=%s",
                self.service,
                host,
                attempt + 1,
                kwargs.get("timeout"),
            )
            response = self._get(url, **kwargs)
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
                return response

            if attempt >= self.settings.max_retries:
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

    def _wait_for_service_slot(self):
        state = self._state_for_service()
        with state["lock"]:
            last_request_at = state["last_request_at"]
            now = self._clock()
            if last_request_at is not None:
                wait = self.settings.request_delay - (now - last_request_at)
                if wait > 0:
                    self._sleep(wait)
                    now = self._clock()
            state["last_request_at"] = now

    def _retry_delay(self, response, attempt: int) -> float:
        retry_after = response.headers.get("Retry-After")
        if retry_after:
            try:
                return max(float(retry_after), 0)
            except ValueError:
                pass
        return min(self.settings.backoff_base * (2 ** attempt), self.settings.backoff_max)

    def _state_for_service(self):
        with self._state_lock:
            state = self._service_state.get(self.service)
            if state is None:
                state = {"lock": threading.Lock(), "last_request_at": None}
                self._service_state[self.service] = state
            return state

    @classmethod
    def reset_service_state(cls):
        with cls._state_lock:
            cls._service_state = {}
