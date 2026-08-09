# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
import os
import inspect
import sys
import threading
import unittest
from unittest.mock import Mock

import requests


ROOT = os.path.dirname(os.path.dirname(__file__))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from meta_standards_converter.helpers.request_helper import (  # noqa: E402
    NCBIApplicationIdentity,
    RateLimitedRequester,
    RequestSettings,
)


class FakeTime:
    def __init__(self):
        self.now = 0
        self.sleeps = []

    def clock(self):
        return self.now

    def sleep(self, seconds):
        self.sleeps.append(seconds)
        self.now += seconds


def response(status_code=200, headers=None):
    item = Mock()
    item.status_code = status_code
    item.headers = headers or {}
    item.raise_for_status = Mock()
    if status_code >= 400:
        item.raise_for_status.side_effect = requests.HTTPError(str(status_code))
    return item


class TestRateLimitedRequester(unittest.TestCase):
    def test_ncbi_application_identity_validates_and_builds_parameters(self):
        identity = NCBIApplicationIdentity(
            tool="fibrosis_atlas",
            email="atlas@example.org",
        )

        self.assertEqual(
            {"tool": "fibrosis_atlas", "email": "atlas@example.org"},
            identity.params(),
        )
        with self.assertRaisesRegex(ValueError, "email"):
            NCBIApplicationIdentity(tool="fibrosis_atlas", email="not-an-email")

    def setUp(self):
        RateLimitedRequester.reset_service_state()

    def test_public_request_boundary_has_explicit_return_types(self):
        self.assertEqual(
            inspect.signature(RateLimitedRequester.get).return_annotation,
            "requests.Response",
        )
        self.assertEqual(
            inspect.signature(RateLimitedRequester.reset_service_state).return_annotation,
            "None",
        )

    def test_get_applies_default_timeout(self):
        get = Mock(return_value=response())
        requester = RateLimitedRequester(
            service="test_timeout",
            settings=RequestSettings(timeout=12, request_delay=0),
            get=get,
        )

        requester.get("https://example.org/data", params={"id": "1"})

        get.assert_called_once_with(
            "https://example.org/data",
            params={"id": "1"},
            timeout=12,
        )

    def test_get_preserves_explicit_timeout(self):
        get = Mock(return_value=response())
        requester = RateLimitedRequester(
            service="test_explicit_timeout",
            settings=RequestSettings(timeout=12, request_delay=0),
            get=get,
        )

        requester.get("https://example.org/data", timeout=3)

        get.assert_called_once_with("https://example.org/data", timeout=3)

    def test_get_sleeps_between_sequential_requests_for_same_service(self):
        fake_time = FakeTime()
        get = Mock(return_value=response())
        settings = RequestSettings(request_delay=0.5)
        first = RateLimitedRequester(
            service="shared_service",
            settings=settings,
            get=get,
            sleep=fake_time.sleep,
            clock=fake_time.clock,
        )
        second = RateLimitedRequester(
            service="shared_service",
            settings=settings,
            get=get,
            sleep=fake_time.sleep,
            clock=fake_time.clock,
        )

        first.get("https://example.org/one")
        second.get("https://example.org/two")

        self.assertEqual([0.5], fake_time.sleeps)
        self.assertEqual(2, get.call_count)

    def test_rate_limit_is_shared_by_host_across_service_labels(self):
        fake_time = FakeTime()
        get = Mock(return_value=response())
        first = RateLimitedRequester(
            service="one",
            settings=RequestSettings(request_delay=0.5),
            get=get,
            sleep=fake_time.sleep,
            clock=fake_time.clock,
        )
        second = RateLimitedRequester(
            service="two",
            settings=RequestSettings(request_delay=0.5),
            get=get,
            sleep=fake_time.sleep,
            clock=fake_time.clock,
        )

        first.get("https://api.example.org/one")
        second.get("https://api.example.org/two")

        self.assertEqual([0.5], fake_time.sleeps)

    def test_rate_limit_is_independent_for_different_hosts(self):
        fake_time = FakeTime()
        get = Mock(return_value=response())
        requester = RateLimitedRequester(
            service="shared",
            settings=RequestSettings(request_delay=0.5),
            get=get,
            sleep=fake_time.sleep,
            clock=fake_time.clock,
        )

        requester.get("https://one.example.org/data")
        requester.get("https://two.example.org/data")

        self.assertEqual([], fake_time.sleeps)

    def test_host_in_flight_limit_bounds_overlapping_requests(self):
        entered = 0
        maximum_entered = 0
        state_lock = threading.Lock()
        release = threading.Event()

        def blocking_get(url, **kwargs):
            nonlocal entered, maximum_entered
            with state_lock:
                entered += 1
                maximum_entered = max(maximum_entered, entered)
            release.wait(timeout=1)
            with state_lock:
                entered -= 1
            return response()

        settings = RequestSettings(request_delay=0, max_in_flight=1)
        requester = RateLimitedRequester(
            service="bounded", settings=settings, get=blocking_get
        )
        first = threading.Thread(
            target=requester.get, args=("https://api.example.org/one",)
        )
        second = threading.Thread(
            target=requester.get, args=("https://api.example.org/two",)
        )
        first.start()
        second.start()
        threading.Event().wait(0.05)
        release.set()
        first.join()
        second.join()

        self.assertEqual(1, maximum_entered)

    def test_request_settings_reject_invalid_host_limits(self):
        with self.assertRaisesRegex(ValueError, "max_in_flight"):
            RequestSettings(max_in_flight=0)

    def test_get_retries_transient_status_using_retry_after(self):
        fake_time = FakeTime()
        get = Mock(side_effect=[
            response(429, headers={"Retry-After": "2"}),
            response(200),
        ])
        requester = RateLimitedRequester(
            service="retry_after_service",
            settings=RequestSettings(request_delay=0.5, max_retries=3),
            get=get,
            sleep=fake_time.sleep,
            clock=fake_time.clock,
        )

        result = requester.get("https://example.org/data")

        self.assertEqual(200, result.status_code)
        self.assertEqual([2.0], fake_time.sleeps)
        self.assertEqual(2, get.call_count)

    def test_get_retries_transient_status_using_exponential_backoff(self):
        fake_time = FakeTime()
        get = Mock(side_effect=[
            response(503),
            response(502),
            response(200),
        ])
        requester = RateLimitedRequester(
            service="backoff_service",
            settings=RequestSettings(request_delay=0, max_retries=3),
            get=get,
            sleep=fake_time.sleep,
            clock=fake_time.clock,
        )

        result = requester.get("https://example.org/data")

        self.assertEqual(200, result.status_code)
        self.assertEqual([0.5, 1.0], fake_time.sleeps)
        self.assertEqual(3, get.call_count)

    def test_get_raises_after_retries_are_exhausted(self):
        fake_time = FakeTime()
        get = Mock(side_effect=[
            response(503),
            response(503),
        ])
        requester = RateLimitedRequester(
            service="exhausted_service",
            settings=RequestSettings(request_delay=0, max_retries=1),
            get=get,
            sleep=fake_time.sleep,
            clock=fake_time.clock,
        )

        with self.assertRaises(requests.HTTPError):
            requester.get("https://example.org/data")

        self.assertEqual([0.5], fake_time.sleeps)
        self.assertEqual(2, get.call_count)

    def test_get_retries_transient_request_exceptions_with_deterministic_backoff(self):
        for exception_type in (
            requests.ConnectionError,
            requests.Timeout,
            requests.exceptions.ChunkedEncodingError,
        ):
            with self.subTest(exception_type=exception_type.__name__):
                RateLimitedRequester.reset_service_state()
                fake_time = FakeTime()
                get = Mock(
                    side_effect=[exception_type("transient"), exception_type("transient"), response(200)]
                )
                requester = RateLimitedRequester(
                    service=f"exception_{exception_type.__name__}",
                    settings=RequestSettings(request_delay=0, max_retries=2),
                    get=get,
                    sleep=fake_time.sleep,
                    clock=fake_time.clock,
                )

                result = requester.get("https://example.org/data")

                self.assertEqual(200, result.status_code)
                self.assertEqual(3, get.call_count)
                self.assertEqual([0.5, 1.0], fake_time.sleeps)

    def test_get_reraises_request_exception_after_exact_attempt_count(self):
        fake_time = FakeTime()
        get = Mock(side_effect=requests.Timeout("still unavailable"))
        requester = RateLimitedRequester(
            service="exception_exhausted",
            settings=RequestSettings(request_delay=0, max_retries=2),
            get=get,
            sleep=fake_time.sleep,
            clock=fake_time.clock,
        )

        with self.assertRaisesRegex(requests.Timeout, "still unavailable"):
            requester.get("https://example.org/data")

        self.assertEqual(3, get.call_count)
        self.assertEqual([0.5, 1.0], fake_time.sleeps)

    def test_get_logs_safe_attempt_status_and_duration_without_parameters(self):
        fake_time = FakeTime()
        get = Mock(return_value=response(200))
        requester = RateLimitedRequester(
            service="safe_service",
            settings=RequestSettings(request_delay=0),
            get=get,
            sleep=fake_time.sleep,
            clock=fake_time.clock,
        )

        with self.assertLogs(
            "meta_standards_converter.helpers.request_helper", level="DEBUG"
        ) as logs:
            requester.get(
                "https://example.org/data?api_key=do-not-log",
                params={"token": "do-not-log"},
            )

        output = "\n".join(logs.output)
        self.assertIn(
            "HTTP request service=safe_service host=example.org attempt=1", output
        )
        self.assertIn("status=200", output)
        self.assertIn("elapsed_seconds=", output)
        self.assertNotIn("do-not-log", output)


if __name__ == "__main__":
    unittest.main()
