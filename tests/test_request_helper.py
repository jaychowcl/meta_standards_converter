# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
import os
import sys
import unittest
from unittest.mock import Mock

import requests


ROOT = os.path.dirname(os.path.dirname(__file__))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from meta_standards_converter.helpers.request_helper import (  # noqa: E402
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
    def setUp(self):
        RateLimitedRequester.reset_service_state()

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
