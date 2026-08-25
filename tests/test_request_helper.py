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
import multiprocessing
from pathlib import Path
import sys
import tempfile
import threading
import time
import unittest
from unittest.mock import Mock, patch

import requests


ROOT = os.path.dirname(os.path.dirname(__file__))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from meta_standards_converter.helpers.request_helper import (  # noqa: E402
    HostRequestCooldownDeferred,
    HostRequestGate,
    NCBIApplicationIdentity,
    RateLimitedRequester,
    RequestSettings,
)


def _take_host_gate_slot(directory, barrier, output):
    gate = HostRequestGate(directory=directory)
    barrier.wait(timeout=5)
    gate.wait("eutils.ncbi.nlm.nih.gov", min_interval_seconds=0.05)
    output.put(time.time())


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
        secret_identity = NCBIApplicationIdentity(
            tool="fibrosis_atlas",
            email="atlas@example.org",
            api_key="do-not-render",
        )
        self.assertEqual(
            {
                "tool": "fibrosis_atlas",
                "email": "atlas@example.org",
                "api_key": "do-not-render",
            },
            secret_identity.params(),
        )
        self.assertNotIn("do-not-render", repr(secret_identity))
        with self.assertLogs(
            "meta_standards_converter.helpers.request_helper",
            level="WARNING",
        ) as logs:
            unidentified = NCBIApplicationIdentity(
                tool="fibrosis_atlas",
                email=None,
            )
            self.assertEqual({"tool": "fibrosis_atlas"}, unidentified.params())
            unidentified.params()
        self.assertEqual(1, sum("contact email" in line for line in logs.output))
        with self.assertRaisesRegex(ValueError, "email"):
            NCBIApplicationIdentity(tool="fibrosis_atlas", email="not-an-email")

    def setUp(self):
        self._gate_directory = tempfile.TemporaryDirectory()
        self._old_gate_directory = os.environ.get("SCIENTIFIC_PROVIDER_GATE_DIR")
        os.environ["SCIENTIFIC_PROVIDER_GATE_DIR"] = self._gate_directory.name
        HostRequestGate.reset_default()
        RateLimitedRequester.reset_service_state()

    def tearDown(self):
        HostRequestGate.reset_default()
        if self._old_gate_directory is None:
            os.environ.pop("SCIENTIFIC_PROVIDER_GATE_DIR", None)
        else:
            os.environ["SCIENTIFIC_PROVIDER_GATE_DIR"] = self._old_gate_directory
        self._gate_directory.cleanup()

    def test_host_gate_paces_separate_processes(self):
        context = multiprocessing.get_context("fork")
        barrier = context.Barrier(2)
        output = context.Queue()
        processes = [
            context.Process(
                target=_take_host_gate_slot,
                args=(self._gate_directory.name, barrier, output),
            )
            for _ in range(2)
        ]

        for process in processes:
            process.start()
        timestamps = sorted(output.get(timeout=5) for _ in processes)
        for process in processes:
            process.join(timeout=5)
            self.assertEqual(0, process.exitcode)

        self.assertGreaterEqual(timestamps[1] - timestamps[0], 0.04)

    def test_host_gate_persists_cooldown_across_instances(self):
        fake_time = FakeTime()
        first = HostRequestGate(directory=self._gate_directory.name)
        second = HostRequestGate(directory=self._gate_directory.name)

        first.defer(
            "api.example.org",
            delay_seconds=12,
            clock=fake_time.clock,
        )
        waited = second.wait(
            "api.example.org",
            min_interval_seconds=0.5,
            sleep=fake_time.sleep,
            clock=fake_time.clock,
        )

        self.assertEqual(12, waited)
        self.assertEqual([12], fake_time.sleeps)

    def test_host_gate_falls_back_when_implicit_runtime_directory_is_read_only(self):
        fallback = Path(self._gate_directory.name)
        with (
            patch.object(
                HostRequestGate,
                "default_directory",
                return_value=Path("/read-only/runtime/gates"),
            ),
            patch.object(HostRequestGate, "fallback_directory", return_value=fallback),
            patch.object(
                HostRequestGate,
                "_validated_directory",
                side_effect=[OSError("read only"), fallback],
            ) as validate,
        ):
            gate = HostRequestGate()

        self.assertEqual(fallback, gate.directory)
        self.assertEqual(2, validate.call_count)

    def test_host_gate_does_not_fall_back_from_an_explicit_directory(self):
        explicit = Path(self._gate_directory.name)
        with patch.object(
            HostRequestGate,
            "_validated_directory",
            side_effect=OSError("read only"),
        ) as validate:
            with self.assertRaisesRegex(OSError, "read only"):
                HostRequestGate(directory=explicit)

        validate.assert_called_once_with(explicit)

    def test_retry_after_http_date_is_honoured_without_eight_second_cap(self):
        fake_time = FakeTime()
        get = Mock(side_effect=[
            response(
                429,
                headers={"Retry-After": "Thu, 01 Jan 1970 00:00:20 GMT"},
            ),
            response(200),
        ])
        requester = RateLimitedRequester(
            service="retry_after_date_service",
            settings=RequestSettings(request_delay=0, max_retries=1),
            get=get,
            sleep=fake_time.sleep,
            clock=fake_time.clock,
        )

        result = requester.get("https://example.org/data")

        self.assertEqual(200, result.status_code)
        self.assertEqual([20.0], fake_time.sleeps)

    def test_retry_after_beyond_inline_budget_defers_without_second_request(self):
        fake_time = FakeTime()
        get = Mock(return_value=response(429, headers={"Retry-After": "60"}))
        requester = RateLimitedRequester(
            service="retry_after_deferred_service",
            settings=RequestSettings(
                request_delay=0,
                max_retries=1,
                max_inline_wait=30,
            ),
            get=get,
            sleep=fake_time.sleep,
            clock=fake_time.clock,
        )

        with self.assertRaises(HostRequestCooldownDeferred) as raised:
            requester.get("https://deferred.example.org/data")

        self.assertEqual(60, raised.exception.retry_at)
        self.assertEqual(1, get.call_count)
        self.assertEqual([], fake_time.sleeps)

    def test_public_request_boundary_has_explicit_return_types(self):
        self.assertEqual(
            inspect.signature(RateLimitedRequester.get).return_annotation,
            "requests.Response",
        )
        self.assertEqual(
            inspect.signature(RateLimitedRequester.reset_service_state).return_annotation,
            "None",
        )

    def test_default_retry_statuses_include_provider_throttling_not_plain_4xx(self):
        statuses = RequestSettings().retry_statuses

        self.assertTrue({403, 408, 425, 429, 500, 502, 503, 504} <= statuses)
        self.assertNotIn(400, statuses)
        self.assertNotIn(404, statuses)

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
            random_value=lambda: 1.0,
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
            random_value=lambda: 1.0,
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
                    random_value=lambda: 1.0,
                )

                host = exception_type.__name__.lower()
                result = requester.get(f"https://{host}.example.org/data")

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
            random_value=lambda: 1.0,
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
