"""The unified boundary retains retry decisions and measures failed requests."""
import json

import pytest
import requests

from meta_standards_converter.converters import Converter
from meta_standards_converter.converters.unified.contracts import InputError
from meta_standards_converter.helpers.request_helper import HostRequestCooldownDeferred
from meta_standards_converter.runtime_contracts import RetryCategory
from meta_standards_converter.sources.contracts import RequestMetrics


class FailingSource:
    kind = "geo_accession"

    def __init__(self, error):
        self.error = error
        self.calls = 0

    def probe(self, source, context):
        return str(source).startswith("GSE")

    def validate(self, spec, context):
        pass

    def load(self, spec, context):
        context.service("counted", lambda: self)
        self.calls += 1
        raise self.error

    def metrics(self):
        return RequestMetrics(self.calls, self.calls, self.calls * 0.25)


def http_error(status):
    response = requests.Response()
    response.status_code = status
    return requests.HTTPError("secret response must not persist", response=response)


@pytest.mark.parametrize("error,category", [
    (requests.Timeout("secret query"), RetryCategory.RETRYABLE),
    (http_error(503), RetryCategory.RETRYABLE),
    (http_error(404), RetryCategory.TERMINAL),
    (InputError("invalid_metadata", "bad input"), RetryCategory.TERMINAL),
])
def test_conversion_failure_retains_structured_retry_decision(error, category):
    converter = Converter(handlers=[FailingSource(error)])
    item = converter.convert("GSE1", out_type="json").items[0]
    assert item.error.retry_category == category
    assert item.error.stage == "load"
    assert "secret" not in json.dumps(item.to_dict())
    assert item.metrics.provider_attempts == 1
    assert converter.metrics() == RequestMetrics(1, 1, 0.25)


def test_cooldown_and_repeated_call_metrics_are_not_lost_or_double_counted():
    source = FailingSource(HostRequestCooldownDeferred("example.test", 12345.0))
    converter = Converter(services={"counted": source}, handlers=[source])
    first = converter.convert("GSE1", out_type="json").items[0]
    second = converter.convert("GSE2", out_type="json").items[0]
    assert first.retry_at == second.retry_at == 12345.0
    assert first.error.retry_category == RetryCategory.RETRYABLE
    assert first.metrics == second.metrics == RequestMetrics(1, 1, 0.25)
    assert converter.metrics() == RequestMetrics(2, 2, 0.5)
