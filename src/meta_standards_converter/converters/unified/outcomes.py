"""Safe failure and request accounting at the whole-conversion boundary."""
from http.client import IncompleteRead
import math

import requests
from urllib3.exceptions import ProtocolError

from meta_standards_converter.runtime_contracts import RetryCategory, SafeErrorEnvelope
from meta_standards_converter.sources.contracts import RequestMetrics, request_metrics


def failure(error, *, stage, provider, item_id):
    response = getattr(error, "response", None)
    status = getattr(response, "status_code", None)
    transient = isinstance(error, (
        TimeoutError, ConnectionError, IncompleteRead, ProtocolError,
        requests.Timeout, requests.ConnectionError,
        requests.exceptions.ChunkedEncodingError,
        requests.exceptions.ContentDecodingError,
    )) or status in {403, 408, 425, 429} or (isinstance(status, int) and status >= 500)
    trace = getattr(error, "request_trace", None) or {}
    transient = transient or any(
        item.get("retry_category") == "retryable"
        for item in trace.get("errors", []) if isinstance(item, dict)
    )
    envelope = SafeErrorEnvelope.from_exception(
        error, provider=provider, stage=stage, item_id=item_id,
        http_status=status,
        retry_category=RetryCategory.RETRYABLE if transient else RetryCategory.TERMINAL,
    )
    retry_at = getattr(error, "retry_at", getattr(response, "_thematic_atlases_retry_at", None))
    if not isinstance(retry_at, (int, float)) or not math.isfinite(retry_at):
        retry_at = None
    return envelope, retry_at


def snapshot(context):
    services = {id(value): value for value in (
        *context.converter.services.values(), *context._services.values(),
    )}
    return {key: request_metrics(value) for key, value in services.items()}


def delta(before, after):
    fields = ("provider_attempts", "retry_count", "rate_wait_seconds")
    return RequestMetrics(*(sum(
        max(0, getattr(value, field) - getattr(before.get(key, RequestMetrics()), field))
        for key, value in after.items()
    ) for field in fields))


def add(left, right):
    return RequestMetrics(
        left.provider_attempts + right.provider_attempts,
        left.retry_count + right.retry_count,
        left.rate_wait_seconds + right.rate_wait_seconds,
    )
