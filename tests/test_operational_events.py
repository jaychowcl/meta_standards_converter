# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
import json
import logging
from unittest.mock import Mock

from meta_standards_converter.helpers.request_helper import RateLimitedRequester, RequestSettings
from meta_standards_converter.operational_events import OperationalEventEmitter


def _response():
    response = Mock(status_code=200, headers={}, content=b"private payload")
    return response


def test_emitter_redacts_jsonl_log_and_exports_both_metrics(tmp_path, caplog):
    emitter = OperationalEventEmitter(
        repository="meta_standards_converter", jsonl_path=tmp_path / "events.jsonl"
    )
    with caplog.at_level(logging.INFO):
        emitter.emit(
            "provider_request",
            component="ena",
            status="success",
            duration_seconds=0.1,
            attributes={"api_key": "secret", "url": "https://x.test/a?token=secret"},
        )
    event = json.loads((tmp_path / "events.jsonl").read_text())
    assert event["attributes"]["api_key"] == "[REDACTED]"
    assert event["attributes"]["url"] == "https://x.test/a"
    assert "secret" not in caplog.text
    assert emitter.export_json(tmp_path / "metrics.json")["events"]["provider_request|success"]["count"] == 1
    assert "operational_events_total" in emitter.export_prometheus(tmp_path / "metrics.prom")


def test_requester_emits_safe_terminal_request_event(tmp_path):
    emitter = OperationalEventEmitter(
        repository="meta_standards_converter", jsonl_path=tmp_path / "events.jsonl"
    )
    RateLimitedRequester(
        "ena_portal",
        settings=RequestSettings(request_delay=0, max_retries=0),
        get=Mock(return_value=_response()),
        event_emitter=emitter,
    ).get("https://example.org/data?token=secret")

    text = (tmp_path / "events.jsonl").read_text()
    assert "private payload" not in text
    assert "secret" not in text
    event = json.loads(text)
    assert event["component"] == "ena_portal"
    assert event["attributes"]["host"] == "example.org"
