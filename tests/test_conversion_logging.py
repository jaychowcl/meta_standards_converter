"""Operation telemetry is independent of display and scientific payloads."""
import json
import logging
from concurrent.futures import ThreadPoolExecutor
from io import StringIO

import pytest

from meta_standards_converter.operational_events import OperationalEventEmitter


def events(caplog):
    return [r.msc_event for r in caplog.records if hasattr(r, "msc_event")]


def test_event_payload_is_attached_and_credentials_are_redacted(caplog):
    with caplog.at_level(logging.DEBUG):
        OperationalEventEmitter(repository="meta_standards_converter").emit(
            "request", component="provider", level="DEBUG",
            attributes={"url": "https://user:secret@example.org/a?token=hidden#private",
                        "message": "retry https://user:secret@example.org/a?token=hidden"})
    event = events(caplog)[0]
    assert event["attributes"]["url"] == "https://example.org/a"
    assert not any(secret in caplog.text for secret in ("secret", "hidden", "private"))


def test_nested_conversion_has_one_lifecycle_and_restores_context(caplog):
    from meta_standards_converter.conversion_logging import logged_conversion, stage, current
    @logged_conversion
    def inner():
        with stage("load") as span:
            span.update(samples=3)
        return []
    @logged_conversion
    def outer():
        return inner()
    with caplog.at_level(logging.DEBUG):
        assert outer() == []
    recorded = events(caplog)
    assert [e["event_type"] for e in recorded].count("conversion_started") == 1
    assert [e["event_type"] for e in recorded].count("conversion_finished") == 1
    loaded = next(e for e in recorded if e["event_type"] == "stage_finished")
    assert loaded["attributes"]["samples"] == 3
    assert loaded["duration_seconds"] >= 0
    assert len({e["run_id"] for e in recorded}) == 1
    assert current() is None


def test_failure_is_reported_once_without_raw_exception(caplog):
    from meta_standards_converter.conversion_logging import logged_conversion, stage, current
    @logged_conversion
    def fail():
        with stage("projection"):
            raise ValueError("secret-payload")
    with caplog.at_level(logging.DEBUG), pytest.raises(ValueError, match="secret-payload"):
        fail()
    assert "secret-payload" not in caplog.text
    failures = [e for e in events(caplog) if e["level"] == "ERROR"]
    assert len(failures) == 1
    assert failures[0]["stage"] == "projection"
    assert current() is None


def test_parallel_operations_and_explicit_worker_context(caplog):
    from meta_standards_converter.conversion_logging import logged_conversion, emit, context_callable
    @logged_conversion
    def convert():
        with ThreadPoolExecutor(1) as pool:
            pool.submit(context_callable(lambda: emit("worker", level="DEBUG"))).result()
    with caplog.at_level(logging.DEBUG), ThreadPoolExecutor(2) as pool:
        list(pool.map(lambda _: convert(), range(2)))
    recorded = events(caplog)
    starts = {e["run_id"] for e in recorded if e["event_type"] == "conversion_started"}
    assert len(starts) == 2
    assert {e["run_id"] for e in recorded if e["event_type"] == "worker"} == starts


def test_heartbeat_is_throttled_with_fake_clock(caplog):
    from meta_standards_converter.conversion_logging import operation, stage
    now = [0.0]
    with caplog.at_level(logging.DEBUG), operation("test", clock=lambda: now[0], heartbeat=False) as run:
        with stage("load"):
            now[0] = 29
            run.tick()
            now[0] = 30
            run.tick()
            run.tick()
            now[0] = 60
            run.tick()
    assert len([e for e in events(caplog) if e["event_type"] == "progress"]) == 2


def test_debug_statistics_collected_without_enabled_logger():
    from meta_standards_converter.conversion_logging import operation, count
    log = logging.getLogger("meta_standards_converter")
    log.setLevel(logging.CRITICAL)
    with operation("test", heartbeat=False) as run:
        count("cache_hits", 2)
        count("provider_attempts")
    assert run.counters == {"cache_hits": 2, "provider_attempts": 1}


@pytest.mark.parametrize("exception", [KeyboardInterrupt, SystemExit])
def test_interruption_cleans_up_context(exception, caplog):
    from meta_standards_converter.conversion_logging import operation, current
    with caplog.at_level(logging.DEBUG), pytest.raises(exception):
        with operation("test", heartbeat=False):
            raise exception()
    assert current() is None
    assert events(caplog)[-1]["status"] == "interrupted"


def test_unified_stage_order_and_final_item_counts(caplog, tmp_path):
    from meta_standards_converter.converters import Converter, InputSpec
    from tests.converters.test_json2tsv import package
    with caplog.at_level(logging.DEBUG):
        result = Converter().convert([InputSpec(package(), id="good"), tmp_path / 'missing.json'],
            out_type="json", enrichment="off", options={"expand_studies": False})
    assert result.status == "partial"
    recorded = events(caplog)
    stages = [e["stage"] for e in recorded if e["event_type"] == "stage_started"]
    assert stages.index("discovery") < stages.index("load") < stages.index("preparation") < stages.index("conversion")
    finished = [e for e in recorded if e["event_type"] == "item_finished"]
    assert [e["status"] for e in finished] == ["complete", "failed"]
    assert recorded[-1]["attributes"]["complete"] == 1
    assert recorded[-1]["attributes"]["failed"] == 1
    assert any(e["level"] == "ERROR" and e["attributes"].get("item_id") == "input_2" for e in recorded)


def test_optional_failure_is_warning_and_reports_retained_output(caplog):
    from meta_standards_converter.converters import Converter
    from tests.converters.test_converter_preparation import native, services
    with caplog.at_level(logging.DEBUG):
        result = Converter(services=services([], failure=True)).convert(native(),
            out_type="json", enrichment="curators", options={"expand_studies":False})
    assert result.items[0].execution == "succeeded"
    assert result.status == "partial"
    recorded = events(caplog)
    assert any(e["level"] == "WARNING" and e["attributes"].get("code") == "source_partial" for e in recorded)
    assert not any(e["level"] == "ERROR" for e in recorded)


def test_fail_fast_finishes_every_input(caplog, tmp_path):
    from meta_standards_converter.converters import Converter, InputSpec
    from tests.converters.test_json2tsv import package
    with caplog.at_level(logging.DEBUG):
        result = Converter().convert([tmp_path/'missing.json', InputSpec(package())],
            out_type="json", options={"fail_fast":True})
    assert [i.status for i in result.items] == ["failed", "skipped"]
    assert [e["status"] for e in events(caplog) if e["event_type"] == "item_finished"] == ["failed", "skipped"]


@pytest.mark.parametrize('provider', ['sra', 'ena'])
def test_native_partial_results_log_failed_study(provider, caplog):
    from tests.test_native_archive_converters import Source, SRA2JSONConverter, ENA2JSONConverter
    cls = SRA2JSONConverter if provider == 'sra' else ENA2JSONConverter
    with caplog.at_level(logging.DEBUG):
        result = cls(source=Source(provider)).convert('SRR11192680')
    assert [s.status for s in result.studies] == ['complete','failed']
    recorded = events(caplog)
    assert len([e for e in recorded if e['event_type']=='conversion_started']) == 1
    assert any(e['level']=='ERROR' for e in recorded)
    assert recorded[-1]['attributes']['studies'] == 2


def test_publication_collision_reports_publication_not_conversion(caplog, tmp_path):
    from meta_standards_converter.converters import Converter
    from tests.converters.test_json2tsv import package
    output = tmp_path/'protected.json'
    output.write_text('original')
    with caplog.at_level(logging.DEBUG):
        result = Converter().convert(package(), outfile=output, out_type='json',
            enrichment='off', options={'expand_studies':False})
    assert result.status == 'failed'
    assert output.read_text() == 'original'
    assert any(e['stage']=='publication' and e['level']=='ERROR' for e in events(caplog))
