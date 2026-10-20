"""Export telemetry and work reuse must preserve the actual MAGE-TAB contract."""
import logging

import pytest

from meta_standards_converter.conversion_logging import operation
from meta_standards_converter.converters.json2ae import JSON2AEConverter
from meta_standards_converter.magetab.constructor import AEConstructor
from meta_standards_converter.miniml import MINiMLCodec
from meta_standards_converter.sources.json import DatasetPackageGroup, SourceLoadResult
from tests.test_archive_export_cleanup import legacy
from tests.test_archive_publications import NoRuns, PubMed
from tests.test_conversion_logging import events


def test_export_stages_report_work_cpu_and_process_peak_without_changing_output(caplog, tmp_path):
    package = MINiMLCodec().decode(legacy()).package
    source = SourceLoadResult((DatasetPackageGroup("study", (package,)),))
    converter = JSON2AEConverter(
        ae_constructor=AEConstructor(pubmed_client=PubMed(), insdc_client=NoRuns())
    )
    expected = converter.convert_loaded(source, enrich=False)
    with caplog.at_level(logging.DEBUG):
        result = converter.convert_loaded(source, enrich=False, out=str(tmp_path))
    assert result == expected
    phases = [e for e in events(caplog)
              if e["event_type"] == "stage_finished" and
              str(e["stage"]).startswith("magetab_")]
    names = {e["stage"] for e in phases}
    assert {
        "magetab_package_resolution", "magetab_model_reconstruction",
        "magetab_normalization", "magetab_protocol_matching",
        "magetab_file_projection", "magetab_path_construction",
        "magetab_column_planning", "magetab_row_rendering",
        "magetab_semantic_overlay", "magetab_validation", "magetab_writing",
    } <= names
    for event in phases:
        metrics = event["attributes"]
        assert metrics["wall_seconds"] >= metrics["exclusive_instrumented_wall_seconds"] >= 0
        assert metrics["process_cpu_seconds"] >= metrics["exclusive_instrumented_cpu_seconds"] >= 0
        assert metrics["peak_rss_scope"] == "process_lifetime_high_water"
        assert metrics["peak_rss_bytes"] is None or metrics["peak_rss_bytes"] > 0
    assert len(list(tmp_path.glob("*.idf.txt"))) == 1
    assert len(list(tmp_path.glob("*.sdrf.txt"))) == 1


def test_export_phase_excludes_nested_instrumented_work_and_restores_failure_context(caplog, monkeypatch):
    from meta_standards_converter.magetab import _performance

    clock = [0.0]
    cpu = [0.0]
    monkeypatch.setattr(_performance.time, "perf_counter", lambda: clock[0])
    monkeypatch.setattr(_performance.time, "process_time", lambda: cpu[0])
    with caplog.at_level(logging.DEBUG), operation("test", heartbeat=False):
        with _performance.phase("outer"):
            clock[0], cpu[0] = 2, 1
            with pytest.raises(ValueError), _performance.phase("inner"):
                clock[0], cpu[0] = 5, 3
                raise ValueError("failed projection")
            clock[0], cpu[0] = 9, 5
        with _performance.phase("next"):
            clock[0], cpu[0] = 10, 6
    completed = {e["stage"]: e for e in events(caplog) if e["event_type"] == "stage_finished"}
    assert completed["inner"]["status"] == "failed"
    assert completed["outer"]["attributes"]["wall_seconds"] == 9
    assert completed["outer"]["attributes"]["exclusive_instrumented_wall_seconds"] == 6
    assert completed["outer"]["attributes"]["exclusive_instrumented_cpu_seconds"] == 3
    assert completed["next"]["attributes"]["exclusive_instrumented_wall_seconds"] == 1
