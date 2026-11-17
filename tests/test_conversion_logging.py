# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Operation telemetry is independent of display and scientific payloads."""
from meta_standards_converter.converters import PreparationPolicy, RuntimeOptions
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


def test_statistics_count_typed_runs_without_serializing_samples(monkeypatch):
    from meta_standards_converter.conversion_logging import statistics
    from meta_standards_converter.miniml.model import MINiMLPackage, Series, SourceInfo, Sample, SRARun
    package=MINiMLPackage(Series(iid='ERP1'),SourceInfo('ENA'),samples=(
        Sample(iid='SAMEA1',sra_runs=(SRARun(run='ERR1'),SRARun(run='ERR2'))),))
    def forbidden(self):raise AssertionError('telemetry must not materialize sample metadata')
    monkeypatch.setattr(Sample,'to_mapping',forbidden)
    assert statistics(package)=={'packages':1,'samples':1,'runs':2}


def test_native_processing_exposes_separate_stage_timings(caplog,tmp_path):
    from meta_standards_converter.conversion_logging import operation
    from meta_standards_converter.miniml.sra_parser import SRAParser
    from meta_standards_converter.metadata.archive_enrichment import merge_archive_metadata
    from meta_standards_converter.metadata.session import PreparationSession
    from meta_standards_converter.converters import Converter
    from tests.test_native_archive_parsers import fixture_records
    with caplog.at_level(logging.DEBUG),operation('stages',heartbeat=False):
        package=SRAParser().parse(fixture_records('sra'))
        merged,_=merge_archive_metadata(package,package)
        with PreparationSession() as session:session.retrieve('native',lambda:merged)
        result=Converter().convert(merged, out_type='json', outdir=tmp_path, preparation=PreparationPolicy(enrichment='off'))
    assert result.items[0].execution=='succeeded'
    timings={e['stage']:e['duration_seconds'] for e in events(caplog) if e['event_type']=='stage_finished'}
    assert {'indexing','entity_construction','supplementary_projection','merging',
            'copying','source_evidence_copying','residual_reconciliation','model_decoding','validation','export'}<=timings.keys()
    assert all(duration>=0 for duration in timings.values())


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
        result = Converter().convert([InputSpec(package(), id='good'), tmp_path / 'missing.json'], out_type='json', preparation=PreparationPolicy(enrichment='off', expand_studies=False))
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
        result = Converter(services=services([], failure=True)).convert(native(), out_type='json', preparation=PreparationPolicy(enrichment='curators', expand_studies=False))
    assert result.items[0].execution == "succeeded"
    assert result.status == "partial"
    recorded = events(caplog)
    assert any(e["level"] == "WARNING" and e["attributes"].get("code") == "source_partial" for e in recorded)
    assert not any(e["level"] == "ERROR" for e in recorded)


def test_fail_fast_finishes_every_input(caplog, tmp_path):
    from meta_standards_converter.converters import Converter, InputSpec
    from tests.converters.test_json2tsv import package
    with caplog.at_level(logging.DEBUG):
        result = Converter().convert([tmp_path / 'missing.json', InputSpec(package())], out_type='json', runtime=RuntimeOptions(fail_fast=True))
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
        result = Converter().convert(package(), outfile=output, out_type='json', preparation=PreparationPolicy(enrichment='off', expand_studies=False))
    assert result.status == 'failed'
    assert output.read_text() == 'original'
    assert any(e['stage']=='publication' and e['level']=='ERROR' for e in events(caplog))


def test_failure_after_first_published_artifact_is_truthful(caplog, tmp_path, monkeypatch):
    from meta_standards_converter.converters import Converter
    from meta_standards_converter.converters import archive_results
    from tests.converters.test_json2tsv import package
    original = archive_results.publish_json
    calls = []
    def fail_second(*args, **kwargs):
        calls.append(args[0])
        if len(calls) == 2:
            raise OSError('private storage detail')
        return original(*args, **kwargs)
    monkeypatch.setattr(archive_results, 'publish_json', fail_second)
    with caplog.at_level(logging.DEBUG):
        result = Converter().convert([package(), package('GSE2', 'GSM2')], outdir=tmp_path, out_type='json', preparation=PreparationPolicy(enrichment='off', expand_studies=False))
    assert result.status=='failed'
    assert len(result.items[0].artifacts)==1
    assert (tmp_path/'GSE1.json').is_file()
    assert not (tmp_path/'GSE2.json').exists()
    assert events(caplog)[-1]['attributes']['artifacts']==1
    assert 'private storage detail' not in caplog.text


def test_context_statistics_do_not_turn_missing_optional_data_into_errors(caplog):
    from meta_standards_converter.converters import Converter
    from tests.converters.test_json2tsv import package
    with caplog.at_level(logging.DEBUG):
        result=Converter().convert(package(), out_type='magetab', preparation=PreparationPolicy(enrichment='off', expand_studies=False))
    assert result.status=='complete'
    assert not [r for r in caplog.records if r.levelno>=logging.WARNING]
    assert any(e['attributes'].get('sdrf_rows')==1 for e in events(caplog))


def test_publication_progress_only_advances_after_commit(caplog, tmp_path):
    from meta_standards_converter.conversion_logging import operation
    from meta_standards_converter.converters.archive_results import publish_json
    with caplog.at_level(logging.DEBUG), operation('test', heartbeat=False) as run:
        publish_json(tmp_path/'one.json', {'safe': True})
        with pytest.raises(FileExistsError):
            publish_json(tmp_path/'one.json', {})
    assert run.counters['artifacts_committed'] == 1
    completed = [e for e in events(caplog) if e['event_type']=='progress' and e['attributes'].get('unit')=='artifacts']
    assert len(completed) == 1
    assert completed[0]['attributes']['completed'] == 1


def test_nested_context_restored_and_semantic_phases_visible_at_info(caplog):
    from meta_standards_converter.conversion_logging import logged_conversion, bind, emit, stage
    @logged_conversion
    def child():
        bind(provider='ena', dataset_id='other')
        with stage('publication'):
            with stage('serialization'):
                pass
    @logged_conversion
    def parent():
        bind(provider='geo', dataset_id='GSE1')
        with stage('conversion'):
            child()
        emit('after_child')
    with caplog.at_level(logging.DEBUG):
        parent()
    recorded=events(caplog)
    assert next(e for e in recorded if e['event_type']=='after_child')['attributes']['provider']=='geo'
    phases=[e for e in recorded if e['event_type']=='stage_started']
    assert next(e for e in phases if e['stage']=='publication')['level']=='INFO'
    assert next(e for e in phases if e['stage']=='serialization')['level']=='DEBUG'


def test_repeated_diagnostics_retained_once_per_result_traversal(caplog):
    from types import SimpleNamespace
    from meta_standards_converter.conversion_logging import operation, result_diagnostics
    from meta_standards_converter.converters.unified.contracts import Diagnostic
    result=SimpleNamespace(diagnostics=[Diagnostic('conflict','same issue','preparation','warning')]*2)
    with caplog.at_level(logging.DEBUG), operation('test', heartbeat=False) as run:
        result_diagnostics(result)
        result_diagnostics(result)
    assert len([e for e in events(caplog) if e['event_type']=='diagnostic'])==2
    assert run.counters['warnings']==2


def test_redaction_removes_authorization_schemes_and_quoted_secrets():
    from meta_standards_converter.operational_events import redact
    value=redact('Authorization: Bearer hidden_bearer Authorization=Basic hidden_basic password="hidden quoted" token=hidden_token')
    assert 'hidden' not in value
    assert '[REDACTED]' in value


def test_output_counts_do_not_claim_cached_samples_were_executed():
    from meta_standards_converter.conversion_logging import statistics
    from meta_standards_converter.expression.catalogue import ConversionResult
    result=ConversionResult('GSE1',sample_h5ads={'GSM1':'/already/committed.h5ad'})
    summary=statistics(result)
    assert summary['sample_outputs']==1
    assert 'samples_converted' not in summary
