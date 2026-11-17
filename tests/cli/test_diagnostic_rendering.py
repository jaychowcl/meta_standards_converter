# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Compact cause summaries retain full evidence and survive failed sinks."""
import io
import json
import logging
import pytest
from types import SimpleNamespace

from meta_standards_converter.cli import convert
from meta_standards_converter.converters import Converter
from tests.converters.test_ae_fallback import services, MAPPINGS
from tests.converters.test_json2tsv import package


def test_fallback_console_groups_one_cause(monkeypatch, capsys, tmp_path):
    deps, _ = services(*MAPPINGS[0])
    monkeypatch.setattr(convert, 'Converter', lambda: Converter(services=deps))
    assert convert.main([MAPPINGS[0][0], '--out-type', 'json', '--out', str(tmp_path),
                         '--enrichment', 'off', '--no-expand-studies', '--progress', 'off']) == 1
    stderr = capsys.readouterr().err
    assert stderr.count('BioStudies reports Study not accessible') == 1
    assert 'partial output written' in stderr
    assert MAPPINGS[0][1] in stderr
    assert 'provider_issue' not in stderr


def test_summary_bounds_groups_without_merging_distinct_causes(tmp_path):
    from meta_standards_converter.cli.common import configure_logging
    from meta_standards_converter.cli.logging_output import cli_entrypoint
    from meta_standards_converter.conversion_logging import operation, item_scope
    from meta_standards_converter.converters.unified.contracts import Diagnostic, ConversionItemResult
    stream = io.StringIO()
    log = tmp_path/'events.jsonl'
    @cli_entrypoint
    def run():
        configure_logging(SimpleNamespace(quiet=False, verbose=0, log_file=str(log), log_format='jsonl', progress='off'), stream=stream)
        item = ConversionItemResult('one', 'GSE1', status='partial', execution='succeeded', completeness='partial',
            artifacts={'json': 'GSE1.json'}, diagnostics=[Diagnostic('source_partial', f'detail-{i}', severity='warning', correlation_id=str(i)) for i in range(8)])
        with operation('test', heartbeat=False), item_scope(item, total=1, index=1):
            pass
        return 1
    assert run() == 1
    assert 'detail-4' in stream.getvalue()
    assert 'detail-5' not in stream.getvalue()
    assert '3 additional diagnostic groups' in stream.getvalue()
    events = [json.loads(line) for line in log.read_text().splitlines()]
    assert len([e for e in events if e['event_type']=='diagnostic']) == 8


def test_log_write_failure_retains_conversion_without_exception_dump(tmp_path, monkeypatch, capsys):
    source = tmp_path/'input.json'
    source.write_text(json.dumps(package()))
    class BrokenStream(io.StringIO):
        def write(self, value):
            raise OSError('private_sink_canary')
    monkeypatch.setattr(logging.FileHandler, '_open', lambda self: BrokenStream())
    code = convert.main([str(source), '--out-type', 'json', '--out', str(tmp_path/'out'),
        '--enrichment', 'off', '--no-expand-studies', '--log-file', str(tmp_path/'events.log')])
    assert code == 1
    assert list((tmp_path/'out').glob('*.json'))
    stderr = capsys.readouterr().err
    assert stderr.count('[reporting_failed]') == 1
    assert 'private_sink_canary' not in stderr and 'Traceback' not in stderr


def test_report_failure_changes_command_summary_but_retains_outputs(tmp_path, monkeypatch, capsys):
    source = tmp_path/'input.json'
    source.write_text(json.dumps(package()))
    def broken(*args, **kwargs):
        raise PermissionError('private_report_canary')
    monkeypatch.setattr(convert, 'atomic_write', broken)
    assert convert.main([str(source), '--out-type', 'json', '--out', str(tmp_path/'out'),
        '--enrichment', 'off', '--no-expand-studies', '--report-json', str(tmp_path/'report.json')]) == 1
    stderr = capsys.readouterr().err
    assert 'Conversion: complete' not in stderr
    assert 'Conversion: partial; reporting failed' in stderr
    assert list((tmp_path/'out').glob('*.json'))
    assert 'private_report_canary' not in stderr


def test_log_close_failure_is_accounted_for_before_final_summary(tmp_path, monkeypatch, capsys):
    source = tmp_path/'input.json'
    source.write_text(json.dumps(package()))
    class BrokenClose(io.StringIO):
        def close(self):
            super().close()
            raise OSError('private_close_canary')
    monkeypatch.setattr(logging.FileHandler, '_open', lambda self: BrokenClose())
    assert convert.main([str(source), '--out-type', 'json', '--out', str(tmp_path/'out'),
        '--enrichment', 'off', '--no-expand-studies', '--log-file', str(tmp_path/'events.log')]) == 1
    stderr = capsys.readouterr().err
    assert stderr.count('[reporting_failed]') == 1
    assert 'Conversion: complete' not in stderr
    assert 'private_close_canary' not in stderr
    assert list((tmp_path/'out').glob('*.json'))


def test_migration_failure_has_final_summary(tmp_path, capsys):
    from meta_standards_converter.cli.miniml_migrate import main
    assert main([str(tmp_path/'missing.json'), str(tmp_path/'out.json')]) == 1
    stderr = capsys.readouterr().err
    assert stderr.count('[missing_file]') == 1
    assert 'Conversion: failed' in stderr


def test_optional_pubmed_failure_is_reported_once(tmp_path, monkeypatch, capsys):
    import requests
    from unittest.mock import Mock
    from meta_standards_converter.metadata.enrichment import MINiMLEnricher
    data = package()
    data['series']['pubmed_id'] = ['123']
    source = tmp_path/'input.json'
    source.write_text(json.dumps(data))
    fetcher = Mock()
    fetcher.pubmed_summary.side_effect = requests.Timeout('private_provider_canary')
    enricher = MINiMLEnricher(pubmed_fetcher=fetcher, insdc_fetcher=Mock())
    monkeypatch.setattr(convert, 'Converter', lambda: Converter(services={'enricher': enricher}))
    assert convert.main([str(source), '--out-type', 'json', '--out', str(tmp_path/'out'),
        '--no-expand-studies']) == 1
    stderr = capsys.readouterr().err
    assert stderr.count('PubMed 123: Timeout') == 1
    assert 'private_provider_canary' not in stderr
    assert fetcher.pubmed_summary.call_count == 1


def test_report_initialization_failure_precedes_conversion(tmp_path, monkeypatch, capsys):
    import tempfile
    from unittest.mock import Mock
    source = tmp_path/'input.json'
    source.write_text(json.dumps(package()))
    conversion = Mock()
    monkeypatch.setattr(convert, 'Converter', conversion)
    def denied(*args, **kwargs):
        raise PermissionError('private_initialization_canary')
    monkeypatch.setattr(tempfile, 'TemporaryFile', denied)
    assert convert.main([str(source), '--out-type', 'json', '--out', str(tmp_path/'out'),
        '--report-json', str(tmp_path/'report.json')]) == 1
    conversion.assert_not_called()
    assert not (tmp_path/'out').exists()
    stderr = capsys.readouterr().err
    assert 'private_initialization_canary' not in stderr
    assert 'reporting_failed' in stderr


def test_interrupted_item_is_not_reported_as_a_provider_failure(tmp_path, monkeypatch, capsys):
    from tests.converters.test_conversion_outcomes import FailingSource
    monkeypatch.setattr(convert, 'Converter', lambda: Converter(handlers=[FailingSource(KeyboardInterrupt())]))
    with pytest.raises(KeyboardInterrupt):
        convert.main(['GSE1', '--out-type', 'json', '--out', str(tmp_path)])
    stderr = capsys.readouterr().err
    assert 'interrupted' in stderr
    assert 'GSE1: failed' not in stderr
    assert 'provider_issue' not in stderr
