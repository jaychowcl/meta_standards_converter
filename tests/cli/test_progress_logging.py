# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Console progress must not alter durable logs or machine-readable results."""
import json
import logging
from io import StringIO
from types import SimpleNamespace

import pytest


def event(kind='stage_started', stage='load', status='running', **attrs):
    return dict(schema_version='1.0', timestamp='now', repository='meta_standards_converter',
                component='test', event_type=kind, level='INFO', run_id='run1',
                stage=stage, status=status, duration_seconds=1, attributes=attrs)


class Terminal(StringIO):
    def isatty(self):
        return True


def test_auto_progress_uses_real_total_and_throttles():
    from meta_standards_converter.cli.logging_output import ProgressRenderer
    stream = Terminal(); now=[0.0]
    renderer = ProgressRenderer(stream, clock=lambda:now[0], width=lambda:90, animate=False)
    renderer.update(event('batch_planned', inputs=10))
    renderer.update(event(item_id='GSE1'))
    before = stream.getvalue()
    renderer.update(event('item_finished', status='failed', finished=1, inputs=10))
    assert stream.getvalue() == before
    now[0] = .25
    renderer.draw()
    assert '1/10 inputs finished' in stream.getvalue()
    assert 'failed=1' in stream.getvalue()
    renderer.close()


def test_unknown_total_plain_and_off_have_no_escape_sequences():
    from meta_standards_converter.cli.logging_output import ProgressRenderer
    stream = StringIO()
    auto = ProgressRenderer(stream, animate=False)
    auto.update(event())
    assert stream.getvalue() == ''
    plain = ProgressRenderer(stream, mode='plain', animate=False)
    plain.update(event())
    plain.update(event('stage_finished', status='partial'))
    assert '[RUN]' in stream.getvalue() and '[PARTIAL]' in stream.getvalue()
    assert '\r' not in stream.getvalue() and '\x1b' not in stream.getvalue()
    assert '%' not in stream.getvalue()
    plain.close()


def test_narrow_terminal_and_log_redraw():
    from meta_standards_converter.cli.logging_output import ProgressRenderer
    stream = Terminal()
    renderer=ProgressRenderer(stream, width=lambda:24, animate=False)
    renderer.update(event(item_id='GSE'+'1'*100))
    renderer.message('WARNING example')
    renderer.close()
    assert 'WARNING example\n' in stream.getvalue()
    assert all(len(line) <= 24 for line in stream.getvalue().split('\r') if line and '\n' not in line)
    assert renderer.closed


@pytest.mark.parametrize('fmt', ['text','jsonl'])
def test_file_keeps_debug_when_console_is_quiet(tmp_path, fmt):
    from meta_standards_converter.cli.common import configure_logging
    from meta_standards_converter.operational_events import OperationalEventEmitter
    path=tmp_path/'events.log'; stream=StringIO()
    args=SimpleNamespace(quiet=True, verbose=0, log_file=str(path), log_format=fmt,
                         log_file_level='DEBUG', progress='plain')
    configure_logging(args, stream=stream)
    logger=logging.getLogger('meta_standards_converter.test')
    logger.debug('detail https://user:secret@example.org/a?token=hidden')
    OperationalEventEmitter(repository='meta_standards_converter').emit('stage_started',component='test',stage='load')
    assert stream.getvalue()==''
    text=path.read_text()
    assert 'detail' in text and 'stage_started' in text
    assert 'secret' not in text and 'hidden' not in text
    if fmt=='jsonl':
        values=[json.loads(line) for line in text.splitlines()]
        assert len(values)==2
        assert all(v['schema_version']=='1.0' for v in values)


def test_cli_quiet_keeps_json_and_no_progress(tmp_path, capsys):
    from meta_standards_converter.cli.convert import main
    from tests.converters.test_json2tsv import package
    path=tmp_path/'source.json'; path.write_text(json.dumps(package()))
    assert main([str(path),'--out-type','json','--outdir',str(tmp_path/'out'),
                 '--enrichment','off','--no-expand-studies','--quiet','--progress','plain','--report-json','-'])==0
    captured=capsys.readouterr()
    assert json.loads(captured.out)['status']=='complete'
    assert captured.err==''


def test_warning_console_groups_but_jsonl_keeps_occurrences(tmp_path):
    from meta_standards_converter.cli.common import configure_logging
    from meta_standards_converter.operational_events import OperationalEventEmitter
    stream=StringIO(); path=tmp_path/'events.jsonl'
    configure_logging(SimpleNamespace(quiet=False,verbose=0,log_file=str(path),log_format='jsonl',progress='off'),stream=stream)
    emitter=OperationalEventEmitter(repository='meta_standards_converter')
    for index in range(4):
        emitter.emit('diagnostic',component='test',level='WARNING',run_id='run',stage='validation',
                     attributes={'item_id':'one','code':'missing_optional','message':f'entity {index}'})
    emitter.emit('conversion_finished',component='test',run_id='run',status='complete')
    assert stream.getvalue().count('entity ')==1
    assert '4 occurrences' in stream.getvalue()
    assert len(path.read_text().splitlines())==5


def test_configure_does_not_remove_caller_handler():
    from meta_standards_converter.cli.common import configure_logging
    logger=logging.getLogger('meta_standards_converter')
    owned=logging.StreamHandler(StringIO()); logger.addHandler(owned)
    configure_logging(SimpleNamespace(quiet=False,verbose=0,log_file=None),stream=StringIO())
    assert owned in logger.handlers


def test_plain_reports_terminal_failure_and_changing_totals():
    from meta_standards_converter.cli.logging_output import ProgressRenderer
    stream = StringIO()
    renderer = ProgressRenderer(stream, mode='plain', animate=False)
    renderer.update(event('batch_planned', inputs=2))
    renderer.update(event('progress', completed=1, total=None, unit='studies'))
    renderer.update(event('batch_planned', inputs=3))
    renderer.update(event('item_finished', status='skipped', item_id='duplicate'))
    renderer.update(event('conversion_finished', status='failed'))
    assert renderer.total == 3
    assert '[SKIP]' in stream.getvalue()
    assert '[FAIL]' in stream.getvalue()
    assert '1 processed' in stream.getvalue()
    assert '\r' not in stream.getvalue()


@pytest.mark.parametrize('mode,quiet', [('off',False), ('auto',True), ('plain',True)])
def test_explicit_off_and_quiet_override_live_terminal(mode, quiet):
    from meta_standards_converter.cli.logging_output import ProgressRenderer
    stream = Terminal()
    renderer = ProgressRenderer(stream, mode=mode, quiet=quiet, animate=False)
    renderer.update(event())
    renderer.update(event('conversion_finished', status='complete'))
    assert stream.getvalue() == ''
    assert renderer.closed


def test_legacy_human_summary_is_stderr(tmp_path, capsys, monkeypatch):
    from meta_standards_converter.cli.geo2json import main
    from meta_standards_converter.converters.geo2json import GEO2JSONConverter
    monkeypatch.setattr(GEO2JSONConverter, 'convert', lambda *a, **kw: [])
    assert main(['GSE1', '--out', str(tmp_path)]) == 0
    captured = capsys.readouterr()
    assert captured.out == ''
    assert 'Conversion: complete' in captured.err


@pytest.mark.parametrize('provider', ['sra', 'ena'])
def test_native_cli_correlates_discovery_and_study_outcomes(provider, tmp_path, capsys, monkeypatch):
    import importlib
    from tests.test_native_archive_converters import Source
    from meta_standards_converter.conversion_logging import emit
    module=importlib.import_module('meta_standards_converter.cli.'+provider+'2json')
    converter_class=getattr(module, provider.upper()+'2JSONConverter')
    class InstrumentedSource(Source):
        def resolve(self, accession):
            emit('resolution_probe', level='DEBUG')
            return super().resolve(accession)
    converter=converter_class(source=InstrumentedSource(provider))
    monkeypatch.setattr(module, provider.upper()+'2JSONConverter', lambda **kwargs: converter)
    path=tmp_path/'events.jsonl'
    assert module.main(['SRR11192680','--out',str(tmp_path/'out'), '--quiet',
                        '--log-file',str(path),'--log-format','jsonl'])==1
    records=[json.loads(line) for line in path.read_text().splitlines()]
    assert len([r for r in records if r['event_type']=='conversion_started'])==1
    assert any(r['event_type']=='resolution_probe' and r['run_id'] for r in records)
    assert any(r['level']=='ERROR' for r in records)
    assert records[-1]['attributes']['partial']==1
    assert capsys.readouterr().out==''


def test_all_display_settings_preserve_artifacts_and_detailed_events(tmp_path, capsys):
    from meta_standards_converter.cli.convert import main
    from tests.converters.test_json2tsv import package
    source=tmp_path/'source.json'; source.write_text(json.dumps(package()))
    contents=[]; milestones=[]
    for verbosity in ([], ['-v'], ['-vv'], ['-q']):
        for mode in ('auto','plain','off'):
            directory=tmp_path/f'run{len(contents)}'
            log=tmp_path/f'run{len(contents)}.jsonl'
            assert main([str(source), '--out-type','json','--out',str(directory),
                         '--enrichment','off','--no-expand-studies','--report-json','-',
                         '--log-file',str(log),'--log-format','jsonl','--progress',mode,*verbosity])==0
            captured=capsys.readouterr()
            assert json.loads(captured.out)['status']=='complete'
            assert '\r' not in captured.err and '\x1b' not in captured.err
            if verbosity==['-q']:
                assert captured.err==''
            contents.append((directory/'GSE1.json').read_bytes())
            records=[json.loads(line) for line in log.read_text().splitlines()]
            milestones.append([(e['event_type'],e['stage'],e['status']) for e in records])
    assert all(value==contents[0] for value in contents)
    assert all(value==milestones[0] for value in milestones)


def test_file_level_is_independent_of_debug_console(tmp_path):
    from meta_standards_converter.cli.common import configure_logging
    stream=StringIO(); path=tmp_path/'filtered.jsonl'
    configure_logging(SimpleNamespace(quiet=False,verbose=2,log_file=str(path),
        log_format='jsonl',log_file_level='ERROR',progress='off'),stream=stream)
    logger=logging.getLogger('meta_standards_converter.test')
    logger.debug('routing detail'); logger.error('failed requested item')
    assert 'routing detail' in stream.getvalue()
    assert [json.loads(line)['level'] for line in path.read_text().splitlines()]==['ERROR']


def test_concurrent_cli_sessions_have_isolated_sinks(tmp_path):
    from concurrent.futures import ThreadPoolExecutor
    from meta_standards_converter.cli.logging_output import cli_entrypoint
    from meta_standards_converter.cli.common import configure_logging
    from meta_standards_converter.conversion_logging import operation, emit
    logger=logging.getLogger('meta_standards_converter')
    original=(logger.level,logger.propagate,list(logger.handlers))
    @cli_entrypoint
    def run(index):
        configure_logging(SimpleNamespace(quiet=True,verbose=0,log_file=str(tmp_path/f'{index}.jsonl'),
            log_format='jsonl',progress='off'),stream=StringIO())
        with operation('test'):
            emit('identity', level='DEBUG', item_id=index)
        return 0
    with ThreadPoolExecutor(2) as pool:
        assert list(pool.map(run, range(2)))==[0,0]
    for index in range(2):
        records=[json.loads(line) for line in (tmp_path/f'{index}.jsonl').read_text().splitlines()]
        assert [r['attributes']['item_id'] for r in records if r['event_type']=='identity']==[index]
        assert len({r['run_id'] for r in records})==1
    assert (logger.level,logger.propagate,list(logger.handlers))==original


def test_interrupted_cli_restores_handlers_and_terminal():
    from meta_standards_converter.cli.logging_output import cli_entrypoint
    from meta_standards_converter.cli.common import configure_logging
    from meta_standards_converter.conversion_logging import operation, stage
    import threading
    stream=Terminal(); logger=logging.getLogger('meta_standards_converter')
    original=(logger.level,logger.propagate,list(logger.handlers))
    @cli_entrypoint
    def run():
        configure_logging(SimpleNamespace(quiet=False,verbose=0,log_file=None,progress='auto'),stream=stream)
        with operation('test'), stage('load'):
            raise KeyboardInterrupt()
    with pytest.raises(KeyboardInterrupt):
        run()
    assert (logger.level,logger.propagate,list(logger.handlers))==original
    assert not any(t.name in {'msc-progress','msc-terminal-progress'} for t in threading.enumerate())
    assert 'interrupted' in stream.getvalue()


def test_usage_failure_during_conversion_does_not_dump_exception(tmp_path, monkeypatch, capsys):
    from meta_standards_converter.cli import convert
    class Broken:
        def convert(self, *args, **kwargs):
            raise ValueError('raw_scientific_payload_canary')
    monkeypatch.setattr(convert,'Converter',Broken)
    with pytest.raises(SystemExit) as raised:
        convert.main(['GSE1','--out-type','json','--out',str(tmp_path)])
    assert raised.value.code==2
    captured=capsys.readouterr()
    assert 'raw_scientific_payload_canary' not in captured.err
    assert 'ValueError' in captured.err
