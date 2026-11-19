"""Verbose curator summaries are bounded; structured locations remain complete."""
from copy import deepcopy
import json

import pytest

from meta_standards_converter.cli.convert import main
from tests.support.conversion import package


@pytest.mark.parametrize('flags,visible,detailed', [([], False, False), (['-v'], True, False),
    (['-vv'], True, True), (['--quiet'], False, False)])
def test_verbose_curation_report_and_scientific_outputs(tmp_path, capsys, flags, visible, detailed):
    data = package()
    data['series']['pubmed_publication'] = [{'title': 'Study paper'}]
    data['series']['variable'] = [{'factor': 'dose'}]
    data['sample'][0]['library_strategy'] = 'RNA-Seq'
    data['sample'][0]['channel'][0]['characteristics'] = [{'name': 'disease', 'value': 'unknown'}]
    data['sample'] = [dict(deepcopy(data['sample'][0]), iid=f'GSM{i}') for i in range(1, 101)]
    for i, sample in enumerate(data['sample'], 1):
        sample['sra_run'] = [{'run': f'SRR{i}'}]
    source = tmp_path/'input.json'; source.write_text(json.dumps(data))
    report = tmp_path/'report.json'; log = tmp_path/'run.jsonl'
    args = [str(source), '--out-type', 'json', '--out', str(tmp_path/'out'),
        '--enrichment', 'off', '--no-expand-studies', '--progress', 'off',
        '--report-json', str(report), '--log-file', str(log), '--log-format', 'jsonl', *flags]
    assert main(args) == 0
    streams = capsys.readouterr()
    assert streams.out == ''
    assert ('[curation_missing]' in streams.err) == visible
    if visible:
        assert 'reported as unknown in 100 samples' in streams.err
        assert 'organism part missing in 100 samples' in streams.err
        if detailed:
            assert '/sample/99/channel/0/characteristics' in streams.err
        else:
            assert 'GSM1' in streams.err and 'GSM2' in streams.err and 'GSM3' in streams.err
            assert 'GSM4 (' not in streams.err and 'GSM100 (' not in streams.err
            assert '--report-json' in streams.err
            assert 'additional diagnostic groups' not in streams.err
    diagnostics = json.loads(report.read_text())['items'][0]['diagnostics']
    gaps = [d for d in diagnostics if d.get('curation', {}).get('field') == 'sample.organism_part']
    assert len(gaps) == 100 and len({d['location'] for d in gaps}) == 100
    if visible:
        for field in {d['curation']['field'] for d in diagnostics if d.get('curation')}:
            assert field.split('.', 1)[-1].replace('_', ' ') in streams.err
        if not detailed:
            examples = [line for line in streams.err.splitlines() if 'Examples: GSM' in line]
            assert examples and all(line.count('(package ') == 3 for line in examples)
    events = [json.loads(line) for line in log.read_text().splitlines()]
    logged = [e for e in events if e['event_type'] == 'diagnostic' and
              e['attributes'].get('curation', {}).get('field') == 'sample.organism_part']
    assert len(logged) == 100
    assert source.read_text() == json.dumps(data)


def test_json_stdout_and_text_log_retain_advisories(tmp_path, capsys):
    source = tmp_path/'input.json'; source.write_text(json.dumps(package()))
    log = tmp_path/'run.log'
    assert main([str(source), '--out-type', 'json', '--out', str(tmp_path/'out'),
        '--enrichment', 'off', '--no-expand-studies', '--report-json', '-', '--quiet',
        '--log-file', str(log)]) == 0
    captured = capsys.readouterr()
    assert captured.err == ''
    assert any(d.get('curation') for d in json.loads(captured.out)['items'][0]['diagnostics'])
    assert 'curation_missing' in log.read_text()


def test_group_counts_records_separately_from_occurrences_and_keeps_lookup_context():
    from meta_standards_converter.cli.summaries import _Groups
    from meta_standards_converter.converters.unified.contracts import CurationContext, Diagnostic
    groups = _Groups(verbose=True)
    for channel, lookup in enumerate(['failed', 'disabled', 'failed', 'incomplete']):
        d = Diagnostic('curation_missing', 'missing', 'curation', 'warning', action='Supply annotation.',
            location=f'/sample/0/channel/{channel}/molecule', curation=CurationContext(
                'GSE1', 'sample', 'GSM1', 'sample.molecule', 'missing', lookup_status=lookup))
        groups.add(d.to_dict()); groups.add(d.to_dict())
    lines = list(groups.lines())
    assert sum('WARNING' in line for line in lines) == 1
    assert 'missing in 1 sample' in lines[0]
    examples = next(line for line in lines if 'Examples:' in line)
    assert examples.count('(package ') == 3
    assert any('failed' in line and 'disabled' in line and 'incomplete' in line for line in lines)


@pytest.mark.parametrize('flags', [[], ['-v'], ['-vv'], ['--quiet']])
def test_failed_lookup_keeps_exit_status_and_existing_visibility(tmp_path, capsys, monkeypatch, flags):
    import requests
    from unittest.mock import Mock
    from meta_standards_converter.sources.pubmed import PubmedWebFetcher
    fetch = Mock(side_effect=requests.Timeout('private-secret'))
    monkeypatch.setattr(PubmedWebFetcher, 'pubmed_summary', fetch)
    data = package(); data['series']['pubmed_id'] = ['123']
    source = tmp_path/'input.json'; source.write_text(json.dumps(data))
    report = tmp_path/'report.json'
    assert main([str(source), '--out-type', 'json', '--out', str(tmp_path/'out'),
        '--enrichment', 'curators', '--no-expand-studies', '--progress', 'off',
        '--report-json', str(report), *flags]) == 1
    captured = capsys.readouterr()
    assert captured.out == ''
    assert ('source_partial' in captured.err) == ('--quiet' not in flags)
    assert 'private-secret' not in captured.err + report.read_text()
    assert json.loads(report.read_text())['status'] == 'partial'
    fetch.assert_called_once()
