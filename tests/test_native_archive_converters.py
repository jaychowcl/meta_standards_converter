from pathlib import Path
import json
import pytest

from meta_standards_converter.converters.sra2json import SRA2JSONConverter
from meta_standards_converter.converters.ena2json import ENA2JSONConverter
from meta_standards_converter.sources.archive_support import Resolution, StudySeed
from tests.test_native_archive_parsers import fixture_records


class Source:
    def __init__(self, provider='sra', partial=False):
        self.provider, self.partial = provider, partial
        self.calls = []
    def resolve(self, accession):
        primary = 'SRP250911' if self.provider == 'sra' else 'PRJNA609050'
        return Resolution([StudySeed('SRP250911', primary), StudySeed('SRP2', 'SRP2')])
    def fetch(self, seed):
        self.calls.append(seed.study)
        if seed.study == 'SRP2':
            raise OSError('failed study')
        records = fixture_records(self.provider)
        if self.partial:
            records.issues.append('sample retrieval incomplete')
        return records


@pytest.mark.parametrize('provider,cls', [('sra', SRA2JSONConverter), ('ena', ENA2JSONConverter)])
def test_outputs_are_study_named_and_failures_independent(tmp_path, provider, cls):
    result = cls(source=Source(provider)).convert('SRR11192680', out=tmp_path, report_path=tmp_path / 'report.json')
    assert len(result.packages) == 1
    # A failed sibling does not discard a successful package.
    name = 'SRP250911' if provider == 'sra' else 'PRJNA609050'
    output = json.loads((tmp_path / (name + '.json')).read_text())
    assert output['series']['iid'] == name
    assert [r.status for r in result.studies] == ['complete', 'failed']
    assert not result.ok
    assert json.loads((tmp_path / 'report.json').read_text())['studies'][1]['status'] == 'failed'
    assert 'failed study' not in json.dumps(output)


def test_partial_package_written_without_diagnostic_columns(tmp_path):
    result = SRA2JSONConverter(source=Source(partial=True)).convert('SRR11192680', out=tmp_path)
    assert result.studies[0].status == 'partial'
    assert 'sample retrieval incomplete' not in (tmp_path / 'SRP250911.json').read_text()


def test_no_overwrite_and_batch_duplicate_resolution(tmp_path):
    source = Source()
    converter = SRA2JSONConverter(source=source)
    seen = set()
    converter.convert('SRR11192680', out=tmp_path, seen_studies=seen)
    converter.convert('SRX7812918', out=tmp_path, seen_studies=seen)
    assert source.calls == ['SRP250911', 'SRP2']
    output = tmp_path / 'SRP250911.json'
    output.write_text('protected')
    result = converter.convert('SRR11192680', out=tmp_path)
    assert result.studies[0].status == 'failed'
    assert output.read_text() == 'protected'


def test_cli_batch_flags_exit_status_and_report(tmp_path):
    from meta_standards_converter.cli.archive import parser_for, run_cli
    class Converter(SRA2JSONConverter):
        def __init__(self, **kwargs):
            super().__init__(source=Source(), **kwargs)
    assert run_cli(parser_for('SRA'), Converter, ['SRR11192680', 'SRX7812918', '--out', str(tmp_path), '--report', str(tmp_path / 'batch.json')]) == 1
    report = json.loads((tmp_path / 'batch.json').read_text())
    assert len(report['imports']) == 2
    assert len(report['imports'][1]['studies']) == 0


def test_primary_filename_collisions_are_qualified():
    from meta_standards_converter.converters.archive_results import output_name
    seeds = [StudySeed('SRP1', 'PRJNA1'), StudySeed('SRP2', 'PRJNA1')]
    assert output_name(seeds[0], seeds) == 'PRJNA1__SRP1.json'
