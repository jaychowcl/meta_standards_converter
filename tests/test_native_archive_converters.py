# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
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


def test_cli_qualifies_shared_primary_across_separate_inputs(tmp_path):
    from meta_standards_converter.cli.archive import parser_for, run_cli
    class MultiSource:
        def resolve(self, accession):
            return Resolution([StudySeed(accession, 'PRJNA1')])
        def fetch(self, seed):
            records = fixture_records('ena'); records.seed = seed
            return records
    class Converter(ENA2JSONConverter):
        def __init__(self, **kwargs):
            super().__init__(source=MultiSource(), **kwargs)
    run_cli(parser_for('ENA'), Converter, ['SRP250911', 'SRP2', '--out', str(tmp_path), '--overwrite'])
    assert (tmp_path / 'PRJNA1__SRP250911.json').exists()
    assert (tmp_path / 'PRJNA1__SRP2.json').exists()


def test_optional_evidence_contains_original_bytes_only_when_requested(tmp_path):
    from meta_standards_converter.sources.archive_support import ArchiveHTTP
    class Response:
        headers = {}
        content = b'<STUDY_SET/>'
        def raise_for_status(self): pass
        def close(self): pass
        def iter_content(self, **kwargs): yield self.content
    class Requester:
        def get(self, *args, **kwargs): return Response()
    http = ArchiveHTTP('ena_portal', requester=Requester())
    http.get('https://www.ebi.ac.uk/ena/browser/api/xml/ERP1')
    assert not list(tmp_path.iterdir())
    http.evidence_dir = tmp_path
    http.get('https://www.ebi.ac.uk/ena/browser/api/xml/ERP1')
    assert next(tmp_path.iterdir()).read_bytes() == b'<STUDY_SET/>'


def test_optional_service_failure_does_not_discard_native_package():
    class Broken:
        def convert(self, *args, **kwargs): raise OSError('unavailable')
        def enrich(self, *args, **kwargs): raise OSError('unavailable')
    converter = SRA2JSONConverter(source=Source(), peer_converter=Broken(), linked_enricher=Broken())
    result = converter.convert('SRR11192680', include_peer=True, enrich_from_geo_ae=True)
    assert result.studies[0].status == 'partial'
    assert result.packages[0].series.iid == 'SRP250911'


def test_cli_evidence_flag_also_covers_accession_resolution(tmp_path):
    from types import SimpleNamespace
    from meta_standards_converter.cli.archive import parser_for, run_cli
    evidence = tmp_path / 'evidence'
    class EvidenceSource(Source):
        http = SimpleNamespace(evidence_dir=None)
        def resolve(self, accession):
            assert self.http.evidence_dir == evidence
            return super().resolve(accession)
    class Converter(SRA2JSONConverter):
        def __init__(self, **kwargs): super().__init__(source=EvidenceSource(), **kwargs)
    run_cli(parser_for('SRA'), Converter, ['SRR11192680', '--out', str(tmp_path), '--evidence-dir', str(evidence)])
    assert (tmp_path / 'SRP250911.json').exists()
