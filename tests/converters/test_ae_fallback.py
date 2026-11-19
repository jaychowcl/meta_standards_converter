# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Fallback is accession loading, independent of optional preparation policy."""
import json
import pytest
from meta_standards_converter.converters import Converter, PreparationPolicy, RuntimeOptions
from meta_standards_converter.sources.errors import MetadataSourceError
from tests.support.ae import MAPPINGS

from tests.support.ae import MAPPINGS, services


@pytest.mark.parametrize('mapping', MAPPINGS)
@pytest.mark.parametrize('preset', ['off', 'standard', 'curators'])
@pytest.mark.parametrize('expand', [False, True])
def test_verified_fallback_is_valid_partial_for_every_policy(mapping, preset, expand, tmp_path):
    deps, original = services(*mapping)
    result = Converter(services=deps).convert(mapping[0], out_type='json', outdir=tmp_path,
        preparation=PreparationPolicy(enrichment=preset, expand_studies=expand))
    item = result.items[0]
    assert (item.execution, item.validation, item.completeness) == ('succeeded', 'valid', 'partial'), result.to_dict()
    assert item.provider == 'ena' and item.source == mapping[0] and item.artifacts
    assert deps['ae_source'].load.call_count == 1  # Linked enrichment reuses the denial.
    deps['ena_source'].resolve_arrayexpress.assert_called_once_with(mapping[0])
    assert any(r['operation'] == 'source_fallback' for r in item.preparation)
    data = item.payload[0].to_mapping()
    assert data['source']['format'] == 'ENA'
    evidence = next(r for r in data['extensions']['insdc']['records'] if r['kind'] == 'source_fallback')
    assert evidence['accession'] == mapping[0] and evidence['metadata']['study'] == mapping[1]
    assert original.get('extensions') is None  # injected package is detached
    saved = next(iter(item.artifacts.values()))
    again = Converter(services=deps).convert(saved, out_type='magetab', outdir=tmp_path / 'tables',
        preparation=PreparationPolicy(enrichment='off', expand_studies=False))
    assert again.items[0].completeness == 'partial' and again.items[0].validation == 'valid', again.to_dict()
    idf = next(path for path in again.items[0].artifacts.values() if path.endswith('.idf.txt'))
    assert 'MSCSourceFallback' in open(idf).read()


@pytest.mark.parametrize('failure', ['source_not_found', 'missing_idf', 'missing_sdrf'])
def test_missing_source_metadata_allows_fallback(failure):
    deps, _ = services(*MAPPINGS[0], failure=failure)
    result = Converter(services=deps).convert(MAPPINGS[0][0], out_type='json', preparation=PreparationPolicy(enrichment='off', expand_studies=False))
    assert result.items[0].execution == 'succeeded'


@pytest.mark.parametrize('failure', ['ambiguous_idf', 'ambiguous_sdrf', 'provider_response_invalid'])
def test_nonavailability_failures_do_not_substitute_source(failure):
    deps, _ = services(*MAPPINGS[0], failure=failure)
    result = Converter(services=deps).convert(MAPPINGS[0][0], out_type='json')
    assert result.status == 'failed'
    assert result.items[0].diagnostics[0].code == failure
    deps['ena_source'].resolve_arrayexpress.assert_not_called()


def test_optional_biosamples_failure_remains_separately_retryable():
    deps, _ = services(*MAPPINGS[0], optional_failure=True)
    result = Converter(services=deps).convert(MAPPINGS[0][0], out_type='json', preparation=PreparationPolicy(enrichment='off', expand_studies=False))
    assert result.items[0].execution == 'succeeded'
    issues = [d.provider_issue for d in result.items[0].diagnostics if d.provider_issue]
    assert {(i.error.provider, i.error.retry_category.value) for i in issues} == {('biosamples', 'retryable'), ('biostudies', 'terminal')}


def test_cli_reports_action_and_no_output(monkeypatch, capsys, tmp_path):
    from meta_standards_converter.cli import convert
    deps, _ = services(*MAPPINGS[0], failure='ambiguous_idf')
    monkeypatch.setattr(convert, 'Converter', lambda: Converter(services=deps))
    assert convert.main([MAPPINGS[0][0], '--out-type', 'json', '--out', str(tmp_path), '--progress', 'off', '--report-json', '-']) == 1
    captured = capsys.readouterr()
    assert json.loads(captured.out)['status'] == 'failed'
    assert 'no output produced' in captured.err and 'ambiguous_idf' in captured.err
    assert 'load' in captured.err and 'Select an IDF' in captured.err


def test_cli_partial_output_and_direct_json_mediated_exports_agree(monkeypatch, capsys, tmp_path):
    from meta_standards_converter.cli import convert
    from pathlib import Path
    deps, _ = services(*MAPPINGS[0])
    monkeypatch.setattr(convert, 'Converter', lambda: Converter(services=deps))
    flags = ['--enrichment', 'off', '--no-expand-studies', '--progress', 'off']
    assert convert.main([MAPPINGS[0][0], '--out-type', 'magetab', '--out', str(tmp_path/'direct'), *flags]) == 1
    assert 'partial output written' in capsys.readouterr().err
    assert convert.main([MAPPINGS[0][0], '--out-type', 'json', '--out', str(tmp_path/'json'), *flags]) == 1
    path = next((tmp_path/'json').glob('*.json'))
    assert convert.main([str(path), '--out-type', 'magetab', '--out', str(tmp_path/'via'), *flags]) == 1
    assert {p.name: p.read_bytes() for p in (tmp_path/'direct').glob('*.txt')} == {p.name: p.read_bytes() for p in (tmp_path/'via').glob('*.txt')}


@pytest.mark.parametrize('error', [MetadataSourceError('fallback_unverified', MAPPINGS[0][0], provider='ena'), TimeoutError('PRIVATE')])
def test_fallback_failure_keeps_both_causes_and_writes_nothing(error, tmp_path):
    deps, _ = services(*MAPPINGS[0])
    deps['ena_source'].resolve_arrayexpress.side_effect = error
    result = Converter(services=deps).convert(MAPPINGS[0][0], out_type='json', outdir=tmp_path)
    assert result.status == 'failed' and not result.items[0].artifacts
    codes = {d.code for d in result.items[0].diagnostics}
    assert {'source_inaccessible', 'fallback_failed'} <= codes
    assert 'PRIVATE' not in json.dumps(result.to_dict())
    assert result.items[0].error.retry_category.value == ('retryable' if isinstance(error, TimeoutError) else 'terminal')


def test_access_denial_cache_survives_default_option_spelling():
    from meta_standards_converter.metadata.session import PreparationSession
    deps, _ = services(*MAPPINGS[0])
    session = PreparationSession()
    for options in ({'sdrf_sources': None}, {}):
        with pytest.raises(MetadataSourceError):
            session.load(deps['ae_source'], MAPPINGS[0][0], **options)
    assert deps['ae_source'].load.call_count == 1


def test_fallback_uses_shared_http_budget_and_resource_profile(tmp_path):
    from tests.support.ae import fetcher, response, INFO
    from meta_standards_converter.sources.metadata import MAGETabMetadataSource, NativeMetadataSource
    from meta_standards_converter.sources.ena import ENASource
    from meta_standards_converter.sources.archive_support import ArchiveHTTP
    ae, calls = fetcher(tmp_path, lambda *args: response({'errorMessage': 'Study not accessible'}, 403))
    ena = NativeMetadataSource('ena', source=ENASource(http=ArchiveHTTP('ena', requester=ae.requester)))
    result = Converter(services={'ae_source': MAGETabMetadataSource(fetcher=ae), 'ena_source': ena}).convert(
        MAPPINGS[0][0], out_type='json', runtime=RuntimeOptions(metadata_max_requests=1))
    assert result.status == 'failed' and len(calls) == 1
    assert 'retrieval_budget_exhausted' in {d.code for d in result.items[0].diagnostics}


@pytest.mark.parametrize('error', [ValueError('bad parse'), RuntimeError('validation')])
def test_parsing_errors_never_fall_back(error):
    deps, _ = services(*MAPPINGS[0])
    deps['ae_source'].load.side_effect = error
    result = Converter(services=deps).convert(MAPPINGS[0][0], out_type='json')
    assert result.status == 'failed'
    deps['ena_source'].resolve_arrayexpress.assert_not_called()


@pytest.mark.parametrize('provider', ['ena', 'sra'])
def test_native_fallback_and_peer_activity_reaches_consumer_metrics(tmp_path, provider):
    from tests.support.ae import fetcher, response
    from meta_standards_converter.sources.metadata import NativeMetadataSource
    from meta_standards_converter.sources.archive_support import ArchiveHTTP
    source, _ = fetcher(tmp_path, lambda *args: response([]))
    native = NativeMetadataSource(provider)
    native.source.http = ArchiveHTTP(provider, requester=source.requester)
    native.source.http.get('https://www.ebi.ac.uk/ena/portal/api/search', fmt='json')
    assert native.metrics().provider_attempts == 1
