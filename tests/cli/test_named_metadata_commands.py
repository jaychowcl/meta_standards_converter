# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Named commands share configuration, outcome and publication contracts."""
import importlib
from pathlib import Path
from unittest.mock import Mock

import pytest

from meta_standards_converter.converters import Converter, PreparationPolicy, RuntimeOptions
from meta_standards_converter.converters.unified.contracts import ConversionBatchResult, ConversionItemResult
from meta_standards_converter.sources.metadata import GEOMetadataSource

COMMANDS = [('geo2json', 'json'), ('ae2json', 'json'), ('json2ae', 'magetab'), ('geo2ae', 'magetab')]


def complete():
    return ConversionBatchResult([ConversionItemResult('one', 'source', status='complete',
        execution='succeeded', completeness='complete', validation='valid')])


@pytest.mark.parametrize('command,target', COMMANDS)
@pytest.mark.parametrize('enrichment', ['standard', 'curators', 'off'])
@pytest.mark.parametrize('expand', [True, False])
def test_named_commands_share_typed_configuration(monkeypatch, tmp_path, command, target, enrichment, expand):
    cli = importlib.import_module('meta_standards_converter.cli.convert')
    converter = Mock()
    converter.convert.return_value = complete()
    monkeypatch.setattr(cli, 'Converter', lambda: converter)
    main = importlib.import_module('meta_standards_converter.cli.' + command).main
    argv = ['input.json' if command == 'json2ae' else 'GSE1', '--out', str(tmp_path),
        '--enrichment', enrichment, '--expand-studies' if expand else '--no-expand-studies',
        '--resource-profile', 'large', '--resource-override', 'max_xml_bytes=4096',
        '--metadata-max-requests', '4', '--metadata-timeout-seconds', '8',
        '--no-peer-provider', '--no-linked-metadata']
    assert main(argv) == 0
    args, kw = converter.convert.call_args
    assert args == ([argv[0]],)
    assert kw['out_type'] == target
    assert kw['preparation'] == PreparationPolicy(enrichment, expand, False, False)
    assert kw['runtime'].resource_profile == 'large'
    assert dict(kw['runtime'].resource_overrides) == {'max_xml_bytes': 4096}
    assert kw['runtime'].metadata_max_requests == 4
    assert kw['runtime'].metadata_timeout_seconds == 8
    if command == 'json2ae':
        assert kw['in_type'] == 'json'


@pytest.mark.parametrize('command,target', COMMANDS)
def test_defaults_batch_order_failure_quiet_and_safe_logging(monkeypatch, tmp_path, capsys, command, target):
    cli = importlib.import_module('meta_standards_converter.cli.convert')
    converter = Mock()
    converter.convert.return_value = complete()
    monkeypatch.setattr(cli, 'Converter', lambda: converter)
    main = importlib.import_module('meta_standards_converter.cli.' + command).main
    assert main(['one', 'two', '--out', str(tmp_path), '--quiet']) == 0
    assert converter.convert.call_args.args == (['one', 'two'],)
    assert converter.convert.call_args.kwargs['preparation'] == PreparationPolicy()
    assert converter.convert.call_args.kwargs['runtime'] == RuntimeOptions()
    assert not capsys.readouterr().out
    converter.convert.side_effect = RuntimeError('secret-token')
    assert main(['one', '--out', str(tmp_path)]) == 1
    captured = capsys.readouterr()
    assert 'secret-token' not in captured.out + captured.err


@pytest.mark.parametrize('command,target', COMMANDS)
@pytest.mark.parametrize('flag', ['--enrich', '--no-enrich', '--related', '--related-series', '--force-in-type'])
def test_obsolete_options_are_rejected(command, target, flag):
    main = importlib.import_module('meta_standards_converter.cli.' + command).main
    with pytest.raises(SystemExit) as error:
        main(['one', flag])
    assert error.value.code == 2


@pytest.mark.parametrize('command,target', COMMANDS)
def test_listing_does_not_construct_converter(monkeypatch, capsys, command, target):
    cli = importlib.import_module('meta_standards_converter.cli.convert')
    monkeypatch.setattr(cli, 'Converter', Mock(side_effect=AssertionError('must not construct')))
    main = importlib.import_module('meta_standards_converter.cli.' + command).main
    assert main(['--list-platform-handlers']) == 0
    assert 'generic' in capsys.readouterr().out


@pytest.mark.parametrize('command,target', [('geo2json','json'), ('geo2ae','magetab')])
@pytest.mark.parametrize('accession', [True, False])
def test_named_cli_and_api_publish_identical_artifacts(monkeypatch, tmp_path, command, target, accession):
    cli = importlib.import_module('meta_standards_converter.cli.convert')
    source = Path(__file__).parents[1] / 'fixtures/studies/GSE100/inputs/geo.xml'
    fetcher = Mock()
    fetcher.fetch_gse_miniml.return_value = source.read_text()
    services = {'geo_source': GEOMetadataSource(fetcher=fetcher)}
    monkeypatch.setattr(cli, 'Converter', lambda: Converter(services=services))
    token = 'GSE100' if accession else str(source)
    policy = PreparationPolicy(enrichment='off', expand_studies=False)
    batch = Converter(services=services).convert(token, out_type=target,
        preparation=policy, outdir=tmp_path/'api')
    assert batch.status == 'complete', batch.to_dict()
    main = importlib.import_module('meta_standards_converter.cli.'+command).main
    assert main([token, '--enrichment', 'off', '--no-expand-studies', '--out', str(tmp_path/'cli')]) == 0
    assert {p.name:p.read_bytes() for p in (tmp_path/'api').iterdir() if p.is_file()} == {
        p.name:p.read_bytes() for p in (tmp_path/'cli').iterdir() if p.is_file()}


@pytest.mark.parametrize('command,accession', [('ae2json', True), ('ae2json', False), ('json2ae', False)])
@pytest.mark.parametrize('enrichment', ['standard', 'curators', 'off'])
@pytest.mark.parametrize('expand', [True, False])
def test_ae_and_saved_json_cli_api_artifact_parity(monkeypatch, tmp_path, command, accession, enrichment, expand):
    from meta_standards_converter.sources.metadata import MAGETabMetadataSource
    cli = importlib.import_module('meta_standards_converter.cli.convert')
    fixtures = Path(__file__).parents[1] / 'fixtures/studies'
    ae_file = fixtures / 'E-MTAB-6486/inputs/E-MTAB-6486.idf.txt'
    ae_source = Mock()
    ae_source.load.return_value = MAGETabMetadataSource().load(str(ae_file))
    enricher = Mock(publication_issues=[], _sra_failures=0)
    enricher.enrich_selected.side_effect = lambda package, **options: package
    class UnavailableRelative:
        def load(self, accession, **options):
            raise OSError('relative absent from the offline fixture')
    services = {'ae_source': ae_source, 'enricher': enricher, 'geo_source': UnavailableRelative()}
    monkeypatch.setattr(cli, 'Converter', lambda: Converter(services=services))
    token = ('E-MTAB-6486' if accession else str(ae_file)) if command == 'ae2json' else str(fixtures / 'GSE328265/inputs/miniml.json')
    target = 'json' if command == 'ae2json' else 'magetab'
    policy = PreparationPolicy(enrichment=enrichment, expand_studies=expand)
    runtime = RuntimeOptions(resource_profile='large', resource_overrides={'max_xml_bytes': 67108864})
    result = Converter(services=services).convert(token, out_type=target,
        preparation=policy, runtime=runtime, outdir=tmp_path/'api')
    # The saved GEO fixture declares relatives not included in this fixture set.
    # Expansion therefore preserves the same usable partial result in both APIs.
    partial = command == 'json2ae' and expand
    assert result.status == ('partial' if partial else 'complete'), result.to_dict()
    main = importlib.import_module('meta_standards_converter.cli.'+command).main
    assert main([token, '--out', str(tmp_path/'cli'), '--enrichment', enrichment,
        '--expand-studies' if expand else '--no-expand-studies',
        '--resource-profile', 'large', '--resource-override', 'max_xml_bytes=67108864']) == (1 if partial else 0)
    assert {p.name:p.read_bytes() for p in (tmp_path/'api').iterdir() if p.is_file()} == {
        p.name:p.read_bytes() for p in (tmp_path/'cli').iterdir() if p.is_file()}
