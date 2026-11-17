# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Every command protects inputs and distinguishes usage from execution failures."""
import pytest


def test_unexpected_value_error_is_execution_failure(tmp_path, monkeypatch, capsys):
    from meta_standards_converter.cli import convert
    class Broken:
        def convert(self, *args, **kwargs):
            raise ValueError('private_payload_canary')
    monkeypatch.setattr(convert, 'Converter', Broken)
    assert convert.main(['GSE1', '--out-type', 'json', '--out', str(tmp_path)]) == 1
    stderr = capsys.readouterr().err
    assert 'conversion_failed' in stderr
    assert 'private_payload_canary' not in stderr


def test_migration_missing_input_is_actionable(tmp_path, capsys):
    from meta_standards_converter.cli.miniml_migrate import main
    destination = tmp_path/'output.json'
    assert main([str(tmp_path/'missing.json'), str(destination)]) == 1
    assert not destination.exists()
    assert 'missing_file' in capsys.readouterr().err


@pytest.mark.parametrize('command', ['json2tsv', 'json2h5ad', 'json2obs', 'miniml_migrate'])
def test_all_file_commands_protect_input_from_log_truncation(command, tmp_path):
    import importlib
    source = tmp_path/'source.json'
    source.write_text('original data')
    args = [str(source)]
    if command == 'miniml_migrate':
        args += [str(tmp_path/'destination.json')]
    else:
        args += ['--outdir', str(tmp_path/'out')]
    args += ['--log-file', str(source)]
    with pytest.raises(SystemExit) as caught:
        importlib.import_module('meta_standards_converter.cli.'+command).main(args)
    assert caught.value.code == 2
    assert source.read_text() == 'original data'
