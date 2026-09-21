# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Unified terminal contracts: routing, reports and safe publication."""
import json
from pathlib import Path
import pytest
from meta_standards_converter.cli import convert
from meta_standards_converter.converters.unified.contracts import ConversionBatchResult, ConversionItemResult
from tests.converters.test_json2tsv import package


@pytest.fixture
def source(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    path = tmp_path / "input.json"
    path.write_text(json.dumps(package()))
    return path


@pytest.mark.parametrize("target,suffix", [("json", ".json"), ("tsv", ".tsv"), ("csv", ".csv"), ("magetab", ".idf.txt")])
def test_real_metadata_outputs(source, target, suffix, capsys):
    assert convert.main([str(source), "--out-type", target, "--enrichment", "off"]) == 0
    assert source.with_name(("E-GEOD-1" if target == "magetab" else "GSE1") + suffix).exists()
    assert "complete" in capsys.readouterr().out


def test_json_stdout_is_machine_readable(source, capsys):
    assert convert.main([str(source), "--out-type", "json", "--report-json", "-", "-v", "--enrichment", "off"]) == 0
    report = json.loads(capsys.readouterr().out)
    assert report["operation"] == "convert"
    assert report["status"] == "complete"
    assert report["items"][0]["artifacts"]


@pytest.mark.parametrize("args", [[], ["x", "--out-type", "bad"], ["x", "--out-type", "json", "--out", "a", "--outfile", "b"], ["x", "--out-type", "json", "--expand-studies", "--no-related"]])
def test_usage_errors(args):
    with pytest.raises(SystemExit) as error:
        convert.main(args)
    assert error.value.code == 2


@pytest.mark.parametrize("flag,expected", [("--list-platform-handlers", "tenx_v3_droplet_single_cell_sequencing"), ("--list-input-types", "geo"), ("--list-output-types", "magetab")])
def test_list_without_conversion(flag, expected, capsys):
    assert convert.main([flag]) == 0
    assert expected in capsys.readouterr().out


def test_omitted_defaults_are_not_forwarded(monkeypatch, tmp_path):
    calls = []
    class Fake:
        def convert(self, *args, **kwargs):
            calls.append((args, kwargs))
            return ConversionBatchResult([ConversionItemResult("GSE1", "GSE1", status="complete")])
    monkeypatch.setattr(convert, "Converter", Fake)
    assert convert.main(["GSE1", "--out-type", "json"]) == 0
    assert calls == [((["GSE1"],), {"out_type": "json", "outdir": ".", "options": {}})]


@pytest.mark.parametrize("flag,value,key,expected", [
    ("--platform-handler", "generic", "platform_handler", "generic"),
    ("--related", None, "expand_studies", True),
    ("--no-expand-studies", None, "expand_studies", False),
    ("--asset", "GSM1=counts.tsv", "asset_specs", ["GSM1=counts.tsv"]),
    ("--sdrf-source", "x.sdrf.txt", "sdrf_sources", ["x.sdrf.txt"]),
    ("--allowed-host", "example.org", "allowed_hosts", ["example.org"]),
    ("--execution-profile", "docker", "profile", "docker"),
    ("--resource-override", "max_xml_bytes=1024", "resource_overrides", {"max_xml_bytes": 1024}),
    ("--no-remove-empty", None, "remove_empty", False),
    ("--allow-processing", None, "allow_processing", True),
])
def test_dedicated_flags_forward_typed_values(monkeypatch, flag, value, key, expected):
    calls = []
    class Fake:
        def convert(self, *args, **kwargs):
            calls.append(kwargs)
            return ConversionBatchResult()
    monkeypatch.setattr(convert, "Converter", Fake)
    argv = ["GSE1", "--out-type", "magetab", flag]
    if value is not None:
        argv.append(value)
    assert convert.main(argv) == 1
    assert calls[0]["options"][key] == expected


def test_exact_file_and_overwrite(source):
    argv = [str(source), "--out-type", "json", "--outfile", "exact.json", "--enrichment", "off"]
    assert convert.main(argv) == 0
    before = Path("exact.json").read_bytes()
    assert convert.main(argv) == 1
    assert Path("exact.json").read_bytes() == before
    assert convert.main(argv + ["--overwrite"]) == 0


def test_partial_fail_fast_and_report(source, capsys):
    assert convert.main(["missing.json", str(source), "--out-type", "json", "--fail-fast", "--report-json", "report.json", "--enrichment", "off"]) == 1
    report = json.loads(Path("report.json").read_text())
    assert [i["status"] for i in report["items"]] == ["failed", "skipped"]
    assert not Path("GSE1.json").exists()
    assert "failed" in capsys.readouterr().out


@pytest.mark.parametrize("report", ["input.json", "GSE1.json"])
def test_report_cannot_overwrite_input_or_conversion_output(source, report):
    before = source.read_bytes()
    argv = [str(source), "--out-type", "json", "--report-json", report, "--overwrite", "--enrichment", "off"]
    if report == "input.json":
        with pytest.raises(SystemExit) as error:
            convert.main(argv)
        assert error.value.code == 2
    else:
        assert convert.main(argv) == 1
        assert any(d["code"] == "output_collision" for d in json.loads(Path(report).read_text())["items"][0]["diagnostics"])
    assert source.read_bytes() == before


def test_existing_report_requires_overwrite_before_conversion(source):
    Path("report.json").write_text("keep")
    with pytest.raises(SystemExit) as error:
        convert.main([str(source), "--out-type", "json", "--report-json", "report.json"])
    assert error.value.code == 2
    assert Path("report.json").read_text() == "keep"
    assert not Path("GSE1.json").exists()


def test_manifest_directory_and_current_output(source):
    Path("manifest.json").write_text(json.dumps({"schema_version": "1.0", "inputs": [
        {"id": "one", "sources": [source.name], "input_options": {"enrichment": "off"}}
    ]}))
    assert convert.main(["--input-manifest", "manifest.json", "--out-type", "csv"]) == 0
    assert Path("GSE1.csv").exists()


def test_irrelevant_flags_fail_before_publication(source):
    with pytest.raises(SystemExit) as error:
        convert.main([str(source), "--out-type", "json", "--platform-handler", "generic"])
    assert error.value.code == 2
    assert not Path("GSE1.json").exists()


@pytest.mark.parametrize("target", ["h5ad", "obs"])
def test_expression_values_survive_cli(source, target, capsys):
    import anndata
    import numpy as np
    import pandas as pd
    matrix = anndata.AnnData(np.array([[1., 2.], [3., 4.]]),
        obs=pd.DataFrame({"sample": ["a", "a"]}, index=["c1", "c2"]),
        var=pd.DataFrame(index=["g1", "g2"]))
    matrix.write_h5ad("source.h5ad")
    args = ["source.h5ad", "--out-type", target, "--out", "out", "--report-json", "-"]
    if target == "obs":
        args += ["--include-var", "--include-uns"]
    assert convert.main(args) == 0
    report = json.loads(capsys.readouterr().out)
    assert report["status"] == "complete"
    artifacts = report["items"][0]["artifacts"]
    if target == "h5ad":
        output = next(p for p in artifacts.values() if p.endswith(".h5ad"))
        np.testing.assert_array_equal(anndata.read_h5ad(output).X, matrix.X)
    else:
        output = next(p for k, p in artifacts.items() if "obs" in k and p.endswith(".csv"))
        obs = pd.read_csv(output, index_col=0)
        assert obs.index.tolist() == ["c1", "c2"]
        assert obs["sample"].tolist() == ["a", "a"]


def test_raw_processing_permission_from_cli(source, capsys):
    raw = Path("reads.fastq")
    raw.write_text("@r\nAC\n+\n!!\n")
    Path("manifest.json").write_text(json.dumps({"schema_version": "1.0", "inputs": [
        {"id": "reads", "sources": [str(raw)], "metadata": str(source)}]}))
    assert convert.main(["--input-manifest", "manifest.json", "--out-type", "h5ad", "--report-json", "-", "--enrichment", "off"]) == 1
    report = json.loads(capsys.readouterr().out)
    assert any(d["code"] == "processing_required" for d in report["items"][0]["diagnostics"])


def test_report_cannot_replace_explicit_asset(source):
    Path("counts.tsv").write_text("gene\tc1\ng1\t1\n")
    with pytest.raises(SystemExit) as error:
        convert.main([str(source), "--out-type", "h5ad", "--asset", "GSM1=counts.tsv", "--report-json", "counts.tsv", "--overwrite"])
    assert error.value.code == 2
    assert Path("counts.tsv").read_text().startswith("gene")


def test_report_cannot_replace_manifest_companion(source):
    Path("counts.tsv").write_text("gene\tc1\ng1\t1\n")
    Path("manifest.json").write_text(json.dumps({"schema_version": "1.0", "inputs": [
        {"id": "matrix", "sources": ["counts.tsv"], "metadata": str(source)}]}))
    with pytest.raises(SystemExit) as error:
        convert.main(["--input-manifest", "manifest.json", "--out-type", "json", "--report-json", str(source), "--overwrite"])
    assert error.value.code == 2


def test_report_excluded_from_recursive_discovery(source, capsys):
    Path("report.json").write_text("{}")
    with pytest.raises(SystemExit) as error:
        convert.main([".", "--out-type", "csv", "--report-json", "report.json", "--overwrite", "--enrichment", "off"])
    assert error.value.code == 2  # An existing file in the explicitly selected input root is protected.


def test_success_survives_later_failed_input(source, capsys):
    assert convert.main([str(source), "missing.json", "--out-type", "json", "--report-json", "-", "--enrichment", "off"]) == 1
    report = json.loads(capsys.readouterr().out)
    assert report["status"] == "partial"
    assert [i["status"] for i in report["items"]] == ["complete", "failed"]
    assert Path("GSE1.json").exists()


@pytest.mark.parametrize('report,outdir', [('output', 'output'), ('output', 'output/nested'), ('.artifact-bundles/report.json', '.')])
def test_auxiliary_files_cannot_be_output_directories_or_bundle_internals(source, report, outdir):
    with pytest.raises(SystemExit) as error:
        convert.main([str(source), '--out-type', 'json', '--out', outdir, '--report-json', report, '--enrichment', 'off'])
    assert error.value.code == 2
    assert not Path('output').exists()


def test_log_cannot_overwrite_input(source):
    original = source.read_bytes()
    with pytest.raises(SystemExit) as error:
        convert.main([str(source), '--out-type', 'json', '--log-file', str(source), '--overwrite'])
    assert error.value.code == 2
    assert source.read_bytes() == original


def test_raw_processing_still_requires_reference(source, capsys):
    Path('reads.fastq').write_text('@r\nAC\n+\n!!\n')
    Path('manifest.json').write_text(json.dumps({'schema_version': '1.0', 'inputs': [
        {'id': 'reads', 'sources': ['reads.fastq'], 'metadata': str(source)}]}))
    assert convert.main(['--input-manifest', 'manifest.json', '--out-type', 'h5ad', '--allow-processing', '--report-json', '-', '--enrichment', 'off']) == 1
    report = json.loads(capsys.readouterr().out)
    assert any(d['code'] == 'reference_required' for d in report['items'][0]['diagnostics'])


@pytest.mark.parametrize('value', ['{}', '[{"scope_id": "GSM1"}]', 'not-json'])
def test_explicit_asset_structure_errors_are_usage_errors(source, value):
    with pytest.raises(SystemExit) as error:
        convert.main([str(source), '--out-type', 'h5ad', '--explicit-assets', value])
    assert error.value.code == 2


@pytest.mark.fake_process
def test_help_and_metadata_cli_do_not_import_scientific_libraries():
    import os
    import subprocess
    import sys
    root = Path(__file__).resolve().parents[2]
    script = '''
import importlib.abc
import sys
class NoScientific(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split('.')[0] in {'anndata', 'scanpy', 'numpy', 'pandas', 'scipy', 'h5py'}:
            raise AssertionError('Scientific import: ' + fullname)
sys.meta_path.insert(0, NoScientific())
from meta_standards_converter.cli.convert import main
assert main(['--list-output-types']) == 0
try:
    main(['--help'])
except SystemExit as exc:
    assert exc.code == 0
'''
    completed = subprocess.run([sys.executable, '-c', script], cwd=root,
        env=os.environ | {'PYTHONPATH': str(root / 'src')}, capture_output=True, text=True)
    assert completed.returncode == 0, completed.stderr
    assert 'msc-convert' in completed.stdout


@pytest.mark.parametrize("kind,target", [("h5ad", "obs"), ("matrix", "h5ad")])
@pytest.mark.parametrize("missing_name", ["anndata", "private_dependency"])
def test_missing_scientific_dependencies_give_install_guidance(tmp_path, monkeypatch, capsys, kind, target, missing_name):
    import builtins

    original = builtins.__import__

    def missing_scientific(name, *args, **kwargs):
        if name.split(".")[0] == "anndata":
            raise ModuleNotFoundError("private environment details", name=missing_name)
        return original(name, *args, **kwargs)

    source = tmp_path / ("input.h5ad" if kind == "h5ad" else "matrix.tsv")
    source.write_bytes(b"\x89HDF\r\n\x1a\n" if kind == "h5ad" else b"gene\tc1\ng1\t1\n")
    monkeypatch.setattr(builtins, "__import__", missing_scientific)
    argv = [str(source), "--in-type", kind, "--out-type", target,
            "--out", str(tmp_path / "out"), "--report-json", "-"]
    if kind == "matrix":
        argv += ["--orientation", "genes-by-observations"]
    assert convert.main(argv) == 1
    report = json.loads(capsys.readouterr().out)
    assert report["status"] == "failed"
    item = report["items"][0]
    assert not item["artifacts"]
    diagnostic = item["diagnostics"][0]
    if missing_name == "anndata":
        assert diagnostic["code"] == "missing_optional_dependency"
        assert "meta-standards-converter[h5ad]" in diagnostic["message"]
        assert "full" in diagnostic["message"]
    else:
        assert diagnostic["code"] == "conversion_failed"
        assert "meta-standards-converter[h5ad]" not in diagnostic["message"]
    assert "private environment details" not in json.dumps(report)
    assert "private_dependency" not in json.dumps(report)
