# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Preflight boundaries must not trigger providers, processing or large loads."""

import json
from pathlib import Path

import pytest

from meta_standards_converter.converters import Converter, InputSpec
from tests.converters.test_json2tsv import package


def test_explicit_missing_path_named_like_accession_is_not_retrieved(
    tmp_path, monkeypatch
):
    monkeypatch.chdir(tmp_path)

    class Never:
        def convert(self, *args, **kwargs):
            pytest.fail("Explicit paths must not reach accession retrieval")

    result = Converter(services={"geo2json": Never()}).convert(
        Path("GSE1"), out_type="json"
    )
    assert result.items[0].diagnostics[0].code == "missing_file"


def test_forcing_accession_cannot_reinterpret_missing_path(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    class Never:
        def convert(self, *args, **kwargs):
            pytest.fail("Missing Path reached provider")

    result = Converter(services={"geo2json": Never()}).convert(
        Path("GSE1"), out_type="json", force_in_type="geo"
    )
    assert result.status == "failed"
    assert result.items[0].diagnostics[0].code == "missing_file"


def test_empty_directory_is_an_input_outcome_even_with_destination(tmp_path):
    source = tmp_path / "empty"
    source.mkdir()
    result = Converter().convert(source, out_type="json", outdir=tmp_path / "out")
    assert result.status == "failed"
    assert result.items[0].source == str(source)


@pytest.mark.parametrize("forced", [None, "matrix"])
def test_root_directory_symlink_is_not_traversed(tmp_path, forced):
    real = tmp_path / "real"
    real.mkdir()
    (real / "input.json").write_text(json.dumps(package()))
    link = tmp_path / "link"
    link.symlink_to(real, target_is_directory=True)
    result = Converter().convert(
        link, out_type="json", force_in_type=forced, runtime_options={"recursive": True}
    )
    assert result.status == "failed"
    assert not result.items[0].dataset_ids
    assert any(d.code == "directory_symlink" for d in result.items[0].diagnostics)


def test_duplicate_source_with_conflicting_settings_is_reported(tmp_path):
    path = tmp_path / "input.json"
    path.write_text(json.dumps(package()))
    result = Converter().convert(
        [InputSpec(path), InputSpec(path, output_options={"aggregate": True})],
        out_type="json",
    )
    assert result.items[1].status == "failed"
    assert any(d.code == "conflicting_input" for d in result.items[1].diagnostics)


def test_irrelevant_companions_and_metadata_are_request_errors():
    with pytest.raises(ValueError, match="companions"):
        Converter().convert(
            InputSpec(package(), companions={"barcodes": "x"}), out_type="json"
        )
    with pytest.raises(ValueError, match="metadata"):
        Converter().convert(InputSpec(package(), metadata=package()), out_type="json")


def test_manifest_binding_removes_unreferenced_companion_from_discovery(tmp_path):
    from shutil import copyfile

    fixture = Path(__file__).parents[2] / "tests/fixtures/studies/E-MTAB-1/inputs"
    source_idf = next(fixture.glob("*.idf.txt"))
    source_sdrf = next(fixture.glob("*.sdrf.txt"))
    idf = tmp_path / "study.idf.txt"
    idf.write_text(
        "\n".join(
            line
            for line in source_idf.read_text().splitlines()
            if not line.startswith("SDRF File")
        )
    )
    copyfile(source_sdrf, tmp_path / "bound.sdrf.txt")
    manifest = tmp_path / "input-manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "schema_version": "1.0",
                "inputs": [
                    {
                        "id": "bound",
                        "sources": ["study.idf.txt"],
                        "companions": {"sdrf": "bound.sdrf.txt"},
                    }
                ],
            }
        )
    )
    result = Converter().convert(tmp_path, out_type="json", input_manifest=manifest)
    assert len(result.items) == 1, result.to_dict()
    assert result.status == "complete", result.to_dict()


def test_standalone_matrix_obeys_memory_admission(tmp_path):
    path = tmp_path / "counts.tsv"
    path.write_text("gene\tc1\ng1\t1\n")
    result = Converter().convert(
        path,
        out_type="h5ad",
        input_options={"orientation": "genes-by-observations"},
        runtime_options={"resource_overrides": {"max_in_memory_matrix_bytes": 1}},
    )
    assert result.status == "failed"
    assert any(d.code == "memory_limit" for d in result.items[0].diagnostics)


def test_raw_missing_reference_does_not_call_injected_runner(tmp_path):
    from meta_standards_converter.converters.json2h5ad import JSON2H5ADConverter

    reads = tmp_path / "reads.fastq"
    reads.write_text("@x\nAC\n+\n!!\n")
    calls = []

    class Runner:
        def process(self, *args, **kwargs):
            calls.append(True)
            raise AssertionError("reference must be checked first")

    result = Converter(
        services={"json2h5ad": JSON2H5ADConverter(pipeline_runner=Runner())}
    ).convert(
        InputSpec(reads, metadata=package()),
        out_type="h5ad",
        outdir=tmp_path / "out",
        runtime_options={"allow_processing": True},
    )
    assert not calls
    assert any(d.code == "reference_required" for d in result.items[0].diagnostics)
    assert not (tmp_path / "out").exists()


def test_raw_metadata_export_does_not_require_processing(tmp_path):
    reads = tmp_path / "reads.fastq"
    reads.write_text("@x\nAC\n+\n!!\n")
    result = Converter().convert(InputSpec(reads, metadata=package()), out_type="json")
    assert result.status == "complete", result.to_dict()


@pytest.mark.parametrize(
    "kind",
    [
        "geo_xml",
        "geo_archive",
        "sra_xml",
        "ena_xml",
        "magetab",
        "h5ad",
        "matrix",
        "fastq",
    ],
)
def test_forced_handlers_reject_wrong_content_without_outputs(kind, tmp_path):
    path = tmp_path / "wrong.json"
    path.write_text(json.dumps(package()))
    result = Converter().convert(
        InputSpec(path, in_type=kind), out_type="h5ad", outdir=tmp_path / "out"
    )
    assert result.status == "failed", result.to_dict()
    assert not result.items[0].artifacts


def test_probe_is_used_for_replaced_registered_handler():
    from meta_standards_converter.converters.unified.handlers import InputHandler

    class Custom(InputHandler):
        def __init__(self):
            super().__init__("miniml")

        def probe(self, source, context):
            return source == 42

        def load(self, spec, context):
            return self._metadata(package())

    result = Converter(handlers=[Custom()]).convert(42, out_type="json")
    assert result.status == "complete", result.to_dict()


def test_standalone_sdrf_finds_unique_idf():
    fixture = Path(__file__).parents[2] / "tests/fixtures/studies/E-MTAB-1/inputs"
    result = Converter().convert(next(fixture.glob("*.sdrf.txt")), out_type="json")
    assert result.status == "complete", result.to_dict()


def test_request_reports_ambiguous_idf_sdrf_association(tmp_path):
    (tmp_path / "one.idf.txt").write_text("SDRF File\tcommon.sdrf.txt\n")
    (tmp_path / "two.idf.txt").write_text("SDRF File\tcommon.sdrf.txt\n")
    (tmp_path / "common.sdrf.txt").write_text("Source Name\tSample Name\nx\tx\n")
    result = Converter().convert(tmp_path, out_type="json")
    assert all(
        any(d.code == "ambiguous_companion" for d in item.diagnostics)
        for item in result.items
    )


def test_routes_are_registered_explicitly():
    from meta_standards_converter.converters.unified.routes import ROUTES

    assert ("metadata", "json") in ROUTES
    assert ("expression", "obs") in ROUTES
    assert ("expression", "json") not in ROUTES


def test_content_identity_is_deterministic_and_sources_are_retained():
    first = Converter().convert(package(), out_type="json").items[0]
    second = (
        Converter()
        .convert(dict(reversed(list(package().items()))), out_type="json")
        .items[0]
    )
    assert first.content_id == second.content_id and len(first.content_id) == 64


def test_manifest_xml_sources_take_precedence_over_directory_grouping(tmp_path):
    from shutil import copyfile

    fixture = Path(__file__).parents[2] / "docs/ena/fixtures/SRX7812918"
    names = [kind + ".xml" for kind in ("study", "sample", "experiment", "run")]
    for name in names:
        copyfile(fixture / name, tmp_path / name)
    manifest = tmp_path / "manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "schema_version": "1.0",
                "inputs": [{"id": "native", "sources": names, "in_type": "ena_xml"}],
            }
        )
    )
    result = Converter().convert(tmp_path, input_manifest=manifest, out_type="json")
    assert len(result.items) == 1, result.to_dict()


def test_metadata_routes_do_not_import_scientific_dependencies(monkeypatch):
    import builtins

    original = builtins.__import__

    def guarded(name, *args, **kwargs):
        if name.split(".")[0] in {
            "anndata",
            "scanpy",
            "pandas",
            "numpy",
            "scipy",
            "h5py",
        }:
            pytest.fail("Metadata-only route imported a scientific dependency")
        return original(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", guarded)
    result = Converter().convert(package(), out_type="tsv")
    assert result.status == "complete"
