# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Observable unified API contracts, including defensive negative routes."""

import copy
import json
from pathlib import Path

import pytest

from meta_standards_converter.converters import (
    Converter,
    InputSpec,
    ConversionBatchResult,
)
from meta_standards_converter.miniml import MINiMLCodec
from tests.converters.test_json2tsv import package, atlas_v1_payload, atlas_v1_dataset


@pytest.mark.parametrize(
    "form", ["mapping", "typed", "list", "file", "curator", "atlas"]
)
def test_metadata_forms_have_identical_json_and_table(form, tmp_path):
    value = package()
    if form == "typed":
        value = MINiMLCodec().decode(value).package
    if form == "list":
        value = [value]
    if form == "curator":
        value = {"miniml_json": value}
    if form == "atlas":
        value = atlas_v1_payload([atlas_v1_dataset("GSE1", value)])
    if form == "file":
        path = tmp_path / "input.json"
        path.write_text(json.dumps(value))
        value = path
    result = Converter().convert(value, out_type="json")
    assert isinstance(result, ConversionBatchResult)
    assert result.status == "complete", result.to_dict()
    assert result.items[0].dataset_ids == ("GSE1",)
    assert result.items[0].payload[0].samples[0].iid == "GSM1"
    table = Converter().convert(value, out_type="tsv")
    assert table.status == "complete", table.to_dict()
    assert table.items[0].payload["GSE1"].rows[0]["msc.sample.accession"] == "GSM1"


def test_does_not_mutate_input_and_writes_exact_json(tmp_path):
    value = package()
    original = copy.deepcopy(value)
    output = tmp_path / "exact.json"
    result = Converter().convert(value, out_type="json", outfile=output)
    assert result.status == "complete", result.to_dict()
    assert json.loads(output.read_text())[0]["sample"][0]["iid"] == "GSM1"
    assert value == original
    before = output.read_bytes()
    rejected = Converter().convert(value, out_type="json", outfile=output)
    assert rejected.status == "failed"
    assert output.read_bytes() == before
    assert (
        Converter()
        .convert(
            value, out_type="json", outfile=output, runtime_options={"overwrite": True}
        )
        .status
        == "complete"
    )


@pytest.mark.parametrize(
    "kwargs",
    [
        {"out_type": "unknown"},
        {"out_type": "json", "force_in_type": "unknown"},
        {"out_type": "json", "outfile": "x", "outdir": "y"},
        {"out_type": "json", "runtime_options": {"overwite": True}},
        {"out_type": "json", "runtime_options": {"overwrite": "false"}},
        {"out_type": "json", "output_options": {"matrix_orientation": "auto"}},
        {"out_type": "json", "input_options": {"related_series": True}},
    ],
)
def test_bad_request_settings_raise(kwargs):
    with pytest.raises((TypeError, ValueError)):
        Converter().convert(package(), **kwargs)


@pytest.mark.parametrize(
    "value,forced",
    [
        ({"arbitrary": "mapping"}, None),
        (package(), "atlas"),
        ({**package(), "miniml_schema_version": "2.0"}, "miniml"),
        ('{"anything": 1}', "json"),
        ('{"anything": 1}', None),
    ],
)
def test_invalid_or_forced_mismatch_is_reported(value, forced):
    result = Converter().convert(value, out_type="json", force_in_type=forced)
    assert result.status == "failed"
    assert result.items[0].diagnostics
    json.dumps(result.to_dict())


def test_inline_json_requires_forcing():
    text = json.dumps(package())
    assert Converter().convert(text, out_type="json").status == "failed"
    assert (
        Converter().convert(text, force_in_type="json", out_type="json").status
        == "complete"
    )


def test_batch_continues_and_fail_fast_records_remaining(tmp_path):
    inputs = [tmp_path / "missing.json", InputSpec(package(), id="good")]
    result = Converter().convert(inputs, out_type="json")
    assert result.status == "partial"
    assert [item.status for item in result.items] == ["failed", "complete"]
    stopped = Converter().convert(
        inputs, out_type="json", runtime_options={"fail_fast": True}
    )
    assert [item.status for item in stopped.items] == ["failed", "skipped"]


def test_directory_discovery_manifest_and_duplicates(tmp_path):
    a = tmp_path / "a.json"
    a.write_text(json.dumps(package()))
    b = tmp_path / "b.json"
    b.write_text(json.dumps(package("GSE2", "GSM2")))
    (tmp_path / "notes.txt").write_text("not a supported input")
    nested = tmp_path / "nested"
    nested.mkdir()
    (nested / "c.json").write_text(json.dumps(package("GSE3", "GSM3")))
    (tmp_path / "loop").symlink_to(tmp_path, target_is_directory=True)
    converter = Converter()
    result = converter.convert(tmp_path, out_type="json")
    assert [i.dataset_ids for i in result.items if i.status == "complete"] == [
        ("GSE1",),
        ("GSE2",),
    ]
    assert any(i.status == "failed" for i in result.items)
    recursive = converter.convert(
        tmp_path, out_type="json", runtime_options={"recursive": True}
    )
    assert ("GSE3",) in [i.dataset_ids for i in recursive.items]
    duplicate = converter.convert([a, a], out_type="json")
    assert [i.status for i in duplicate.items] == ["complete", "skipped"]
    manifest = tmp_path / "manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "schema_version": "1.0",
                "inputs": [
                    {"id": "chosen", "sources": ["b.json"], "in_type": "miniml"}
                ],
            }
        )
    )
    selected = converter.convert(None, out_type="json", input_manifest=manifest)
    assert selected.items[0].id == "chosen"
    assert selected.items[0].dataset_ids == ("GSE2",)


def test_output_cardinality_and_safe_names(tmp_path):
    converter = Converter()
    with pytest.raises(ValueError):
        converter.convert(
            [InputSpec(package(), id="a"), InputSpec(package("GSE2", "GSM2"), id="b")],
            out_type="json",
            outfile=tmp_path / "one.json",
        )
    result = converter.convert(
        [package(), package("GSE2", "GSM2")],
        out_type="json",
        outfile=tmp_path / "one.json",
    )
    assert result.status == "failed"
    aggregated = converter.convert(
        [package(), package("GSE2", "GSM2")],
        out_type="json",
        outfile=tmp_path / "one.json",
        output_options={"aggregate": True},
    )
    assert aggregated.status == "complete"
    assert len(json.loads((tmp_path / "one.json").read_text())) == 2
    unsafe = converter.convert(
        package("../escape", "GSM1"), out_type="json", outdir=tmp_path / "out"
    )
    assert unsafe.status == "failed"
    assert not (tmp_path / "escape.json").exists()


@pytest.mark.parametrize("target,delimiter", [("csv", ","), ("tsv", "\t")])
def test_table_file_contains_sample_rows(tmp_path, target, delimiter):
    import csv

    output = tmp_path / ("samples." + target)
    result = Converter().convert(package(), out_type=target, outfile=output)
    assert result.status == "complete", result.to_dict()
    with output.open() as stream:
        rows = list(csv.DictReader(stream, delimiter=delimiter))
    assert len(rows) == 1 and rows[0]["msc.sample.accession"] == "GSM1"


def test_raw_processing_requires_explicit_permission(tmp_path):
    fastq = tmp_path / "sample.fastq"
    fastq.write_text("@x\nAC\n+\n!!\n")
    result = Converter().convert(
        InputSpec(fastq, metadata=package()), out_type="h5ad", outdir=tmp_path / "out"
    )
    assert result.status == "failed"
    assert any(d.code == "processing_required" for d in result.items[0].diagnostics)
    assert not (tmp_path / "out").exists()


def test_invalid_loaded_metadata_is_not_published(tmp_path):
    value = package()
    value["series"]["sample_ref"] = ["MISSING_SAMPLE"]
    result = Converter().convert(
        value, out_type="json", outfile=tmp_path / "invalid.json"
    )
    assert result.status == "failed"
    assert not (tmp_path / "invalid.json").exists()
    assert any(d.code == "invalid_metadata" for d in result.items[0].diagnostics)


def test_manifest_companions_are_not_discovered_twice(tmp_path):
    from shutil import copyfile

    root = Path(__file__).parents[2] / "tests/fixtures/studies/E-MTAB-1/inputs"
    for source in root.iterdir():
        copyfile(source, tmp_path / source.name)
    idf = next(tmp_path.glob("*.idf.txt"))
    sdrf = next(tmp_path.glob("*.sdrf.txt"))
    manifest = {
        "schema_version": "1.0",
        "inputs": [
            {"id": "bound", "sources": [str(idf)], "companions": {"sdrf": str(sdrf)}}
        ],
    }
    result = Converter().convert(tmp_path, out_type="json", input_manifest=manifest)
    assert len(result.items) == 1, result.to_dict()


def test_no_destination_does_not_use_current_directory(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    result = Converter().convert(package(), out_type="csv")
    assert result.status == "complete"
    assert not list(tmp_path.iterdir())


def test_projector_invalid_permission_is_separate_from_partial_publication(tmp_path):
    from meta_standards_converter.converters import JSON2TSVConverter
    from meta_standards_converter.metadata.projection.tabular import (
        TabularMetadataProjection,
    )

    class Projector:
        def project_sample(self, **kwargs):
            return TabularMetadataProjection(
                {"x": "value"}, errors=("invalid scientific annotation",)
            )

    converter = Converter(
        services={"json2tsv": JSON2TSVConverter(metadata_projectors=[Projector()])}
    )
    rejected = converter.convert(
        package(), out_type="csv", outfile=tmp_path / "invalid.csv"
    )
    assert rejected.items[0].validation == "invalid"
    assert not (tmp_path / "invalid.csv").exists()
    permitted = converter.convert(
        package(),
        out_type="csv",
        outfile=tmp_path / "permitted.csv",
        output_options={"allow_invalid": True},
    )
    assert permitted.status == "partial"
    assert permitted.items[0].execution == "succeeded"
    assert permitted.items[0].validation == "invalid"
    assert permitted.items[0].completeness == "complete"


def test_failed_detection_keeps_its_own_stage_and_origins(tmp_path):
    result = Converter().convert(
        [InputSpec(package()), tmp_path / "missing.json"],
        out_type="json",
        outdir=tmp_path / "out",
    )
    failed = result.items[1]
    assert failed.diagnostics[-1].stage == "input"
    assert failed.origins == (str(tmp_path / "missing.json"),)


def test_skipped_atlas_datasets_keep_partial_completeness(tmp_path):
    source = (
        Path(__file__).parents[2]
        / "tests/fixtures/edge_cases/atlas-groups/inputs/atlas.json"
    )
    result = Converter().convert(source, out_type="csv", outdir=tmp_path)
    assert result.status == "partial", result.to_dict()
    assert result.items[0].completeness == "partial"
    assert result.items[0].validation == "valid"
    assert (tmp_path / "GSE100.csv").is_file()
