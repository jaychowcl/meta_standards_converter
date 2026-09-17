# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Source routing preserves provider identity and native parser semantics."""

import io
import json
import tarfile
from pathlib import Path
from types import SimpleNamespace

import pytest
from meta_standards_converter.converters import Converter, InputSpec
from meta_standards_converter.converters.archive_results import (
    ArchiveImportResult,
    StudyImportOutcome,
)
from meta_standards_converter.miniml import MINiMLCodec
from tests.converters.test_json2tsv import package

ROOT = Path(__file__).parents[2]


@pytest.mark.parametrize("target", ["json", "magetab"])
def test_geo_accession_and_local_xml_produce_equivalent_artifact_bytes(
    target, tmp_path
):
    from meta_standards_converter.converters.geo2ae import GEO2AEConverter
    from meta_standards_converter.converters.geo2json import GEO2JSONConverter
    from meta_standards_converter.converters.json2ae import JSON2AEConverter
    from meta_standards_converter.magetab.constructor import AEConstructor
    from meta_standards_converter.miniml.geo_parser import GEOParser

    source = ROOT / "tests/fixtures/studies/GSE100/inputs/geo.xml"
    text = source.read_text()

    class Fetcher:
        def fetch_gse_miniml(self, **kwargs):
            return text

    class Parser:
        def parse(self, **kwargs):
            return GEOParser().parse(
                kwargs["miniml"], remove_empty=kwargs.get("remove_empty", True)
            )

    class Enricher:
        def enrich(self, data):
            return data

    constructor = AEConstructor(
        evidence_resolver=SimpleNamespace(
            publications=lambda data: [], sample_runs=lambda handler, kind: {}
        )
    )
    services = {
        "geo2ae": GEO2AEConverter(
            geo_fetcher=Fetcher(),
            parser=Parser(),
            enricher=Enricher(),
            ae_constructor=constructor,
        ),
        "geo2json": GEO2JSONConverter(
            geo_fetcher=Fetcher(), parser=Parser(), enricher=Enricher()
        ),
        "json2ae": JSON2AEConverter(ae_constructor=constructor),
    }
    accession = Converter(services=services).convert(
        "GSE100", out_type=target, outdir=tmp_path / "accession"
    )
    local = Converter(services=services).convert(
        source,
        out_type=target,
        outdir=tmp_path / "local",
        output_options={"enrich": False} if target == "magetab" else None,
    )
    assert accession.status == local.status == "complete", (
        accession.to_dict(),
        local.to_dict(),
    )
    expected = {
        p.name: p.read_bytes()
        for p in (tmp_path / "accession").iterdir()
        if p.suffix in {".json", ".txt"}
    }
    actual = {
        p.name: p.read_bytes()
        for p in (tmp_path / "local").iterdir()
        if p.suffix in {".json", ".txt"}
    }
    assert expected == actual
    assert len(expected) == (1 if target == "json" else 2)


def test_enrichment_is_applied_at_the_route_stage_that_owns_it():
    calls = []

    class Exporter:
        def convert_loaded(self, loaded, **options):
            calls.append(options["enrich"])
            return [
                [
                    ["Investigation Accession", "GSE1"],
                    ["SDRF File", [["Sample Name"], ["GSM1"]]],
                ]
            ]

    class Importer:
        def convert(self, *args, **kwargs):
            return [MINiMLCodec().decode(package()).package]

    converter = Converter(services={"json2ae": Exporter(), "geo2json": Importer()})
    local = converter.convert(
        ROOT / "tests/fixtures/studies/GSE100/inputs/geo.xml", out_type="magetab"
    )
    assert local.status == "complete"
    assert calls == [True]
    from meta_standards_converter.converters.unified.handlers import InputHandler
    from meta_standards_converter.converters.unified.engine import ExecutionContext
    from meta_standards_converter.converters.unified.options import RUNTIME_DEFAULTS
    from contextlib import ExitStack

    with ExitStack() as stack:
        context = ExecutionContext(
            converter, RUNTIME_DEFAULTS.copy(), "json", None, None, stack
        )
        loaded = InputHandler("geo_accession").load(InputSpec("GSE1"), context)
        assert loaded.enrichment_applied is True


@pytest.mark.parametrize(
    "accession,provider",
    [
        ("SRR1", "sra"),
        ("ERR1", "ena"),
        ("DRR1", "ena"),
        ("PRJNA1", "sra"),
        ("PRJEB1", "ena"),
        ("SAMN1", "sra"),
    ],
)
def test_accession_defaults_and_partial_propagation(accession, provider, tmp_path):
    calls = []

    class Archive:
        def __init__(self, name):
            self.name = name

        def convert(self, value, **options):
            calls.append((self.name, value, options))
            return ArchiveImportResult(
                value,
                [
                    StudyImportOutcome(
                        "SRP1",
                        "SRP1",
                        "partial",
                        MINiMLCodec().decode(package()).package,
                        issues=["optional source unavailable"],
                    )
                ],
            )

    converter = Converter(
        services={"sra2json": Archive("sra"), "ena2json": Archive("ena")}
    )
    result = converter.convert(accession, out_type="csv", outdir=tmp_path)
    assert result.status == "partial", result.to_dict()
    assert calls[0][0] == provider and result.items[0].provider == provider
    assert result.items[0].artifacts
    forced = converter.convert(accession, out_type="json", force_in_type="sra")
    assert forced.items[0].provider == "sra"


def test_missing_paths_never_reach_provider(tmp_path):
    class Never:
        def convert(self, *a, **kw):
            raise AssertionError("must not retrieve")

    result = Converter(services={"ae2json": Never()}).convert(
        tmp_path / "E-MTAB-1.idf.txt", out_type="json"
    )
    assert result.status == "failed"
    assert result.items[0].diagnostics[0].code == "missing_file"


def test_geo_xml_and_archive_have_same_packages(tmp_path):
    source = ROOT / "tests/fixtures/studies/GSE100/inputs/geo.xml"
    converter = Converter()
    xml = converter.convert(source, out_type="json")
    assert xml.status == "complete", xml.to_dict()
    archive = tmp_path / "family.tgz"
    with tarfile.open(archive, "w:gz") as tar:
        tar.add(source, arcname="GSE100_family.xml")
    packed = converter.convert(archive, out_type="json")
    assert packed.status == "complete", packed.to_dict()
    assert packed.items[0].payload == xml.items[0].payload
    with tarfile.open(archive, "w:gz") as tar:
        member = tarfile.TarInfo("../escape.xml")
        member.size = 1
        tar.addfile(member, io.BytesIO(b"x"))
    assert converter.convert(archive, out_type="json").status == "failed"
    assert not (tmp_path.parent / "escape.xml").exists()


@pytest.mark.parametrize("provider", ["sra", "ena"])
def test_native_record_bundle(provider):
    directory = ROOT / "docs" / provider / "fixtures/SRX7812918"
    paths = (
        [directory / "sra-efetch.xml"]
        if provider == "sra"
        else [
            directory / (kind + ".xml")
            for kind in ("study", "sample", "experiment", "run")
        ]
    )
    result = Converter().convert(
        InputSpec(paths, in_type=provider + "_xml"), out_type="json"
    )
    assert result.status in {"complete", "partial"}, result.to_dict()
    sample = result.items[0].payload[0].to_mapping()["sample"][0]
    assert sample["sra_run"][0]["run"] == "SRR11192680"


def test_magetab_directory_consumes_sdrf_once():
    directory = ROOT / "tests/fixtures/studies/E-MTAB-1/inputs"
    result = Converter().convert(directory, out_type="json")
    assert len(result.items) == 1, result.to_dict()
    assert result.status == "complete", result.to_dict()
    assert result.items[0].payload[0].samples


def test_untrusted_remote_url_is_rejected_without_network():
    result = Converter().convert("http://127.0.0.1/private.json", out_type="json")
    assert result.status == "failed"


def test_manifest_rejects_unknown_fields_and_conflicting_ids():
    invalid = {
        "schema_version": "1.0",
        "inputs": [{"id": "x", "sources": ["GSE1"], "typo": True}],
    }
    with pytest.raises(ValueError):
        Converter().convert(out_type="json", input_manifest=invalid)


def test_local_xml_force_does_not_bypass_root_validation():
    result = Converter().convert("<ROOT/>", force_in_type="geo_xml", out_type="json")
    assert result.status == "failed"


def test_discovery_groups_ena_companions_by_directory(tmp_path):
    import shutil

    original = ROOT / "docs/ena/fixtures/SRX7812918"
    for kind in ("study", "sample", "experiment", "run"):
        shutil.copyfile(original / (kind + ".xml"), tmp_path / (kind + ".xml"))
    result = Converter().convert(tmp_path, out_type="json")
    assert len(result.items) == 1, result.to_dict()
    assert result.status in {"complete", "partial"}, result.to_dict()
    assert result.items[0].dataset_ids == ("PRJNA609050",)


def test_native_conflicting_duplicate_records_are_rejected(tmp_path):
    original = ROOT / "docs/ena/fixtures/SRX7812918"
    paths = [
        original / (kind + ".xml") for kind in ("study", "sample", "experiment", "run")
    ]
    other = tmp_path / "other.xml"
    other.write_text(
        paths[0].read_text().replace("Microbiota analysis", "Conflicting analysis")
    )
    result = Converter().convert(
        InputSpec([*paths, other], in_type="ena_xml"), out_type="json"
    )
    assert result.status == "failed"
    assert any(d.code == "conflicting_records" for d in result.items[0].diagnostics)


def test_same_destination_from_different_inputs_is_never_overwritten(tmp_path):
    result = Converter().convert(
        [InputSpec(package(), id="a"), InputSpec(package(), id="b")],
        out_type="json",
        outdir=tmp_path,
        runtime_options={"overwrite": True},
    )
    assert [i.status for i in result.items] == ["complete", "failed"]
    assert any(d.code == "output_collision" for d in result.items[1].diagnostics)


@pytest.mark.parametrize("provider", ["sra", "ena"])
def test_native_records_object_is_defensive_and_not_mutated(provider):
    from copy import deepcopy
    from xml.etree import ElementTree as ET
    from meta_standards_converter.sources.archive_support import StudyRecords, StudySeed

    directory = ROOT / "docs" / provider / "fixtures/SRX7812918"
    paths = (
        [directory / "sra-efetch.xml"]
        if provider == "sra"
        else [
            directory / (kind + ".xml")
            for kind in ("study", "sample", "experiment", "run")
        ]
    )
    roots = [ET.fromstring(path.read_text()) for path in paths]
    study = next(
        root.find(".//EXPERIMENT/STUDY_REF").get("accession")
        for root in roots
        if root.find(".//EXPERIMENT/STUDY_REF") is not None
    )
    records = StudyRecords(StudySeed(study, study), xml=roots)
    before = deepcopy(records)
    result = Converter().convert(
        records, out_type="json", force_in_type=provider + "_records"
    )
    assert result.status in {"complete", "partial"}, result.to_dict()
    assert records.issues == before.issues
    assert [ET.tostring(root) for root in records.xml] == [
        ET.tostring(root) for root in before.xml
    ]
    rejected = Converter().convert(
        StudyRecords(StudySeed(study, study)),
        out_type="json",
        force_in_type=provider + "_records",
    )
    assert rejected.status == "failed"
