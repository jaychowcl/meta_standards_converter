# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================

"""The supported canonical boundary is v3, including fresh legacy ingestion."""
import copy
import json

import pytest

from meta_standards_converter.miniml import MINiMLCodec, MINiMLPackage, MINiMLV1Migrator
from meta_standards_converter.sources.json import JSONPackageSource


def obsolete_package():
    return {"miniml_schema_version": "2.0", "source": {"format": "test"}, "series": {"iid": "GSE1"}}


@pytest.mark.parametrize("decode", [MINiMLCodec().decode, MINiMLPackage.from_mapping, MINiMLV1Migrator().migrate])
def test_v2_is_rejected_at_every_decode_boundary(decode):
    with pytest.raises(ValueError, match="3.0"):
        decode(obsolete_package())


def test_v2_migration_is_not_a_public_api():
    import meta_standards_converter.miniml as miniml
    assert not hasattr(miniml, "MINiMLV2Migrator")
    assert not hasattr(MINiMLCodec, "migrate_v2")


def test_v2_migration_cli_does_not_write_output(tmp_path):
    from meta_standards_converter.cli.miniml_migrate import main
    source, destination = tmp_path / "old.json", tmp_path / "new.json"
    source.write_text(json.dumps(obsolete_package()))
    with pytest.raises(ValueError, match="3.0"):
        main([str(source), str(destination)])
    assert not destination.exists()


@pytest.mark.parametrize("envelope", [False, True])
def test_json_source_rejects_v2(envelope, tmp_path):
    payload = obsolete_package()
    if envelope:
        payload = {"miniml_json": payload}
    source = tmp_path / "input.json"
    source.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="3.0"):
        JSONPackageSource().load(source)


@pytest.mark.parametrize("kind", ["ae", "tsv", "h5ad", "obs"])
@pytest.mark.parametrize("atlas", [False, True])
def test_converters_reject_v2_without_destination_artifacts(kind, atlas, tmp_path):
    from pathlib import Path
    from meta_standards_converter.converters import (
        JSON2AEConverter, JSON2TSVConverter, JSON2H5ADConverter, JSON2OBSConverter,
    )
    payload = obsolete_package()
    if atlas:
        document = json.loads((Path(__file__).parent / "fixtures/edge_cases/atlas-groups/inputs/atlas.json").read_text())
        document["datasets"] = document["datasets"][:1]
        document["datasets"][0]["metadata"] = payload
        document["summary"].update(dataset_count=1, failed_dataset_count=0)
        payload = document
    source, destination = tmp_path / "old.json", tmp_path / "destination"
    source.write_text(json.dumps(payload))
    conversions = {
        "ae": lambda: JSON2AEConverter().convert(str(source), enrich=False, platform_handler="generic"),
        "tsv": lambda: JSON2TSVConverter().convert_source(source, destination),
        "h5ad": lambda: JSON2H5ADConverter().convert(str(source), out=str(destination)),
        "obs": lambda: JSON2OBSConverter().convert(source, outdir=destination),
    }
    with pytest.raises(ValueError, match="3.0"):
        conversions[kind]()
    assert not destination.exists()


def test_legacy_ingestion_constructs_v3_without_annotation_intermediates(monkeypatch):
    original = MINiMLPackage.from_mapping.__func__
    seen = []

    def capture(cls, value):
        seen.append(copy.deepcopy(value))
        assert value["miniml_schema_version"] == "3.0"
        assert '"annotations"' not in json.dumps(value)
        return original(cls, value)

    monkeypatch.setattr(MINiMLPackage, "from_mapping", classmethod(capture))
    legacy = {"series": {"iid": "GSE1"}, "sample": [{"iid": "GSM1", "channel": [{
        "characteristics": [{"tag": "disease", "value": "raw case"},
                            {"tag": "hz_disease", "value": "fibrosis"},
                            {"tag": "hz_disease_id", "value": "MONDO:1"},
                            {"tag": "hz_disease_hierarchy_depth", "value": 0},
                            {"tag": "hz_disease(1)", "value": "lung disease"},
                            {"tag": "hz_disease_id(1)", "value": "MONDO:2"}]}]}]}
    before = copy.deepcopy(legacy)
    package = MINiMLV1Migrator().migrate(legacy).package.to_mapping()
    rows = package["sample"][0]["channel"][0]["characteristics"]
    assert {row["name"]: row["value"] for row in rows} == {
        "disease": "raw case", "hz_disease": "fibrosis", "hz_disease_id": "MONDO:1",
        "hz_disease_hierarchy_depth": 0, "hz_disease(1)": "lung disease", "hz_disease_id(1)": "MONDO:2"}
    assert seen
    assert legacy == before


def test_breaking_distribution_version():
    import tomllib
    from pathlib import Path
    config = tomllib.loads((Path(__file__).parents[1] / "pyproject.toml").read_text())
    assert config["project"]["version"] == "8.0.0"


def test_legacy_channel_groups_preserve_explicit_indexes():
    legacy = {"series": {"iid": "GSE1"}, "sample": [{"iid": "GSM1", "channel": [{
        "hz_disease(2)": "fibrosis", "hz_disease_id(2)": "MONDO:1",
        "hz_disease_hierarchy_depth(2)": 0,
    }]}]}
    package = MINiMLV1Migrator().migrate(legacy).package.to_mapping()
    rows = package["sample"][0]["channel"][0]["characteristics"]
    assert {row["name"]: row["value"] for row in rows} == {
        "hz_disease(2)": "fibrosis", "hz_disease_id(2)": "MONDO:1",
        "hz_disease_hierarchy_depth(2)": 0,
    }


@pytest.mark.parametrize("command", ["json2ae", "json2tsv", "json2h5ad", "json2obs"])
def test_cli_v2_failure_does_not_publish_files(command, tmp_path):
    import importlib
    source = tmp_path / "old.json"
    source.write_text(json.dumps(obsolete_package()))
    output = tmp_path / "output"
    args = [str(source), "--outdir" if command == "json2obs" else "--out", str(output)]
    assert importlib.import_module("meta_standards_converter.cli." + command).main(args) == 1
    assert not list(output.rglob("*"))


def test_legacy_assay_value_and_unit_construct_v3_directly(monkeypatch):
    original = MINiMLPackage.from_mapping.__func__

    def capture(cls, value):
        assert value["miniml_schema_version"] == "3.0"
        assert '"annotations"' not in json.dumps(value)
        return original(cls, value)

    monkeypatch.setattr(MINiMLPackage, "from_mapping", classmethod(capture))
    legacy = {"series": {"iid": "GSE1"}, "mage_tab": {"model": {
        "protocols": [{"name": "treatment"}],
        "assay_paths": [{"steps": [
            {"kind": "protocol", "value": "treatment"},
            {"kind": "attribute", "attribute_type": "parameter value", "name": "duration",
             "value": "30", "hz_value": "thirty", "unit": "min", "hz_unit": "minute",
             "hz_unit_id": "UO:0000031", "hz_unit_onto": "UO", "hz_unit_hierarchy_depth": 0},
        ]}],
    }}}
    result = MINiMLV1Migrator().migrate(legacy).package.to_mapping()
    parameter = result["series"]["assay_paths"][0]["steps"][0]["parameter_values"][0]
    assert parameter["value"] == "30"
    assert parameter["hz_value"] == "thirty"
    assert parameter["unit"] == {
        "value": "min", "hz_unit": "minute", "hz_unit_id": "UO:0000031",
        "hz_unit_onto": "UO", "hz_unit_hierarchy_depth": 0,
    }
