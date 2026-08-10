# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from __future__ import annotations

import json
from pathlib import Path

from meta_standards_converter.cli.miniml_migrate import main
import meta_standards_converter.miniml as miniml
from meta_standards_converter.miniml import MINiMLCodec
from tests.test_msc_miniml_v2 import package_v2


def legacy_package() -> dict:
    return {
        "miniml_schema_version": "1.0",
        "series": {"iid": "E-MTAB-1"},
        "mage_tab": {"version": "1.1"},
    }


def test_json_schema_is_not_public_or_packaged() -> None:
    package_dir = Path(miniml.__file__).resolve().parent

    assert not hasattr(miniml, "miniml_schema_path")
    assert not list(package_dir.glob("*.schema.json"))


def test_python_model_accepts_migrated_and_complete_typed_packages() -> None:
    migrated = MINiMLCodec().migrate_v1(legacy_package()).package.to_mapping()

    assert MINiMLCodec().decode(migrated, strict=True).package.to_mapping() == migrated
    assert MINiMLCodec().decode(package_v2(), strict=True).package.to_mapping() == package_v2()


def test_miniml_migrate_cli_writes_v2_and_diagnostics(tmp_path, capsys) -> None:
    source = tmp_path / "legacy.json"
    destination = tmp_path / "v2.json"
    source.write_text(json.dumps(legacy_package()), encoding="utf-8")

    assert main([str(source), str(destination)]) == 0

    payload = json.loads(destination.read_text(encoding="utf-8"))
    assert payload["miniml_schema_version"] == "2.0"
    assert payload["source"]["format"] == "MAGE-TAB"
    assert json.loads(capsys.readouterr().out)["packages_migrated"] == 1


def test_v1_migrator_folds_private_harmonization_fields_into_annotations() -> None:
    legacy = legacy_package()
    legacy["sample"] = [{
        "iid": "GSM1",
        "channel": [{
            "pre_hz_label": "Homo sapiens",
            "hz_species_name": {
                "value": "Homo sapiens",
                "id": "NCBITaxon:9606",
                "onto": "NCBITaxon",
            },
            "characteristics": [
                {"tag": "disease", "value": "raw case"},
                {"tag": "hz_disease_category", "value": "lung carcinoma"},
                {"tag": "hz_disease_category_id", "value": "MONDO:0008903"},
                {"tag": "hz_disease_category_onto", "value": "MONDO"},
            ],
        }],
    }]

    migrated = MINiMLCodec().migrate_v1(legacy).package.to_mapping()
    channel = migrated["sample"][0]["channel"][0]
    characteristics = {item["name"]: item for item in channel["characteristics"]}

    assert "hz_" not in json.dumps(migrated)
    assert characteristics["organism"]["annotations"][0]["field"] == "species_name"
    assert characteristics["organism"]["annotations"][0]["term_accession_number"] == "NCBITaxon:9606"
    assert characteristics["disease"]["annotations"][0]["field"] == "disease_category"
    assert characteristics["disease"]["annotations"][0]["term_accession_number"] == "MONDO:0008903"
