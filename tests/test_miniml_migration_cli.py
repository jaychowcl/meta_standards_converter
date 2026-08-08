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

from jsonschema import Draft202012Validator

from meta_standards_converter.cli.miniml_migrate import main
from meta_standards_converter.miniml import MINiMLCodec, miniml_schema_path
from tests.test_msc_miniml_v2 import package_v2


def legacy_package() -> dict:
    return {
        "miniml_schema_version": "1.0",
        "series": {"iid": "E-MTAB-1"},
        "mage_tab": {"version": "1.1"},
    }


def test_v2_schema_accepts_migrated_package() -> None:
    schema = json.loads(miniml_schema_path().read_text(encoding="utf-8"))
    migrated = MINiMLCodec().migrate_v1(legacy_package()).package.to_mapping()

    Draft202012Validator(schema).validate(migrated)
    assert miniml_schema_path().name == "miniml-package-v2.schema.json"


def test_v2_schema_accepts_complete_typed_package() -> None:
    schema = json.loads(miniml_schema_path().read_text(encoding="utf-8"))

    Draft202012Validator(schema).validate(package_v2())


def test_miniml_migrate_cli_writes_v2_and_diagnostics(tmp_path, capsys) -> None:
    source = tmp_path / "legacy.json"
    destination = tmp_path / "v2.json"
    source.write_text(json.dumps(legacy_package()), encoding="utf-8")

    assert main([str(source), str(destination)]) == 0

    payload = json.loads(destination.read_text(encoding="utf-8"))
    assert payload["miniml_schema_version"] == "2.0"
    assert payload["source"]["format"] == "MAGE-TAB"
    assert json.loads(capsys.readouterr().out)["packages_migrated"] == 1
