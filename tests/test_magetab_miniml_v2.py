from __future__ import annotations

from meta_standards_converter.ae_handlers.ae_constructor import AEConstructor
from meta_standards_converter.ae_handlers.ae_parser import AEParser

from tests.test_ae2json import resolved_input


def test_magetab_parser_folds_semantics_into_msc_miniml_v2() -> None:
    package = AEParser().parse(resolved_input())
    payload = package.to_mapping()

    assert payload["miniml_schema_version"] == "2.0"
    assert payload["source"]["format"] == "MAGE-TAB"
    assert {item["kind"] for item in payload["source"]["documents"]} == {"idf", "sdrf"}
    assert [item["name"] for item in payload["series"]["protocols"]] == [
        "P-collect",
        "P-extract",
    ]
    assert payload["series"]["assay_paths"]
    assert "mage_tab" not in payload


def test_msc_miniml_v2_renders_semantic_magetab() -> None:
    package = AEParser().parse(resolved_input())

    rows = AEConstructor().miniml2magetab(package)

    by_label = {row[0]: row[1:] for row in rows if row}
    assert by_label["Investigation Accession"][0] == "E-MTAB-1"
    assert by_label["Protocol Name"][:2] == ["P-collect", "P-extract"]
    assert isinstance(by_label["SDRF File"][0], list)
