# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from __future__ import annotations

from pathlib import Path

from meta_standards_converter.magetab.writer import MAGETabWriter
from meta_standards_converter.magetab.constructor import AEConstructor
from meta_standards_converter.magetab.parser import AEParser

from tests.test_ae2json import resolved_input


def test_magetab_parser_folds_semantics_into_msc_miniml_v2() -> None:
    package = AEParser().parse(resolved_input())
    payload = package.to_mapping()

    assert payload["miniml_schema_version"] == "3.0"
    assert payload["source"]["format"] == "MAGE-TAB"
    assert {item["kind"] for item in payload["source"]["documents"]} == {"idf", "sdrf"}
    assert [item["name"] for item in payload["series"]["protocols"]] == [
        "P-collect",
        "P-extract",
    ]
    assert payload["series"]["protocols"][0]["type"] == {
        "value": "sample collection protocol",
        "term_source_ref": "EFO",
        "term_accession_number": "EFO:0005518",
    }
    assert payload["series"]["pubmed_publication"][0]["status_term_source_ref"] == "EFO"
    assert payload["series"]["pubmed_publication"][0]["status_term_accession_number"] == "EFO:0000001"
    assert payload["series"]["assay_paths"]
    assert "mage_tab" not in payload


def test_msc_miniml_v2_renders_semantic_magetab() -> None:
    package = AEParser().parse(resolved_input())

    rows = AEConstructor().miniml2magetab(package)

    by_label = {row[0]: row[1:] for row in rows if row}
    assert by_label["Investigation Accession"][0] == "E-MTAB-1"
    assert by_label["Protocol Name"][:2] == ["P-collect", "P-extract"]
    assert by_label["Publication Status Term Source REF"] == ["EFO"]
    assert by_label["Publication Status Term Accession Number"] == ["EFO:0000001"]
    assert by_label["Protocol Term Source REF"][:2] == ["EFO", "EFO"]
    assert by_label["Protocol Term Accession Number"][:2] == [
        "EFO:0005518",
        "EFO:0002944",
    ]
    assert not {
        "Status Term Source Ref",
        "Status Term Accession Number",
        "Protocol Type Term Source REF",
        "Protocol Type Term Accession Number",
    }.intersection(by_label)
    assert isinstance(by_label["SDRF File"][0], list)


def test_canonical_idf_companion_labels_parse_without_unmapped_warnings() -> None:
    canonical_idf = (
        resolved_input().idf.text
        .replace("Status Term Source Ref", "Publication Status Term Source REF")
        .replace("Status Term Accession Number", "Publication Status Term Accession Number")
        .replace("Protocol Type Term Source REF", "Protocol Term Source REF")
        .replace("Protocol Type Term Accession Number", "Protocol Term Accession Number")
    )
    parser = AEParser()

    package = parser.parse(resolved_input(idf=canonical_idf))

    protocol_type = package.to_mapping()["series"]["protocols"][0]["type"]
    assert protocol_type["term_source_ref"] == "EFO"
    assert protocol_type["term_accession_number"] == "EFO:0005518"
    assert not [
        warning
        for warning in parser.warnings
        if "Protocol Term Source REF" in warning
        or "Protocol Term Accession Number" in warning
    ]


def test_canonical_idf_companion_labels_take_precedence_over_legacy_aliases() -> None:
    mixed_idf = (
        resolved_input().idf.text
        .replace(
            "Status Term Source Ref\tEFO\n",
            "Publication Status Term Source REF\tCANONICAL\n"
            "Status Term Source Ref\tLEGACY\n",
        )
        .replace(
            "Status Term Accession Number\tEFO:0000001\n",
            "Publication Status Term Accession Number\tCANONICAL:1\n"
            "Status Term Accession Number\tLEGACY:1\n",
        )
        .replace(
            "Protocol Type Term Source REF\tEFO\tEFO\n",
            "Protocol Term Source REF\tCANONICAL\tCANONICAL\n"
            "Protocol Type Term Source REF\tLEGACY\tLEGACY\n",
        )
        .replace(
            "Protocol Type Term Accession Number\tEFO:0005518\tEFO:0002944\n",
            "Protocol Term Accession Number\tCANONICAL:2\tCANONICAL:3\n"
            "Protocol Type Term Accession Number\tLEGACY:2\tLEGACY:3\n",
        )
    )

    payload = AEParser().parse(resolved_input(idf=mixed_idf)).to_mapping()

    publication = payload["series"]["pubmed_publication"][0]
    protocol_type = payload["series"]["protocols"][0]["type"]
    assert publication["status_term_source_ref"] == "CANONICAL"
    assert publication["status_term_accession_number"] == "CANONICAL:1"
    assert protocol_type["term_source_ref"] == "CANONICAL"
    assert protocol_type["term_accession_number"] == "CANONICAL:2"


def test_semantic_magetab_writer_publishes_only_canonical_idf_labels(tmp_path) -> None:
    constructor = AEConstructor()
    rows = constructor.miniml2magetab(AEParser().parse(resolved_input()))

    idf_path = Path(MAGETabWriter().write(rows, out=str(tmp_path)))
    idf_text = idf_path.read_text(encoding="utf-8")
    idf_labels = {line.split("\t", 1)[0] for line in idf_text.splitlines()}

    for label in (
        "Publication Status Term Source REF",
        "Publication Status Term Accession Number",
        "Protocol Term Source REF",
        "Protocol Term Accession Number",
    ):
        assert label in idf_labels
    for label in (
        "Status Term Source Ref",
        "Status Term Accession Number",
        "Protocol Type Term Source REF",
        "Protocol Type Term Accession Number",
    ):
        assert label not in idf_labels
