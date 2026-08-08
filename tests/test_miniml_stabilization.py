from __future__ import annotations

import json

import pytest
from jsonschema import Draft202012Validator

from meta_standards_converter.ae_handlers.ae_constructor import AEConstructor
from meta_standards_converter.ae_handlers.ae_parser import AEParser
from meta_standards_converter.miniml import (
    AssayNode,
    MINiMLCodec,
    MINiMLModelError,
    MINiMLPackage,
    Series,
    SourceInfo,
    miniml_schema_path,
)
from tests.test_ae2json import IDF, resolved_input


def _parse_sdrf(header: list[str], row: list[str]):
    text = "\n".join(("\t".join(header), "\t".join(row))) + "\n"
    return AEParser().parse(resolved_input(idf=IDF, sdrfs=[text]))


def _render_sdrf(package):
    rows = AEConstructor().miniml2magetab(package)
    return next(row[1] for row in rows if row and row[0] == "SDRF File")


def test_schema_and_model_accept_the_same_supported_series_and_assay_fields():
    payload = {
        "miniml_schema_version": "2.0",
        "source": {"format": "MAGE-TAB"},
        "database": [], "organization": [], "contributor": [], "platform": [], "sample": [],
        "series": {
            "iid": "E-TEST-1",
            "contributor_ref": [{"ref": "external-person"}],
            "assay_paths": [{"steps": [{"kind": "derived_array_data_file", "name": "result.txt"}]}],
        },
    }
    canonical = MINiMLPackage.from_mapping(payload).to_mapping()
    Draft202012Validator(json.loads(miniml_schema_path().read_text())).validate(canonical)


def test_strict_codec_revalidates_public_dataclass_construction():
    invalid = MINiMLPackage(series=Series(), source=SourceInfo(""))
    with pytest.raises(MINiMLModelError):
        MINiMLCodec().decode(invalid, strict=True)
    with pytest.raises(MINiMLModelError):
        MINiMLCodec().encode(invalid)
    with pytest.raises(MINiMLModelError):
        AssayNode("bogus", "")


def test_blank_and_external_protocol_references_are_compatible():
    blank = _parse_sdrf(["Source Name", "Protocol REF", "Sample Name"], ["s1", "", "x1"])
    assert [step["kind"] for step in blank["series"]["assay_paths"][0]["steps"]] == ["source", "sample"]

    external = _parse_sdrf(["Source Name", "Protocol REF", "Sample Name"], ["s1", "EXTERNAL:P1", "x1"])
    application = external["series"]["assay_paths"][0]["steps"][1]
    assert application == {"kind": "protocol_application", "protocol_ref": "EXTERNAL:P1"}
    assert any(issue.code == "external_protocol_reference" for issue in external.validate())


def test_repeated_sample_layers_remain_an_ordered_assay_path():
    package = _parse_sdrf(
        ["Source Name", "Sample Name", "Protocol REF", "Sample Name"],
        ["s1", "parent", "P-collect", "child"],
    )
    steps = package["series"]["assay_paths"][0]["steps"]
    assert [(step["kind"], step.get("name")) for step in steps] == [
        ("source", "s1"), ("sample", "parent"),
        ("protocol_application", None), ("sample", "child"),
    ]


def test_protocol_application_metadata_and_ontology_companions_round_trip():
    package = _parse_sdrf(
        [
            "Source Name", "Protocol REF", "Performer", "Date", "Parameter Value[temp]",
            "Unit[temperature unit]", "Term Source REF", "Term Accession Number",
            "Comment[note]", "Assay Name", "Technology Type", "Term Source REF",
            "Term Accession Number",
        ],
        [
            "s1", "P-collect", "Jane Doe", "2025-01-02", "37", "degree Celsius", "UO",
            "UO:0000027", "critical", "a1", "RNA sequencing", "EFO", "EFO:0002770",
        ],
    )
    steps = package["series"]["assay_paths"][0]["steps"]
    application = next(step for step in steps if step["kind"] == "protocol_application")
    assert application["performer"] == "Jane Doe"
    assert application["date"] == "2025-01-02"
    assert application["parameter_values"][0]["unit_type"] == "temperature unit"
    assert application["parameter_values"][0]["comments"] == [{"name": "note", "value": "critical"}]
    assay = next(step for step in steps if step["kind"] == "assay")
    assert assay["technology_type"]["term_accession_number"] == "EFO:0002770"
    rendered = _render_sdrf(package)
    assert rendered[0] == [
        "Source Name", "Protocol REF", "Performer", "Date", "Parameter Value[temp]",
        "Unit[temperature unit]", "Term Source REF", "Term Accession Number",
        "Comment[note]", "Assay Name", "Technology Type", "Term Source REF",
        "Term Accession Number",
    ]


def test_factor_qualifier_and_idf_semantics_are_typed():
    idf = IDF.replace(
        "Experimental Design\ttranscription profiling by high throughput sequencing\n",
        "Experimental Design\ttranscription profiling by high throughput sequencing\n"
        "Experimental Design Term Source REF\tEFO\n"
        "Experimental Design Term Accession Number\tEFO:0002768\n"
        "Person Roles\tinvestigator\nPerson Roles Term Source Ref\tEFO\n"
        "Person Roles Term Accession Number\tEFO:0000001\nComment[Goal]\tTest mechanism\n",
    )
    text = "Source Name\tFactor Value[age] (time)\tSample Name\ns1\t5\tx1\n"
    package = AEParser().parse(resolved_input(idf=idf, sdrfs=[text]))
    assert package.series.types[0].term_accession_number == "EFO:0002768"
    assert package.contributors[0].roles[0].value == "investigator"
    assert package.series.comments[0].name == "Goal"
    assert package.series.experiment_date == "2024-01-02"
    factor = next(step for step in package["series"]["assay_paths"][0]["steps"] if step["kind"] == "source")["factor_values"][0]
    assert factor["name"] == "age" and factor["qualifier"] == "time"


def test_multiple_documents_render_without_cross_document_reordering():
    paths = [
        {"document": "a.sdrf", "steps": [{"kind": "source", "name": "s1"}, {"kind": "protocol_application", "protocol_ref": "P"}, {"kind": "sample", "name": "x1"}]},
        {"document": "b.sdrf", "steps": [{"kind": "source", "name": "s2"}, {"kind": "sample", "name": "x2"}, {"kind": "protocol_application", "protocol_ref": "P"}, {"kind": "assay", "name": "a2"}]},
    ]
    from meta_standards_converter.ae_handlers.ae_model import render_miniml_assay_documents
    documents = render_miniml_assay_documents(paths)
    assert documents["a.sdrf"][0] == ["Source Name", "Protocol REF", "Sample Name"]
    assert documents["b.sdrf"][0] == ["Source Name", "Sample Name", "Protocol REF", "Assay Name"]
