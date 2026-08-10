# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from __future__ import annotations

from meta_standards_converter.converters.json2h5ad import JSON2H5ADConverter
from meta_standards_converter.converters.miniml_metadata import MINiMLMetadataService


def _sample(iid: str) -> dict:
    return {"iid": iid, "channel": [{}]}


def test_typed_protocol_projection_is_sample_bound_and_preserves_declared_ontology() -> None:
    first = _sample("S1")
    second = _sample("S2")
    package = {
        "series": {
            "protocols": [
                {
                    "name": "P-collect",
                    "type": {
                        "value": "sample collection protocol",
                        "term_source_ref": "EFO",
                        "term_accession_number": "EFO:0005518",
                    },
                },
                {
                    "name": "P-extract",
                    "type": {
                        "value": "nucleic acid extraction protocol",
                        "term_source_ref": "OBI",
                        "term_accession_number": "OBI:0000257",
                    },
                },
            ],
            "assay_paths": [
                {
                    "steps": [
                        {"kind": "sample", "name": "S1", "sample_ref": "S1"},
                        {"kind": "protocol_application", "protocol_ref": "P-collect"},
                    ]
                },
                {
                    "steps": [
                        {"kind": "sample", "name": "S2", "sample_ref": "S2"},
                        {"kind": "protocol_application", "protocol_ref": "P-extract"},
                    ]
                },
            ],
        },
        "database": [],
    }

    first_values = MINiMLMetadataService().sample_metadata_values(first, package)
    second_values = MINiMLMetadataService().sample_metadata_values(second, package)

    assert first_values["protocol_types"] == ("sample collection protocol",)
    assert first_values["protocol_term_source_refs"] == ("EFO",)
    assert first_values["protocol_term_accession_numbers"] == ("EFO:0005518",)
    assert second_values["protocol_types"] == ("nucleic acid extraction protocol",)
    assert second_values["protocol_term_source_refs"] == ("OBI",)
    assert second_values["protocol_term_accession_numbers"] == ("OBI:0000257",)


def test_typed_protocol_projection_falls_back_to_all_declared_without_sample_binding() -> None:
    package = {
        "series": {
            "protocols": [
                {"name": "P1", "type": {"value": "collection", "term_source_ref": "SRC", "term_accession_number": "SRC:1"}},
                {"name": "P2", "type": {"value": "extraction", "term_source_ref": "SRC", "term_accession_number": "SRC:2"}},
            ]
        },
        "database": [],
    }

    values = MINiMLMetadataService().sample_metadata_values(_sample("S1"), package)

    assert values["protocol_types"] == ("collection", "extraction")
    assert values["protocol_term_accession_numbers"] == ("SRC:1", "SRC:2")


def test_legacy_protocol_inference_is_used_only_without_typed_protocols() -> None:
    sample = {"iid": "S1", "channel": [{"extract_protocol": "legacy details"}]}

    values = MINiMLMetadataService().sample_metadata_values(
        sample, {"series": {}, "database": []}
    )

    assert values["protocol_types"] == ("nucleic acid extraction protocol",)
    assert values["protocol_term_source_refs"] == ("EFO",)


def test_material_type_prefers_sample_bound_assay_nodes_then_typed_characteristic() -> None:
    sample = {
        "iid": "S1",
        "channel": [
            {
                "source": {"value": "lung"},
                "characteristics": [
                    {"name": "material type", "value": "whole organism"}
                ],
            }
        ],
    }
    package = {
        "series": {
            "assay_paths": [
                {
                    "steps": [
                        {
                            "kind": "sample",
                            "name": "S1",
                            "sample_ref": "S1",
                            "material_type": {"value": "tissue specimen"},
                        }
                    ]
                }
            ]
        },
        "database": [],
    }

    service = MINiMLMetadataService()
    assert service.sample_metadata_values(sample, package)["material_type"] == (
        "tissue specimen",
    )

    package["series"]["assay_paths"] = []
    assert service.sample_metadata_values(sample, package)["material_type"] == (
        "whole organism",
    )


def test_h5ad_projection_uses_the_shared_typed_protocol_and_material_semantics() -> None:
    sample = _sample("S1")
    package = {
        "series": {
            "protocols": [
                {
                    "name": "P1",
                    "type": {
                        "value": "bespoke protocol",
                        "term_source_ref": "TEST",
                        "term_accession_number": "TEST:1",
                    },
                }
            ],
            "assay_paths": [
                {
                    "steps": [
                        {
                            "kind": "sample",
                            "sample_ref": "S1",
                            "material_type": {"value": "fresh specimen"},
                        },
                        {"kind": "protocol_application", "protocol_ref": "P1"},
                    ]
                }
            ],
        },
        "database": [],
    }

    metadata = JSON2H5ADConverter()._sample_metadata_values(sample, package)

    assert metadata["material_type"] == ("fresh specimen",)
    assert metadata["protocol_types"] == ("bespoke protocol",)
    assert metadata["protocol_term_source_refs"] == ("TEST",)
    assert metadata["protocol_term_accession_numbers"] == ("TEST:1",)
