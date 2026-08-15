# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from __future__ import annotations

from copy import deepcopy

import pytest

from meta_standards_converter.miniml import (
    HarmonizedValue,
    MINiMLCodec,
    MINiMLHarmonizationPatch,
    MINiMLModelError,
    apply_miniml_harmonization_patch,
    harmonization_provenance_index,
    iter_harmonization_operations,
    iter_harmonization_patches,
    miniml_source_fingerprint,
)


def _package(accession: str = "GSE1", source: str = "Human ADPKD kidney") -> dict:
    return {
        "miniml_schema_version": "3.0",
        "source": {"format": "test"},
        "sample": [
            {
                "iid": f"{accession}-S1",
                "channel": [{"source": {"value": source}}],
            }
        ],
        "series": {"iid": accession},
    }


def _operation(*, path: str = "/sample/0/channel/0/source", label: str = "ADPKD") -> dict:
    return {
        "path": path,
        "harmonized_value": {
            "field": "sample_disease_name",
            "value": "autosomal dominant polycystic kidney disease",
            "term_source_ref": "DOID",
            "term_accession_number": "DOID:898",
            "hierarchy_depth": 0,
        },
        "source_evidence": {
            "target_id": "target-1",
            "source_field": "sample_disease_name",
            "source_label": label,
            "source_field_path": path,
            "source_value_path": f"{path}/value",
            "derivation": "split",
            "match_kind": "exact_span",
        },
    }


def test_patch_31_applies_hz_values_and_retains_immutable_provenance() -> None:
    source = _package()
    patch = MINiMLHarmonizationPatch(
        base_sha256=miniml_source_fingerprint(source),
        adds=(_operation(),),
    )

    result = apply_miniml_harmonization_patch(source, patch)
    occurrence = result["sample"][0]["channel"][0]["source"]
    assert occurrence["hz_sample_disease_name"] == (
        "autosomal dominant polycystic kidney disease"
    )
    assert occurrence["hz_sample_disease_name_hierarchy_depth"] == 0

    fragments = iter_harmonization_patches(result)
    assert len(fragments) == 1
    fragment = fragments[0]
    assert fragment["schema_version"] == "1.0"
    assert fragment["patch_schema_version"] == "3.1"
    assert fragment["package_index"] == 0
    assert fragment["status"] == "applied"
    assert fragment["operations"][0]["source_evidence"]["source_label"] == "ADPKD"

    # Applying the identical document-level patch again is a no-op even though
    # the biological content now contains the values added by the first pass.
    assert apply_miniml_harmonization_patch(result, patch) == result
    assert source == _package()


def test_provenance_index_canonicalizes_a_package_only_once(monkeypatch) -> None:
    source = _package()
    patch = MINiMLHarmonizationPatch(
        base_sha256=miniml_source_fingerprint(source),
        adds=(_operation(),),
    )
    result = apply_miniml_harmonization_patch(source, patch)
    decode_calls = 0
    original_decode = MINiMLCodec.decode

    def counted_decode(self, *args, **kwargs):
        nonlocal decode_calls
        decode_calls += 1
        return original_decode(self, *args, **kwargs)

    monkeypatch.setattr(MINiMLCodec, "decode", counted_decode)

    index = harmonization_provenance_index(result)

    assert index["GSE1-S1"]["sample_disease_name"][0]["source_label"] == "ADPKD"
    assert decode_calls == 1


def test_patch_30_is_readable_and_has_no_source_provenance() -> None:
    source = _package()
    old = MINiMLHarmonizationPatch.from_mapping(
        {
            "schema_version": "3.0",
            "miniml_schema_version": "3.0",
            "base_sha256": miniml_source_fingerprint(source),
            "adds": [
                {
                    "path": "/sample/0/channel/0/source",
                    "harmonized_value": {
                        "field": "tissue_name",
                        "value": "kidney",
                    },
                }
            ],
        }
    )
    assert old.schema_version == "3.0"
    result = apply_miniml_harmonization_patch(source, old)
    operation = next(iter(iter_harmonization_operations(result)))
    assert operation["harmonized_value"]["value"] == "kidney"
    assert "source_evidence" not in operation


@pytest.mark.parametrize(
    ("source", "label", "kind"),
    [
        ("Human ADPKD kidney", "PKD", "exact_span"),
        ("Human ADPKD kidney", "ADPKD", "exact_value"),
    ],
)
def test_invalid_exact_source_claim_rejects_the_patch(
    source: str, label: str, kind: str
) -> None:
    document = _package(source=source)
    operation = _operation(label=label)
    operation["source_evidence"]["match_kind"] = kind
    patch = MINiMLHarmonizationPatch(
        base_sha256=miniml_source_fingerprint(document), adds=(operation,)
    )
    with pytest.raises(MINiMLModelError, match="source evidence"):
        apply_miniml_harmonization_patch(document, patch)


def test_interpreted_evidence_is_auditable_but_not_an_exact_claim() -> None:
    document = _package(source="ICM")
    operation = _operation(label="ischemic cardiomyopathy")
    operation["source_evidence"]["match_kind"] = "interpreted"
    patch = MINiMLHarmonizationPatch(
        base_sha256=miniml_source_fingerprint(document), adds=(operation,)
    )
    result = apply_miniml_harmonization_patch(document, patch)
    index = harmonization_provenance_index(result)
    assert index["GSE1-S1"]["sample_disease_name"][0]["match_kind"] == "interpreted"


def test_multi_package_patch_partitions_and_rebases_paths() -> None:
    document = [_package("GSE1"), _package("GSE2", "Human control kidney")]
    operation_one = _operation(path="/0/sample/0/channel/0/source")
    operation_two = _operation(
        path="/1/sample/0/channel/0/source", label="control"
    )
    operation_two["harmonized_value"] = HarmonizedValue(
        field="sample_disease_status",
        value="normal",
        term_source_ref="PATO",
        term_accession_number="PATO:0000461",
        hierarchy_depth=0,
    ).to_annotation_mapping()
    operation_two["harmonized_value"].pop("index")
    operation_two["source_evidence"]["source_field_path"] = (
        "/1/sample/0/channel/0/source"
    )
    operation_two["source_evidence"]["source_value_path"] = (
        "/1/sample/0/channel/0/source/value"
    )
    patch = MINiMLHarmonizationPatch(
        base_sha256=miniml_source_fingerprint(document),
        adds=(operation_one, operation_two),
    )

    result = apply_miniml_harmonization_patch(document, patch)
    for package in result:
        fragment = iter_harmonization_patches(package)[0]
        assert len(fragment["operations"]) == 1
        assert fragment["operations"][0]["path"] == "/sample/0/channel/0/source"
    assert result[1]["sample"][0]["channel"][0]["source"]["hz_sample_disease_status"] == "normal"


def test_fingerprint_ignores_only_retained_harmonization_fragments() -> None:
    source = _package()
    baseline = miniml_source_fingerprint(source)
    with_fragment = deepcopy(source)
    with_fragment["extensions"] = {
        "vendor": {"retained": True},
        "msc_harmonization": {"schema_version": "1.0", "patches": []},
    }
    without_fragment = deepcopy(with_fragment)
    without_fragment["extensions"].pop("msc_harmonization")
    assert miniml_source_fingerprint(with_fragment) == miniml_source_fingerprint(
        without_fragment
    )
    assert miniml_source_fingerprint(with_fragment) != baseline
