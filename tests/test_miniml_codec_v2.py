# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from __future__ import annotations

from dataclasses import FrozenInstanceError

import pytest

from meta_standards_converter.miniml import (
    FASTQFile,
    MINiMLBatchDecodeResult,
    MINiMLCodec,
    MINiMLCompatibilityError,
    MINiMLDecodeResult,
    MINiMLPackage,
    PubMedPublication,
    SRARun,
    Sample,
    Series,
    SourceInfo,
)


def package_payload() -> dict:
    return {
        "miniml_schema_version": "2.0",
        "source": {"format": "test"},
        "series": {
            "iid": "GSE1",
            "pubmed_publication": [
                {
                    "pubmed_id": "123",
                    "doi": "10.1/example",
                    "title": "Study",
                }
            ],
            "sample_ref": [{"ref": "GSM1"}],
        },
        "sample": [
            {
                "iid": "GSM1",
                "sra_accession": ["SRX1"],
                "ena_accession": ["ERP1"],
                "sra_run": [
                    {
                        "run": "SRR1",
                        "study": "ERP1",
                        "fastq_files": [
                            {
                                "filename": "reads.fastq.gz",
                                "uri": "https://example/reads.fastq.gz",
                                "md5": "0123456789abcdef0123456789abcdef",
                            }
                        ],
                    }
                ],
            }
        ],
    }


def test_codec_decodes_typed_enrichment_from_v2_input() -> None:
    result = MINiMLCodec().decode(package_payload())

    assert isinstance(result, MINiMLDecodeResult)
    assert isinstance(result.package, MINiMLPackage)
    assert isinstance(result.package.series.pubmed_publications[0], PubMedPublication)
    assert isinstance(result.package.samples[0].sra_runs[0], SRARun)
    assert isinstance(result.package.samples[0].sra_runs[0].fastq_files[0], FASTQFile)
    assert MINiMLCodec().encode(result.package)["miniml_schema_version"] == "2.0"


def test_codec_returns_warnings_and_strict_mode_promotes_them() -> None:
    payload = package_payload()
    payload["series"]["sample_ref"].append({"ref": "GSM404"})

    compatible = MINiMLCodec().decode(payload)

    assert [issue.code for issue in compatible.diagnostics] == ["unresolved_reference"]
    with pytest.raises(MINiMLCompatibilityError, match="unresolved_reference"):
        MINiMLCodec().decode(payload, strict=True)


def test_codec_decodes_one_or_many_packages() -> None:
    single = MINiMLCodec().decode_many(package_payload())
    multiple = MINiMLCodec().decode_many([package_payload(), package_payload()])

    assert isinstance(single, MINiMLBatchDecodeResult)
    assert len(single.packages) == 1
    assert len(multiple.packages) == 2


def test_typed_package_is_immutable_and_preserves_typed_annotations() -> None:
    payload = package_payload()
    payload["sample"][0]["channel"] = [
        {
            "source": "lung",
            "characteristics": [{
                "name": "organism",
                "value": "Homo sapiens",
                "annotations": [
                    {
                        "field": "organism",
                        "value": "Homo sapiens",
                        "term_source_ref": "NCBITaxon",
                        "term_accession_number": "NCBITaxon:9606",
                    }
                ],
            }],
        }
    ]
    package = MINiMLCodec().decode(payload).package

    annotation = package.samples[0].channels[0].characteristics[0].annotations[0]
    assert annotation.value == "Homo sapiens"
    assert annotation.term_accession_number == "NCBITaxon:9606"
    with pytest.raises(FrozenInstanceError):
        package.version = "changed"  # type: ignore[misc]
    with pytest.raises(FrozenInstanceError):
        package.samples[0].channels[0].characteristics[0].annotations[0].value = "changed"  # type: ignore[misc]


def test_direct_model_construction_deep_freezes_collections_and_extensions() -> None:
    root_extensions = {"nested": {"values": []}}
    series_extensions = {"vendor": []}
    samples = [Sample(iid="GSM1")]
    package = MINiMLPackage(
        series=Series(iid="GSE1", extras=series_extensions),
        source=SourceInfo("test"),
        samples=samples,
        extensions=root_extensions,
    )

    root_extensions["nested"]["values"].append("mutated")
    series_extensions["vendor"].append("mutated")
    samples.append(Sample(iid="GSM2"))

    assert package.to_mapping()["extensions"] == {"nested": {"values": []}}
    assert package.series.to_mapping()["extensions"] == {"vendor": []}
    assert tuple(item.iid for item in package.samples) == ("GSM1",)
    with pytest.raises(TypeError):
        package.extensions["new"] = "value"  # type: ignore[index]


def test_codec_dump_is_deterministic_and_replaces_existing_file(tmp_path) -> None:
    codec = MINiMLCodec()
    package = codec.decode(package_payload()).package
    destination = tmp_path / "package.json"
    destination.write_text("old", encoding="utf-8")

    codec.dump(package, destination)
    first = destination.read_bytes()
    codec.dump(package, destination)

    assert destination.read_bytes() == first
    assert not list(tmp_path.glob(".*.tmp"))


def test_channel_annotations_cover_harmonized_channel_scalars() -> None:
    payload = package_payload()
    payload["sample"][0]["channel"] = [{
        "source": "lung",
        "annotations": [{
            "field": "tissue_name",
            "value": "lung",
            "term_source_ref": "UBERON",
            "term_accession_number": "UBERON:0002048",
        }],
    }]

    package = MINiMLCodec().decode(payload).package

    assert package.samples[0].channels[0].annotations[0].field == "tissue_name"
    assert package.to_mapping()["sample"][0]["channel"][0]["annotations"][0][
        "term_accession_number"
    ] == "UBERON:0002048"
