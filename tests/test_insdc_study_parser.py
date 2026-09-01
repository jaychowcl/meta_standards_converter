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

from meta_standards_converter.insdc_handlers.study_models import (
    ProviderDocument,
    StudyFetchResult,
)
from meta_standards_converter.insdc_handlers.study_parser import INSDCStudyParser
from meta_standards_converter.miniml import MINiMLCodec


ROOT = Path(__file__).resolve().parents[1]


def provider_document(path: Path, *, kind: str, uri: str | None = None) -> ProviderDocument:
    return ProviderDocument(
        kind=kind,
        name=path.name,
        uri=uri or f"https://example.test/{path.name}",
        media_type="application/json" if path.suffix == ".json" else "application/xml",
        content=path.read_bytes(),
    )


def encoded(package):
    return MINiMLCodec().encode(package)


def warning_codes(package) -> list[str]:
    return [item["code"] for item in encoded(package)["extensions"]["insdc"]["warnings"]]


def test_sra_fixture_projects_source_faithful_miniml() -> None:
    fixture = ROOT / "docs" / "sra" / "fixtures" / "SRX017289"
    result = StudyFetchResult(
        provider="sra",
        requested_accession="SRX017289",
        study_accession="SRP002056",
        documents=(
            provider_document(fixture / "sra-efetch.xml", kind="sra_efetch"),
            provider_document(fixture / "biosample-efetch.xml", kind="biosample_efetch"),
            provider_document(fixture / "bioproject-efetch.xml", kind="bioproject_efetch"),
            provider_document(fixture / "pubmed-esummary.xml", kind="pubmed_esummary"),
        ),
    )

    package = INSDCStudyParser().parse(result)
    data = encoded(package)

    assert data["miniml_schema_version"] == "3.0"
    assert data["source"]["format"] == "NCBI SRA"
    assert {item["value"] for item in data["series"]["accession"]} >= {
        "SRP002056",
        "PRJNA123835",
        "GSE18729",
    }
    assert data["series"]["pubmed_id"] == ["20133686"]
    assert data["series"]["pubmed_publication"][0]["title"]
    assert data["sample"][0]["iid"] == "SRS011830"
    assert data["sample"][0]["channel"][0]["organism"][0]["taxid"] == "6239"
    assert {item["name"]: item["value"] for item in data["sample"][0]["channel"][0]["characteristics"]}["strain"] == "Wild type, fog-2(q71)"
    run = data["sample"][0]["sra_run"][0]
    assert run["run"] == "SRR037073"
    assert run["experiment"] == "SRX017289"
    assert run["library_strategy"] == "RNA-Seq"
    assert len(data["source"]["documents"]) == 4
    assert all(len(item["sha256"]) == 64 for item in data["source"]["documents"])
    assert "linked_geo_enrichment_disabled" in warning_codes(package)
    assert data["extensions"]["insdc"]["submission_origin"] == "geo_brokered"
    MINiMLCodec().decode(data, strict=True)


def test_ena_fixture_aligns_files_and_preserves_geography() -> None:
    fixture = ROOT / "docs" / "ena" / "fixtures" / "SRX7812918"
    result = StudyFetchResult(
        provider="ena",
        requested_accession="SRX7812918",
        study_accession="SRP250911",
        documents=(
            provider_document(fixture / "study.xml", kind="ena_study"),
            provider_document(fixture / "sample.xml", kind="ena_sample"),
            provider_document(fixture / "experiment.xml", kind="ena_experiment"),
            provider_document(fixture / "run.xml", kind="ena_run"),
            provider_document(fixture / "file-report.json", kind="ena_file_report"),
            provider_document(fixture / "taxonomy.json", kind="ena_taxonomy"),
            provider_document(fixture / "xref.json", kind="ena_xref"),
        ),
    )

    data = encoded(INSDCStudyParser().parse(result))

    sample = data["sample"][0]
    characteristics = sample["channel"][0]["characteristics"]
    assert [(item["name"], item["value"]) for item in characteristics if item["name"] in {"geo_loc_name", "lat_lon"}] == [
        ("geo_loc_name", "Canada: London"),
        ("lat_lon", "43.01 N 81.27 W"),
    ]
    files = sample["sra_run"][0]["fastq_files"]
    assert [item["filename"] for item in files] == [
        "SRR11192680_1.fastq.gz",
        "SRR11192680_2.fastq.gz",
    ]
    assert [item["bytes"] for item in files] == ["4216429", "4680355"]
    assert data["extensions"]["insdc"]["submission_origin"] == "ncbi_direct"
    assert data["series"].get("pubmed_id", []) == []


def test_duplicate_attributes_and_missing_terms_are_not_collapsed() -> None:
    xml = b"""<EXPERIMENT_PACKAGE_SET><EXPERIMENT_PACKAGE>
      <EXPERIMENT accession="SRX1"><STUDY_REF accession="SRP1"/><DESIGN>
        <SAMPLE_DESCRIPTOR accession="SRS1"/><LIBRARY_DESCRIPTOR>
          <LIBRARY_STRATEGY>OTHER</LIBRARY_STRATEGY><LIBRARY_SOURCE>OTHER</LIBRARY_SOURCE>
          <LIBRARY_SELECTION>other</LIBRARY_SELECTION><LIBRARY_LAYOUT><SINGLE/></LIBRARY_LAYOUT>
        </LIBRARY_DESCRIPTOR></DESIGN><PLATFORM><ILLUMINA><INSTRUMENT_MODEL>Illumina NovaSeq 6000</INSTRUMENT_MODEL></ILLUMINA></PLATFORM></EXPERIMENT>
      <STUDY accession="SRP1"><DESCRIPTOR><STUDY_TITLE>Study</STUDY_TITLE></DESCRIPTOR></STUDY>
      <SAMPLE accession="SRS1"><TITLE>Sample</TITLE><SAMPLE_NAME><TAXON_ID>9606</TAXON_ID><SCIENTIFIC_NAME>Homo sapiens</SCIENTIFIC_NAME></SAMPLE_NAME>
        <SAMPLE_ATTRIBUTES>
          <SAMPLE_ATTRIBUTE><TAG>tissue</TAG><VALUE>blood</VALUE></SAMPLE_ATTRIBUTE>
          <SAMPLE_ATTRIBUTE><TAG>tissue</TAG><VALUE>plasma</VALUE></SAMPLE_ATTRIBUTE>
          <SAMPLE_ATTRIBUTE><TAG>collection_date</TAG><VALUE>missing: not collected</VALUE></SAMPLE_ATTRIBUTE>
        </SAMPLE_ATTRIBUTES></SAMPLE>
      <RUN_SET><RUN accession="SRR1"><EXPERIMENT_REF accession="SRX1"/></RUN></RUN_SET>
    </EXPERIMENT_PACKAGE></EXPERIMENT_PACKAGE_SET>"""
    result = StudyFetchResult(
        provider="sra",
        requested_accession="SRP1",
        study_accession="SRP1",
        documents=(ProviderDocument("sra_efetch", "study.xml", "https://example.test/study.xml", "application/xml", xml),),
    )

    data = encoded(INSDCStudyParser().parse(result))
    values = [(item["name"], item["value"]) for item in data["sample"][0]["channel"][0]["characteristics"]]

    assert values == [
        ("tissue", "blood"),
        ("tissue", "plasma"),
        ("collection_date", "missing: not collected"),
    ]


def test_mixed_experiment_values_remain_run_local_and_warn() -> None:
    xml = b"""<EXPERIMENT_PACKAGE_SET>
      <EXPERIMENT_PACKAGE><EXPERIMENT accession="DRX1"><STUDY_REF accession="DRP1"/><DESIGN><SAMPLE_DESCRIPTOR accession="DRS1"/><LIBRARY_DESCRIPTOR><LIBRARY_STRATEGY>OTHER</LIBRARY_STRATEGY><LIBRARY_SOURCE>GENOMIC</LIBRARY_SOURCE><LIBRARY_SELECTION>other</LIBRARY_SELECTION><LIBRARY_LAYOUT><PAIRED/></LIBRARY_LAYOUT></LIBRARY_DESCRIPTOR></DESIGN></EXPERIMENT><STUDY accession="DRP1"><DESCRIPTOR><STUDY_TITLE>Study</STUDY_TITLE></DESCRIPTOR></STUDY><SAMPLE accession="DRS1"><TITLE>Sample</TITLE></SAMPLE><RUN_SET><RUN accession="DRR1"><EXPERIMENT_REF accession="DRX1"/></RUN></RUN_SET></EXPERIMENT_PACKAGE>
      <EXPERIMENT_PACKAGE><EXPERIMENT accession="DRX2"><STUDY_REF accession="DRP1"/><DESIGN><SAMPLE_DESCRIPTOR accession="DRS1"/><LIBRARY_DESCRIPTOR><LIBRARY_STRATEGY>OTHER</LIBRARY_STRATEGY><LIBRARY_SOURCE>GENOMIC</LIBRARY_SOURCE><LIBRARY_SELECTION>other</LIBRARY_SELECTION><LIBRARY_LAYOUT><SINGLE/></LIBRARY_LAYOUT></LIBRARY_DESCRIPTOR></DESIGN></EXPERIMENT><STUDY accession="DRP1"><DESCRIPTOR><STUDY_TITLE>Study</STUDY_TITLE></DESCRIPTOR></STUDY><SAMPLE accession="DRS1"><TITLE>Sample</TITLE></SAMPLE><RUN_SET><RUN accession="DRR2"><EXPERIMENT_REF accession="DRX2"/></RUN></RUN_SET></EXPERIMENT_PACKAGE>
    </EXPERIMENT_PACKAGE_SET>"""
    result = StudyFetchResult(
        provider="sra",
        requested_accession="DRP1",
        study_accession="DRP1",
        documents=(ProviderDocument("sra_efetch", "study.xml", None, "application/xml", xml),),
    )

    package = INSDCStudyParser().parse(result)
    sample = encoded(package)["sample"][0]

    assert "library_layout" not in sample
    assert [item["library_layout"] for item in sample["sra_run"]] == ["PAIRED", "SINGLE"]
    assert "non_unanimous_experiment_value" in warning_codes(package)


def test_valid_ena_study_without_hierarchy_is_degraded_not_rejected() -> None:
    xml = b"""<STUDY_SET><STUDY accession="DRP000158" broker_name="DDBJ"><IDENTIFIERS><PRIMARY_ID>DRP000158</PRIMARY_ID><SECONDARY_ID>PRJDA39855</SECONDARY_ID></IDENTIFIERS><DESCRIPTOR><STUDY_TITLE>Wild rice Genome</STUDY_TITLE><STUDY_ABSTRACT>Sequencing of wild rice.</STUDY_ABSTRACT></DESCRIPTOR></STUDY></STUDY_SET>"""
    result = StudyFetchResult(
        provider="ena",
        requested_accession="DRP000158",
        study_accession="DRP000158",
        documents=(ProviderDocument("ena_study", "DRP000158.xml", None, "application/xml", xml),),
    )

    package = INSDCStudyParser().parse(result)
    data = encoded(package)

    assert data["sample"] == []
    assert data["series"]["iid"] == "DRP000158"
    assert "metadata_only_study" in warning_codes(package)


def test_provider_document_digest_is_over_downloaded_bytes() -> None:
    document = ProviderDocument("test", "value.json", None, "application/json", json.dumps({"a": 1}).encode())
    assert document.sha256 == "f9d86028c6e0d64e225186f96acb69338b2c59764df79162107f5c4bb34d1310"


def test_same_accession_has_common_typed_core_and_provider_specific_provenance() -> None:
    sra_fixture = ROOT / "docs" / "sra" / "fixtures" / "SRX017289"
    ena_fixture = ROOT / "docs" / "ena" / "fixtures" / "SRX017289"
    sra = StudyFetchResult(
        "sra",
        "SRX017289",
        "SRP002056",
        (provider_document(sra_fixture / "sra-efetch.xml", kind="sra_efetch"),),
    )
    ena = StudyFetchResult(
        "ena",
        "SRX017289",
        "SRP002056",
        tuple(
            provider_document(ena_fixture / name, kind=kind)
            for kind, name in (
                ("ena_study", "study.xml"),
                ("ena_sample", "sample.xml"),
                ("ena_experiment", "experiment.xml"),
                ("ena_run", "run.xml"),
                ("ena_file_report", "file-report.json"),
            )
        ),
    )

    sra_data = encoded(INSDCStudyParser().parse(sra))
    ena_data = encoded(INSDCStudyParser().parse(ena))

    for data in (sra_data, ena_data):
        assert data["series"]["iid"] == "SRP002056"
        assert data["sample"][0]["iid"] == "SRS011830"
        assert data["sample"][0]["sra_run"][0]["experiment"] == "SRX017289"
        assert data["sample"][0]["sra_run"][0]["run"] == "SRR037073"
        assert data["sample"][0]["sra_run"][0]["library_strategy"] == "RNA-Seq"
    assert sra_data["extensions"]["insdc"]["provider"] == "sra"
    assert ena_data["extensions"]["insdc"]["provider"] == "ena"
    assert ena_data["sample"][0]["sra_run"][0]["fastq_files"]
    assert "fastq_files" not in sra_data["sample"][0]["sra_run"][0]


def test_shorter_ena_file_columns_do_not_shift_later_files() -> None:
    fixture = ROOT / "docs" / "ena" / "fixtures" / "SRX7812918"
    file_report = json.dumps(
        [
            {
                "run_accession": "SRR11192680",
                "fastq_ftp": "host/one.fastq.gz;host/two.fastq.gz;host/three.fastq.gz",
                "fastq_md5": "first-md5",
                "fastq_bytes": "10;20",
            }
        ]
    ).encode()
    documents = [
        provider_document(fixture / name, kind=kind)
        for kind, name in (
            ("ena_study", "study.xml"),
            ("ena_sample", "sample.xml"),
            ("ena_experiment", "experiment.xml"),
            ("ena_run", "run.xml"),
        )
    ]
    documents.append(
        ProviderDocument(
            "ena_file_report",
            "file-report.json",
            None,
            "application/json",
            file_report,
        )
    )
    result = StudyFetchResult(
        "ena", "SRX7812918", "SRP250911", tuple(documents)
    )

    files = encoded(INSDCStudyParser().parse(result))["sample"][0]["sra_run"][0]["fastq_files"]

    assert files == [
        {"uri": "ftp://host/one.fastq.gz", "filename": "one.fastq.gz", "md5": "first-md5", "bytes": "10"},
        {"uri": "ftp://host/two.fastq.gz", "filename": "two.fastq.gz", "bytes": "20"},
        {"uri": "ftp://host/three.fastq.gz", "filename": "three.fastq.gz"},
    ]


def test_unmapped_archive_attributes_and_target_provenance_are_ordered() -> None:
    fixture = ROOT / "docs" / "sra" / "fixtures" / "SRX017289"
    result = StudyFetchResult(
        "sra",
        "SRX017289",
        "SRP002056",
        (provider_document(fixture / "sra-efetch.xml", kind="sra_efetch"),),
    )

    extension = encoded(INSDCStudyParser().parse(result))["extensions"]["insdc"]

    assert extension["provenance"][0]["target"] == "sample/SRS011830"
    assert extension["provenance"][-1]["target"] == "series"
    assert any(
        item["source_path"].endswith("parent_bioproject")
        and item["value"] == "PRJNA121635"
        for item in extension["unmapped"]
    )


def test_analysis_and_assembly_documents_remain_ordered_unmapped_records() -> None:
    study = b'<STUDY_SET><STUDY accession="ERP1"><DESCRIPTOR><STUDY_TITLE>Study</STUDY_TITLE></DESCRIPTOR></STUDY></STUDY_SET>'
    analysis = b'<ANALYSIS_SET><ANALYSIS accession="ERZ1"><TITLE>Analysis title</TITLE><ANALYSIS_TYPE><SEQUENCE_ANNOTATION/></ANALYSIS_TYPE></ANALYSIS></ANALYSIS_SET>'
    assembly = b'<ASSEMBLY_SET><ASSEMBLY accession="GCA_1"><NAME>Assembly name</NAME><ASSEMBLY_LEVEL>chromosome</ASSEMBLY_LEVEL></ASSEMBLY></ASSEMBLY_SET>'
    result = StudyFetchResult(
        "ena",
        "ERP1",
        "ERP1",
        (
            ProviderDocument("ena_study", "study.xml", None, "application/xml", study),
            ProviderDocument("ena_analysis", "analysis.xml", None, "application/xml", analysis),
            ProviderDocument("ena_assembly", "assembly.xml", None, "application/xml", assembly),
        ),
    )

    unmapped = encoded(INSDCStudyParser().parse(result))["extensions"]["insdc"]["unmapped"]

    assert any(item["source_path"] == "ANALYSIS[ERZ1]/TITLE" and item["value"] == "Analysis title" for item in unmapped)
    assert any(item["source_path"] == "ASSEMBLY[GCA_1]/NAME" and item["value"] == "Assembly name" for item in unmapped)


def test_sra_organizations_are_deduplicated_and_contacts_are_typed() -> None:
    fixture = ROOT / "docs" / "sra" / "fixtures" / "SRX7812918"
    content = (fixture / "sra-efetch.xml").read_bytes()
    result = StudyFetchResult(
        "sra",
        "SRX7812918",
        "SRP250911",
        (
            ProviderDocument("sra_efetch", "resolve.xml", None, "application/xml", content),
            ProviderDocument("sra_efetch", "study.xml", None, "application/xml", content),
        ),
    )

    data = encoded(INSDCStudyParser().parse(result))

    assert data["organization"] == [
        {
            "iid": "organization-1",
            "name": "University of Western Ontario",
            "address": {
                "line": ["Microbiology and Immunology", "University of Western Ontario", "1151 Richmond Street"],
                "city": "London",
                "province": "Ontario",
                "postal_code": "N6A3K7",
                "country": "Canada",
            },
        }
    ]
    assert data["contributor"][0]["person"] == {"first": "Brendan", "last": "Daisley"}
    assert data["contributor"][0]["email"] == "bdaisley@uwo.ca"
    assert data["series"]["contact_ref"] == [{"ref": "contact-1"}]
