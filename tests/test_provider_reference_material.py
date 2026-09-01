# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
import csv
import hashlib
import json
import xml.etree.ElementTree as ET
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
SRA = DOCS / "sra"
ENA = DOCS / "ena"

SRA_SCHEMAS = {
    "SRA.analysis.xsd",
    "SRA.common.xsd",
    "SRA.experiment.xsd",
    "SRA.receipt.xsd",
    "SRA.run.xsd",
    "SRA.sample.xsd",
    "SRA.study.xsd",
    "SRA.submission.xsd",
}
ENA_SCHEMAS = {
    "ENA.assembly.xsd",
    "ENA.checklist.xsd",
    "ENA.embl.xsd",
    "ENA.project.xsd",
    "ENA.root.xsd",
    "ENA.sample_group.xsd",
    "ENA.taxonomy.xsd",
    "ENA.webin.xsd",
}
ENA_RESULT_TYPES = {
    "study",
    "read_study",
    "sample",
    "read_experiment",
    "read_run",
    "analysis",
    "analysis_study",
    "assembly",
    "taxon",
}
ENA_CONTROLLED_FIELDS = {
    "analysis_type",
    "assembly_level",
    "broker_name",
    "category",
    "checklist",
    "datahub",
    "genome_representation",
    "host_sex",
    "instrument_model",
    "instrument_platform",
    "library_layout",
    "library_selection",
    "library_source",
    "library_strategy",
    "ncbi_reporting_standard",
    "sex",
    "tag",
    "tax_division",
}
FIXTURES = {
    "SRX017289": {
        "study": "SRP002056",
        "sample": "SRS011830",
        "experiment": "SRX017289",
        "run": "SRR037073",
        "bioproject": "PRJNA123835",
        "biosample": "SAMN00009557",
        "pubmed": "20133686",
        "organism": "Caenorhabditis elegans",
        "strategy": "RNA-Seq",
        "source": "TRANSCRIPTOMIC",
        "selection": "size fractionation",
        "layout": "SINGLE",
        "instrument": "Illumina Genome Analyzer II",
    },
    "SRX7812918": {
        "study": "SRP250911",
        "sample": "SRS6225446",
        "experiment": "SRX7812918",
        "run": "SRR11192680",
        "bioproject": "PRJNA609050",
        "biosample": "SAMN14218700",
        "pubmed": None,
        "organism": "human gut metagenome",
        "strategy": "AMPLICON",
        "source": "GENOMIC",
        "selection": "PCR",
        "layout": "PAIRED",
        "instrument": "Illumina MiSeq",
    },
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _manifest(root: Path) -> dict:
    return json.loads((root / "manifest.json").read_text(encoding="utf-8"))


def _tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _find_external_id(node: ET.Element, namespace: str) -> str | None:
    return next(
        (
            item.text
            for item in node.findall(".//IDENTIFIERS/EXTERNAL_ID")
            if (item.get("namespace") or "").lower() == namespace.lower()
        ),
        None,
    )


def _sra_package(fixture: str) -> ET.Element:
    root = ET.parse(SRA / "fixtures" / fixture / "sra-efetch.xml").getroot()
    package = root.find(".//EXPERIMENT_PACKAGE")
    assert package is not None
    return package


def test_required_schema_snapshots_are_present_and_well_formed():
    assert {path.name for path in (SRA / "schemas").glob("SRA.*.xsd")} == SRA_SCHEMAS
    assert {path.name for path in (ENA / "schemas").glob("ENA.*.xsd")} == ENA_SCHEMAS
    assert (SRA / "schemas" / "ncbi" / "biosample.xsd").is_file()
    assert (SRA / "schemas" / "ncbi" / "bioproject.xsd").is_file()

    for root in (SRA / "schemas", ENA / "schemas"):
        for path in root.rglob("*.xsd"):
            parsed = ET.parse(path).getroot()
            assert parsed.tag == "{http://www.w3.org/2001/XMLSchema}schema"


def test_ena_openapi_snapshots_cover_the_expected_services_and_paths():
    expected = {
        "portal-api.openapi.json": "/search",
        "browser-api.openapi.json": "/xml/{accession}",
        "xref-api.openapi.json": "/json/search",
        "taxonomy-api.openapi.json": "/tax-id/**",
    }
    for filename, path in expected.items():
        document = json.loads((ENA / "openapi" / filename).read_text(encoding="utf-8"))
        assert document["openapi"].startswith("3.")
        assert path in document["paths"]


def test_ena_field_and_value_catalogues_are_complete():
    results = _tsv_rows(ENA / "catalogues" / "results.tsv")
    assert ENA_RESULT_TYPES <= {row["resultId"] for row in results}

    for result in ENA_RESULT_TYPES:
        for folder in ("return-fields", "search-fields"):
            rows = _tsv_rows(ENA / "catalogues" / folder / f"{result}.tsv")
            assert rows
            assert {"columnId", "description", "type"} <= set(rows[0])

    controlled_paths = ENA / "catalogues" / "controlled-vocab"
    assert {path.stem for path in controlled_paths.glob("*.json")} == ENA_CONTROLLED_FIELDS
    for path in controlled_paths.glob("*.json"):
        values = json.loads(path.read_text(encoding="utf-8"))
        assert isinstance(values, list)
        assert all({"value", "description"} <= set(item) for item in values)

    assert {item["value"] for item in json.loads(
        (controlled_paths / "library_layout.json").read_text(encoding="utf-8")
    )} == {"SINGLE", "PAIRED"}

    checklist_catalogue = {
        item["value"]
        for item in json.loads(
            (controlled_paths / "checklist.json").read_text(encoding="utf-8")
        )
    }
    checklist_paths = list((ENA / "catalogues" / "checklists").glob("ERC*.xml"))
    checklist_xml_ids = {path.stem for path in checklist_paths}
    unavailable_xml_ids = {
        "ERC000004",
        "ERC000042",
        "ERC000043",
        "ERC000044",
        "ERC000045",
        "ERC000047",
        "ERC000048",
        "ERC000049",
        "ERC000050",
        "ERC000051",
        "ERC000052",
        "ERC000053",
        "ERC000055",
        "ERC000056",
        "ERC000058",
        "ERC100001",
    }
    assert len(checklist_catalogue) == 47
    assert len(checklist_paths) == 31
    assert checklist_xml_ids | unavailable_xml_ids == checklist_catalogue
    assert checklist_xml_ids.isdisjoint(unavailable_xml_ids)
    for path in checklist_paths:
        root = ET.parse(path).getroot()
        assert root.find(".//CHECKLIST") is not None or root.tag == "CHECKLIST"

    availability = (ENA / "checklist-availability.md").read_text(encoding="utf-8")
    for accession in unavailable_xml_ids:
        assert accession in availability


def test_ncbi_catalogues_are_present_and_well_formed():
    for filename in ("einfo-sra.xml", "biosample-packages.xml", "biosample-attributes.xml"):
        ET.parse(SRA / "catalogues" / filename)


def test_manifests_cover_every_vendored_artifact_and_match_hashes():
    for root, provider in ((SRA, "sra"), (ENA, "ena")):
        manifest = _manifest(root)
        assert manifest["snapshot_schema_version"] == "1.0"
        assert manifest["provider"] == provider
        assert manifest["retrieved_at"]

        entries = manifest["artifacts"]
        paths = {entry["path"] for entry in entries}
        expected_paths = {
            path.relative_to(root).as_posix()
            for path in root.rglob("*")
            if path.is_file() and path.suffix != ".md" and path.name != "manifest.json"
        }
        assert paths == expected_paths

        for entry in entries:
            path = root / entry["path"]
            assert entry["url"].startswith("https://")
            assert entry["content_type"]
            assert entry["sha256"] == _sha256(path)
            assert "etag" in entry
            assert "last_modified" in entry
            assert "fixture_accession" in entry


def test_fixture_chains_preserve_expected_sra_metadata():
    for fixture, expected in FIXTURES.items():
        package = _sra_package(fixture)
        study = package.find("./STUDY")
        sample = package.find("./SAMPLE")
        experiment = package.find("./EXPERIMENT")
        run = package.find(".//RUN_SET/RUN")
        assert study is not None and sample is not None and experiment is not None and run is not None

        assert study.get("accession") == expected["study"]
        assert sample.get("accession") == expected["sample"]
        assert experiment.get("accession") == expected["experiment"]
        assert run.get("accession") == expected["run"]
        assert _find_external_id(study, "BioProject") == expected["bioproject"]
        assert _find_external_id(sample, "BioSample") == expected["biosample"]
        assert sample.findtext("./SAMPLE_NAME/SCIENTIFIC_NAME") == expected["organism"]
        assert experiment.findtext(".//LIBRARY_STRATEGY") == expected["strategy"]
        assert experiment.findtext(".//LIBRARY_SOURCE") == expected["source"]
        assert experiment.findtext(".//LIBRARY_SELECTION") == expected["selection"]
        layout = experiment.find(".//LIBRARY_LAYOUT")
        assert layout is not None and list(layout)[0].tag == expected["layout"]
        assert experiment.findtext(".//INSTRUMENT_MODEL") == expected["instrument"]

        pubmed_ids = [
            link.findtext("ID")
            for link in study.findall(".//XREF_LINK")
            if (link.findtext("DB") or "").lower() == "pubmed"
        ]
        if expected["pubmed"]:
            assert expected["pubmed"] in pubmed_ids
            ET.parse(SRA / "fixtures" / fixture / "pubmed-esummary.xml")
        else:
            assert not pubmed_ids
            assert not (SRA / "fixtures" / fixture / "pubmed-esummary.xml").exists()

        ET.parse(SRA / "fixtures" / fixture / "biosample-efetch.xml")
        ET.parse(SRA / "fixtures" / fixture / "bioproject-efetch.xml")


def test_fixture_chains_preserve_expected_ena_metadata_and_reports():
    for fixture, expected in FIXTURES.items():
        fixture_root = ENA / "fixtures" / fixture
        study = ET.parse(fixture_root / "study.xml").getroot().find(".//STUDY")
        sample = ET.parse(fixture_root / "sample.xml").getroot().find(".//SAMPLE")
        experiment = ET.parse(fixture_root / "experiment.xml").getroot().find(".//EXPERIMENT")
        run = ET.parse(fixture_root / "run.xml").getroot().find(".//RUN")
        assert study is not None and sample is not None and experiment is not None and run is not None

        assert study.get("accession") == expected["study"]
        assert sample.get("accession") == expected["sample"]
        assert experiment.get("accession") == expected["experiment"]
        assert run.get("accession") == expected["run"]
        assert sample.findtext("./SAMPLE_NAME/SCIENTIFIC_NAME") == expected["organism"]
        assert experiment.findtext(".//LIBRARY_STRATEGY") == expected["strategy"]
        assert experiment.findtext(".//LIBRARY_SOURCE") == expected["source"]
        assert experiment.findtext(".//LIBRARY_SELECTION") == expected["selection"]
        layout = experiment.find(".//LIBRARY_LAYOUT")
        assert layout is not None and list(layout)[0].tag == expected["layout"]
        assert experiment.findtext(".//INSTRUMENT_MODEL") == expected["instrument"]

        report = json.loads((fixture_root / "file-report.json").read_text(encoding="utf-8"))
        assert report and report[0]["run_accession"] == expected["run"]
        assert report[0].get("fastq_ftp")
        for name in ("portal-study.json", "portal-sample.json", "portal-read_experiment.json", "portal-read_run.json", "xref.json", "taxonomy.json"):
            json.loads((fixture_root / name).read_text(encoding="utf-8"))


def test_provider_reference_docs_are_routed_to_implemented_converters():
    for provider in ("sra", "ena"):
        assert (DOCS / provider / "README.md").is_file()
        assert (DOCS / provider / "expected-fields.md").is_file()

    codebase = (DOCS / "codebase.md").read_text(encoding="utf-8")
    index = (DOCS / "index.md").read_text(encoding="utf-8")
    assert '<a id="insdc-provider-reference-material"></a>' in codebase
    assert "codebase.md#insdc-provider-reference-material" in index
    assert "`SRAStudyFetcher` and `ENAStudyFetcher` now retrieve complete study graphs" in codebase
    assert "codebase.md#workflow-insdc2json" in index
