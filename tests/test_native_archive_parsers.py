"""Native metadata contracts using the independently vendored provider responses."""
from pathlib import Path
import copy
import json
import xml.etree.ElementTree as ET

import pytest

from meta_standards_converter.sources.archive_support import StudyRecords, StudySeed
from meta_standards_converter.miniml.sra_parser import SRAParser
from meta_standards_converter.miniml.ena_parser import ENAParser
from meta_standards_converter.miniml import MINiMLCodec

DOCS = Path(__file__).resolve().parents[1] / "docs"


def fixture_records(provider, experiment="SRX7812918"):
    study, project = (("SRP250911", "PRJNA609050") if experiment == "SRX7812918"
                      else ("SRP002056", "PRJNA123835"))
    records = StudyRecords(StudySeed(study, project if provider == "ena" else study))
    directory = DOCS / provider / "fixtures" / experiment
    if provider == "sra":
        records.xml.append(ET.parse(directory / "sra-efetch.xml").getroot())
        for kind in ("biosample", "bioproject"):
            records.xml.append(ET.parse(directory / f"{kind}-efetch.xml").getroot())
    else:
        records.xml.extend(ET.parse(directory / f"{kind}.xml").getroot()
                           for kind in ("study", "sample", "experiment", "run"))
        records.indexed["read_run"] = json.loads((directory / "file-report.json").read_text())
    return records


@pytest.mark.parametrize("provider,parser", [("sra", SRAParser), ("ena", ENAParser)])
def test_native_identity_biology_and_relationships(provider, parser):
    package = parser().parse(fixture_records(provider))
    data = package.to_mapping()
    assert data["series"]["iid"] == ("SRP250911" if provider == "sra" else "PRJNA609050")
    assert len(data["sample"]) == 1
    sample = data["sample"][0]
    assert sample["iid"] == ("SRS6225446" if provider == "sra" else "SAMN14218700")
    assert sample["channel"][0]["organism"][0]["taxid"] == "408170"
    attrs = sample["channel"][0]["characteristics"]
    assert any(a["name"] == "host" and a["value"] == "Homo sapiens" for a in attrs)
    assert sample["sra_run"][0]["run"] == "SRR11192680"
    assert sample["sra_run"][0]["library_layout"] == "PAIRED"
    assert data["series"]["assay_paths"]
    assert data["extensions"]["insdc"]["version"] == "1.0"
    assert data["extensions"]["insdc"]["records"]
    assert "completeness" not in json.dumps(data)
    MINiMLCodec().decode(data, strict=True)


def test_ena_parallel_files_preserve_duplicates_and_short_hash_columns():
    records = fixture_records("ena")
    records.indexed["read_run"] = [{"run_accession": "SRR11192680",
        "fastq_ftp": "ftp.sra.ebi.ac.uk/a.fastq.gz;ftp.sra.ebi.ac.uk/a.fastq.gz",
        "fastq_bytes": "10;20", "fastq_md5": "a" * 32}]
    files = ENAParser().parse(records).to_mapping()["sample"][0]["sra_run"][0]["fastq_files"]
    assert len(files) == 2
    assert files[0]["uri"] == files[1]["uri"]
    assert files[1].get("md5") is None
    assert files[1]["bytes"] == "20"


def test_repeating_attributes_units_and_missing_literals_survive():
    records = fixture_records("sra")
    attrs = records.xml[0].find(".//SAMPLE/SAMPLE_ATTRIBUTES")
    for value in ("missing: restricted access", "5"):
        attr = ET.SubElement(attrs, "SAMPLE_ATTRIBUTE")
        for tag, text in (("TAG", "dose"), ("VALUE", value), ("UNITS", "mg")):
            ET.SubElement(attr, tag).text = text
    values = SRAParser().parse(records).to_mapping()["sample"][0]["channel"][0]["characteristics"]
    assert [(a["value"], a["unit"]["value"]) for a in values if a["name"] == "dose"] == [
        ("missing: restricted access", "mg"), ("5", "mg")]


def test_multiple_libraries_do_not_create_samples_or_choose_scalar():
    records = fixture_records("sra")
    original = records.xml[0].find("EXPERIMENT_PACKAGE")
    second = copy.deepcopy(original)
    experiment = second.find("EXPERIMENT")
    experiment.set("accession", "SRX999999")
    experiment.find("IDENTIFIERS/PRIMARY_ID").text = "SRX999999"
    experiment.find("DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_STRATEGY").text = "WGS"
    run = second.find("RUN_SET/RUN")
    run.set("accession", "SRR999999")
    run.find("EXPERIMENT_REF").set("accession", "SRX999999")
    records.xml[0].append(second)
    sample = SRAParser().parse(records).to_mapping()["sample"][0]
    assert len(sample["sra_run"]) == 2
    assert "library_strategy" not in sample
    assert {r["library_strategy"] for r in sample["sra_run"]} == {"AMPLICON", "WGS"}


def test_sra_archives_are_not_fastq_and_statistics_remain_scoped():
    data = SRAParser().parse(fixture_records("sra")).to_mapping()
    run = data["sample"][0]["sra_run"][0]
    assert run["read_lengths"]
    assert all(not f["filename"].endswith("sralite") for f in run.get("fastq_files", []))
    assert any("sralite" in str(r) for r in data["extensions"]["insdc"]["records"])


def test_linked_biosample_attributes_and_project_publications_are_mapped():
    records = fixture_records('sra', 'SRX017289')
    biosample = records.xml[1].find('BioSample')
    attr = ET.SubElement(biosample.find('Attributes'), 'Attribute', attribute_name='tissue')
    attr.text = 'explicit tissue'
    records.xml.append(ET.fromstring('''<PubmedArticleSet><PubmedArticle><MedlineCitation><PMID>20133686</PMID><Article><ArticleTitle>Publication title</ArticleTitle><AuthorList><Author><LastName>Example</LastName><ForeName>Alice</ForeName></Author></AuthorList></Article></MedlineCitation><PubmedData><ArticleIdList><ArticleId IdType="doi">10.1/example</ArticleId></ArticleIdList></PubmedData></PubmedArticle></PubmedArticleSet>'''))
    data = SRAParser().parse(records).to_mapping()
    assert any(c['name'] == 'tissue' for c in data['sample'][0]['channel'][0]['characteristics'])
    assert data['series']['pubmed_publication'][0]['doi'] == '10.1/example'


def test_analysis_files_are_sample_scoped_only_with_explicit_references():
    records = fixture_records('ena')
    records.xml.append(ET.fromstring('''<ANALYSIS_SET><ANALYSIS accession="ERZ1"><STUDY_REF accession="SRP250911"/><SAMPLE_REF accession="SRS6225446"/><TITLE>assembly</TITLE><DESCRIPTION>Assembly method</DESCRIPTION><ANALYSIS_TYPE><SEQUENCE_ASSEMBLY><PROGRAM>assembler</PROGRAM><COVERAGE>30</COVERAGE></SEQUENCE_ASSEMBLY></ANALYSIS_TYPE><FILES><FILE filename="https://example.org/assembly.fasta" filetype="fasta" checksum_method="SHA256" checksum="abc"/></FILES></ANALYSIS></ANALYSIS_SET>'''))
    data = ENAParser().parse(records).to_mapping()
    link = data['sample'][0]['supplementary_data'][0]
    assert link['value'].endswith('assembly.fasta')
    assert 'checksum' not in link
    assert data['series']['protocols'][0]['software'] == ['assembler']
