# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Independent checks linking reviewed expectations to provider/source evidence."""
import csv
import hashlib
import json
import xml.etree.ElementTree as ET
from pathlib import Path
from importlib.metadata import version

import numpy as np
import pandas as pd
from tests.support.contracts import GEO, PBMC, AE


def test_geo_expectation_preserves_provider_facts():
    source = ET.parse(GEO / "inputs/geo.xml").getroot()
    ns = {"m": source.tag.split("}")[0][1:]}
    sample = source.find("m:Sample", ns)
    series = source.find("m:Series", ns)
    package = json.loads((GEO / "expected/geo2json/GSE328265.json").read_text())[0]
    assert len(package["sample"]) == 1
    expected_sample = package["sample"][0]
    assert expected_sample["iid"] == sample.attrib["iid"] == "GSM9651991"
    assert expected_sample["title"] == sample.findtext("m:Title", namespaces=ns)
    assert package["series"]["title"] == series.findtext("m:Title", namespaces=ns)
    assert package["series"]["sample_ref"] == [{"ref": "GSM9651991"}]
    channel = expected_sample["channel"][0]
    assert channel["organism"] == [{"value": "Homo sapiens", "taxid": "9606"}]
    assert channel["source"]["value"] == sample.findtext("m:Channel/m:Source", namespaces=ns)
    assert channel["molecule"]["value"] == sample.findtext("m:Channel/m:Molecule", namespaces=ns)
    assert channel["extract_protocol"] == " ".join(sample.findtext("m:Channel/m:Extract-Protocol", namespaces=ns).split())
    assert expected_sample["relation"] == [dict(n.attrib) for n in sample.findall("m:Relation", ns)]
    pubmed = ET.parse(GEO / "provider_responses/pubmed.xml")
    publication = package["series"]["pubmed_publication"][0]
    assert publication["doi"] == pubmed.findtext(".//Item[@Name='DOI']")
    assert publication["title"] == pubmed.findtext(".//Item[@Name='Title']")
    assert publication["author_list"] == ", ".join(n.text for n in pubmed.findall(".//Item[@Name='AuthorList']/Item"))
    sra = ET.parse(GEO / "provider_responses/sra.xml")
    run = expected_sample["sra_run"][0]
    assert run["run"] == sra.find(".//RUN").attrib["accession"]
    assert run["experiment"] == sra.find(".//EXPERIMENT").attrib["accession"]
    assert run["instrument_model"] == sra.findtext(".//INSTRUMENT_MODEL")
    assert run["library_strategy"] == sra.findtext(".//LIBRARY_STRATEGY")
    assert "fastq_files" not in run  # Frozen ENA report has no available FASTQ URLs.


def test_magetab_expectations_preserve_protocol_bindings_and_evidence():
    for converter in ("geo2ae", "json2ae"):
        expected = GEO / "expected" / converter
        with (expected / "E-GEOD-328265.idf.txt").open() as handle:
            idf = {row[0]: row[1:] for row in csv.reader(handle, delimiter="\t")}
        with (expected / "E-GEOD-328265.sdrf.txt").open() as handle:
            rows = list(csv.reader(handle, delimiter="\t"))
        assert len(rows) == 2
        header, sample = rows
        assert len(header) == len(sample)
        assert sample[header.index("Source Name")] == "GSM9651991"
        assert sample[header.index("Comment[ENA_RUN]")] == "SRR37976431"
        assert sample[header.index("Characteristics[tissue]")] == "CSF"
        assert sample[header.index("Characteristics[organism]")] == "Homo sapiens"
        for index, column in enumerate(header):
            if column == "Protocol REF": assert sample[index] in idf["Protocol Name"]
        assert idf["PubMed ID"] == ["42129775"]
        assert idf["Publication DOI"] == ["10.1186/s12974-026-03861-9"]
        assert idf["Public Release Date"] == ["2026-05-17"]
        assert idf["Protocol Type"] == ["normalization data transformation protocol", "nucleic acid extraction protocol", "nucleic acid library construction protocol", "sample collection protocol", "nucleic acid sequencing protocol"]


def test_pbmc_expectation_preserves_actual_counts_and_identity():
    source = pd.read_csv(PBMC / "inputs/counts.tsv", sep="\t", index_col=0)
    provenance = json.loads((PBMC / "provenance.json").read_text())
    expected = json.loads((PBMC / "expected/json2h5ad/PBMC3K_SAMPLE.h5ad.json").read_text())
    counts = [[0, 0, 0, 0, 1, 0], [0, 2, 0, 1, 0, 1], [1, 0, 1, 0, 1, 0], [9, 0, 0, 1, 0, 3]]
    assert source.T.values.tolist() == provenance["expected_counts_observations_by_genes"] == counts
    assert expected["X"] == {"kind": "csr", "dtype": "float64", "shape": [4, 6], "values": counts}
    assert expected["obs"]["index"]["values"] == [x + "-PBMC3K_SAMPLE" for x in source.columns]
    assert expected["var"]["index"]["values"] == list(source.index)
    retained = expected["uns"]["msc_miniml"]["packages_json"]
    source_package = json.loads((PBMC / "inputs/miniml.json").read_text())
    assert json.loads(retained) == [{**source_package, "database": [], "organization": [], "contributor": [], "platform": []}]
    assert expected["uns"]["meta_standards_converter"]["converter_version"] == version("meta-standards-converter") == "6.0.0"


def test_fixture_checksum_inventory():
    root = Path(__file__).resolve().parents[1] / "fixtures"
    inventory = json.loads((root / "contracts.sha256.json").read_text())
    for path, digest in inventory.items():
        assert hashlib.sha256((root / path).read_bytes()).hexdigest() == digest, path


def test_public_magetab_expectation_preserves_sample_grouping_and_fastqs():
    with (AE / "inputs/E-MTAB-6486.sdrf.txt").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    package = json.loads((AE / "expected/ae2json/E-MTAB-6486.json").read_text())[0]
    assert len(rows) == 12
    assert [s["iid"] for s in package["sample"]] == list(dict.fromkeys(r["Source Name"] for r in rows))
    assert len(package["sample"]) == 6
    for sample in package["sample"]:
        source_rows = [r for r in rows if r["Source Name"] == sample["iid"]]
        assert sample["library_layout"] == source_rows[0]["Comment[LIBRARY_LAYOUT]"]
        assert sample["channel"][0]["organism"][0]["value"] == source_rows[0]["Characteristics[organism]"]
        assert len(sample["sra_run"]) == 1
        assert sample["sra_run"][0]["run"] == source_rows[0]["Comment[ENA_RUN]"]
        assert [f["uri"] for f in sample["sra_run"][0]["fastq_files"]] == [r["Comment[FASTQ_URI]"] for r in source_rows]
