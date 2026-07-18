# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
import json
import copy
import os
import sys
import tempfile
import unittest
from unittest.mock import MagicMock


ROOT = os.path.dirname(os.path.dirname(__file__))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from meta_standards_converter.ae_handlers.ae_constructor import AEConstructor  # noqa: E402
from meta_standards_converter.ae_handlers.ae_idf_handlers import IDFConstructor  # noqa: E402
from meta_standards_converter.ae_handlers.ae_sdrf_handlers import SDRFConstructor  # noqa: E402
from meta_standards_converter.ae_handlers.ae_webfetcher import (  # noqa: E402
    MAGETabInput,
    TextResource,
)
from meta_standards_converter.converters.ae2json import ae2json  # noqa: E402


IDF = """MAGE-TAB Version\t1.1
Investigation Title\tExample study
Investigation Accession\tE-MTAB-1
Comment[SecondaryAccession]\tGSE123
Comment[SecondaryAccessionTermSourceRef]\tGEO
Experimental Design\tRNA-seq
Experimental Factor Name\tdisease
Experimental Factor Type\tdisease
Person Last Name\tDoe
Person First Name\tJane
Person Email\tjane@example.org
Person Affiliation\tExample Institute
Public Release Date\t2025-01-02
PubMed ID\t12345
Publication DOI\t10.1/example
Publication Author List\tDoe J
Publication Title\tExample paper
Publication Status\tpublished
Status Term Source Ref\tEFO
Status Term Accession Number\tEFO:0000001
Protocol Name\tP-collect\tP-extract
Protocol Type\tsample collection protocol\tnucleic acid extraction protocol
Protocol Description\tCollect samples\tExtract material
Protocol Type Term Source REF\tEFO\tEFO
Protocol Type Term Accession Number\tEFO:0005518\tEFO:0002944
SDRF File\tstudy.sdrf.txt
Term Source Name\tEFO
Term Source File\thttps://www.ebi.ac.uk/efo/
Mystery Row\tkeep me
"""

SDRF_HEADER = [
    "Source Name",
    "Comment[Sample_title]",
    "Characteristics[organism]",
    "Term Source REF",
    "Characteristics[disease]",
    "Factor Value[disease]",
    "Protocol REF",
    "Extract Name",
    "Material Type",
    "Array Design REF",
    "Technology Type",
    "Comment[LIBRARY_SOURCE]",
    "Comment[LIBRARY_STRATEGY]",
    "Comment[ENA_RUN]",
    "Comment[FASTQ_URI]",
    "Comment[MD5]",
    "Mystery Column",
]


def sdrf(rows=None):
    rows = rows or [
        ["GSM1", "Sample one", "Homo sapiens", "EFO", "case", "case", "P-extract", "GSM1", "RNA", "A-TEST-1", "sequencing assay", "TRANSCRIPTOMIC", "RNA-SEQ", "ERR1", "https://example/1.fastq.gz", "aaa", "x"],
        ["GSM1", "Sample one", "Homo sapiens", "EFO", "case", "case", "P-extract", "GSM1", "RNA", "A-TEST-1", "sequencing assay", "TRANSCRIPTOMIC", "RNA-SEQ", "ERR1", "https://example/2.fastq.gz", "bbb", "y"],
    ]
    return "\n".join("\t".join(row) for row in [SDRF_HEADER, *rows]) + "\n"


def resolved_input(idf=IDF, sdrfs=None):
    return MAGETabInput(
        idf=TextResource("study.idf.txt", idf, "memory:idf"),
        sdrfs=tuple(
            TextResource(f"study{index}.sdrf.txt", text, f"memory:sdrf:{index}")
            for index, text in enumerate(sdrfs or [sdrf()], start=1)
        ),
        source="E-MTAB-1",
        source_kind="accession",
    )


class TestAE2JSONConverter(unittest.TestCase):
    def test_maps_idf_and_sdrf_to_miniml_compatible_package(self):
        fetcher = MagicMock()
        fetcher.resolve.return_value = resolved_input()

        packages = ae2json(fetcher=fetcher).convert("E-MTAB-1")

        self.assertEqual(1, len(packages))
        package = packages[0]
        self.assertIsNone(package["version"])
        self.assertEqual("Example study", package["series"]["title"])
        self.assertEqual(
            ["GSE123", "E-MTAB-1"],
            [item["value"] for item in package["series"]["accession"]],
        )
        self.assertEqual([{"factor": "disease", "type": "disease"}], package["series"]["variable"])
        self.assertEqual("12345", package["series"]["pubmed_publication"][0]["pubmed_id"])
        self.assertEqual("Doe", package["contributor"][0]["person"]["last"])
        self.assertEqual("EFO", package["database"][0]["iid"])

        sample = package["sample"][0]
        self.assertEqual("GSM1", sample["iid"])
        self.assertEqual("Sample one", sample["title"])
        self.assertEqual({"ref": "A-TEST-1"}, sample["platform_ref"])
        self.assertEqual("TRANSCRIPTOMIC", sample["library_source"])
        self.assertEqual("Extract material", sample["channel"][0]["extract_protocol"])
        self.assertIn(
            {"tag": "disease", "value": "case"},
            sample["channel"][0]["characteristics"],
        )
        self.assertEqual(1, len(sample["sra_run"]))
        self.assertEqual(
            ["https://example/1.fastq.gz", "https://example/2.fastq.gz"],
            [item["uri"] for item in sample["sra_run"][0]["fastq_files"]],
        )
        self.assertEqual([{"ref": "GSM1"}], package["series"]["sample_ref"])

    def test_preserves_unmapped_rows_and_columns_and_logs_warnings(self):
        fetcher = MagicMock()
        fetcher.resolve.return_value = resolved_input()

        with self.assertLogs("meta_standards_converter", level="WARNING") as logs:
            package = ae2json(fetcher=fetcher).convert("E-MTAB-1")[0]

        extension = package["mage_tab"]
        self.assertEqual("1.1", extension["version"])
        self.assertEqual("Mystery Row", extension["unmapped_idf_rows"][0]["label"])
        self.assertEqual("Mystery Column", extension["unmapped_sdrf_columns"][0]["header"])
        self.assertEqual(["x", "y"], extension["unmapped_sdrf_columns"][0]["values"])
        self.assertTrue(extension["warnings"])
        self.assertIn("unmapped", "\n".join(logs.output).lower())

    def test_merges_samples_across_multiple_sdrfs_in_first_seen_order(self):
        second = sdrf(rows=[
            ["GSM2", "Sample two", "Mus musculus", "EFO", "control", "control", "P-extract", "GSM2", "RNA", "A-TEST-2", "sequencing assay", "TRANSCRIPTOMIC", "RNA-SEQ", "ERR2", "https://example/3.fastq.gz", "ccc", "z"],
        ])
        fetcher = MagicMock()
        fetcher.resolve.return_value = resolved_input(sdrfs=[sdrf(), second])

        package = ae2json(fetcher=fetcher).convert("E-MTAB-1")[0]

        self.assertEqual(["GSM1", "GSM2"], [sample["iid"] for sample in package["sample"]])
        self.assertEqual(["A-TEST-1", "A-TEST-2"], [item["iid"] for item in package["platform"]])

    def test_conflicting_scalar_keeps_first_value_and_warns(self):
        rows = [
            ["GSM1", "First", "Homo sapiens", "EFO", "case", "case", "P-extract", "GSM1", "RNA", "A-TEST-1", "sequencing assay", "TRANSCRIPTOMIC", "RNA-SEQ", "ERR1", "https://example/1.fastq.gz", "aaa", "x"],
            ["GSM1", "Second", "Homo sapiens", "EFO", "case", "case", "P-extract", "GSM1", "RNA", "A-TEST-1", "sequencing assay", "TRANSCRIPTOMIC", "RNA-SEQ", "ERR1", "https://example/2.fastq.gz", "bbb", "y"],
        ]
        fetcher = MagicMock()
        fetcher.resolve.return_value = resolved_input(sdrfs=[sdrf(rows)])

        package = ae2json(fetcher=fetcher).convert("E-MTAB-1")[0]

        self.assertEqual("First", package["sample"][0]["title"])
        self.assertTrue(any("conflicting" in warning for warning in package["mage_tab"]["warnings"]))

    def test_library_protocol_does_not_replace_extraction_protocol(self):
        idf = IDF.replace(
            "Protocol Name\tP-collect\tP-extract\n"
            "Protocol Type\tsample collection protocol\tnucleic acid extraction protocol\n"
            "Protocol Description\tCollect samples\tExtract material\n",
            "Protocol Name\tP-collect\tP-library\tP-extract\n"
            "Protocol Type\tsample collection protocol\tnucleic acid library construction protocol\tnucleic acid extraction protocol\n"
            "Protocol Description\tCollect samples\tBuild library\tExtract material\n",
        )
        header = list(SDRF_HEADER)
        protocol_index = header.index("Protocol REF")
        header.insert(protocol_index + 1, "Protocol REF")
        row = [
            "GSM1", "Sample one", "Homo sapiens", "EFO", "case", "case",
            "P-library", "P-extract", "GSM1", "RNA", "A-TEST-1",
            "sequencing assay", "TRANSCRIPTOMIC", "RNA-SEQ", "ERR1",
            "https://example/1.fastq.gz", "aaa", "x",
        ]
        text = "\n".join("\t".join(values) for values in [header, row]) + "\n"
        fetcher = MagicMock()
        fetcher.resolve.return_value = resolved_input(idf=idf, sdrfs=[text])

        package = ae2json(fetcher=fetcher).convert("E-MTAB-1")[0]

        self.assertEqual("Extract material", package["sample"][0]["channel"][0]["extract_protocol"])

    def test_rejects_malformed_sdrf_row_width(self):
        fetcher = MagicMock()
        fetcher.resolve.return_value = resolved_input(sdrfs=["Source Name\tSample Name\nS1\n"])

        with self.assertRaisesRegex(ValueError, "columns"):
            ae2json(fetcher=fetcher).convert("E-MTAB-1")

    def test_writes_package_list_using_investigation_accession(self):
        fetcher = MagicMock()
        fetcher.resolve.return_value = resolved_input()

        with tempfile.TemporaryDirectory() as tmpdir:
            packages = ae2json(fetcher=fetcher).convert("E-MTAB-1", out=tmpdir)
            path = os.path.join(tmpdir, "E-MTAB-1.json")
            with open(path, encoding="utf-8") as handle:
                written = json.load(handle)

        self.assertEqual(packages, written)

    def test_semantic_round_trip_through_json2ae(self):
        fetcher = MagicMock()
        fetcher.resolve.return_value = resolved_input()
        package = ae2json(fetcher=fetcher).convert("E-MTAB-1")[0]
        pubmed = MagicMock()
        insdc = MagicMock()
        constructor = AEConstructor(
            idf_constructor=IDFConstructor(pubmed_fetcher=pubmed),
            sdrf_constructor=SDRFConstructor(insdc_fetcher=insdc),
        )

        magetab = constructor.miniml2magetab(package)
        rows = {row[0]: row for row in magetab}
        rendered_sdrf = rows["SDRF File"][1]

        self.assertEqual("Example study", rows["Investigation Title"][1])
        self.assertEqual("GSM1", rendered_sdrf[1][rendered_sdrf[0].index("Source Name")])
        self.assertIn("Characteristics[disease]", rendered_sdrf[0])
        self.assertIn("https://example/1.fastq.gz", str(rendered_sdrf))

    def test_records_versioned_lossless_roundtrip_sidecar(self):
        fetcher = MagicMock()
        fetcher.resolve.return_value = resolved_input()

        package = ae2json(fetcher=fetcher).convert("E-MTAB-1")[0]
        roundtrip = package["mage_tab"]["roundtrip"]

        self.assertEqual(1, roundtrip["schema_version"])
        self.assertEqual(64, len(roundtrip["semantic_sha256"]))
        self.assertEqual("Mystery Row", roundtrip["idf_rows"][-1][0])
        self.assertEqual("study1.sdrf.txt", roundtrip["sdrfs"][0]["name"])
        self.assertEqual("Mystery Column", roundtrip["sdrfs"][0]["rows"][0][-1])

    def test_unchanged_package_reuses_original_source_tables(self):
        fetcher = MagicMock()
        fetcher.resolve.return_value = resolved_input()
        package = ae2json(fetcher=fetcher).convert("E-MTAB-1")[0]

        magetab = AEConstructor().miniml2magetab(package)
        rows = {row[0]: row for row in magetab}

        self.assertEqual(["Investigation Title", "Example study"], rows["Investigation Title"])
        self.assertEqual(["Mystery Row", "keep me"], rows["Mystery Row"])
        self.assertEqual("Mystery Column", rows["SDRF File"][1][0][-1])
        self.assertEqual("x", rows["SDRF File"][1][1][-1])

    def test_edited_json_wins_while_unmapped_metadata_is_restored(self):
        fetcher = MagicMock()
        fetcher.resolve.return_value = resolved_input()
        package = ae2json(fetcher=fetcher).convert("E-MTAB-1")[0]
        edited = copy.deepcopy(package)
        edited["series"]["title"] = "Edited title"

        magetab = AEConstructor().miniml2magetab(edited)
        rows = {row[0]: row for row in magetab}

        self.assertEqual(["Investigation Title", "Edited title"], rows["Investigation Title"])
        self.assertEqual(["Mystery Row", "keep me"], rows["Mystery Row"])
        self.assertIn("Mystery Column", rows["SDRF File"][1][0])


if __name__ == "__main__":
    unittest.main()
