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
from meta_standards_converter.ae_handlers.ae_model import overlay_core, render_model  # noqa: E402
from meta_standards_converter.ae_handlers.ae_parser import AEParser  # noqa: E402
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
    def test_typed_model_preserves_ragged_and_label_only_idf_rows(self):
        idf = IDF.replace(
            "Protocol Description\tCollect samples\tExtract material\n",
            "Protocol Description\tCollect samples\tExtract material\n"
            "Protocol Parameters\ttemperature\n"
            "Quality Control Type\n",
        )
        fetcher = MagicMock()
        fetcher.resolve.return_value = resolved_input(idf=idf)

        package = ae2json(fetcher=fetcher).convert("E-MTAB-1")[0]
        rendered = render_model(package["mage_tab"]["model"])

        self.assertIn(["Protocol Parameters", "temperature"], rendered)
        self.assertIn(["Quality Control Type"], rendered)

    def test_builds_editable_typed_magetab_model(self):
        idf = IDF.replace(
            "Protocol Description\tCollect samples\tExtract material\n",
            "Protocol Description\tCollect samples\tExtract material\n"
            "Protocol Hardware\tfreezer\tcentrifuge\n"
            "Protocol Software\tLIMS 2\tExtractSoft\n"
            "Protocol Parameters\ttemperature\tspeed\n"
            "Protocol Contact\tJane Doe\tJohn Doe\n"
            "Quality Control Type\tbiological replicate\n"
            "Quality Control Term Source REF\tEFO\n"
            "Quality Control Term Accession Number\tEFO:0000001\n"
            "Replicate Type\ttechnical replicate\n"
            "Normalization Type\tquantile normalization\n",
        )
        header = [
            "Source Name", "Sample Name", "Protocol REF", "Extract Name",
            "Assay Name", "Hybridization Name", "Scan Name",
            "Characteristics[age]", "Unit", "Term Source REF",
            "Term Accession Number", "Comment[cell barcode size]",
        ]
        rows = [
            ["source-1", "sample-1", "P-extract", "extract-1", "assay-1", "hyb-1", "scan-1", "5", "year", "UO", "UO:0000036", "16"],
            ["source-1", "sample-1", "P-extract", "extract-1", "assay-2", "hyb-2", "scan-2", "5", "year", "UO", "UO:0000036", "16"],
        ]
        text = "\n".join("\t".join(values) for values in [header, *rows]) + "\n"
        fetcher = MagicMock()
        fetcher.resolve.return_value = resolved_input(idf=idf, sdrfs=[text])

        package = ae2json(fetcher=fetcher).convert("E-MTAB-1")[0]
        model = package["mage_tab"]["model"]

        self.assertEqual(1, model["schema_version"])
        self.assertEqual("centrifuge", model["protocols"][1]["hardware"])
        self.assertEqual("ExtractSoft", model["protocols"][1]["software"])
        self.assertEqual("speed", model["protocols"][1]["parameters"])
        self.assertEqual("John Doe", model["protocols"][1]["contact"])
        self.assertEqual("biological replicate", model["declarations"]["quality_control"][0]["value"])
        self.assertEqual("EFO:0000001", model["declarations"]["quality_control"][0]["term_accession_number"])
        self.assertEqual(2, len(model["assay_paths"]))
        self.assertEqual(["assay-1", "assay-2"], [path["binding"]["assay_name"] for path in model["assay_paths"]])
        age = next(step for step in model["assay_paths"][0]["steps"] if step["kind"] == "attribute")
        self.assertEqual("year", age["unit"])
        self.assertEqual("UO:0000036", age["term_accession_number"])
        barcode = next(step for step in model["assay_paths"][0]["steps"] if step.get("name") == "cell barcode size")
        self.assertEqual("16", barcode["value"])

    def test_model_edits_render_without_merging_into_miniml_fields(self):
        fetcher = MagicMock()
        fetcher.resolve.return_value = resolved_input()
        package = ae2json(fetcher=fetcher).convert("E-MTAB-1")[0]
        model = package["mage_tab"]["model"]
        model["protocols"][1]["hardware"] = "edited centrifuge"

        magetab = AEConstructor().miniml2magetab(package)
        rows = {row[0]: row for row in magetab}

        self.assertEqual("edited centrifuge", rows["Protocol Hardware"][2])
        self.assertNotIn("edited centrifuge", str(package["sample"]))

    def test_core_edits_win_while_typed_model_preserves_path_multiplicity(self):
        header = ["Source Name", "Sample Name", "Assay Name", "Hybridization Name", "Scan Name"]
        rows = [
            ["source-1", "sample-1", "assay-1", "hyb-1", "scan-1"],
            ["source-1", "sample-1", "assay-2", "hyb-2", "scan-2"],
        ]
        text = "\n".join("\t".join(values) for values in [header, *rows]) + "\n"
        fetcher = MagicMock()
        fetcher.resolve.return_value = resolved_input(sdrfs=[text])
        package = ae2json(fetcher=fetcher).convert("E-MTAB-1")[0]
        package["series"]["title"] = "Edited core title"
        package["mage_tab"]["model"]["assay_paths"][1]["steps"][-1]["value"] = "edited-scan-2"

        magetab = AEConstructor().miniml2magetab(package)
        rows_by_label = {row[0]: row for row in magetab}
        rendered_sdrf = rows_by_label["SDRF File"][1]

        self.assertEqual("Edited core title", rows_by_label["Investigation Title"][1])
        self.assertEqual(2, len(rendered_sdrf) - 1)
        self.assertEqual(["assay-1", "assay-2"], [row[2] for row in rendered_sdrf[1:]])
        self.assertEqual("edited-scan-2", rendered_sdrf[2][4])

    def test_core_harmonization_columns_are_unioned_while_model_preserves_multiplicity(self):
        fetcher = MagicMock()
        fetcher.resolve.return_value = resolved_input()
        package = ae2json(fetcher=fetcher).convert("E-MTAB-1")[0]
        package["sample"][0]["channel"][0]["characteristics"].extend(
            [
                {"tag": "hz_cell_type", "value": "regulatory T cell"},
                {"tag": "hz_cell_type_id", "value": "CL:0000815"},
                {"tag": "hz_cell_type_onto", "value": "cl"},
            ]
        )

        magetab = AEConstructor().miniml2magetab(package)
        rendered_sdrf = next(row[1] for row in magetab if row[0] == "SDRF File")

        for label, value in (
            ("Characteristics[hz_cell_type]", "regulatory T cell"),
            ("Characteristics[hz_cell_type_id]", "CL:0000815"),
            ("Characteristics[hz_cell_type_onto]", "cl"),
        ):
            with self.subTest(label=label):
                index = rendered_sdrf[0].index(label)
                self.assertEqual([value, value], [row[index] for row in rendered_sdrf[1:]])
        self.assertEqual(2, len(rendered_sdrf) - 1)
        self.assertIn("Mystery Column", rendered_sdrf[0])

    def test_overlay_unions_allowlisted_idf_rows_and_nonstructural_sdrf_columns(self):
        model_sdrf = [
            ["Source Name", "Characteristics[disease]", "Protocol REF", "Assay Name", "Mystery Column"],
            ["sample-1", "old", "P-1", "assay-1", "keep-1"],
            ["sample-1", "old", "P-1", "assay-2", "keep-2"],
            ["sample-2", "old", "P-1", "assay-3", "keep-3"],
        ]
        core_sdrf = [
            [
                "Source Name",
                "Characteristics[disease]",
                "Characteristics[hz_cell_type]",
                "Characteristics[hz_cell_type_id]",
                "Characteristics[hz_cell_type_onto]",
                "Extract Name",
                "Protocol REF",
                "Assay Name",
            ],
            [
                "sample-1",
                "edited",
                "regulatory T cell",
                "CL:0000815",
                "cl",
                "extract-1",
                "P-generated",
                "generated-assay",
            ],
        ]
        model = [
            ["MAGE-TAB Version", "1.1"],
            ["Investigation Accession", "E-MTAB-1"],
            ["Mystery Row", "keep me"],
            ["SDRF File", model_sdrf],
        ]
        core = [
            ["Investigation Accession", "E-MTAB-1"],
            ["Experiment Description", "Edited description"],
            ["Generated Custom Row", "drop me"],
            ["SDRF File", core_sdrf],
        ]

        result = overlay_core(model, core)
        labels = [row[0] for row in result]
        rendered = next(row[1] for row in result if row[0] == "SDRF File")

        self.assertEqual(
            [
                "MAGE-TAB Version",
                "Investigation Accession",
                "Mystery Row",
                "Experiment Description",
                "SDRF File",
            ],
            labels,
        )
        self.assertNotIn("Generated Custom Row", labels)
        self.assertEqual(
            [
                "Source Name",
                "Characteristics[disease]",
                "Characteristics[hz_cell_type]",
                "Characteristics[hz_cell_type_id]",
                "Characteristics[hz_cell_type_onto]",
                "Protocol REF",
                "Assay Name",
                "Mystery Column",
            ],
            rendered[0],
        )
        self.assertNotIn("Extract Name", rendered[0])
        self.assertEqual("edited", rendered[1][1])
        self.assertEqual("edited", rendered[2][1])
        self.assertEqual("", rendered[3][1])
        self.assertEqual("regulatory T cell", rendered[1][2])
        self.assertEqual("regulatory T cell", rendered[2][2])
        self.assertEqual("", rendered[3][2])
        self.assertEqual(["assay-1", "assay-2", "assay-3"], [row[6] for row in rendered[1:]])
        self.assertEqual(["keep-1", "keep-2", "keep-3"], [row[7] for row in rendered[1:]])

    def test_overlay_matches_repeated_headers_by_occurrence_and_inserts_companion_group(self):
        model = [
            [
                "Source Name",
                "Characteristics[age]",
                "Unit",
                "Characteristics[age]",
                "Unit",
                "Protocol REF",
            ],
            ["sample-1", "old-1", "old-unit-1", "old-2", "old-unit-2", "P-1"],
        ]
        core = [
            [
                "Source Name",
                "Characteristics[age]",
                "Unit",
                "Characteristics[age]",
                "Unit",
                "Characteristics[hz_age]",
                "Term Source REF",
                "Term Accession Number",
                "Protocol REF",
            ],
            ["sample-1", "10", "year", "20", "month", "adult", "EFO", "EFO:0001272", "P-2"],
        ]

        result = overlay_core([["SDRF File", model]], [["SDRF File", core]])
        rendered = result[0][1]

        self.assertEqual(
            [
                "Source Name",
                "Characteristics[age]",
                "Unit",
                "Characteristics[age]",
                "Unit",
                "Characteristics[hz_age]",
                "Term Source REF",
                "Term Accession Number",
                "Protocol REF",
            ],
            rendered[0],
        )
        self.assertEqual(
            ["sample-1", "10", "year", "20", "month", "adult", "EFO", "EFO:0001272", "P-1"],
            rendered[1],
        )

    def test_recognizes_all_generated_single_cell_comment_headers(self):
        generated_headers = (
            "Comment[cdna read]",
            "Comment[cdna read offset]",
            "Comment[cdna read size]",
            "Comment[cell barcode offset]",
            "Comment[cell barcode read]",
            "Comment[cell barcode size]",
            "Comment[end bias]",
            "Comment[input molecule]",
            "Comment[library construction]",
            "Comment[primer]",
            "Comment[LIBRARY_STRAND]",
            "Comment[sample barcode offset]",
            "Comment[sample barcode read]",
            "Comment[sample barcode size]",
            "Comment[single cell isolation]",
            "Comment[spike in]",
            "Comment[umi barcode offset]",
            "Comment[umi barcode read]",
            "Comment[umi barcode size]",
        )

        self.assertTrue(all(AEParser()._known_sdrf_header(label) for label in generated_headers))

    def test_maps_idf_and_sdrf_to_miniml_compatible_package(self):
        fetcher = MagicMock()
        fetcher.resolve.return_value = resolved_input()

        packages = ae2json(fetcher=fetcher).convert("E-MTAB-1")

        self.assertEqual(1, len(packages))
        package = packages[0]
        self.assertEqual("magetabv1.1", package["version"])
        self.assertEqual(
            "https://www.ebi.ac.uk/biostudies/misc/MAGE-TABv1.1_2011_07_28.pdf",
            package["schema_location"],
        )
        self.assertEqual("E-MTAB-1", package["series"]["iid"])
        self.assertEqual("1.1", package["mage_tab"]["version"])
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

    def test_series_iid_prefers_explicit_arrayexpress_accession(self):
        idf = IDF.replace(
            "Investigation Accession\tE-MTAB-1\n",
            "Investigation Accession\tE-MTAB-1\n"
            "Comment[ArrayExpressAccession]\tE-MTAB-999\n",
        )
        fetcher = MagicMock()
        fetcher.resolve.return_value = resolved_input(idf=idf)

        package = ae2json(fetcher=fetcher).convert("E-MTAB-1")[0]

        self.assertEqual("E-MTAB-999", package["series"]["iid"])

    def test_series_iid_uses_investigation_fallback_not_geo_secondary(self):
        idf = IDF.replace(
            "Investigation Accession\tE-MTAB-1\n",
            "Investigation Accession\tLOCAL-STUDY-1\n",
        )
        fetcher = MagicMock()
        fetcher.resolve.return_value = resolved_input(idf=idf)

        package = ae2json(fetcher=fetcher).convert("LOCAL-STUDY-1")[0]

        self.assertEqual("LOCAL-STUDY-1", package["series"]["iid"])
        self.assertNotEqual("GSE123", package["series"]["iid"])

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
        self.assertEqual(64, len(roundtrip["model_sha256"]))
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
