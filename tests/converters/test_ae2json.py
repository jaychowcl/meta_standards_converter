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


ROOT = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from meta_standards_converter.magetab.constructor import AEConstructor  # noqa: E402
from meta_standards_converter.magetab.idf import IDFConstructor  # noqa: E402
from meta_standards_converter.magetab.semantics import (  # noqa: E402
    MAGETabModelError,
    overlay_core,
    render_model,
    validate_model,
)
from meta_standards_converter.magetab.parser import AEParser  # noqa: E402
from meta_standards_converter.magetab.sdrf.constructor import SDRFConstructor
from meta_standards_converter.sources.magetab import (  # noqa: E402
    MAGETabInput,
    TextResource,
)
from meta_standards_converter.converters.ae2json import AE2JSONConverter  # noqa: E402
from meta_standards_converter.miniml import MINiMLCodec, MINiMLPackage  # noqa: E402
from meta_standards_converter.runtime_contracts import get_resource_profile  # noqa: E402


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
    def test_resource_profile_and_explicit_source_hosts_reach_default_fetcher(self):
        profile = get_resource_profile(
            "standard", overrides={"max_xml_bytes": 4096}
        )

        converter = AE2JSONConverter(
            resource_profile=profile,
            source_hosts=("metadata.example.org",),
        )

        self.assertIs(profile, converter.fetcher.resource_profile)
        self.assertEqual(
            frozenset({"metadata.example.org"}),
            converter.fetcher.retrieval_policy.allowed_hosts,
        )

    def test_model_validation_rejects_unsupported_or_malformed_models(self):
        with self.assertRaisesRegex(MAGETabModelError, "schema_version"):
            validate_model({"schema_version": 2})
        with self.assertRaisesRegex(MAGETabModelError, "assay_paths"):
            validate_model({
                "schema_version": 1,
                "idf_layout": [],
                "protocols": [],
                "declarations": {},
                "assay_paths": "invalid",
                "sdrfs": [],
                "investigation_fields": [],
            })

    def test_harmonized_model_annotations_render_and_parse_additively(self):
        header = [
            "Source Name", "Sample Name", "Protocol REF",
            "Parameter Value[duration]", "Unit", "Assay Name",
        ]
        rows = [["source-1", "sample-1", "P-extract", "30", "minutes", "assay-1"]]
        text = "\n".join("\t".join(values) for values in [header, *rows]) + "\n"
        fetcher = MagicMock()
        fetcher.resolve.return_value = resolved_input(sdrfs=[text])
        package = AE2JSONConverter(fetcher=fetcher).convert("E-MTAB-1")[0]
        payload = package.to_mapping()
        step = next(
            item for item in payload["series"]["assay_paths"][0]["steps"]
            if item.get("kind") == "protocol_application"
        )
        unit = step["parameter_values"][0]["unit"]
        unit.update({
            "hz_unit": "minute",
            "hz_unit_onto": "uo",
            "hz_unit_id": "UO:0000031",
        })
        package = MINiMLCodec().decode(payload).package

        rendered = AEConstructor().miniml2magetab(package)
        table = next(row[1] for row in rendered if row[0] == "SDRF File")
        parameter = table[0].index("Parameter Value[duration]")
        self.assertEqual("30", table[1][parameter])
        self.assertEqual("minutes", table[1][parameter + 1])
        self.assertIn("Comment[hz_unit]", table[0])

        idf_rows = []
        for row in rendered:
            if row[0] == "SDRF File":
                idf_rows.append(["SDRF File", "roundtrip.sdrf.txt"])
            else:
                idf_rows.append(row)
        reparsed = AEParser().parse(resolved_input(
            idf="\n".join("\t".join(str(value) for value in row) for row in idf_rows) + "\n",
            sdrfs=["\n".join("\t".join(str(value) for value in row) for row in table) + "\n"],
        ))
        reparsed_step = next(
            item for item in reparsed["series"]["assay_paths"][0]["steps"]
            if item.get("kind") == "protocol_application"
        )["parameter_values"][0]
        self.assertEqual("30", reparsed_step["value"])
        self.assertEqual("minutes", reparsed_step["unit"]["value"])
    def test_typed_model_preserves_ragged_and_label_only_idf_rows(self):
        idf = IDF.replace(
            "Protocol Description\tCollect samples\tExtract material\n",
            "Protocol Description\tCollect samples\tExtract material\n"
            "Protocol Parameters\ttemperature\n"
            "Quality Control Type\n",
        )
        fetcher = MagicMock()
        fetcher.resolve.return_value = resolved_input(idf=idf)

        package = AE2JSONConverter(fetcher=fetcher).convert("E-MTAB-1")[0]
        self.assertEqual(["temperature"], package["series"]["protocols"][0]["parameters"])
        self.assertEqual([{"value": ""}], package["series"]["quality_controls"])

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

        package = AE2JSONConverter(fetcher=fetcher).convert("E-MTAB-1")[0]
        series = package["series"]
        self.assertEqual(["centrifuge"], series["protocols"][1]["hardware"])
        self.assertEqual(["ExtractSoft"], series["protocols"][1]["software"])
        self.assertEqual(["speed"], series["protocols"][1]["parameters"])
        self.assertEqual(["John Doe"], series["protocols"][1]["contacts"])
        self.assertEqual("biological replicate", series["quality_controls"][0]["value"])
        self.assertEqual("EFO:0000001", series["quality_controls"][0]["term_accession_number"])
        self.assertEqual(2, len(series["assay_paths"]))
        assay_names = [
            next(step["name"] for step in path["steps"] if step["kind"] == "assay")
            for path in series["assay_paths"]
        ]
        self.assertEqual(["assay-1", "assay-2"], assay_names)
        age = next(
            characteristic
            for step in series["assay_paths"][0]["steps"]
            for characteristic in step.get("characteristics", [])
            if characteristic["name"] == "age"
        )
        self.assertEqual("year", age["unit"]["value"])
        self.assertEqual("UO:0000036", age["unit"]["term_accession_number"])
        barcode = next(comment for comment in age["comments"] if comment["name"] == "cell barcode size")
        self.assertEqual("16", barcode["value"])

    def test_model_edits_render_without_merging_into_miniml_fields(self):
        fetcher = MagicMock()
        fetcher.resolve.return_value = resolved_input()
        package = AE2JSONConverter(fetcher=fetcher).convert("E-MTAB-1")[0]
        payload = package.to_mapping()
        payload["series"]["protocols"][1]["hardware"] = ["edited centrifuge"]
        package = MINiMLCodec().decode(payload).package

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
        package = AE2JSONConverter(fetcher=fetcher).convert("E-MTAB-1")[0]
        payload = package.to_mapping()
        payload["series"]["title"] = "Edited core title"
        payload["series"]["assay_paths"][1]["steps"][-1]["name"] = "edited-scan-2"
        package = MINiMLCodec().decode(payload).package

        magetab = AEConstructor().miniml2magetab(package)
        rows_by_label = {row[0]: row for row in magetab}
        rendered_sdrf = rows_by_label["SDRF File"][1]

        self.assertEqual("Edited core title", rows_by_label["Investigation Title"][1])
        self.assertEqual(2, len(rendered_sdrf) - 1)
        assay_index = rendered_sdrf[0].index("Assay Name")
        scan_index = rendered_sdrf[0].index("Scan Name")
        self.assertEqual(["assay-1", "assay-2"], [row[assay_index] for row in rendered_sdrf[1:]])
        self.assertEqual("edited-scan-2", rendered_sdrf[2][scan_index])

    def test_harmonized_groups_export_while_paths_preserve_multiplicity(self):
        fetcher = MagicMock()
        fetcher.resolve.return_value = resolved_input()
        package = AE2JSONConverter(fetcher=fetcher).convert("E-MTAB-1")[0]
        payload = package.to_mapping()
        characteristics = payload["sample"][0]["channel"][0]["characteristics"]
        characteristics.extend([
            {"name": "hz_disease", "value": "disease"},
            {"name": "hz_disease_onto", "value": "mondo"},
            {"name": "hz_disease_id", "value": "MONDO:0000001"},
        ])
        package = MINiMLCodec().decode(payload).package

        magetab = AEConstructor().miniml2magetab(package)
        rendered_sdrf = next(row[1] for row in magetab if row[0] == "SDRF File")

        self.assertEqual(2, len(rendered_sdrf) - 1)
        self.assertIn("Characteristics[hz_disease]", rendered_sdrf[0])
        rendered = package["sample"][0]["channel"][0]["characteristics"]
        self.assertIn({"name": "hz_disease_id", "value": "MONDO:0000001"}, rendered)

    def test_ecto_and_pcl_characteristics_parse_and_render_additively(self):
        header = [
            "Source Name",
            "Sample Name",
            "Characteristics[exposure]", "Term Source REF", "Term Accession Number",
            "Characteristics[cell state]", "Term Source REF", "Term Accession Number",
            "Assay Name",
        ]
        row = [
            "source-1",
            "sample-1",
            "exposure to bleomycin via injection", "ecto", "ECTO:0900222",
            "Fbl_24",
            "pcl", "PCL:0015251",
            "assay-1",
        ]
        text = "\n".join("\t".join(values) for values in [header, row]) + "\n"
        fetcher = MagicMock()
        fetcher.resolve.return_value = resolved_input(sdrfs=[text])

        package = AE2JSONConverter(fetcher=fetcher).convert("E-MTAB-1")[0]
        characteristics = {
            item["name"]: item
            for item in package["sample"][0]["channel"][0]["characteristics"]
        }
        rendered = AEConstructor().miniml2magetab(package)
        sdrf = next(item[1] for item in rendered if item[0] == "SDRF File")

        self.assertEqual("ECTO:0900222", characteristics["exposure"]["term_accession_number"])
        self.assertEqual("PCL:0015251", characteristics["cell state"]["term_accession_number"])
        self.assertIn("Characteristics[exposure]", sdrf[0])
        self.assertIn("Characteristics[cell state]", sdrf[0])

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
        self.assertEqual("old", rendered[3][1])
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

    def test_overlay_leaves_ambiguous_inserted_values_blank(self):
        model = [
            ["Source Name", "Protocol REF"],
            ["sample-1", "P-1"],
        ]
        core = [
            ["Source Name", "Characteristics[hz_disease]", "Protocol REF"],
            ["sample-1", "disease-a", "P-generated"],
            ["sample-1", "disease-b", "P-generated"],
        ]

        result = overlay_core([["SDRF File", model]], [["SDRF File", core]])
        rendered = result[0][1]

        self.assertEqual(
            ["Source Name", "Characteristics[hz_disease]", "Protocol REF"],
            rendered[0],
        )
        self.assertEqual(["sample-1", "", "P-1"], rendered[1])

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

        packages = AE2JSONConverter(fetcher=fetcher).convert("E-MTAB-1")

        self.assertEqual(1, len(packages))
        package = packages[0]
        self.assertEqual("3.0", package["miniml_schema_version"])
        self.assertEqual(
            "https://www.ebi.ac.uk/biostudies/misc/MAGE-TABv1.1_2011_07_28.pdf",
            package["source"]["schema_location"],
        )
        self.assertEqual("E-MTAB-1", package["series"]["iid"])
        self.assertEqual("MAGE-TAB", package["source"]["format"])
        self.assertEqual("Example study", package["series"]["title"])
        self.assertEqual(
            ["GSE123", "E-MTAB-1"],
            [item["value"] for item in package["series"]["accession"]],
        )
        self.assertEqual(
            [{"factor": "disease state", "name": "disease", "type": {"value": "disease"}}],
            package["series"]["variable"],
        )
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
            {"name": "disease", "value": "case"},
            sample["channel"][0]["characteristics"],
        )
        self.assertEqual(1, len(sample["sra_run"]))
        self.assertEqual(
            ["https://example/1.fastq.gz", "https://example/2.fastq.gz"],
            [item["uri"] for item in sample["sra_run"][0]["fastq_files"]],
        )
        self.assertEqual([{"ref": "GSM1"}], package["series"]["sample_ref"])

    def test_parser_emits_canonical_versioned_package(self):
        fetcher = MagicMock()
        fetcher.resolve.return_value = resolved_input()

        package = AE2JSONConverter(fetcher=fetcher).convert("E-MTAB-1")[0]

        self.assertIsInstance(package, MINiMLPackage)
        self.assertEqual("3.0", package.miniml_schema_version)
        self.assertEqual(package, MINiMLPackage.from_mapping(package.to_mapping()))

    def test_e_mtab_6486_values_are_losslessly_normalized_for_strict_miniml(self):
        fixture_dir = os.path.join(ROOT, "tests", "fixtures", "arrayexpress")
        with open(
            os.path.join(fixture_dir, "E-MTAB-6486.idf.txt"), encoding="utf-8"
        ) as handle:
            idf = handle.read()
        with open(
            os.path.join(fixture_dir, "E-MTAB-6486.sdrf.txt"), encoding="utf-8"
        ) as handle:
            frozen_sdrf = handle.read()
        fetcher = MagicMock()
        fetcher.resolve.return_value = resolved_input(
            idf=idf,
            sdrfs=[frozen_sdrf],
        )

        package = AE2JSONConverter(fetcher=fetcher).convert("E-MTAB-6486")[0]
        mapping = MINiMLCodec().decode(
            package.to_mapping(), strict=True
        ).package.to_mapping()

        self.assertEqual(
            "total RNA",
            mapping["sample"][0]["channel"][0]["molecule"]["value"],
        )
        self.assertIn(
            {"name": "material type", "value": "cell"},
            mapping["sample"][0]["channel"][0]["characteristics"],
        )
        self.assertEqual(
            [{"factor": "agent", "name": "compound", "type": {"value": "compound"}}],
            mapping["series"]["variable"],
        )
        self.assertGreaterEqual(
            {item["iid"] for item in mapping["database"]},
            {"ArrayExpress", "ENA"},
        )

    def test_series_iid_prefers_explicit_arrayexpress_accession(self):
        idf = IDF.replace(
            "Investigation Accession\tE-MTAB-1\n",
            "Investigation Accession\tE-MTAB-1\n"
            "Comment[ArrayExpressAccession]\tE-MTAB-999\n",
        )
        fetcher = MagicMock()
        fetcher.resolve.return_value = resolved_input(idf=idf)

        package = AE2JSONConverter(fetcher=fetcher).convert("E-MTAB-1")[0]

        self.assertEqual("E-MTAB-999", package["series"]["iid"])

    def test_series_iid_uses_investigation_fallback_not_geo_secondary(self):
        idf = IDF.replace(
            "Investigation Accession\tE-MTAB-1\n",
            "Investigation Accession\tLOCAL-STUDY-1\n",
        )
        fetcher = MagicMock()
        fetcher.resolve.return_value = resolved_input(idf=idf)

        package = AE2JSONConverter(fetcher=fetcher).convert("LOCAL-STUDY-1")[0]

        self.assertEqual("LOCAL-STUDY-1", package["series"]["iid"])
        self.assertNotEqual("GSE123", package["series"]["iid"])

    def test_preserves_unmapped_rows_and_columns_and_logs_warnings(self):
        fetcher = MagicMock()
        fetcher.resolve.return_value = resolved_input()

        with self.assertLogs("meta_standards_converter", level="WARNING") as logs:
            package = AE2JSONConverter(fetcher=fetcher).convert("E-MTAB-1")[0]

        self.assertNotIn("mage_tab", package)
        self.assertEqual(["idf", "sdrf"], [item["kind"] for item in package["source"]["documents"]])
        comments = [comment["name"] for path in package["series"]["assay_paths"] for step in path["steps"] for comment in step.get("comments", [])]
        self.assertIn("Mystery Column", comments)
        self.assertIn("unmapped", "\n".join(logs.output).lower())

    def test_merges_samples_across_multiple_sdrfs_in_first_seen_order(self):
        second = sdrf(rows=[
            ["GSM2", "Sample two", "Mus musculus", "EFO", "control", "control", "P-extract", "GSM2", "RNA", "A-TEST-2", "sequencing assay", "TRANSCRIPTOMIC", "RNA-SEQ", "ERR2", "https://example/3.fastq.gz", "ccc", "z"],
        ])
        fetcher = MagicMock()
        fetcher.resolve.return_value = resolved_input(sdrfs=[sdrf(), second])

        package = AE2JSONConverter(fetcher=fetcher).convert("E-MTAB-1")[0]

        self.assertEqual(["GSM1", "GSM2"], [sample["iid"] for sample in package["sample"]])
        self.assertEqual(["A-TEST-1", "A-TEST-2"], [item["iid"] for item in package["platform"]])

    def test_conflicting_scalar_keeps_first_value_and_warns(self):
        rows = [
            ["GSM1", "First", "Homo sapiens", "EFO", "case", "case", "P-extract", "GSM1", "RNA", "A-TEST-1", "sequencing assay", "TRANSCRIPTOMIC", "RNA-SEQ", "ERR1", "https://example/1.fastq.gz", "aaa", "x"],
            ["GSM1", "Second", "Homo sapiens", "EFO", "case", "case", "P-extract", "GSM1", "RNA", "A-TEST-1", "sequencing assay", "TRANSCRIPTOMIC", "RNA-SEQ", "ERR1", "https://example/2.fastq.gz", "bbb", "y"],
        ]
        fetcher = MagicMock()
        fetcher.resolve.return_value = resolved_input(sdrfs=[sdrf(rows)])

        with self.assertLogs("meta_standards_converter", level="WARNING") as logs:
            package = AE2JSONConverter(fetcher=fetcher).convert("E-MTAB-1")[0]

        self.assertEqual("First", package["sample"][0]["title"])
        self.assertIn("conflicting title", "\n".join(logs.output))

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

        package = AE2JSONConverter(fetcher=fetcher).convert("E-MTAB-1")[0]

        self.assertEqual("Extract material", package["sample"][0]["channel"][0]["extract_protocol"])

    def test_rejects_malformed_sdrf_row_width(self):
        fetcher = MagicMock()
        fetcher.resolve.return_value = resolved_input(sdrfs=["Source Name\tSample Name\nS1\n"])

        with self.assertRaisesRegex(ValueError, "columns"):
            AE2JSONConverter(fetcher=fetcher).convert("E-MTAB-1")

    def test_writes_package_list_using_investigation_accession(self):
        fetcher = MagicMock()
        fetcher.resolve.return_value = resolved_input()

        with tempfile.TemporaryDirectory() as tmpdir:
            packages = AE2JSONConverter(fetcher=fetcher).convert("E-MTAB-1", out=tmpdir)
            path = os.path.join(tmpdir, "E-MTAB-1.json")
            with open(path, encoding="utf-8") as handle:
                written = json.load(handle)

        self.assertEqual(MINiMLCodec().encode_many(packages), written)

    def test_semantic_round_trip_through_json2ae(self):
        fetcher = MagicMock()
        fetcher.resolve.return_value = resolved_input()
        package = AE2JSONConverter(fetcher=fetcher).convert("E-MTAB-1")[0]
        pubmed = MagicMock()
        insdc = MagicMock()
        constructor = AEConstructor(
            pubmed_client=pubmed,
            insdc_client=insdc,
        )

        magetab = constructor.miniml2magetab(package)
        rows = {row[0]: row for row in magetab}
        rendered_sdrf = rows["SDRF File"][1]

        self.assertEqual("Example study", rows["Investigation Title"][1])
        self.assertEqual("GSM1", rendered_sdrf[1][rendered_sdrf[0].index("Source Name")])
        self.assertIn("Characteristics[disease]", rendered_sdrf[0])
        self.assertIn("https://example/1.fastq.gz", str(rendered_sdrf))

    def test_records_source_documents_without_raw_roundtrip_sidecar(self):
        fetcher = MagicMock()
        fetcher.resolve.return_value = resolved_input()

        package = AE2JSONConverter(fetcher=fetcher).convert("E-MTAB-1")[0]
        self.assertNotIn("mage_tab", package)
        self.assertEqual(
            ["study.idf.txt", "study1.sdrf.txt"],
            [item["name"] for item in package["source"]["documents"]],
        )

    def test_unchanged_package_renders_semantic_source_content(self):
        fetcher = MagicMock()
        fetcher.resolve.return_value = resolved_input()
        package = AE2JSONConverter(fetcher=fetcher).convert("E-MTAB-1")[0]

        magetab = AEConstructor().miniml2magetab(package)
        rows = {row[0]: row for row in magetab}

        self.assertEqual(["Investigation Title", "Example study"], rows["Investigation Title"])
        self.assertNotIn("Mystery Row", rows)
        self.assertEqual("Comment[Mystery Column]", rows["SDRF File"][1][0][-1])
        self.assertEqual("x", rows["SDRF File"][1][1][-1])

    def test_edited_json_wins_while_semantic_sdrf_comments_are_retained(self):
        fetcher = MagicMock()
        fetcher.resolve.return_value = resolved_input()
        package = AE2JSONConverter(fetcher=fetcher).convert("E-MTAB-1")[0]
        edited_payload = package.to_mapping()
        edited_payload["series"]["title"] = "Edited title"
        edited = MINiMLCodec().decode(edited_payload).package

        magetab = AEConstructor().miniml2magetab(edited)
        rows = {row[0]: row for row in magetab}

        self.assertEqual(["Investigation Title", "Edited title"], rows["Investigation Title"])
        self.assertNotIn("Mystery Row", rows)
        self.assertIn("Comment[Mystery Column]", rows["SDRF File"][1][0])


if __name__ == "__main__":
    unittest.main()
