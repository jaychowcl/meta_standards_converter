# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
import os
import re
import sys
import unittest
from unittest.mock import Mock, patch


ROOT = os.path.dirname(os.path.dirname(__file__))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from meta_standards_converter.ae_handlers.ae_constructor import ProtocolRegistry  # noqa: E402
from meta_standards_converter.ae_handlers.ae_sdrf_handlers import (  # noqa: E402
    SDRFConstructor,
    _BulkSequencingSDRFHandler,
    _DropletSingleCellSequencingSDRFHandler,
    _PlateSingleCellSequencingSDRFHandler,
    _SequencingSDRFHandler,
    _TenXV2DropletSingleCellSequencingSDRFHandler,
    _TenXV3DropletSingleCellSequencingSDRFHandler,
    classify_file,
)
from meta_standards_converter.insdc_handlers.insdc_webfetcher import INSDCWebfetcher  # noqa: E402


class EmptyINSDCFetcher(INSDCWebfetcher):
    def fetch_sra_runs(self, accession: str) -> list:
        return []


class Parent(SDRFConstructor):
    def __init__(self):
        super().__init__(insdc_fetcher=EmptyINSDCFetcher())


class FixedTechParent(Parent):
    def __init__(self, tech_type):
        super().__init__()
        self.tech_type = tech_type

    def _detect_sdrf_technology(self, data: dict) -> str:
        return self.tech_type


def base_data(sample, platform_technology="high-throughput sequencing"):
    return {
        "series": {
            "accession": [{"value": "GSE1"}],
            "sample_ref": [{"ref": sample["iid"]}],
            "title": "example",
        },
        "platform": [
            {
                "iid": "GPL1",
                "title": "Example platform",
                "technology": platform_technology,
                "accession": [{"value": "GPL1"}],
            }
        ],
        "sample": [sample],
    }


class TestSDRFGraphHandlers(unittest.TestCase):
    def cell(self, sdrf, row_index, label):
        return sdrf[row_index][sdrf[0].index(label)]

    def test_renders_blank_strings_not_none(self):
        sample = {
            "iid": "GSM1",
            "title": "sample 1",
            "accession": [{"value": "GSM1"}],
            "platform_ref": {"ref": "GPL1"},
            "channel": [{"source": "source 1", "molecule": "total RNA"}],
            "library_strategy": "RNA-Seq",
            "library_source": "transcriptomic",
            "library_selection": "cDNA",
        }

        sdrf = Parent()._miniml2sdrf(base_data(sample))

        self.assertNotIn(None, [cell for row in sdrf for cell in row])
        self.assertIn("Protocol REF", sdrf[0])
        protocol_values = [
            sdrf[1][index]
            for index, label in enumerate(sdrf[0])
            if label == "Protocol REF"
        ]
        self.assertIn("", protocol_values)

    def test_source_name_uses_gsm_accession_and_preserves_geo_source(self):
        sample = {
            "iid": "GSM1",
            "title": "sample 1",
            "accession": [{"value": "GSM1"}],
            "platform_ref": {"ref": "GPL1"},
            "channel": [{"source": "whole embryo"}],
        }

        sdrf = Parent()._miniml2sdrf(base_data(sample, platform_technology="other"))

        self.assertEqual("GSM1", self.cell(sdrf, 1, "Source Name"))
        self.assertEqual("whole embryo", self.cell(sdrf, 1, "Comment[Sample_source_name]"))

    def test_organism_part_is_required_and_can_be_blank(self):
        sample = {
            "iid": "GSM1",
            "title": "sample 1",
            "accession": [{"value": "GSM1"}],
            "platform_ref": {"ref": "GPL1"},
            "channel": [{}],
        }

        sdrf = Parent()._miniml2sdrf(base_data(sample, platform_technology="other"))

        self.assertIn("Characteristics[organism part]", sdrf[0])
        self.assertEqual("", self.cell(sdrf, 1, "Characteristics[organism part]"))

    def test_requested_characteristics_are_required_and_can_be_blank(self):
        sample = {
            "iid": "GSM1",
            "title": "sample 1",
            "accession": [{"value": "GSM1"}],
            "platform_ref": {"ref": "GPL1"},
            "channel": [{}],
        }

        sdrf = Parent()._miniml2sdrf(base_data(sample, platform_technology="other"))

        for label in (
            "Characteristics[organism]",
            "Characteristics[developmental stage]",
            "Characteristics[disease]",
            "Characteristics[genotype]",
        ):
            self.assertIn(label, sdrf[0])
            self.assertEqual("", self.cell(sdrf, 1, label))

    def test_requested_characteristics_use_json_values(self):
        sample = {
            "iid": "GSM1",
            "title": "sample 1",
            "accession": [{"value": "GSM1"}],
            "platform_ref": {"ref": "GPL1"},
            "channel": [
                {
                    "organism": [{"name": "Homo sapiens"}],
                    "characteristics": [
                        {"tag": "developmental stage", "value": "adult"},
                        {"tag": "disease", "value": "normal"},
                        {"tag": "genotype", "value": "wild type"},
                    ],
                }
            ],
        }

        sdrf = Parent()._miniml2sdrf(base_data(sample, platform_technology="other"))

        self.assertEqual("Homo sapiens", self.cell(sdrf, 1, "Characteristics[organism]"))
        self.assertEqual("adult", self.cell(sdrf, 1, "Characteristics[developmental stage]"))
        self.assertEqual("normal", self.cell(sdrf, 1, "Characteristics[disease]"))
        self.assertEqual("wild type", self.cell(sdrf, 1, "Characteristics[genotype]"))

    def test_organism_part_uses_explicit_value_before_fallback(self):
        sample = {
            "iid": "GSM1",
            "title": "sample 1",
            "accession": [{"value": "GSM1"}],
            "platform_ref": {"ref": "GPL1"},
            "channel": [
                {
                    "source": "whole embryo",
                    "characteristics": [
                        {"tag": "organism part", "value": "wing disc"},
                        {"tag": "tissue", "value": "embryo"},
                    ],
                }
            ],
        }

        sdrf = Parent()._miniml2sdrf(base_data(sample, platform_technology="other"))

        self.assertEqual("wing disc", self.cell(sdrf, 1, "Characteristics[organism part]"))

    def test_greedy_geo_fallback_comments_are_not_emitted(self):
        sample = {
            "iid": "GSM1",
            "title": "sample 1",
            "accession": [{"value": "GSM1"}],
            "status": [{"release_date": "2024-01-01"}],
            "platform_ref": {"ref": "GPL1"},
            "supplementary_data": [
                {
                    "value": "counts.tsv",
                    "type": "TSV",
                    "checksum": "abc",
                    "extra_attributes": {"custom_attr": "sample extra"},
                }
            ],
            "channel": [
                {
                    "source": "source 1",
                    "extra_fields": [{"tag": "custom_field", "text": "channel extra"}],
                }
            ],
        }
        data = base_data(sample, platform_technology="other")
        data["platform"][0]["manufacturer"] = "Acme"

        sdrf = Parent()._miniml2sdrf(data)

        self.assertEqual("sample 1", self.cell(sdrf, 1, "Comment[Sample_title]"))
        self.assertEqual("source 1", self.cell(sdrf, 1, "Comment[Sample_source_name]"))
        self.assertFalse(any(header.startswith("Comment[GEO_") for header in sdrf[0]))

    def test_repeated_characteristics_render_as_repeated_columns(self):
        sample = {
            "iid": "GSM1",
            "title": "sample 1",
            "accession": [{"value": "GSM1"}],
            "platform_ref": {"ref": "GPL1"},
            "channel": [
                {
                    "source": "source 1",
                    "characteristics": [
                        {"tag": "disease", "value": "case"},
                        {"tag": "disease", "value": "treated"},
                    ],
                }
            ],
        }

        sdrf = Parent()._miniml2sdrf(base_data(sample, platform_technology="other"))

        self.assertEqual(2, sdrf[0].count("Characteristics[disease]"))
        indices = [i for i, label in enumerate(sdrf[0]) if label == "Characteristics[disease]"]
        self.assertEqual(["case", "treated"], [sdrf[1][i] for i in indices])

    def test_greedy_geo_fallback_does_not_emit_channel_metadata(self):
        sample = {
            "iid": "GSM1",
            "title": "sample 1",
            "accession": [{"value": "GSM1"}],
            "platform_ref": {"ref": "GPL1"},
            "channel": [
                {
                    "source": "source 1",
                    "extract_protocol": "extract protocol",
                    "characteristics": [
                        {"tag": "measurement", "value": "42"},
                    ],
                    "extra_fields": [{"tag": "custom_field", "text": "channel extra"}],
                }
            ],
        }

        sdrf = Parent()._miniml2sdrf(base_data(sample, platform_technology="other"))

        self.assertEqual("42", self.cell(sdrf, 1, "Characteristics[measurement]"))
        self.assertEqual("source 1", self.cell(sdrf, 1, "Comment[Sample_source_name]"))
        self.assertFalse(any(header.startswith("Comment[GEO_") for header in sdrf[0]))
        self.assertNotIn("Comment[GEO_sample_channel_source]", sdrf[0])
        self.assertNotIn("Comment[GEO_sample_channel_extract_protocol]", sdrf[0])
        self.assertNotIn("Comment[GEO_sample_channel_characteristics_tag]", sdrf[0])
        self.assertNotIn("Comment[GEO_sample_channel_characteristics_value]", sdrf[0])
        self.assertNotIn("Comment[GEO_sample_channel_extra_fields_text]", sdrf[0])

    def test_two_channel_array_emits_two_paths_and_labels(self):
        sample = {
            "iid": "GSM1",
            "title": "sample 1",
            "accession": [{"value": "GSM1"}],
            "platform_ref": {"ref": "GPL1"},
            "hybridization_protocol": "hyb",
            "channel": [
                {
                    "source": "source A",
                    "label": "Cy3",
                    "label_protocol": "label A",
                    "extract_protocol": "extract A",
                    "characteristics": [{"tag": "condition", "value": "A"}],
                },
                {
                    "source": "source B",
                    "label": "Cy5",
                    "label_protocol": "label B",
                    "extract_protocol": "extract B",
                    "characteristics": [{"tag": "condition", "value": "B"}],
                },
            ],
            "supplementary_data": [{"value": "file.CEL"}],
        }

        sdrf = Parent()._miniml2sdrf(base_data(sample, platform_technology="expression array"))

        self.assertEqual(3, len(sdrf))
        label_index = sdrf[0].index("Label")
        self.assertEqual(["Cy3", "Cy5"], [sdrf[1][label_index], sdrf[2][label_index]])

    def test_array_handler_omits_derived_files_but_keeps_raw_array_files(self):
        sample = {
            "iid": "GSM1",
            "title": "sample 1",
            "accession": [{"value": "GSM1"}],
            "platform_ref": {"ref": "GPL1"},
            "hybridization_protocol": "hyb",
            "channel": [{"source": "source 1", "extract_protocol": "extract"}],
            "supplementary_data": [
                {"value": "raw.CEL"},
                {"value": "counts.tsv.gz"},
                {"value": "archive.zip"},
            ],
        }

        sdrf = Parent()._miniml2sdrf(base_data(sample, platform_technology="expression array"))

        self.assertIn("Array Data File", sdrf[0])
        self.assertEqual("raw.CEL", self.cell(sdrf, 1, "Array Data File"))
        self.assertNotIn("Derived Array Data Matrix File", sdrf[0])
        self.assertNotIn("Derived Array Data File", sdrf[0])
        self.assertNotIn("counts.tsv.gz", sdrf[1])
        self.assertNotIn("archive.zip", sdrf[1])

    def test_protocol_registry_reuses_identical_text(self):
        registry = ProtocolRegistry(series_accession="GSE1")

        self.assertEqual("P-GSE1-1", registry.get_ref(kind="extract", text="same"))
        self.assertEqual("P-GSE1-1", registry.get_ref(kind="extract", text="same"))
        self.assertEqual("P-GSE1-2", registry.get_ref(kind="extract", text="different"))
        self.assertEqual("P-GSE1-3", registry.get_ref(kind="label", text="same"))
        self.assertIsNone(registry.get_ref(kind="extract", text=""))
        self.assertEqual(
            [
                {
                    "ref": "P-GSE1-1",
                    "kind": "extract",
                    "label": "Extract-Protocol",
                    "text": "same",
                },
                {
                    "ref": "P-GSE1-2",
                    "kind": "extract",
                    "label": "Extract-Protocol",
                    "text": "different",
                },
                {
                    "ref": "P-GSE1-3",
                    "kind": "label",
                    "label": "Label-Protocol",
                    "text": "same",
                },
            ],
            registry.records(),
        )

    def test_protocol_registry_creates_and_reuses_required_placeholders(self):
        registry = ProtocolRegistry(series_accession="GSE1")

        self.assertIsNone(registry.get_ref(kind="sample collection", text=""))
        self.assertEqual(
            "P-GSE1-1",
            registry.ensure_required(
                kind="sample collection",
                label="Sample-Collection-Protocol",
            ),
        )
        self.assertEqual(
            "P-GSE1-1",
            registry.ensure_required(
                kind="sample collection",
                label="Sample-Collection-Protocol",
            ),
        )
        self.assertEqual(
            [
                {
                    "ref": "P-GSE1-1",
                    "kind": "sample collection",
                    "label": "Sample-Collection-Protocol",
                    "text": "",
                    "required": True,
                },
            ],
            registry.records(),
        )

    def test_more_than_two_fastqs_are_preserved_and_derived_files_are_omitted(self):
        sample = {
            "iid": "GSM1",
            "title": "sample 1",
            "accession": [{"value": "GSM1"}],
            "platform_ref": {"ref": "GPL1"},
            "channel": [{"source": "source 1", "extract_protocol": "extract"}],
            "library_strategy": "RNA-Seq",
            "library_source": "transcriptomic",
            "library_selection": "cDNA",
            "supplementary_data": [
                {"value": "counts1.tsv.gz"},
                {"value": "counts2.mtx.gz"},
                {"value": "other.h5ad"},
            ],
        }
        handler = _SequencingSDRFHandler(parent=Parent(), data=base_data(sample))
        handler.sra_runs = lambda sample: [
            {
                "geo_sample": "GSM1",
                "experiment": "SRX1",
                "run": "SRR1",
                "scan_name": "SRR1",
                "fastq_files": [
                    {"filename": "r1.fastq.gz", "uri": "ftp://example/r1.fastq.gz"},
                    {"filename": "r2.fastq.gz", "uri": "ftp://example/r2.fastq.gz"},
                    {"filename": "i1.fastq.gz", "uri": "ftp://example/i1.fastq.gz"},
                ],
            }
        ]

        sdrf = handler.build()

        self.assertNotIn("Array Data File", sdrf[0])
        self.assertNotIn("Derived Array Data Matrix File", sdrf[0])
        self.assertNotIn("Derived Array Data File", sdrf[0])
        self.assertIn("Comment[read1 file]", sdrf[0])
        self.assertIn("Comment[read2 file]", sdrf[0])
        self.assertIn("Comment[read3 file]", sdrf[0])
        self.assertNotIn("Comment[derived data file]", sdrf[0])
        self.assertIn("r1.fastq.gz", sdrf[1])
        self.assertIn("r2.fastq.gz", sdrf[1])
        self.assertIn("i1.fastq.gz", sdrf[1])
        self.assertNotIn("counts1.tsv.gz", sdrf[1])
        self.assertNotIn("counts2.mtx.gz", sdrf[1])
        self.assertNotIn("other.h5ad", sdrf[1])

    def test_sdrf_uses_enriched_sra_runs_without_fetching(self):
        fetcher = Mock()
        sample = {
            "iid": "GSM1",
            "title": "sample 1",
            "accession": [{"value": "GSM1"}],
            "platform_ref": {"ref": "GPL1"},
            "relation": [{"type": "SRA", "target": "SRX1"}],
            "channel": [{"source": "source 1", "extract_protocol": "extract"}],
            "sra_run": [
                {
                    "geo_sample": "GSM1",
                    "experiment": "SRX1",
                    "run": "SRR1",
                    "scan_name": "SRR1",
                    "library_layout": "PAIRED",
                    "library_strategy": "RNA-Seq",
                    "library_source": "TRANSCRIPTOMIC",
                    "library_selection": "cDNA",
                    "fastq_files": [{"filename": "r1.fastq.gz", "uri": "ftp://example/r1.fastq.gz"}],
                }
            ],
        }

        sdrf = SDRFConstructor(insdc_fetcher=fetcher)._miniml2sdrf(base_data(sample))

        self.assertEqual("SRR1", self.cell(sdrf, 1, "Comment[ENA_RUN]"))
        self.assertEqual("ftp://example/r1.fastq.gz", self.cell(sdrf, 1, "Comment[FASTQ_URI]"))
        fetcher._extract_sra.assert_not_called()
        fetcher.fetch_sra_runs.assert_not_called()

    def test_bulk_sequencing_handler_emits_one_row_per_fastq(self):
        sample = {
            "iid": "GSM1",
            "title": "bulk RNA-seq sample",
            "accession": [{"value": "GSM1"}],
            "platform_ref": {"ref": "GPL1"},
            "channel": [
                {
                    "source": "source 1",
                    "extract_protocol": "extract",
                    "characteristics": [{"tag": "condition", "value": "treated"}],
                }
            ],
            "library_strategy": "RNA-Seq",
            "library_source": "transcriptomic",
            "library_selection": "cDNA",
            "supplementary_data": [
                {"value": "counts1.tsv.gz"},
                {"value": "counts2.mtx.gz"},
            ],
        }
        data = base_data(sample)
        data["series"]["variable"] = [{"name": "condition"}]
        sra_runs = lambda sample: [
            {
                "geo_sample": "GSM1",
                "experiment": "SRX1",
                "run": "SRR1",
                "scan_name": "SRR1",
                "library_layout": "PAIRED",
                "library_source": "transcriptomic",
                "library_strategy": "RNA-Seq",
                "library_selection": "cDNA",
                "fastq_files": [
                    {
                        "filename": "r1.fastq.gz",
                        "uri": "ftp://example/r1.fastq.gz",
                        "md5": "md5-r1",
                    },
                    {
                        "filename": "r2.fastq.gz",
                        "md5": "md5-r2",
                    },
                    {
                        "filename": "i1.fastq.gz",
                        "uri": "ftp://example/i1.fastq.gz",
                        "md5": "md5-i1",
                    },
                ],
            }
        ]
        bulk_handler = _BulkSequencingSDRFHandler(parent=Parent(), data=data)
        bulk_handler.sra_runs = sra_runs

        sdrf = bulk_handler.build()

        self.assertEqual(4, len(sdrf))
        self.assertIn("Comment[FASTQ_URI]", sdrf[0])
        self.assertIn("Comment[MD5]", sdrf[0])
        self.assertNotIn("Comment[read1 file]", sdrf[0])
        self.assertNotIn("Comment[read2 file]", sdrf[0])
        self.assertNotIn("Comment[read3 file]", sdrf[0])
        self.assertEqual(
            [
                "ftp://example/r1.fastq.gz",
                "r2.fastq.gz",
                "ftp://example/i1.fastq.gz",
            ],
            [self.cell(sdrf, index, "Comment[FASTQ_URI]") for index in range(1, 4)],
        )
        self.assertEqual(
            ["md5-r1", "md5-r2", "md5-i1"],
            [self.cell(sdrf, index, "Comment[MD5]") for index in range(1, 4)],
        )
        self.assertEqual(["GSM1", "GSM1", "GSM1"], [self.cell(sdrf, index, "Source Name") for index in range(1, 4)])
        self.assertEqual(["SRR1", "SRR1", "SRR1"], [self.cell(sdrf, index, "Comment[ENA_RUN]") for index in range(1, 4)])
        self.assertEqual(["treated", "treated", "treated"], [self.cell(sdrf, index, "Factor Value[condition]") for index in range(1, 4)])
        self.assertNotIn("Comment[derived data file]", sdrf[0])
        for index in range(1, 4):
            self.assertNotIn("counts1.tsv.gz", sdrf[index])
            self.assertNotIn("counts2.mtx.gz", sdrf[index])

    def test_bulk_sequencing_inheritance_tree(self):
        self.assertTrue(issubclass(_BulkSequencingSDRFHandler, _SequencingSDRFHandler))
        self.assertTrue(issubclass(_PlateSingleCellSequencingSDRFHandler, _BulkSequencingSDRFHandler))

    def test_tenx_version_handlers_inherit_droplet_single_cell_handler(self):
        self.assertTrue(issubclass(
            _TenXV2DropletSingleCellSequencingSDRFHandler,
            _DropletSingleCellSequencingSDRFHandler,
        ))
        self.assertTrue(issubclass(
            _TenXV3DropletSingleCellSequencingSDRFHandler,
            _DropletSingleCellSequencingSDRFHandler,
        ))

    def test_tenx_version_handler_dispatch(self):
        sample = {
            "iid": "GSM1",
            "title": "single-cell sample",
            "description": "10x Chromium",
            "accession": [{"value": "GSM1"}],
            "platform_ref": {"ref": "GPL1"},
            "channel": [{"source": "source 1", "extract_protocol": "extract"}],
        }

        for tech_type, handler_class in (
            (
                "tenx_v2_droplet_single_cell_sequencing",
                _TenXV2DropletSingleCellSequencingSDRFHandler,
            ),
            (
                "tenx_v3_droplet_single_cell_sequencing",
                _TenXV3DropletSingleCellSequencingSDRFHandler,
            ),
        ):
            with self.subTest(tech_type=tech_type):
                parent = FixedTechParent(tech_type)
                with patch.object(
                    handler_class,
                    "build",
                    return_value=[["Source Name"], ["GSM1"]],
                ) as build:
                    sdrf = parent._miniml2sdrf(base_data(sample))

                build.assert_called_once_with()
                self.assertEqual([["Source Name"], ["GSM1"]], sdrf)

    def test_tenx_v2_handler_emits_fixed_library_attributes(self):
        sample = {
            "iid": "GSM1",
            "title": "10x v2 sample",
            "description": "10x Chromium single cell v2",
            "accession": [{"value": "GSM1"}],
            "platform_ref": {"ref": "GPL1"},
            "channel": [{"source": "source 1", "extract_protocol": "extract"}],
        }

        sdrf = FixedTechParent("tenx_v2_droplet_single_cell_sequencing")._miniml2sdrf(base_data(sample))

        expected = {
            "Comment[cdna read]": "read2",
            "Comment[cdna read offset]": "0",
            "Comment[cdna read size]": "98",
            "Comment[cell barcode offset]": "0",
            "Comment[cell barcode read]": "read1",
            "Comment[cell barcode size]": "16",
            "Comment[end bias]": "3 prime tag",
            "Comment[input molecule]": "polyA RNA",
            "Comment[library construction]": "10xV2",
            "Comment[primer]": "oligo-dT",
            "Comment[LIBRARY_STRAND]": "not applicable",
            "Comment[sample barcode offset]": "0",
            "Comment[sample barcode read]": "index1",
            "Comment[sample barcode size]": "8",
            "Comment[single cell isolation]": "10x technology",
            "Comment[spike in]": "",
            "Comment[umi barcode offset]": "16",
            "Comment[umi barcode read]": "read1",
            "Comment[umi barcode size]": "10",
        }
        for label, value in expected.items():
            with self.subTest(label=label):
                self.assertEqual(value, self.cell(sdrf, 1, label))

    def test_tenx_v3_handler_emits_fixed_library_attributes(self):
        sample = {
            "iid": "GSM1",
            "title": "10x v3 sample",
            "description": "10x Chromium single cell v3",
            "accession": [{"value": "GSM1"}],
            "platform_ref": {"ref": "GPL1"},
            "channel": [{"source": "source 1", "extract_protocol": "extract"}],
        }

        sdrf = FixedTechParent("tenx_v3_droplet_single_cell_sequencing")._miniml2sdrf(base_data(sample))

        expected = {
            "Comment[cdna read]": "read2",
            "Comment[cdna read offset]": "0",
            "Comment[cdna read size]": "91",
            "Comment[cell barcode offset]": "0",
            "Comment[cell barcode read]": "read1",
            "Comment[cell barcode size]": "16",
            "Comment[end bias]": "3 prime tag",
            "Comment[input molecule]": "polyA RNA",
            "Comment[library construction]": "10xV3",
            "Comment[primer]": "oligo-dT",
            "Comment[LIBRARY_STRAND]": "not applicable",
            "Comment[sample barcode offset]": "0",
            "Comment[sample barcode read]": "index1",
            "Comment[sample barcode size]": "8",
            "Comment[single cell isolation]": "10x technology",
            "Comment[spike in]": "",
            "Comment[umi barcode offset]": "16",
            "Comment[umi barcode read]": "read1",
            "Comment[umi barcode size]": "12",
        }
        for label, value in expected.items():
            with self.subTest(label=label):
                self.assertEqual(value, self.cell(sdrf, 1, label))

    def test_sequencing_subclasses_inherit_file_comments(self):
        sample = {
            "iid": "GSM1",
            "title": "single-cell sample",
            "description": "10x Chromium",
            "accession": [{"value": "GSM1"}],
            "platform_ref": {"ref": "GPL1"},
            "channel": [{"source": "source 1", "extract_protocol": "extract"}],
            "library_strategy": "RNA-Seq",
            "raw_data": [{"value": "ftp://example/read1.fastq.gz"}],
            "supplementary_data": [{"value": "counts.tsv.gz"}],
        }

        for tech_type in (
            "single_cell_sequencing",
            "droplet_single_cell_sequencing",
            "spatial_sequencing",
        ):
            sdrf = FixedTechParent(tech_type)._miniml2sdrf(base_data(sample))
            self.assertNotIn("Array Data File", sdrf[0])
            self.assertNotIn("Derived Array Data Matrix File", sdrf[0])
            self.assertEqual("ftp://example/read1.fastq.gz", self.cell(sdrf, 1, "Comment[read1 file]"))
            self.assertNotIn("Comment[derived data file]", sdrf[0])
            self.assertNotIn("counts.tsv.gz", sdrf[1])

        plate_sdrf = FixedTechParent("plate_single_cell_sequencing")._miniml2sdrf(base_data(sample))
        self.assertNotIn("Array Data File", plate_sdrf[0])
        self.assertNotIn("Derived Array Data Matrix File", plate_sdrf[0])
        self.assertNotIn("Comment[read1 file]", plate_sdrf[0])
        self.assertEqual("ftp://example/read1.fastq.gz", self.cell(plate_sdrf, 1, "Comment[FASTQ_URI]"))
        self.assertNotIn("Comment[derived data file]", plate_sdrf[0])
        self.assertNotIn("counts.tsv.gz", plate_sdrf[1])

    def test_greedy_sra_fallback_comments_are_not_emitted(self):
        sample = {
            "iid": "GSM1",
            "title": "sample 1",
            "accession": [{"value": "GSM1"}],
            "platform_ref": {"ref": "GPL1"},
            "channel": [{"source": "source 1", "extract_protocol": "extract"}],
            "library_strategy": "RNA-Seq",
        }
        handler = _SequencingSDRFHandler(parent=Parent(), data=base_data(sample))
        handler.sra_runs = lambda sample: [
            {
                "geo_sample": "GSM1",
                "experiment": "SRX1",
                "run": "SRR1",
                "scan_name": "SRR1",
                "submitted_file_name": "submitted.bam",
                "md5": "run-md5",
                "library_adapter": "custom adapter",
                "library_layout": "PAIRED",
                "run_alias": "lane 1",
                "fastq_files": [
                    {
                        "filename": "r1.fastq.gz",
                        "uri": "ftp://example/r1.fastq.gz",
                        "md5": "fastq-md5",
                        "lane": "1",
                    }
                ],
            }
        ]

        sdrf = handler.build()

        self.assertEqual("PAIRED", self.cell(sdrf, 1, "Comment[LIBRARY_LAYOUT]"))
        self.assertEqual("SRX1", self.cell(sdrf, 1, "Comment[ENA_EXPERIMENT]"))
        self.assertEqual("SRR1", self.cell(sdrf, 1, "Comment[ENA_RUN]"))
        self.assertEqual("submitted.bam", self.cell(sdrf, 1, "Comment[SUBMITTED_FILE_NAME]"))
        self.assertIn("run-md5", sdrf[1])
        self.assertEqual("r1.fastq.gz", self.cell(sdrf, 1, "Comment[read1 file]"))
        self.assertEqual("ftp://example/r1.fastq.gz", self.cell(sdrf, 1, "Comment[FASTQ_URI]"))
        self.assertIn("fastq-md5", sdrf[1])
        self.assertFalse(any(header.startswith("Comment[SRA_library_") for header in sdrf[0]))
        self.assertFalse(any(header.startswith("Comment[SRA_run_") for header in sdrf[0]))
        self.assertFalse(any(header.startswith("Comment[SRA_fastq_") for header in sdrf[0]))

    def test_geo_values_win_over_conflicting_sra_values(self):
        sample = {
            "iid": "GSM1",
            "title": "sample 1",
            "accession": [{"value": "GSM1"}],
            "platform_ref": {"ref": "GPL1"},
            "channel": [{"source": "source 1", "extract_protocol": "extract"}],
            "library_strategy": "RNA-Seq",
            "library_source": "transcriptomic",
            "instrument_model": {"predefined": "GeoSeq 1"},
        }
        handler = _SequencingSDRFHandler(parent=Parent(), data=base_data(sample))
        handler.sra_runs = lambda sample: [
            {
                "geo_sample": "GSM2",
                "experiment": "SRX1",
                "run": "SRR1",
                "scan_name": "SRR1",
                "library_strategy": "ChIP-Seq",
                "library_source": "genomic",
                "instrument_model": "SraSeq 2",
                "fastq_files": [],
            }
        ]

        sdrf = handler.build()

        self.assertEqual("RNA-Seq", self.cell(sdrf, 1, "Comment[LIBRARY_STRATEGY]"))
        self.assertEqual("TRANSCRIPTOMIC", self.cell(sdrf, 1, "Comment[LIBRARY_SOURCE]"))
        self.assertEqual("GeoSeq 1", self.cell(sdrf, 1, "Comment[INSTRUMENT_MODEL]"))
        self.assertEqual("GSM1", self.cell(sdrf, 1, "Assay Name"))
        self.assertTrue(any("using GEO value" in warning for warning in handler.audit.warnings))
        self.assertTrue(any("using GEO accession" in warning for warning in handler.audit.warnings))

    def test_sra_library_source_fallback_is_uppercase(self):
        sample = {
            "iid": "GSM1",
            "title": "sample 1",
            "accession": [{"value": "GSM1"}],
            "platform_ref": {"ref": "GPL1"},
            "channel": [{"source": "source 1", "extract_protocol": "extract"}],
            "library_strategy": "RNA-Seq",
        }
        handler = _SequencingSDRFHandler(parent=Parent(), data=base_data(sample))
        handler.sra_runs = lambda sample: [
            {
                "geo_sample": "GSM1",
                "experiment": "SRX1",
                "run": "SRR1",
                "scan_name": "SRR1",
                "library_source": "metagenomic",
                "fastq_files": [],
            }
        ]

        sdrf = handler.build()

        self.assertEqual("METAGENOMIC", self.cell(sdrf, 1, "Comment[LIBRARY_SOURCE]"))

    def test_generated_sdrf_comment_labels_use_normalized_bracket_style(self):
        sample = {
            "iid": "GSM1",
            "title": "sample 1",
            "description": "description 1",
            "accession": [{"value": "GSM1"}],
            "platform_ref": {"ref": "GPL1"},
            "channel": [{"source": "source 1"}],
            "supplementary_data": [
                {"type": "CEL", "text": "ftp://example/raw.CEL.gz"},
                {"type": "TXT", "text": "ftp://example/counts.txt.gz"},
            ],
            "sra_run": [
                {
                    "sample": "SRS1",
                    "experiment": "SRX1",
                    "run": "SRR1",
                    "instrument_model": "GeoSeq 1",
                    "fastq_files": [{"filename": "r1.fastq.gz", "uri": "ftp://example/r1.fastq.gz"}],
                },
            ],
        }

        sdrf = Parent()._miniml2sdrf(base_data(sample))

        for label in sdrf[0]:
            if label.startswith("Comment"):
                self.assertIsNone(re.match(r"^Comment\s+\[", label))
                self.assertRegex(label, r"^Comment\[[^]]+\]$")

    def test_direct_sdrf_rendering_strips_quotes_from_values(self):
        sample = {
            "iid": "GSM'1",
            "title": "\"sample\" 1",
            "accession": [{"value": "GSM'1"}],
            "platform_ref": {"ref": "GPL1"},
            "channel": [
                {
                    "source": "John's \"sample\"",
                    "characteristics": [{"tag": "condition", "value": "\"treated\""}],
                }
            ],
        }

        sdrf = Parent()._miniml2sdrf(base_data(sample, platform_technology="other"))

        self.assertEqual("GSM1", self.cell(sdrf, 1, "Source Name"))
        self.assertEqual("Johns sample", self.cell(sdrf, 1, "Comment[Sample_source_name]"))
        self.assertEqual("Johns sample", self.cell(sdrf, 1, "Characteristics[organism part]"))
        self.assertEqual("treated", self.cell(sdrf, 1, "Characteristics[condition]"))

    def test_expanded_file_classification(self):
        self.assertEqual("sequencing_raw", classify_file("reads.cram"))
        self.assertEqual("array_raw", classify_file("scan.tiff.gz"))
        self.assertEqual("array_raw", classify_file("raw.exp"))
        self.assertEqual("array_raw", classify_file("raw.rpt"))
        self.assertEqual("array_raw", classify_file("raw.cab"))

    def test_chip_assay_uses_standard_sequencing_handler(self):
        sample = {
            "iid": "GSM1",
            "title": "chip-seq sample",
            "accession": [{"value": "GSM1"}],
            "platform_ref": {"ref": "GPL1"},
            "channel": [
                {
                    "source": "source 1",
                    "extract_protocol": "extract",
                    "characteristics": [{"tag": "chip antibody", "value": "H3K27ac"}],
                }
            ],
            "library_strategy": "ChIP-Seq",
            "library_source": "genomic",
            "library_selection": "ChIP",
        }

        sdrf = Parent()._miniml2sdrf(base_data(sample))

        self.assertIn("Assay Name", sdrf[0])
        self.assertEqual("sequencing assay", self.cell(sdrf, 1, "Technology Type"))
        self.assertNotIn("Factor Value[chip antibody]", sdrf[0])

    def test_assay_terms_do_not_create_sequencing_tech_types(self):
        for title in (
            "ChIP-Seq sample",
            "methylation bisulfite sample",
            "small RNA miRNA non coding sample",
            "multiome sample",
            "scATAC sample",
            "sc-ATAC sample",
            "ATAC-seq sample",
        ):
            sample = {
                "iid": "GSM1",
                "title": title,
                "accession": [{"value": "GSM1"}],
                "platform_ref": {"ref": "GPL1"},
                "channel": [{"source": "source 1", "extract_protocol": "extract"}],
                "library_strategy": title,
            }

            self.assertEqual("bulk_sequencing", Parent()._detect_sdrf_technology(base_data(sample)))

    def test_assay_terms_do_not_create_array_tech_types(self):
        for title in (
            "ChIP-chip array sample",
            "genotyping SNP array sample",
            "comparative genomic hybridization CGH copy number sample",
            "methylation bisulfite array sample",
        ):
            sample = {
                "iid": "GSM1",
                "title": title,
                "accession": [{"value": "GSM1"}],
                "platform_ref": {"ref": "GPL1"},
                "channel": [{"source": "source 1", "extract_protocol": "extract"}],
            }

            self.assertEqual("array", Parent()._detect_sdrf_technology(base_data(sample, platform_technology="expression array")))

    def test_droplet_single_cell_detection_and_output(self):
        for text in ("10x Chromium", "droplet single-cell", "Chromium single cell"):
            sample = {
                "iid": "GSM1",
                "title": "single-cell sample",
                "description": text,
                "accession": [{"value": "GSM1"}],
                "platform_ref": {"ref": "GPL1"},
                "channel": [{"source": "source 1", "extract_protocol": "extract"}],
            }
            self.assertEqual("droplet_single_cell_sequencing", Parent()._detect_sdrf_technology(base_data(sample)))

        sample = {
            "iid": "GSM1",
            "title": "10x v3 sample",
            "description": "10x Chromium single cell v3",
            "accession": [{"value": "GSM1"}],
            "platform_ref": {"ref": "GPL1"},
            "channel": [{"source": "source 1", "extract_protocol": "extract"}],
        }

        sdrf = Parent()._miniml2sdrf(base_data(sample))

        self.assertIn("Comment[cell barcode read]", sdrf[0])
        self.assertEqual("read1", self.cell(sdrf, 1, "Comment[cell barcode read]"))
        self.assertEqual("10xV3", self.cell(sdrf, 1, "Comment[library construction]"))
        self.assertEqual("10x technology", self.cell(sdrf, 1, "Comment[single cell isolation]"))

    def test_generic_single_cell_does_not_emit_droplet_read_geometry(self):
        sample = {
            "iid": "GSM1",
            "title": "10x sample",
            "description": "10x Chromium",
            "accession": [{"value": "GSM1"}],
            "platform_ref": {"ref": "GPL1"},
            "channel": [{"source": "source 1", "extract_protocol": "extract"}],
        }

        sdrf = FixedTechParent("single_cell_sequencing")._miniml2sdrf(base_data(sample))

        self.assertEqual("sequencing assay", self.cell(sdrf, 1, "Technology Type"))
        self.assertNotIn("Comment[cell barcode read]", sdrf[0])
        self.assertNotIn("Comment[umi barcode read]", sdrf[0])

    def test_spatial_and_plate_single_cell_detection_precedence(self):
        spatial_sample = {
            "iid": "GSM1",
            "title": "Visium single cell sample",
            "description": "10x Visium spatial assay",
            "accession": [{"value": "GSM1"}],
            "platform_ref": {"ref": "GPL1"},
            "channel": [{"source": "source 1", "extract_protocol": "extract"}],
        }
        plate_sample = {
            "iid": "GSM2",
            "title": "single-cell sample",
            "description": "single-cell plate assay",
            "accession": [{"value": "GSM2"}],
            "platform_ref": {"ref": "GPL1"},
            "channel": [{"source": "source 1", "extract_protocol": "extract"}],
        }

        self.assertEqual("spatial_sequencing", Parent()._detect_sdrf_technology(base_data(spatial_sample)))
        self.assertEqual("plate_single_cell_sequencing", Parent()._detect_sdrf_technology(base_data(plate_sample)))

    def test_no_global_empty_term_source_or_unit_columns(self):
        sample = {
            "iid": "GSM1",
            "title": "sample 1",
            "accession": [{"value": "GSM1"}],
            "platform_ref": {"ref": "GPL1"},
            "channel": [
                {
                    "source": "source 1",
                    "organism": [{"name": "Homo sapiens"}],
                    "characteristics": [{"tag": "tissue", "value": "blood"}],
                }
            ],
        }

        sdrf = Parent()._miniml2sdrf(base_data(sample, platform_technology="other"))

        self.assertNotIn("Term Accession Number", sdrf[0])
        self.assertNotIn("Unit[dimensionless unit]", sdrf[0])


if __name__ == "__main__":
    unittest.main()
