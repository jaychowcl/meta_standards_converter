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
import tempfile
import unittest
from unittest.mock import ANY, Mock, patch


ROOT = os.path.dirname(os.path.dirname(__file__))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from meta_standards_converter.ae_handlers.ae_constructor import AEConstructor, ProtocolRegistry  # noqa: E402
from meta_standards_converter.ae_handlers.ae_idf_handlers import (  # noqa: E402
    IDFConstructor,
    _ArrayPlatformIDFHandler,
    _BasePlatformIDFHandler,
    _BulkSequencingPlatformIDFHandler,
    _DropletSingleCellSequencingPlatformIDFHandler,
    _GenericPlatformIDFHandler,
    _PlateSingleCellSequencingPlatformIDFHandler,
    _SequencingPlatformIDFHandler,
    _SingleCellSequencingPlatformIDFHandler,
    _SpatialSequencingPlatformIDFHandler,
)
from meta_standards_converter.ae_handlers.ae_sdrf_handlers import SDRFConstructor  # noqa: E402
from meta_standards_converter.harmonizers.harmonizers import Harmonizer  # noqa: E402
from meta_standards_converter.geo_handlers.geo_parser import GEOParser  # noqa: E402


class TestIDFConstructor(unittest.TestCase):
    def row(self, rows, label):
        return next(row for row in rows if row[0] == label)

    def test_geoprotocols2efo_returns_original_label_for_unknown_protocol(self):
        self.assertEqual(
            ["Custom-Protocol", None, None],
            Harmonizer().geoprotocols2efo(protocol_type="Custom-Protocol"),
        )

    def test_geoprotocols2efo_still_rejects_blank_protocol(self):
        with self.assertRaises(ValueError):
            Harmonizer().geoprotocols2efo(protocol_type="")

    def test_geoprotocols2efo_maps_required_protocols(self):
        self.assertEqual(
            ["sample collection protocol", "EFO", "EFO_0005518"],
            Harmonizer().geoprotocols2efo(protocol_type="Sample-Collection-Protocol"),
        )
        self.assertEqual(
            ["nucleic acid sequencing protocol", "EFO", "EFO_0004170"],
            Harmonizer().geoprotocols2efo(protocol_type="Nucleic-Acid-Sequencing-Protocol"),
        )

    def test_idf_protocols_from_registry_uses_harmonizer_fallback(self):
        registry = ProtocolRegistry(series_accession="GSE1")
        registry.get_ref(kind="custom", text="custom protocol", label="Custom-Protocol")

        rows = IDFConstructor()._idf_protocols(data={}, protocol_registry=registry)

        self.assertEqual(["Protocol Name", "P-GSE1-1", "P-GSE1-2"], self.row(rows, "Protocol Name"))
        self.assertEqual(
            ["Protocol Type", "Custom-Protocol", "sample collection protocol"],
            self.row(rows, "Protocol Type"),
        )
        self.assertEqual(
            ["Protocol Type Term Source REF", None, "EFO"],
            self.row(rows, "Protocol Type Term Source REF"),
        )
        self.assertEqual(
            ["Protocol Type Term Accession Number", None, "EFO_0005518"],
            self.row(rows, "Protocol Type Term Accession Number"),
        )
        self.assertEqual(["Protocol Description", "custom protocol", None], self.row(rows, "Protocol Description"))
        self.assertEqual(["Protocol Hardware", None, None], self.row(rows, "Protocol Hardware"))
        self.assertEqual(["Protocol Software", None, None], self.row(rows, "Protocol Software"))
        self.assertFalse(any(row[0] == "Protocol Parameters" for row in rows))
        self.assertFalse(any(row[0] == "Protocol Contact" for row in rows))

    def test_idf_protocols_adds_required_sample_collection_for_unknown_technology(self):
        rows = IDFConstructor()._idf_protocols(
            data={"series": {"accession": [{"value": "GSE1"}]}},
            technology_type="unknown",
        )

        self.assertEqual(["Protocol Name", "P-GSE1-1"], self.row(rows, "Protocol Name"))
        self.assertEqual(["Protocol Type", "sample collection protocol"], self.row(rows, "Protocol Type"))
        self.assertEqual(["Protocol Type Term Source REF", "EFO"], self.row(rows, "Protocol Type Term Source REF"))
        self.assertEqual(
            ["Protocol Type Term Accession Number", "EFO_0005518"],
            self.row(rows, "Protocol Type Term Accession Number"),
        )
        self.assertEqual(["Protocol Description", None], self.row(rows, "Protocol Description"))
        self.assertEqual(["Protocol Hardware", None], self.row(rows, "Protocol Hardware"))
        self.assertEqual(["Protocol Software", None], self.row(rows, "Protocol Software"))
        self.assertFalse(any(row[0] == "Protocol Parameters" for row in rows))
        self.assertFalse(any(row[0] == "Protocol Contact" for row in rows))

    def test_idf_protocols_adds_required_sequencing_protocol_when_missing(self):
        rows = IDFConstructor()._idf_protocols(
            data={"series": {"accession": [{"value": "GSE1"}]}},
            technology_type="bulk_sequencing",
        )

        self.assertEqual(
            ["Protocol Name", "P-GSE1-1", "P-GSE1-2"],
            self.row(rows, "Protocol Name"),
        )
        self.assertEqual(
            ["Protocol Type", "sample collection protocol", "nucleic acid sequencing protocol"],
            self.row(rows, "Protocol Type"),
        )
        self.assertEqual(
            ["Protocol Type Term Source REF", "EFO", "EFO"],
            self.row(rows, "Protocol Type Term Source REF"),
        )
        self.assertEqual(
            ["Protocol Type Term Accession Number", "EFO_0005518", "EFO_0004170"],
            self.row(rows, "Protocol Type Term Accession Number"),
        )
        self.assertEqual(["Protocol Hardware", None, None], self.row(rows, "Protocol Hardware"))
        self.assertEqual(["Protocol Software", None, None], self.row(rows, "Protocol Software"))
        self.assertFalse(any(row[0] == "Protocol Parameters" for row in rows))
        self.assertFalse(any(row[0] == "Protocol Contact" for row in rows))

    def test_idf_protocols_treats_tenx_versions_as_sequencing(self):
        for technology_type in (
            "tenx_v2_droplet_single_cell_sequencing",
            "tenx_v3_droplet_single_cell_sequencing",
        ):
            with self.subTest(technology_type=technology_type):
                rows = IDFConstructor()._idf_protocols(
                    data={"series": {"accession": [{"value": "GSE1"}]}},
                    technology_type=technology_type,
                )

                self.assertIn(
                    "nucleic acid sequencing protocol",
                    self.row(rows, "Protocol Type"),
                )

    def test_idf_protocols_does_not_duplicate_existing_required_protocol_types(self):
        registry = ProtocolRegistry(series_accession="GSE1")
        registry.get_ref(
            kind="submitted sample collection",
            text="actual sample collection",
            label="Sample-Collection-Protocol",
        )
        registry.get_ref(
            kind="submitted sequencing",
            text="actual sequencing",
            label="Nucleic-Acid-Sequencing-Protocol",
        )

        rows = IDFConstructor()._idf_protocols(
            data={},
            protocol_registry=registry,
            technology_type="sequencing",
        )

        self.assertEqual(["Protocol Name", "P-GSE1-1", "P-GSE1-2"], self.row(rows, "Protocol Name"))
        self.assertEqual(
            ["Protocol Type", "sample collection protocol", "nucleic acid sequencing protocol"],
            self.row(rows, "Protocol Type"),
        )
        self.assertEqual(
            ["Protocol Description", "actual sample collection", "actual sequencing"],
            self.row(rows, "Protocol Description"),
        )

    def test_append_required_protocol_rows_skips_missing_optional_rows(self):
        rows = [
            ["Protocol Name"],
            ["Protocol Type"],
            ["Protocol Type Term Source REF"],
            ["Protocol Type Term Accession Number"],
            ["Protocol Description"],
            ["Protocol Hardware"],
            ["Protocol Software"],
        ]

        IDFConstructor()._append_required_protocol_rows(
            rows=rows,
            series_accession="GSE1",
            technology_type="bulk_sequencing",
        )

        self.assertEqual(["Protocol Name", "P-GSE1-1", "P-GSE1-2"], self.row(rows, "Protocol Name"))
        self.assertEqual(
            ["Protocol Type", "sample collection protocol", "nucleic acid sequencing protocol"],
            self.row(rows, "Protocol Type"),
        )
        self.assertEqual(["Protocol Hardware", None, None], self.row(rows, "Protocol Hardware"))
        self.assertEqual(["Protocol Software", None, None], self.row(rows, "Protocol Software"))
        self.assertFalse(any(row[0] == "Protocol Parameters" for row in rows))
        self.assertFalse(any(row[0] == "Protocol Contact" for row in rows))

    def test_miniml2idf_starts_with_magetab_version_row(self):
        constructor = IDFConstructor()

        with patch.object(constructor, "_idf_investigations", return_value=[]), \
             patch.object(constructor, "_idf_experimental", return_value=[]), \
             patch.object(constructor, "_idf_platform_specific", return_value=[]), \
             patch.object(constructor, "_idf_persons", return_value=[]), \
             patch.object(constructor, "_idf_qc_rep_norm", return_value=[]), \
             patch.object(constructor, "_idf_dates", return_value=[]), \
             patch.object(constructor, "_idf_publications", return_value=[]), \
             patch.object(constructor, "_idf_experiments", return_value=[]), \
             patch.object(constructor, "_idf_protocols", return_value=[]), \
             patch.object(constructor, "_idf_term_source", return_value=[]):
            rows = constructor.miniml2idf(data={})

        self.assertEqual(["MAGE-TAB Version", "1.1"], rows[0])

    def test_miniml2idf_accepts_missing_technology_type_for_direct_calls(self):
        rows = IDFConstructor().miniml2idf(data={"series": [], "sample": [], "contributor": []})

        self.assertEqual(["MAGE-TAB Version", "1.1"], rows[0])
        self.assertIn(["SDRF File"], rows)

    def test_idf_investigations_emits_related_experiment_for_superseries_relation(self):
        rows = IDFConstructor()._idf_investigations(
            data={
                "series": {
                    "title": "Example",
                    "accession": [{"value": "GSE1", "database": "GEO"}],
                    "relation": [{"type": "SuperSeries of", "target": "GSE2"}],
                },
            }
        )

        self.assertEqual(["Comment[RelatedExperiment]", "GSE2"], self.row(rows, "Comment[RelatedExperiment]"))

    def test_idf_investigations_emits_related_experiment_for_subseries_relation(self):
        rows = IDFConstructor()._idf_investigations(
            data={
                "series": {
                    "title": "Example",
                    "accession": [{"value": "GSE2", "database": "GEO"}],
                    "relation": [{"type": "SubSeries of", "target": "GSE1"}],
                },
            }
        )

        self.assertEqual(["Comment[RelatedExperiment]", "GSE1"], self.row(rows, "Comment[RelatedExperiment]"))

    def test_idf_investigations_dedupes_related_experiments_in_encounter_order(self):
        rows = IDFConstructor()._idf_investigations(
            data={
                "series": {
                    "title": "Example",
                    "accession": [{"value": "GSE1", "database": "GEO"}],
                    "relation": [
                        {
                            "type": "SuperSeries of GSE3",
                            "target": "gse2",
                            "comment": "also includes GSE3",
                        },
                        {
                            "type": "SubSeries of",
                            "target": "GSE4",
                            "comment": "linked to gSe2",
                        },
                    ],
                },
            }
        )

        self.assertEqual(
            ["Comment[RelatedExperiment]", "GSE3", "GSE2", "GSE4"],
            self.row(rows, "Comment[RelatedExperiment]"),
        )

    def test_idf_investigations_omits_related_experiment_for_unrelated_relations(self):
        rows = IDFConstructor()._idf_investigations(
            data={
                "series": {
                    "title": "Example",
                    "accession": [{"value": "GSE1", "database": "GEO"}],
                    "relation": [
                        {"type": "SRA", "target": "SRX1"},
                        {"type": "BioSample", "target": "SAMN1"},
                    ],
                },
            }
        )

        self.assertFalse(any(row[0] == "Comment[RelatedExperiment]" for row in rows))

    def test_idf_investigations_omits_related_experiment_when_relation_missing(self):
        rows = IDFConstructor()._idf_investigations(
            data={
                "series": {
                    "title": "Example",
                    "accession": [{"value": "GSE1", "database": "GEO"}],
                },
            }
        )

        self.assertEqual(
            [
                ["Investigation Title", "Example"],
                ["Comment[SecondaryAccession]", "GSE1"],
                ["Comment[SecondaryAccessionTermSourceRef]", "GEO"],
                ["Comment[ArrayExpressAccession]", "E-GEOD-1"],
            ],
            rows,
        )

    def test_idf_investigations_merges_secondary_accessions_with_aligned_sources(self):
        rows = IDFConstructor()._idf_investigations(
            data={
                "series": {
                    "title": "Example",
                    "accession": [{"value": "GSE1", "database": "GEO"}],
                },
                "sample": [
                    {"ena_accession": ["ERP137216", "SRP999"]},
                    {"ena_accession": ["DRP123"]},
                ],
            }
        )

        self.assertEqual(
            ["Comment[SecondaryAccession]", "GSE1", "ERP137216", "SRP999", "DRP123"],
            self.row(rows, "Comment[SecondaryAccession]"),
        )
        self.assertEqual(
            ["Comment[SecondaryAccessionTermSourceRef]", "GEO", "ENA", "SRA", "DRA"],
            self.row(rows, "Comment[SecondaryAccessionTermSourceRef]"),
        )
        self.assertEqual(
            1,
            sum(row[0] == "Comment[SecondaryAccession]" for row in rows),
        )

    def test_idf_investigations_dedupes_secondary_accessions_and_preserves_source_alignment(self):
        rows = IDFConstructor()._idf_investigations(
            data={
                "series": {
                    "title": "Example",
                    "accession": [
                        {"value": "GSE1", "database": "GEO"},
                        {"value": "custom1", "database": "CustomDB"},
                    ],
                },
                "sample": [
                    {"ena_accession": ["erp1", "ERP1", "XYZ1"]},
                    {"ena_accession": ["SRP2", "erp1"]},
                ],
            }
        )

        self.assertEqual(
            ["Comment[SecondaryAccession]", "GSE1", "custom1", "erp1", "XYZ1", "SRP2"],
            self.row(rows, "Comment[SecondaryAccession]"),
        )
        self.assertEqual(
            ["Comment[SecondaryAccessionTermSourceRef]", "GEO", "CustomDB", "ENA", None, "SRA"],
            self.row(rows, "Comment[SecondaryAccessionTermSourceRef]"),
        )

    def test_miniml2idf_emits_one_merged_secondary_accession_row(self):
        rows = IDFConstructor().miniml2idf(
            data={
                "series": {
                    "title": "Example",
                    "accession": [{"value": "GSE1", "database": "GEO"}],
                },
                "sample": [
                    {"ena_accession": ["ERP137216", "SRP999"]},
                    {"ena_accession": ["DRP123", "ERP137216"]},
                ],
                "contributor": [],
            },
            technology_type="bulk_sequencing",
        )

        secondary_rows = [row for row in rows if row[0] == "Comment[SecondaryAccession]"]
        self.assertEqual(
            [["Comment[SecondaryAccession]", "GSE1", "ERP137216", "SRP999", "DRP123"]],
            secondary_rows,
        )
        self.assertEqual(
            ["Comment[SecondaryAccessionTermSourceRef]", "GEO", "ENA", "SRA", "DRA"],
            self.row(rows, "Comment[SecondaryAccessionTermSourceRef]"),
        )

    def test_generated_idf_comment_labels_use_normalized_bracket_style(self):
        rows = IDFConstructor().miniml2idf(
            data={
                "series": {
                    "title": "Example",
                    "accession": [{"value": "GSE1", "database": "GEO"}],
                    "relation": [{"type": "SuperSeries of", "target": "GSE2"}],
                    "status": [{"release_date": "2026-05-18"}],
                },
                "sample": [{"sra_run": [{"run": "ERR5385036"}], "ena_accession": ["ERP137216"]}],
                "contributor": [],
            },
            technology_type="bulk_sequencing",
        )

        for row in rows:
            label = row[0]
            if isinstance(label, str) and label.startswith("Comment"):
                self.assertIsNone(re.match(r"^Comment\s+\[", label))
                self.assertRegex(label, r"^Comment\[[^]]+\]$")

    def test_miniml2idf_moves_all_comment_rows_to_bottom(self):
        constructor = IDFConstructor()
        with patch.object(constructor, "_current_idf_date", return_value="2026-05-26"):
            rows = constructor.miniml2idf(
                data={
                    "series": {
                        "title": "Example",
                        "accession": [{"value": "GSE1", "database": "GEO"}],
                        "relation": [{"type": "SuperSeries of", "target": "GSE2"}],
                        "status": [
                            {
                                "submission_date": "2026-05-18",
                                "release_date": "2026-06-05",
                                "last_update_date": "2026-06-06",
                            }
                        ],
                    },
                    "sample": [],
                    "contributor": [],
                }
            )

        comment_indexes = [
            index
            for index, row in enumerate(rows)
            if constructor._is_comment_row(row)
        ]
        self.assertEqual(
            list(range(min(comment_indexes), len(rows))),
            comment_indexes,
        )
        self.assertLess(rows.index(["Investigation Title", "Example"]), comment_indexes[0])
        self.assertLess(rows.index(["Date of Experiment", "2026-05-18"]), comment_indexes[0])
        self.assertLess(rows.index(["Public Release Date", "2026-06-05"]), comment_indexes[0])
        self.assertLess(rows.index(["SDRF File"]), comment_indexes[0])
        term_source_index = next(index for index, row in enumerate(rows) if row[0] == "Term Source Name")
        self.assertLess(term_source_index, comment_indexes[0])
        self.assertEqual(
            [
                ["Comment[SecondaryAccession]", "GSE1"],
                ["Comment[SecondaryAccessionTermSourceRef]", "GEO"],
                ["Comment[ArrayExpressAccession]", "E-GEOD-1"],
                ["Comment[RelatedExperiment]", "GSE2"],
                ["Comment[GEOReleaseDate]", "2026-06-05"],
                ["Comment[GEOLastUpdateDate]", "2026-06-06"],
                ["Comment[ArrayExpressSubmissionDate]", "2026-05-26"],
            ],
            rows[comment_indexes[0]:],
        )

    def test_miniml2idf_moves_platform_comment_rows_to_bottom(self):
        rows = IDFConstructor().miniml2idf(
            data={
                "series": {"title": "Example", "accession": [{"value": "GSE1", "database": "GEO"}]},
                "sample": [{"sra_run": [{"run": "ERR5385036"}]}],
                "contributor": [],
            },
            technology_type="droplet_single_cell_sequencing",
        )

        comment_rows = [row for row in rows if IDFConstructor()._is_comment_row(row)]
        self.assertEqual(comment_rows, rows[-len(comment_rows):])
        self.assertIn(["Comment[AEExperiment]"], comment_rows)
        self.assertIn(["Comment[SequenceDataURI]", "http://www.ebi.ac.uk/ena/data/view/ERR5385036"], comment_rows)
        self.assertIn(["Comment[AEExpectedClusters]"], comment_rows)
        term_source_index = next(index for index, row in enumerate(rows) if row[0] == "Term Source Name")
        self.assertLess(term_source_index, rows.index(["Comment[AEExperiment]"]))

    def test_miniml2idf_moves_experiment_description_after_investigation_title(self):
        rows = IDFConstructor().miniml2idf(
            data={
                "series": {
                    "title": "Example",
                    "summary": "Summary text",
                    "overall_design": "Design text",
                    "accession": [{"value": "GSE1", "database": "GEO"}],
                },
                "sample": [],
                "contributor": [],
            }
        )

        title_index = rows.index(["Investigation Title", "Example"])
        self.assertEqual(
            ["Experiment Description", "Summary text. Design text"],
            rows[title_index + 1],
        )
        comment_rows = [row for row in rows if IDFConstructor()._is_comment_row(row)]
        self.assertEqual(comment_rows, rows[-len(comment_rows):])

    def test_move_experiment_description_after_title_preserves_rows_when_title_missing(self):
        rows = [
            ["MAGE-TAB Version", "1.1"],
            ["Experiment Description", "Summary text"],
            ["Protocol Name", "P-1"],
        ]

        self.assertEqual(
            rows,
            IDFConstructor()._move_experiment_description_after_title(rows=list(rows)),
        )

    def test_move_experiment_description_after_title_preserves_rows_when_description_missing(self):
        rows = [
            ["MAGE-TAB Version", "1.1"],
            ["Investigation Title", "Example"],
            ["Protocol Name", "P-1"],
        ]

        self.assertEqual(
            rows,
            IDFConstructor()._move_experiment_description_after_title(rows=list(rows)),
        )

    def test_is_comment_row_detects_only_top_level_comment_rows(self):
        constructor = IDFConstructor()

        self.assertTrue(constructor._is_comment_row(["Comment[ArrayExpressAccession]", "E-GEOD-1"]))
        self.assertTrue(constructor._is_comment_row(["  comment[lowercase]  ", "value"]))
        self.assertFalse(constructor._is_comment_row([]))
        self.assertFalse(constructor._is_comment_row(["Investigation Title", "Example"]))
        self.assertFalse(constructor._is_comment_row([None, "value"]))
        self.assertFalse(constructor._is_comment_row("Comment[ArrayExpressAccession]"))

    def test_idf_platform_specific_dispatches_to_handler_build(self):
        class RecordingPlatformIDFHandler(_BasePlatformIDFHandler):
            def build(self) -> list:
                self.parent.seen_platform_handler = {
                    "data": self.data,
                    "technology_type": self.technology_type,
                }
                return [["Comment[PlatformSpecific]", self.technology_type]]

        constructor = IDFConstructor()
        data = {"series": [], "sample": [], "contributor": []}
        constructor._platform_idf_handler_class = Mock(return_value=RecordingPlatformIDFHandler)

        with patch.object(constructor, "_idf_publications", return_value=[]), \
             patch.object(constructor, "_idf_protocols", return_value=[]), \
             patch.object(constructor, "_idf_term_source", return_value=[]):
            rows = constructor.miniml2idf(
                data=data,
                technology_type="bulk_sequencing",
            )

        constructor._platform_idf_handler_class.assert_called_once_with(technology_type="bulk_sequencing")
        self.assertEqual(
            {"data": data, "technology_type": "bulk_sequencing"},
            constructor.seen_platform_handler,
        )
        self.assertIn(["Comment[PlatformSpecific]", "bulk_sequencing"], rows)

    def test_platform_idf_handler_inheritance_tree(self):
        self.assertTrue(issubclass(_SequencingPlatformIDFHandler, _BasePlatformIDFHandler))
        self.assertTrue(issubclass(_BulkSequencingPlatformIDFHandler, _SequencingPlatformIDFHandler))
        self.assertTrue(issubclass(_PlateSingleCellSequencingPlatformIDFHandler, _BulkSequencingPlatformIDFHandler))
        self.assertTrue(issubclass(_SingleCellSequencingPlatformIDFHandler, _SequencingPlatformIDFHandler))
        self.assertTrue(issubclass(_DropletSingleCellSequencingPlatformIDFHandler, _SingleCellSequencingPlatformIDFHandler))
        self.assertTrue(issubclass(_SpatialSequencingPlatformIDFHandler, _SingleCellSequencingPlatformIDFHandler))
        self.assertTrue(issubclass(_ArrayPlatformIDFHandler, _BasePlatformIDFHandler))
        self.assertTrue(issubclass(_GenericPlatformIDFHandler, _BasePlatformIDFHandler))

    def test_platform_idf_handler_dispatch_map_matches_technology_keys(self):
        constructor = IDFConstructor()

        expected = {
            "plate_single_cell_sequencing": _PlateSingleCellSequencingPlatformIDFHandler,
            "droplet_single_cell_sequencing": _DropletSingleCellSequencingPlatformIDFHandler,
            "tenx_v2_droplet_single_cell_sequencing": _DropletSingleCellSequencingPlatformIDFHandler,
            "tenx_v3_droplet_single_cell_sequencing": _DropletSingleCellSequencingPlatformIDFHandler,
            "single_cell_sequencing": _SingleCellSequencingPlatformIDFHandler,
            "spatial_sequencing": _SpatialSequencingPlatformIDFHandler,
            "bulk_sequencing": _BulkSequencingPlatformIDFHandler,
            "sequencing": _SequencingPlatformIDFHandler,
            "array": _ArrayPlatformIDFHandler,
            "unknown": _GenericPlatformIDFHandler,
            None: _GenericPlatformIDFHandler,
        }

        for technology_type, handler_class in expected.items():
            with self.subTest(technology_type=technology_type):
                self.assertIs(
                    handler_class,
                    constructor._platform_idf_handler_class(technology_type=technology_type),
                )

    def test_non_sequencing_platform_idf_handlers_do_not_change_rendered_rows(self):
        constructor = IDFConstructor()
        data = {"series": [], "sample": [], "contributor": []}

        rows_without_technology = constructor.miniml2idf(data=data)
        rows_with_array_technology = constructor.miniml2idf(data=data, technology_type="array")

        self.assertEqual(rows_without_technology, rows_with_array_technology)

    def test_sequencing_platform_idf_handler_emits_single_sequence_data_uri(self):
        rows = IDFConstructor()._idf_platform_specific(
            data={"sample": [{"ena_accession": ["ERP137216"], "sra_run": [{"run": "ERR5385036"}]}]},
            technology_type="sequencing",
        )

        self.assertEqual(
            [
                ["Comment[AEExperiment]"],
                ["Comment[AEExperimentType]"],
                ["Comment[AECurator]"],
                ["Comment[SequenceDataURI]", "http://www.ebi.ac.uk/ena/data/view/ERR5385036"],
            ],
            rows,
        )

    def test_sequencing_platform_idf_handler_does_not_emit_secondary_accessions(self):
        rows = IDFConstructor()._idf_platform_specific(
            data={
                "sample": [
                    {"ena_accession": ["ERP137216", "ERP137216"]},
                    {"ena_accession": ["SRP999", "ERP137216"]},
                ]
            },
            technology_type="sequencing",
        )

        self.assertEqual(
            [
                ["Comment[AEExperiment]"],
                ["Comment[AEExperimentType]"],
                ["Comment[AECurator]"],
            ],
            rows,
        )

    def test_sequencing_platform_idf_handler_emits_min_max_sequence_data_uri(self):
        rows = IDFConstructor()._idf_platform_specific(
            data={
                "sample": [
                    {"sra_run": [{"run": "ERR5385041"}, {"run": "ERR5385036"}, {"run": "ERR5385038"}]},
                ]
            },
            technology_type="sequencing",
        )

        self.assertEqual(
            [
                ["Comment[AEExperiment]"],
                ["Comment[AEExperimentType]"],
                ["Comment[AECurator]"],
                ["Comment[SequenceDataURI]", "http://www.ebi.ac.uk/ena/data/view/ERR5385036-ERR5385041"],
            ],
            rows,
        )

    def test_sequencing_platform_idf_handler_deduplicates_runs_across_samples(self):
        rows = IDFConstructor()._idf_platform_specific(
            data={
                "sample": [
                    {"sra_run": [{"run": "ERR5385036"}, {"run": "ERR5385036"}]},
                    {"sra_run": [{"run": "ERR5385041"}]},
                ]
            },
            technology_type="sequencing",
        )

        self.assertEqual(
            [
                ["Comment[AEExperiment]"],
                ["Comment[AEExperimentType]"],
                ["Comment[AECurator]"],
                ["Comment[SequenceDataURI]", "http://www.ebi.ac.uk/ena/data/view/ERR5385036-ERR5385041"],
            ],
            rows,
        )

    def test_sequencing_platform_idf_handler_groups_mixed_prefixes(self):
        rows = IDFConstructor()._idf_platform_specific(
            data={
                "sample": [
                    {
                        "sra_run": [
                            {"run": "ERR5385041"},
                            {"run": "SRR100"},
                            {"run": "ERR5385036"},
                            {"run": "SRR105"},
                        ]
                    },
                ]
            },
            technology_type="sequencing",
        )

        self.assertEqual(
            [
                ["Comment[AEExperiment]"],
                ["Comment[AEExperimentType]"],
                ["Comment[AECurator]"],
                [
                    "Comment[SequenceDataURI]",
                    "http://www.ebi.ac.uk/ena/data/view/ERR5385036-ERR5385041",
                    "http://www.ebi.ac.uk/ena/data/view/SRR100-SRR105",
                ]
            ],
            rows,
        )

    def test_sequencing_platform_idf_handler_ignores_missing_and_malformed_runs(self):
        with self.assertLogs("meta_standards_converter.ae_handlers.ae_idf_handlers", level="WARNING") as logs:
            rows = IDFConstructor()._idf_platform_specific(
                data={
                    "sample": [
                        {
                            "sra_run": [
                                {"run": None},
                                {"run": "ERRABC"},
                                {"run": "not a run"},
                                {"run": "err5385036"},
                            ]
                        },
                    ]
                },
                technology_type="sequencing",
            )

        self.assertEqual(
            [
                ["Comment[AEExperiment]"],
                ["Comment[AEExperimentType]"],
                ["Comment[AECurator]"],
                ["Comment[SequenceDataURI]", "http://www.ebi.ac.uk/ena/data/view/ERR5385036"],
            ],
            rows,
        )
        self.assertTrue(any("SRA run accession missing" in message for message in logs.output))
        self.assertTrue(any("ERRABC" in message for message in logs.output))
        self.assertTrue(any("NOT A RUN" in message for message in logs.output))

    def test_sequencing_platform_idf_handler_returns_empty_comment_rows_without_valid_runs(self):
        with self.assertLogs("meta_standards_converter.ae_handlers.ae_idf_handlers", level="WARNING") as logs:
            rows = IDFConstructor()._idf_platform_specific(
                data={"sample": [{"sra_run": [{"run": None}, {"run": "ERRABC"}]}]},
                technology_type="sequencing",
            )

        self.assertEqual(
            [
                ["Comment[AEExperiment]"],
                ["Comment[AEExperimentType]"],
                ["Comment[AECurator]"],
            ],
            rows,
        )
        self.assertTrue(any("SRA run accession missing" in message for message in logs.output))
        self.assertTrue(any("ERRABC" in message for message in logs.output))
        self.assertTrue(any("SequenceDataURI row skipped" in message for message in logs.output))

    def test_sequencing_platform_idf_handler_emits_empty_comment_rows_without_sra_runs(self):
        with self.assertLogs("meta_standards_converter.ae_handlers.ae_idf_handlers", level="WARNING") as logs:
            rows = IDFConstructor()._idf_platform_specific(
                data={"sample": []},
                technology_type="sequencing",
            )

        self.assertEqual(
            [
                ["Comment[AEExperiment]"],
                ["Comment[AEExperimentType]"],
                ["Comment[AECurator]"],
            ],
            rows,
        )
        self.assertTrue(any("SequenceDataURI row skipped" in message for message in logs.output))

    def test_sequencing_platform_idf_subclasses_emit_sequence_data_uri(self):
        data = {"sample": [{"sra_run": [{"run": "ERR5385036"}, {"run": "ERR5385041"}]}]}
        technology_types = [
            "bulk_sequencing",
            "spatial_sequencing",
            "sequencing",
        ]

        for technology_type in technology_types:
            with self.subTest(technology_type=technology_type):
                rows = IDFConstructor()._idf_platform_specific(data=data, technology_type=technology_type)
                self.assertEqual(
                    [
                        ["Comment[AEExperiment]"],
                        ["Comment[AEExperimentType]"],
                        ["Comment[AECurator]"],
                        ["Comment[SequenceDataURI]", "http://www.ebi.ac.uk/ena/data/view/ERR5385036-ERR5385041"],
                    ],
                    rows,
                )

    def test_single_cell_platform_idf_handlers_emit_ae_experiment_type(self):
        data = {"sample": [{"sra_run": [{"run": "ERR5385036"}, {"run": "ERR5385041"}]}]}
        technology_types = [
            "single_cell_sequencing",
            "droplet_single_cell_sequencing",
            "plate_single_cell_sequencing",
        ]

        for technology_type in technology_types:
            with self.subTest(technology_type=technology_type):
                rows = IDFConstructor()._idf_platform_specific(data=data, technology_type=technology_type)
                self.assertIn(
                    ["Comment[AEExperimentType]", "RNA-seq of coding RNA from single cells"],
                    rows,
                )

    def test_spatial_platform_idf_handler_keeps_blank_ae_experiment_type(self):
        rows = IDFConstructor()._idf_platform_specific(
            data={"sample": [{"sra_run": [{"run": "ERR5385036"}]}]},
            technology_type="spatial_sequencing",
        )

        self.assertIn(["Comment[AEExperimentType]"], rows)
        self.assertNotIn(
            ["Comment[AEExperimentType]", "RNA-seq of coding RNA from single cells"],
            rows,
        )

    def test_droplet_single_cell_platform_idf_handler_emits_droplet_empty_comment_rows(self):
        rows = IDFConstructor()._idf_platform_specific(
            data={"sample": [{"sra_run": [{"run": "ERR5385036"}, {"run": "ERR5385041"}]}]},
            technology_type="droplet_single_cell_sequencing",
        )

        self.assertEqual(
            [
                ["Comment[AEExperiment]"],
                ["Comment[AEExperimentType]", "RNA-seq of coding RNA from single cells"],
                ["Comment[AECurator]"],
                ["Comment[SequenceDataURI]", "http://www.ebi.ac.uk/ena/data/view/ERR5385036-ERR5385041"],
                ["Comment[AEExpectedClusters]"],
                ["Comment[AEAdditionalAttributes]"],
                ["Comment[AEBatchEffect]"],
            ],
            rows,
        )

    def test_droplet_single_cell_platform_idf_handler_emits_droplet_rows_without_runs(self):
        rows = IDFConstructor()._idf_platform_specific(
            data={"sample": []},
            technology_type="droplet_single_cell_sequencing",
        )

        self.assertEqual(
            [
                ["Comment[AEExperiment]"],
                ["Comment[AEExperimentType]", "RNA-seq of coding RNA from single cells"],
                ["Comment[AECurator]"],
                ["Comment[AEExpectedClusters]"],
                ["Comment[AEAdditionalAttributes]"],
                ["Comment[AEBatchEffect]"],
            ],
            rows,
        )

    def test_droplet_single_cell_platform_idf_rows_do_not_leak_to_other_handlers(self):
        data = {"sample": [{"sra_run": [{"run": "ERR5385036"}]}]}
        droplet_labels = {
            "Comment[AEExpectedClusters]",
            "Comment[AEAdditionalAttributes]",
            "Comment[AEBatchEffect]",
        }
        technology_types = [
            "single_cell_sequencing",
            "spatial_sequencing",
            "bulk_sequencing",
            "plate_single_cell_sequencing",
            "sequencing",
            "array",
            "unknown",
        ]

        for technology_type in technology_types:
            with self.subTest(technology_type=technology_type):
                rows = IDFConstructor()._idf_platform_specific(data=data, technology_type=technology_type)
                self.assertTrue(droplet_labels.isdisjoint(row[0] for row in rows))

    def test_non_sequencing_platform_idf_handlers_do_not_emit_sequence_data_uri(self):
        data = {"sample": [{"sra_run": [{"run": "ERR5385036"}, {"run": "ERR5385041"}]}]}
        technology_types = ["array", "generic", "unknown", None]

        for technology_type in technology_types:
            with self.subTest(technology_type=technology_type):
                rows = IDFConstructor()._idf_platform_specific(data=data, technology_type=technology_type)
                self.assertEqual([], rows)

    def test_idf_publications_prefers_enriched_pubmed_records(self):
        fetcher = Mock()
        rows = IDFConstructor(pubmed_fetcher=fetcher)._idf_publications(
            {
                "series": {
                    "pubmed_id": ["123"],
                    "pubmed_publication": [
                        {
                            "pubmed_id": "123",
                            "doi": "10.1000/example",
                            "author_list": "Ada Lovelace",
                            "title": "Example paper",
                            "status": "published",
                            "status_term_source_ref": "EFO",
                            "status_term_accession_number": "EFO_0001796",
                        }
                    ],
                }
            }
        )

        self.assertEqual(["PubMed ID", "123"], self.row(rows, "PubMed ID"))
        self.assertEqual(["Publication DOI", "10.1000/example"], self.row(rows, "Publication DOI"))
        self.assertEqual(["Publication Author List", "Ada Lovelace"], self.row(rows, "Publication Author List"))
        self.assertEqual(["Publication Title", "Example paper"], self.row(rows, "Publication Title"))
        self.assertEqual(["Publication Status", "published"], self.row(rows, "Publication Status"))
        self.assertEqual(["Status Term Source Ref", "EFO"], self.row(rows, "Status Term Source Ref"))
        self.assertEqual(
            ["Status Term Accession Number", "EFO_0001796"],
            self.row(rows, "Status Term Accession Number"),
        )
        fetcher.pubmed_summary.assert_not_called()

    def test_idf_publications_warns_when_pubmed_lookup_has_no_pubmed_id(self):
        with self.assertLogs("meta_standards_converter.ae_handlers.ae_idf_handlers", level="WARNING") as logs:
            rows = IDFConstructor()._idf_publications({"series": {}})

        self.assertEqual(["PubMed ID"], self.row(rows, "PubMed ID"))
        self.assertTrue(any("no PubMed ID found" in message for message in logs.output))

    def test_person_address_starts_with_organization_only(self):
        rows = IDFConstructor()._idf_persons(
            {
                "contributor": [
                    {
                        "person": {"first": "Ada", "last": "Lovelace"},
                        "affiliation": "Analytical Engine Lab",
                        "organization": "Royal Society",
                        "address": {
                            "line": ["12 Computation Way"],
                            "city": "London",
                            "country": "United Kingdom",
                        },
                    },
                ],
            }
        )

        self.assertEqual(
            [
                "Person Address",
                "Royal Society, 12 Computation Way, London, United Kingdom",
            ],
            self.row(rows, "Person Address"),
        )
        self.assertEqual(
            ["Person Affiliation", "Royal Society"],
            self.row(rows, "Person Affiliation"),
        )

    def test_person_address_preserves_address_only_contributors(self):
        rows = IDFConstructor()._idf_persons(
            {
                "contributor": [
                    {
                        "person": {"first": "Grace", "last": "Hopper"},
                        "address": {
                            "line": ["1 Compiler Plaza"],
                            "state": "New York",
                            "country": "USA",
                        },
                    },
                ],
            }
        )

        self.assertEqual(
            ["Person Address", "1 Compiler Plaza, New York, USA"],
            self.row(rows, "Person Address"),
        )
        self.assertEqual(["Person Affiliation", None], self.row(rows, "Person Affiliation"))

    def test_idf_experimental_uses_tags_with_multiple_labels_as_factors(self):
        rows = IDFConstructor()._idf_experimental(
            {
                "sample": [
                    {
                        "channel": [
                            {
                                "characteristics": [
                                    {"tag": "sex", "value": "male"},
                                    {"tag": "tissue", "value": "brain"},
                                ],
                            },
                        ],
                    },
                    {
                        "channel": [
                            {
                                "characteristics": [
                                    {"tag": "sex", "value": "female"},
                                    {"tag": "tissue", "value": "brain"},
                                ],
                            },
                        ],
                    },
                ],
            }
        )

        self.assertEqual(["Experimental Factor Name", "sex"], self.row(rows, "Experimental Factor Name"))
        self.assertEqual(["Experimental Factor Type", "sex"], self.row(rows, "Experimental Factor Type"))
        self.assertEqual(
            ["Experimental Factor Term Source REF", None],
            self.row(rows, "Experimental Factor Term Source REF"),
        )
        self.assertEqual(
            ["Experimental Factor Term Accession Number", None],
            self.row(rows, "Experimental Factor Term Accession Number"),
        )

    def test_idf_experimental_ignores_untagged_and_blank_characteristics(self):
        rows = IDFConstructor()._idf_experimental(
            {
                "sample": [
                    {
                        "channel": [
                            {
                                "characteristics": [
                                    "untagged characteristic",
                                    {"tag": "", "value": "male"},
                                    {"tag": "sex", "value": ""},
                                    {"tag": "condition", "value": "treated"},
                                ],
                            },
                        ],
                    },
                    {
                        "channel": [
                            {
                                "characteristics": [
                                    {"tag": "condition", "value": "treated"},
                                ],
                            },
                        ],
                    },
                ],
            }
        )

        self.assertEqual(["Experimental Factor Name"], self.row(rows, "Experimental Factor Name"))
        self.assertEqual(["Experimental Factor Type"], self.row(rows, "Experimental Factor Type"))

    def test_idf_experimental_preserves_factor_order_and_normalizes_whitespace(self):
        rows = IDFConstructor()._idf_experimental(
            {
                "sample": [
                    {
                        "channel": [
                            {
                                "characteristics": [
                                    {"tag": "dose level", "value": "low"},
                                    {"tag": " time\tpoint ", "value": " day\n1 "},
                                ],
                            },
                            {
                                "characteristics": [
                                    {"tag": "Dose   Level", "value": "high"},
                                    {"tag": "time point", "value": "day  2"},
                                ],
                            },
                        ],
                    },
                ],
            }
        )

        self.assertEqual(
            ["Experimental Factor Name", "dose level", "time point"],
            self.row(rows, "Experimental Factor Name"),
        )
        self.assertEqual(
            ["Experimental Factor Type", "dose level", "time point"],
            self.row(rows, "Experimental Factor Type"),
        )
        self.assertEqual(
            ["Experimental Factor Term Accession Number", None, None],
            self.row(rows, "Experimental Factor Term Accession Number"),
        )

    def test_idf_qc_rep_norm_returns_label_only_rows(self):
        rows = IDFConstructor()._idf_qc_rep_norm(data={})

        self.assertEqual(
            [
                ["Quality Control Type"],
                ["Quality Control Term Source REF"],
                ["Quality Control Term Accession Number"],
                ["Replicate Type"],
                ["Replicate Term Source REF"],
                ["Replicate Term Accession Number"],
                ["Normalization Type"],
                ["Normalization Term Source REF"],
                ["Normalization Term Accession Number"],
            ],
            rows,
        )

    def test_miniml2idf_omits_qc_rep_norm_rows(self):
        constructor = IDFConstructor()

        with patch.object(constructor, "_idf_publications", return_value=[]), \
             patch.object(constructor, "_idf_protocols", return_value=[]), \
             patch.object(constructor, "_idf_term_source", return_value=[]):
            rows = constructor.miniml2idf(
                data={
                    "series": {"title": "Example"},
                    "contributor": [],
                    "sample": [],
                }
            )

        self.assertNotIn(["Quality Control Type"], rows)
        self.assertNotIn(["Quality Control Term Source REF"], rows)
        self.assertNotIn(["Replicate Type"], rows)
        self.assertNotIn(["Normalization Type"], rows)

    def test_idf_dates_normalize_parseable_dates(self):
        constructor = IDFConstructor()
        with patch.object(constructor, "_current_idf_date", return_value="2026-05-26"):
            rows = constructor._idf_dates(
                {
                    "series": {
                        "status": [
                            {
                                "submission_date": "2026-05-18",
                                "release_date": "05/06/2026",
                                "last_update_date": "5 June 2026",
                            },
                            {
                                "submission_date": "not-a-date",
                                "release_date": "18-05-2026",
                                "last_update_date": "",
                            },
                        ],
                    },
                }
            )

        self.assertEqual(
            ["Date of Experiment", "2026-05-18", "not-a-date"],
            self.row(rows, "Date of Experiment"),
        )
        self.assertEqual(
            ["Public Release Date", "2026-05-18"],
            self.row(rows, "Public Release Date"),
        )
        self.assertEqual(
            ["Comment[GEOReleaseDate]", "2026-06-05", "2026-05-18"],
            self.row(rows, "Comment[GEOReleaseDate]"),
        )
        self.assertEqual(
            ["Comment[GEOLastUpdateDate]", "2026-06-05", ""],
            self.row(rows, "Comment[GEOLastUpdateDate]"),
        )
        self.assertEqual(
            ["Comment[ArrayExpressSubmissionDate]", "2026-05-26"],
            self.row(rows, "Comment[ArrayExpressSubmissionDate]"),
        )

    def test_idf_dates_arrayexpress_submission_date_is_one_current_date_for_multiple_series(self):
        constructor = IDFConstructor()
        with patch.object(constructor, "_current_idf_date", return_value="2026-05-26"):
            rows = constructor._idf_dates(
                {
                    "series": [
                        {"status": [{"submission_date": "2026-05-18"}]},
                        {"status": [{"submission_date": "2026-05-19"}]},
                    ],
                }
            )

        self.assertEqual(
            ["Comment[ArrayExpressSubmissionDate]", "2026-05-26"],
            self.row(rows, "Comment[ArrayExpressSubmissionDate]"),
        )

    def test_idf_dates_public_release_date_is_empty_without_release_dates(self):
        rows = IDFConstructor()._idf_dates({"series": {"status": [{"submission_date": "2026-05-18"}]}})

        self.assertEqual(["Public Release Date", None], self.row(rows, "Public Release Date"))
        self.assertEqual(["Comment[GEOReleaseDate]", None], self.row(rows, "Comment[GEOReleaseDate]"))

    def test_idf_dates_public_release_date_ignores_unparseable_release_dates(self):
        rows = IDFConstructor()._idf_dates(
            {
                "series": {
                    "status": [
                        {"release_date": "not-a-date"},
                        {"release_date": "2026-07-01"},
                    ],
                },
            }
        )

        self.assertEqual(["Public Release Date", "2026-07-01"], self.row(rows, "Public Release Date"))
        self.assertEqual(
            ["Comment[GEOReleaseDate]", "not-a-date", "2026-07-01"],
            self.row(rows, "Comment[GEOReleaseDate]"),
        )

    def test_idf_dates_public_release_date_is_empty_when_all_release_dates_are_unparseable(self):
        rows = IDFConstructor()._idf_dates(
            {"series": {"status": [{"release_date": "not-a-date"}, {"release_date": ""}]}}
        )

        self.assertEqual(["Public Release Date", None], self.row(rows, "Public Release Date"))
        self.assertEqual(
            ["Comment[GEOReleaseDate]", "not-a-date", ""],
            self.row(rows, "Comment[GEOReleaseDate]"),
        )

    def test_miniml2idf_places_empty_sdrf_row_before_term_sources(self):
        constructor = IDFConstructor()

        with patch.object(constructor, "_idf_investigations", return_value=[]), \
             patch.object(constructor, "_idf_experimental", return_value=[]), \
             patch.object(constructor, "_idf_platform_specific", return_value=[]), \
             patch.object(constructor, "_idf_persons", return_value=[]), \
             patch.object(constructor, "_idf_qc_rep_norm", return_value=[]), \
             patch.object(constructor, "_idf_dates", return_value=[]), \
             patch.object(constructor, "_idf_publications", return_value=[]), \
             patch.object(constructor, "_idf_experiments", return_value=[]), \
             patch.object(constructor, "_idf_protocols", return_value=[["Protocol Name", "P-1"]]), \
             patch.object(constructor, "_idf_term_source", return_value=[["Term Source Name"]]):
            rows = constructor.miniml2idf(data={})

        self.assertEqual(["Protocol Name", "P-1"], rows[-3])
        self.assertEqual(["SDRF File"], rows[-2])
        self.assertEqual(["Term Source Name"], rows[-1])


class TestAEConstructor(unittest.TestCase):
    def row(self, rows, label):
        return next(row for row in rows if row[0] == label)

    def protocol_description_by_name(self, rows):
        names = self.row(rows, "Protocol Name")[1:]
        descriptions = self.row(rows, "Protocol Description")[1:]
        return dict(zip(names, descriptions))

    def non_empty_sdrf_protocol_refs(self, sdrf):
        return [
            row[index]
            for index, label in enumerate(sdrf[0])
            if label == "Protocol REF"
            for row in sdrf[1:]
            if row[index]
        ]

    def detection_data(
        self,
        sample=None,
        platform_technology="high-throughput sequencing",
        series=None,
    ):
        return {
            "series": series or {
                "accession": [{"value": "GSE1"}],
                "sample_ref": [{"ref": "GSM1"}],
                "title": "Example series",
            },
            "platform": [
                {
                    "iid": "GPL1",
                    "technology": platform_technology,
                    "title": "Example platform",
                    "accession": [{"value": "GPL1"}],
                }
            ],
            "contributor": [],
            "sample": [
                sample or {
                    "iid": "GSM1",
                    "accession": [{"value": "GSM1"}],
                    "platform_ref": {"ref": "GPL1"},
                    "channel": [{"source": "source 1"}],
                }
            ],
        }

    def test_miniml2magetab_builds_sdrf_then_idf_and_inserts_sdrf(self):
        data = {"series": []}
        sdrf = [["Source Name"], ["sample 1"]]
        idf = [["Investigation Title", "Example"], ["SDRF File"], ["Term Source Name"]]
        idf_constructor = Mock()
        idf_constructor.miniml2idf.return_value = idf
        sdrf_constructor = Mock()
        sdrf_constructor._miniml2sdrf.return_value = sdrf

        result = AEConstructor(
            idf_constructor=idf_constructor,
            sdrf_constructor=sdrf_constructor,
        ).miniml2magetab(data=data)

        sdrf_constructor._miniml2sdrf.assert_called_once_with(
            data=data,
            protocol_registry=ANY,
            technology_type="generic",
        )
        protocol_registry = sdrf_constructor._miniml2sdrf.call_args.kwargs["protocol_registry"]
        idf_constructor.miniml2idf.assert_called_once_with(
            data=data,
            protocol_registry=protocol_registry,
            technology_type="generic",
        )
        sdrf_constructor._add_sdrf_to_idf.assert_not_called()
        self.assertEqual(
            [["Investigation Title", "Example"], ["SDRF File", sdrf], ["Term Source Name"]],
            result,
        )

    def test_miniml2magetab_replaces_existing_sdrf_cell_and_preserves_extra_cells(self):
        data = {"series": []}
        sdrf = [["Source Name"], ["sample 1"]]
        idf_constructor = Mock()
        idf_constructor.miniml2idf.return_value = [
            ["Investigation Title", "Example"],
            ["SDRF File", "old.sdrf.txt", "curator note"],
        ]
        sdrf_constructor = Mock()
        sdrf_constructor._miniml2sdrf.return_value = sdrf

        result = AEConstructor(
            idf_constructor=idf_constructor,
            sdrf_constructor=sdrf_constructor,
        ).miniml2magetab(data=data)

        sdrf_constructor._miniml2sdrf.assert_called_once_with(
            data=data,
            protocol_registry=ANY,
            technology_type="generic",
        )
        protocol_registry = sdrf_constructor._miniml2sdrf.call_args.kwargs["protocol_registry"]
        idf_constructor.miniml2idf.assert_called_once_with(
            data=data,
            protocol_registry=protocol_registry,
            technology_type="generic",
        )
        self.assertEqual(["SDRF File", sdrf, "curator note"], self.row(result, "SDRF File"))

    def test_miniml2magetab_passes_detected_technology_type_to_constructors(self):
        data = self.detection_data(
            sample={
                "iid": "GSM1",
                "accession": [{"value": "GSM1"}],
                "platform_ref": {"ref": "GPL1"},
                "library_strategy": "RNA-Seq",
                "channel": [{"source": "source 1"}],
            }
        )
        idf_constructor = Mock()
        idf_constructor.miniml2idf.return_value = [["SDRF File"]]
        sdrf_constructor = Mock()
        sdrf_constructor._miniml2sdrf.return_value = [["Source Name"], ["GSM1"]]

        AEConstructor(
            idf_constructor=idf_constructor,
            sdrf_constructor=sdrf_constructor,
        ).miniml2magetab(data=data)

        sdrf_constructor._miniml2sdrf.assert_called_once_with(
            data=data,
            protocol_registry=ANY,
            technology_type="bulk_sequencing",
        )
        protocol_registry = sdrf_constructor._miniml2sdrf.call_args.kwargs["protocol_registry"]
        idf_constructor.miniml2idf.assert_called_once_with(
            data=data,
            protocol_registry=protocol_registry,
            technology_type="bulk_sequencing",
        )

    def test_detect_ae_technology_returns_bulk_sequencing_for_regular_sra(self):
        sample = {
            "iid": "GSM1",
            "accession": [{"value": "GSM1"}],
            "platform_ref": {"ref": "GPL1"},
            "library_strategy": "RNA-Seq",
            "channel": [{"source": "source 1"}],
        }

        self.assertEqual(
            "bulk_sequencing",
            AEConstructor()._detect_ae_technology(self.detection_data(sample=sample)),
        )

    def test_detect_ae_technology_returns_array_for_array_platform(self):
        self.assertEqual(
            "array",
            AEConstructor()._detect_ae_technology(
                self.detection_data(platform_technology="expression array")
            ),
        )

    def test_detect_ae_technology_returns_droplet_single_cell(self):
        sample = {
            "iid": "GSM1",
            "accession": [{"value": "GSM1"}],
            "platform_ref": {"ref": "GPL1"},
            "library_source": "single cell",
            "description": "10x Chromium RNA-seq",
            "channel": [{"source": "source 1"}],
        }

        self.assertEqual(
            "droplet_single_cell_sequencing",
            AEConstructor()._detect_ae_technology(self.detection_data(sample=sample)),
        )

    def test_detect_ae_technology_returns_tenx_v2_droplet_single_cell(self):
        sample = {
            "iid": "GSM1",
            "accession": [{"value": "GSM1"}],
            "platform_ref": {"ref": "GPL1"},
            "library_source": "single cell",
            "description": "10x Chromium single cell v2",
            "channel": [{"source": "source 1"}],
        }

        self.assertEqual(
            "tenx_v2_droplet_single_cell_sequencing",
            AEConstructor()._detect_ae_technology(self.detection_data(sample=sample)),
        )

    def test_detect_ae_technology_returns_tenx_v3_droplet_single_cell(self):
        sample = {
            "iid": "GSM1",
            "accession": [{"value": "GSM1"}],
            "platform_ref": {"ref": "GPL1"},
            "library_source": "single cell",
            "description": "10x Chromium single cell v3",
            "channel": [{"source": "source 1"}],
        }

        self.assertEqual(
            "tenx_v3_droplet_single_cell_sequencing",
            AEConstructor()._detect_ae_technology(self.detection_data(sample=sample)),
        )

    def test_detect_ae_technology_returns_spatial_sequencing(self):
        sample = {
            "iid": "GSM1",
            "accession": [{"value": "GSM1"}],
            "platform_ref": {"ref": "GPL1"},
            "library_source": "single cell",
            "description": "10x Visium spatial transcriptomics",
            "channel": [{"source": "source 1"}],
        }

        self.assertEqual(
            "spatial_sequencing",
            AEConstructor()._detect_ae_technology(self.detection_data(sample=sample)),
        )

    def test_detect_ae_technology_keeps_spatial_precedence_over_tenx_version(self):
        sample = {
            "iid": "GSM1",
            "accession": [{"value": "GSM1"}],
            "platform_ref": {"ref": "GPL1"},
            "library_source": "single cell",
            "description": "10x Visium spatial transcriptomics v3",
            "channel": [{"source": "source 1"}],
        }

        self.assertEqual(
            "spatial_sequencing",
            AEConstructor()._detect_ae_technology(self.detection_data(sample=sample)),
        )

    def test_detect_ae_technology_returns_plate_single_cell(self):
        sample = {
            "iid": "GSM1",
            "accession": [{"value": "GSM1"}],
            "platform_ref": {"ref": "GPL1"},
            "library_source": "single cell",
            "description": "plate-based single-cell RNA-seq",
            "channel": [{"source": "source 1"}],
        }

        self.assertEqual(
            "plate_single_cell_sequencing",
            AEConstructor()._detect_ae_technology(self.detection_data(sample=sample)),
        )

    def test_detect_ae_technology_returns_generic_without_platform_signal(self):
        self.assertEqual(
            "generic",
            AEConstructor()._detect_ae_technology(
                self.detection_data(platform_technology="other")
            ),
        )

    def test_miniml2magetab_raises_when_idf_has_no_sdrf_row(self):
        idf_constructor = Mock()
        idf_constructor.miniml2idf.return_value = [["Investigation Title", "Example"]]
        sdrf_constructor = Mock()
        sdrf_constructor._miniml2sdrf.return_value = [["Source Name"], ["sample 1"]]

        with self.assertRaisesRegex(ValueError, "IDF does not contain an SDRF File row"):
            AEConstructor(
                idf_constructor=idf_constructor,
                sdrf_constructor=sdrf_constructor,
            ).miniml2magetab(data={"series": []})

    def test_normalize_magetab_rows_still_accepts_legacy_mixed_payloads(self):
        sdrf = [["Source Name"], ["sample 1"]]

        rows = AEConstructor()._normalize_magetab_rows(
            ["MAGE-TAB Version, 1.1", "SDRF file", sdrf]
        )

        self.assertEqual(
            [
                ["MAGE-TAB Version", "1.1"],
                ["SDRF File", sdrf],
            ],
            rows,
        )

    def test_parser_package_can_be_converted_to_magetab(self):
        miniml = """\
<MINiML>
  <Platform iid="GPL1">
    <Accession database="GEO">GPL1</Accession>
    <Title>Example platform</Title>
    <Technology>other</Technology>
  </Platform>
  <Sample iid="GSM1">
    <Accession database="GEO">GSM1</Accession>
    <Title>Example sample</Title>
    <Platform-Ref ref="GPL1" />
    <Channel>
      <Source>whole embryo</Source>
      <Characteristics tag="organism part">embryo</Characteristics>
    </Channel>
  </Sample>
  <Series iid="GSE1">
    <Accession database="GEO">GSE1</Accession>
    <Title>Example series</Title>
    <Summary>Example summary</Summary>
    <Overall-Design>Example design</Overall-Design>
    <Sample-Ref ref="GSM1" />
  </Series>
</MINiML>
"""
        parsed = GEOParser().parse(miniml=miniml)[0]

        magetab = AEConstructor().miniml2magetab(data=parsed)

        self.assertEqual(["Investigation Title", "Example series"], self.row(magetab, "Investigation Title"))
        sdrf = self.row(magetab, "SDRF File")[1]
        self.assertEqual("GSM1", sdrf[1][sdrf[0].index("Source Name")])

    def test_miniml2magetab_strips_quotes_from_idf_and_nested_sdrf_values(self):
        data = {
            "series": {
                "accession": [{"value": "GSE1", "database": "GEO"}],
                "sample_ref": [{"ref": "GSM1"}],
                "title": "\"Quoted\" study",
                "summary": "John's \"summary\"",
                "overall_design": "\"design\"",
                "pubmed_publication": [
                    {
                        "pubmed_id": "PMID'1",
                        "doi": "10.\"quoted\"",
                        "author_list": "O'Neil A",
                        "title": "\"Publication\"",
                        "status": "published",
                        "status_term_source_ref": "EFO",
                        "status_term_accession_number": "EFO_0001796",
                    }
                ],
            },
            "platform": [
                {
                    "iid": "GPL1",
                    "technology": "other",
                    "accession": [{"value": "GPL1"}],
                },
            ],
            "contributor": [
                {
                    "person": {"first": "Ann'e", "last": "\"Smith\""},
                    "email": "a\"b@example.org",
                    "organization": "Org's",
                }
            ],
            "sample": [
                {
                    "iid": "GSM1",
                    "accession": [{"value": "GSM1"}],
                    "platform_ref": {"ref": "GPL1"},
                    "channel": [
                        {
                            "source": "John's \"sample\"",
                            "characteristics": [{"tag": "condition", "value": "\"treated\""}],
                        }
                    ],
                },
            ],
        }

        magetab = AEConstructor().miniml2magetab(data=data)
        sdrf = self.row(magetab, "SDRF File")[1]

        self.assertEqual(["Investigation Title", "Quoted study"], self.row(magetab, "Investigation Title"))
        self.assertEqual(
            ["Experiment Description", "Johns summary. design"],
            self.row(magetab, "Experiment Description"),
        )
        self.assertEqual(["Person Last Name", "Smith"], self.row(magetab, "Person Last Name"))
        self.assertEqual(["Person First Name", "Anne"], self.row(magetab, "Person First Name"))
        self.assertEqual(["Person Email", "ab@example.org"], self.row(magetab, "Person Email"))
        self.assertEqual(["Publication DOI", "10.quoted"], self.row(magetab, "Publication DOI"))
        self.assertEqual("Johns sample", sdrf[1][sdrf[0].index("Comment[Sample_source_name]")])
        self.assertEqual("treated", sdrf[1][sdrf[0].index("Characteristics[condition]")])

    def test_magetab2file_strips_quotes_from_written_idf_and_sdrf_values(self):
        magetab = [
            ["MAGE-TAB Version", "1.1"],
            ["Investigation Title", "\"Quoted\" study"],
            ["Comment[ArrayExpressAccession]", "E-GEOD-1"],
            ["SDRF File", [["Source Name", "Comment[Sample_title]"], ["GSM'1", "\"sample\" title"]]],
        ]

        with tempfile.TemporaryDirectory() as tmpdir:
            idf_path = AEConstructor().magetab2file(magetab=magetab, out=tmpdir)
            sdrf_path = os.path.join(tmpdir, "E-GEOD-1.sdrf.txt")

            with open(idf_path, encoding="utf-8") as handle:
                idf_text = handle.read()
            with open(sdrf_path, encoding="utf-8") as handle:
                sdrf_text = handle.read()

        self.assertIn("Investigation Title\tQuoted study\n", idf_text)
        self.assertIn("GSM1\tsample title\n", sdrf_text)
        self.assertNotIn('"', idf_text + sdrf_text)
        self.assertNotIn("'", idf_text + sdrf_text)

    def test_magetab_protocol_refs_match_idf_protocol_names_for_array(self):
        data = {
            "series": {
                "accession": [{"value": "GSE1"}],
                "sample_ref": [{"ref": "GSM1"}, {"ref": "GSM2"}],
                "title": "example",
            },
            "platform": [
                {
                    "iid": "GPL1",
                    "technology": "expression array",
                    "accession": [{"value": "GPL1"}],
                },
            ],
            "contributor": [],
            "sample": [
                {
                    "iid": "GSM1",
                    "accession": [{"value": "GSM1"}],
                    "platform_ref": {"ref": "GPL1"},
                    "hybridization_protocol": "hyb protocol",
                    "scan_protocol": "scan protocol",
                    "channel": [
                        {
                            "source": "source 1",
                            "extract_protocol": "shared extract",
                            "label": "Cy3",
                            "label_protocol": "label protocol",
                        },
                    ],
                },
                {
                    "iid": "GSM2",
                    "accession": [{"value": "GSM2"}],
                    "platform_ref": {"ref": "GPL1"},
                    "hybridization_protocol": "hyb protocol",
                    "scan_protocol": "scan protocol",
                    "channel": [
                        {
                            "source": "source 2",
                            "extract_protocol": "shared extract",
                            "label": "Cy3",
                            "label_protocol": "label protocol",
                        },
                    ],
                },
            ],
        }

        magetab = AEConstructor().miniml2magetab(data=data)
        sdrf = self.row(magetab, "SDRF File")[1]
        protocol_names = set(self.row(magetab, "Protocol Name")[1:])
        refs = self.non_empty_sdrf_protocol_refs(sdrf=sdrf)
        descriptions = self.protocol_description_by_name(rows=magetab)

        self.assertTrue(refs)
        self.assertTrue(set(refs).issubset(protocol_names))
        self.assertEqual(
            {
                "shared extract",
                "label protocol",
                "hyb protocol",
                "scan protocol",
                "",
            },
            set(descriptions.values()),
        )
        extract_refs = [
            ref
            for ref, description in descriptions.items()
            if description == "shared extract"
        ]
        self.assertEqual(1, len(extract_refs))
        self.assertEqual(2, refs.count(extract_refs[0]))

    def test_magetab_sequencing_library_construction_uses_sra_text_in_idf(self):
        sra_fetcher = Mock()
        sra_fetcher._extract_sra.return_value = ["SRX1"]
        sra_fetcher.fetch_sra_runs.return_value = [
            {
                "geo_sample": "GSM1",
                "experiment": "SRX1",
                "run": "SRR1",
                "scan_name": "SRR1",
                "library_layout": "PAIRED",
                "library_strategy": "RNA-Seq",
                "library_source": "TRANSCRIPTOMIC",
                "library_selection": "cDNA",
                "fastq_files": [],
            }
        ]

        data = {
            "series": {
                "accession": [{"value": "GSE1"}],
                "sample_ref": [{"ref": "GSM1"}],
                "title": "example",
            },
            "platform": [
                {
                    "iid": "GPL1",
                    "technology": "high-throughput sequencing",
                    "accession": [{"value": "GPL1"}],
                },
            ],
            "contributor": [],
            "sample": [
                {
                    "iid": "GSM1",
                    "accession": [{"value": "GSM1"}],
                    "platform_ref": {"ref": "GPL1"},
                    "library_strategy": "RNA-Seq",
                    "library_source": "transcriptomic",
                    "library_selection": "cDNA",
                    "relation": [{"type": "SRA", "target": "SRX1"}],
                    "channel": [{"source": "source 1", "extract_protocol": "extract"}],
                },
            ],
        }

        magetab = AEConstructor(
            sdrf_constructor=SDRFConstructor(insdc_fetcher=sra_fetcher)
        ).miniml2magetab(data=data)
        sdrf = self.row(magetab, "SDRF File")[1]
        protocol_names = set(self.row(magetab, "Protocol Name")[1:])
        refs = self.non_empty_sdrf_protocol_refs(sdrf=sdrf)
        descriptions = self.protocol_description_by_name(rows=magetab)

        self.assertTrue(set(refs).issubset(protocol_names))
        self.assertIn(
            "extract | RNA-Seq | transcriptomic | cDNA | PAIRED | RNA-Seq | TRANSCRIPTOMIC | cDNA",
            set(descriptions.values()),
        )

    def test_magetab_sequencing_includes_required_protocol_defs_and_refs(self):
        data = {
            "series": {
                "accession": [{"value": "GSE1"}],
                "sample_ref": [{"ref": "GSM1"}],
                "title": "example",
            },
            "platform": [
                {
                    "iid": "GPL1",
                    "technology": "high-throughput sequencing",
                    "accession": [{"value": "GPL1"}],
                },
            ],
            "contributor": [],
            "sample": [
                {
                    "iid": "GSM1",
                    "accession": [{"value": "GSM1"}],
                    "platform_ref": {"ref": "GPL1"},
                    "library_strategy": "RNA-Seq",
                    "library_source": "transcriptomic",
                    "library_selection": "cDNA",
                    "channel": [{"source": "source 1", "extract_protocol": "extract"}],
                },
            ],
        }

        magetab = AEConstructor().miniml2magetab(data=data)
        sdrf = self.row(magetab, "SDRF File")[1]
        names_by_type = dict(zip(self.row(magetab, "Protocol Type")[1:], self.row(magetab, "Protocol Name")[1:]))

        self.assertEqual("P-GSE1-3", names_by_type["sample collection protocol"])
        self.assertEqual("P-GSE1-4", names_by_type["nucleic acid sequencing protocol"])
        self.assertIn("P-GSE1-3", self.non_empty_sdrf_protocol_refs(sdrf=sdrf))
        self.assertIn("P-GSE1-4", self.non_empty_sdrf_protocol_refs(sdrf=sdrf))
        self.assertIn("EFO", self.row(magetab, "Term Source Name"))

    def test_magetab_generic_includes_sample_collection_protocol_ref_only(self):
        data = {
            "series": {
                "accession": [{"value": "GSE1"}],
                "sample_ref": [{"ref": "GSM1"}],
                "title": "example",
            },
            "platform": [
                {
                    "iid": "GPL1",
                    "technology": "other",
                    "accession": [{"value": "GPL1"}],
                },
            ],
            "contributor": [],
            "sample": [
                {
                    "iid": "GSM1",
                    "accession": [{"value": "GSM1"}],
                    "platform_ref": {"ref": "GPL1"},
                    "channel": [{"source": "source 1"}],
                },
            ],
        }

        magetab = AEConstructor().miniml2magetab(data=data)
        sdrf = self.row(magetab, "SDRF File")[1]

        self.assertEqual(["Protocol Type", "sample collection protocol"], self.row(magetab, "Protocol Type"))
        self.assertNotIn("nucleic acid sequencing protocol", self.row(magetab, "Protocol Type"))
        self.assertIn("P-GSE1-1", self.non_empty_sdrf_protocol_refs(sdrf=sdrf))

    def test_blank_required_protocol_refs_warn_while_sample_collection_is_defined(self):
        data = {
            "series": {
                "accession": [{"value": "GSE1"}],
                "sample_ref": [{"ref": "GSM1"}],
                "title": "example",
            },
            "platform": [
                {
                    "iid": "GPL1",
                    "technology": "expression array",
                    "accession": [{"value": "GPL1"}],
                },
            ],
            "contributor": [],
            "sample": [
                {
                    "iid": "GSM1",
                    "accession": [{"value": "GSM1"}],
                    "platform_ref": {"ref": "GPL1"},
                    "channel": [{"source": "source 1"}],
                },
            ],
        }
        constructor = AEConstructor()

        magetab = constructor.miniml2magetab(data=data)
        sdrf = self.row(magetab, "SDRF File")[1]

        self.assertIn("Protocol REF", sdrf[0])
        self.assertEqual(["Protocol Name", "P-GSE1-1"], self.row(magetab, "Protocol Name"))
        self.assertEqual(["Protocol Type", "sample collection protocol"], self.row(magetab, "Protocol Type"))
        self.assertIn("P-GSE1-1", self.non_empty_sdrf_protocol_refs(sdrf=sdrf))
        self.assertTrue(
            any("Protocol text missing" in warning for warning in constructor.sdrf_constructor.last_sdrf_audit.warnings)
        )


if __name__ == "__main__":
    unittest.main()
