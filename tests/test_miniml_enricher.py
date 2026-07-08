# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
import os
import sys
import unittest
import xml.etree.ElementTree as ET
from unittest.mock import Mock, call

import requests


ROOT = os.path.dirname(os.path.dirname(__file__))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from meta_standards_converter.enrichers.miniml_enricher import MINiMLEnricher  # noqa: E402


class TestMINiMLEnricher(unittest.TestCase):
    def test_enrich_adds_pubmed_publications(self):
        pubmed_fetcher = Mock()
        pubmed_fetcher.pubmed_summary.side_effect = [
            ("doi-1", "authors 1", "title 1", "published", "EFO", "EFO_0001796"),
            ("doi-2", "authors 2", "title 2", "published", "EFO", "EFO_0001796"),
        ]
        data = {"series": {"pubmed_id": ["123", "123", "456"]}}

        enriched = MINiMLEnricher(pubmed_fetcher=pubmed_fetcher, insdc_fetcher=Mock()).enrich(data=data)

        self.assertEqual(
            [
                {
                    "pubmed_id": "123",
                    "doi": "doi-1",
                    "author_list": "authors 1",
                    "title": "title 1",
                    "status": "published",
                    "status_term_source_ref": "EFO",
                    "status_term_accession_number": "EFO_0001796",
                },
                {
                    "pubmed_id": "456",
                    "doi": "doi-2",
                    "author_list": "authors 2",
                    "title": "title 2",
                    "status": "published",
                    "status_term_source_ref": "EFO",
                    "status_term_accession_number": "EFO_0001796",
                },
            ],
            enriched["series"]["pubmed_publication"],
        )
        self.assertEqual(["123", "123", "456"], enriched["series"]["pubmed_id"])
        self.assertEqual(
            [call(pubmed_id="123"), call(pubmed_id="456")],
            pubmed_fetcher.pubmed_summary.call_args_list,
        )

    def test_enrich_adds_sra_accessions_and_runs(self):
        insdc_fetcher = Mock()
        insdc_fetcher._extract_sra.side_effect = [["SRX1", "SRX1"], ["ERR2"]]
        insdc_fetcher.fetch_sra_runs.side_effect = [
            [{"run": "SRR1", "study": "ERP137216"}],
            [{"run": "ERR2", "study": "ERP137216"}, {"run": "ERR3", "study": "SRP999"}],
        ]
        data = {
            "sample": [
                {
                    "iid": "GSM1",
                    "relation": [
                        {"type": "SRA", "target": "https://example/SRX1"},
                        {"type": "sra", "target": "ERR2"},
                    ],
                }
            ]
        }

        enriched = MINiMLEnricher(pubmed_fetcher=Mock(), insdc_fetcher=insdc_fetcher).enrich(data=data)

        self.assertEqual(["SRX1", "ERR2"], enriched["sample"][0]["sra_accession"])
        self.assertEqual(
            [
                {"run": "SRR1", "study": "ERP137216"},
                {"run": "ERR2", "study": "ERP137216"},
                {"run": "ERR3", "study": "SRP999"},
            ],
            enriched["sample"][0]["sra_run"],
        )
        self.assertEqual(["ERP137216", "SRP999"], enriched["sample"][0]["ena_accession"])
        self.assertEqual(
            [call(accession="SRX1"), call(accession="ERR2")],
            insdc_fetcher.fetch_sra_runs.call_args_list,
        )

    def test_fetch_errors_create_empty_enrichment_without_raising(self):
        pubmed_fetcher = Mock()
        pubmed_fetcher.pubmed_summary.side_effect = requests.RequestException("no pubmed")
        insdc_fetcher = Mock()
        insdc_fetcher._extract_sra.return_value = ["SRX1"]
        insdc_fetcher.fetch_sra_runs.side_effect = ET.ParseError("bad xml")
        data = {
            "series": {"pubmed_id": ["123"]},
            "sample": [{"iid": "GSM1", "relation": [{"type": "SRA", "target": "SRX1"}]}],
        }

        enriched = MINiMLEnricher(pubmed_fetcher=pubmed_fetcher, insdc_fetcher=insdc_fetcher).enrich(data=data)

        self.assertEqual(
            {
                "pubmed_id": "123",
                "doi": None,
                "author_list": None,
                "title": None,
                "status": None,
                "status_term_source_ref": None,
                "status_term_accession_number": None,
            },
            enriched["series"]["pubmed_publication"][0],
        )
        self.assertEqual(["SRX1"], enriched["sample"][0]["sra_accession"])
        self.assertEqual([], enriched["sample"][0]["sra_run"])
        self.assertNotIn("ena_accession", enriched["sample"][0])

    def test_enrich_sra_does_not_add_ena_accession_without_study_values(self):
        insdc_fetcher = Mock()
        insdc_fetcher._extract_sra.return_value = ["SRX1"]
        insdc_fetcher.fetch_sra_runs.return_value = [{"run": "SRR1", "study": None}]
        data = {"sample": [{"iid": "GSM1", "relation": [{"type": "SRA", "target": "SRX1"}]}]}

        enriched = MINiMLEnricher(pubmed_fetcher=Mock(), insdc_fetcher=insdc_fetcher).enrich(data=data)

        self.assertEqual([{"run": "SRR1", "study": None}], enriched["sample"][0]["sra_run"])
        self.assertNotIn("ena_accession", enriched["sample"][0])


if __name__ == "__main__":
    unittest.main()
