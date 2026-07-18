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
import tempfile
import unittest
from unittest.mock import Mock, call


ROOT = os.path.dirname(os.path.dirname(__file__))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from meta_standards_converter.ae_handlers.ae_webfetcher import AEWebFetcher  # noqa: E402


IDF = "MAGE-TAB Version\t1.1\nInvestigation Accession\tE-MTAB-1\nSDRF File\tstudy.sdrf.txt\n"
SDRF = "Source Name\tCharacteristics[organism]\nS1\tHomo sapiens\n"


def response(*, text=None, payload=None):
    item = Mock()
    item.text = text
    item.content = text.encode("utf-8") if text is not None else b""
    item.json.return_value = payload
    item.raise_for_status = Mock()
    return item


class TestAEWebFetcher(unittest.TestCase):
    def test_local_idf_resolves_relative_sdrf(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            idf_path = os.path.join(tmpdir, "study.idf.txt")
            sdrf_path = os.path.join(tmpdir, "study.sdrf.txt")
            with open(idf_path, "w", encoding="utf-8") as handle:
                handle.write(IDF)
            with open(sdrf_path, "w", encoding="utf-8") as handle:
                handle.write(SDRF)

            resolved = AEWebFetcher().resolve(idf_path)

        self.assertEqual(IDF, resolved.idf.text)
        self.assertEqual("study.idf.txt", resolved.idf.name)
        self.assertEqual(["study.sdrf.txt"], [item.name for item in resolved.sdrfs])
        self.assertEqual(SDRF, resolved.sdrfs[0].text)

    def test_http_idf_resolves_relative_sdrf_in_memory(self):
        requester = Mock()
        requester.get.side_effect = [response(text=IDF), response(text=SDRF)]

        resolved = AEWebFetcher(requester=requester).resolve(
            "https://example.org/studies/study.idf.txt"
        )

        self.assertEqual(
            [
                call("https://example.org/studies/study.idf.txt"),
                call("https://example.org/studies/study.sdrf.txt"),
            ],
            requester.get.call_args_list,
        )
        self.assertEqual(SDRF, resolved.sdrfs[0].text)

    def test_explicit_sdrf_override_replaces_idf_references(self):
        requester = Mock()
        requester.get.side_effect = [response(text=IDF), response(text=SDRF)]

        resolved = AEWebFetcher(requester=requester).resolve(
            "https://example.org/study.idf.txt",
            sdrf_sources=["https://override.example/study.sdrf.txt"],
        )

        self.assertEqual(
            call("https://override.example/study.sdrf.txt"),
            requester.get.call_args_list[-1],
        )
        self.assertEqual("study.sdrf.txt", resolved.sdrfs[0].name)

    def test_accession_discovers_and_downloads_idf_and_all_sdrfs(self):
        requester = Mock()
        base = "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/001/E-MTAB-1"
        files = {
            "data": [
                {"Name": "E-MTAB-1.idf.txt", "Type": "IDF File", "path": "E-MTAB-1.idf.txt"},
                {"Name": "part1.sdrf.txt", "Type": "SDRF File", "path": "metadata/part1.sdrf.txt"},
                {"Name": "part2.sdrf.txt", "Type": "SDRF File", "path": "part2.sdrf.txt"},
            ]
        }

        def get(url):
            if url.endswith("/api/v1/files/E-MTAB-1"):
                return response(payload=files)
            if url.endswith("/api/v1/studies/E-MTAB-1/info"):
                return response(payload={"httpLink": base})
            if url.endswith("/Files/E-MTAB-1.idf.txt"):
                return response(text=IDF)
            if url.endswith("/Files/metadata/part1.sdrf.txt"):
                return response(text=SDRF)
            if url.endswith("/Files/part2.sdrf.txt"):
                return response(text=SDRF.replace("S1", "S2"))
            raise AssertionError(url)

        requester.get.side_effect = get
        resolved = AEWebFetcher(requester=requester).resolve("E-MTAB-1")

        self.assertEqual("accession", resolved.source_kind)
        self.assertEqual("E-MTAB-1", resolved.source)
        self.assertEqual(["part1.sdrf.txt", "part2.sdrf.txt"], [x.name for x in resolved.sdrfs])
        self.assertEqual(5, requester.get.call_count)

    def test_accession_rejects_ambiguous_idf_discovery(self):
        requester = Mock()
        requester.get.return_value = response(payload={
            "data": [
                {"Name": "one.idf.txt", "Type": "IDF File", "path": "one.idf.txt"},
                {"Name": "two.idf.txt", "Type": "IDF File", "path": "two.idf.txt"},
                {"Name": "study.sdrf.txt", "Type": "SDRF File", "path": "study.sdrf.txt"},
            ]
        })

        with self.assertRaisesRegex(ValueError, "exactly one IDF"):
            AEWebFetcher(requester=requester).resolve("E-MTAB-1")


if __name__ == "__main__":
    unittest.main()
