# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
import os
import json
import sys
import tempfile
import unittest
from unittest.mock import Mock, call


ROOT = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from meta_standards_converter.sources.magetab import AEWebFetcher  # noqa: E402
from meta_standards_converter.retrieval import (  # noqa: E402
    RetrievalPolicy,
    RetrievalSecurityError,
)
from meta_standards_converter.runtime_contracts import get_resource_profile  # noqa: E402
from meta_standards_converter.xml_safety import XMLSizeLimitError  # noqa: E402


IDF = "MAGE-TAB Version\t1.1\nInvestigation Accession\tE-MTAB-1\nSDRF File\tstudy.sdrf.txt\n"
SDRF = "Source Name\tCharacteristics[organism]\nS1\tHomo sapiens\n"


def response(*, text=None, payload=None, headers=None, status_code=200):
    item = Mock()
    body = (
        text.encode("utf-8")
        if text is not None
        else json.dumps(payload).encode("utf-8")
        if payload is not None
        else b""
    )
    item.text = text
    item.content = body
    item.headers = headers or {}
    item.status_code = status_code
    item.iter_content.side_effect = lambda chunk_size: [body]
    item.json.return_value = payload
    item.raise_for_status = Mock()
    item.close = Mock()
    return item


def public_policy(profile=None):
    profile = profile or get_resource_profile("standard")
    return RetrievalPolicy(
        resource_profile=profile,
        allowed_hosts=frozenset({"example.org", "override.example"}),
        allowed_host_suffixes=frozenset({"ebi.ac.uk"}),
        allowed_schemes=frozenset({"https"}),
        resolver=lambda host, port, **kwargs: [
            (None, None, None, None, ("93.184.216.34", port))
        ],
    )


class TestAEWebFetcher(unittest.TestCase):
    def test_typed_resource_profile_configures_requester_and_is_preserved(self):
        profile = get_resource_profile(
            "standard",
            overrides={"max_xml_bytes": 4096, "network_workers": 3},
        )

        fetcher = AEWebFetcher(resource_profile=profile)

        self.assertIs(profile, fetcher.resource_profile)
        self.assertEqual((10, 60), fetcher.requester.settings.timeout)
        self.assertEqual(3, fetcher.requester.settings.max_in_flight)

    def test_local_metadata_is_rejected_before_reading_past_profile_limit(self):
        profile = get_resource_profile(
            "standard",
            overrides={"max_xml_bytes": len(IDF.encode("utf-8")) - 1},
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            idf_path = os.path.join(tmpdir, "study.idf.txt")
            with open(idf_path, "w", encoding="utf-8") as handle:
                handle.write(IDF)

            with self.assertRaises(XMLSizeLimitError):
                AEWebFetcher(resource_profile=profile).resolve(idf_path)

    def test_remote_metadata_is_streamed_under_profile_limit(self):
        profile = get_resource_profile(
            "standard",
            overrides={"max_xml_bytes": len(IDF.encode("utf-8")) - 1},
        )
        requester = Mock()
        requester.get.return_value = response(
            text=IDF,
            headers={"Content-Length": str(len(IDF.encode("utf-8")))},
        )

        with self.assertRaises(XMLSizeLimitError):
            AEWebFetcher(
                requester=requester,
                resource_profile=profile,
                retrieval_policy=public_policy(profile),
            ).resolve("https://example.org/study.idf.txt")

        requester.get.assert_called_once_with(
            "https://example.org/study.idf.txt",
            stream=True,
            allow_redirects=False,
        )

    def test_http_metadata_rejects_unapproved_egress_before_request(self):
        profile = get_resource_profile("standard")
        requester = Mock()
        fetcher = AEWebFetcher(
            requester=requester,
            resource_profile=profile,
            retrieval_policy=public_policy(profile),
        )

        with self.assertRaises(RetrievalSecurityError):
            fetcher.resolve("https://attacker.invalid/study.idf.txt")

        requester.get.assert_not_called()

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

        resolved = AEWebFetcher(
            requester=requester,
            retrieval_policy=public_policy(),
        ).resolve(
            "https://example.org/studies/study.idf.txt"
        )

        self.assertEqual(
            [
                call(
                    "https://example.org/studies/study.idf.txt",
                    stream=True,
                    allow_redirects=False,
                ),
                call(
                    "https://example.org/studies/study.sdrf.txt",
                    stream=True,
                    allow_redirects=False,
                ),
            ],
            requester.get.call_args_list,
        )
        self.assertEqual(SDRF, resolved.sdrfs[0].text)

    def test_explicit_sdrf_override_replaces_idf_references(self):
        requester = Mock()
        requester.get.side_effect = [response(text=IDF), response(text=SDRF)]

        resolved = AEWebFetcher(
            requester=requester,
            retrieval_policy=public_policy(),
        ).resolve(
            "https://example.org/study.idf.txt",
            sdrf_sources=["https://override.example/study.sdrf.txt"],
        )

        self.assertEqual(
            call(
                "https://override.example/study.sdrf.txt",
                stream=True,
                allow_redirects=False,
            ),
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

        def get(url, **kwargs):
            if url.endswith("/api/v1/files/E-MTAB-1"):
                self.assertEqual({"start": 0, "length": 100}, kwargs["params"])
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
        resolved = AEWebFetcher(
            requester=requester,
            retrieval_policy=public_policy(),
        ).resolve("E-MTAB-1")

        self.assertEqual("accession", resolved.source_kind)
        self.assertEqual("E-MTAB-1", resolved.source)
        self.assertEqual(["part1.sdrf.txt", "part2.sdrf.txt"], [x.name for x in resolved.sdrfs])
        self.assertEqual(5, requester.get.call_count)

    def test_accession_paginates_biostudies_file_discovery(self):
        requester = Mock()
        base = "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/001/E-MTAB-1"
        pages = {
            0: {
                "recordsFiltered": 2,
                "data": [
                    {
                        "Name": "E-MTAB-1.idf.txt",
                        "Type": "IDF File",
                        "path": "E-MTAB-1.idf.txt",
                    }
                ],
            },
            1: {
                "recordsFiltered": 2,
                "data": [
                    {
                        "Name": "E-MTAB-1.sdrf.txt",
                        "Type": "SDRF File",
                        "path": "E-MTAB-1.sdrf.txt",
                    }
                ],
            },
        }

        def get(url, **kwargs):
            if url.endswith("/api/v1/files/E-MTAB-1"):
                return response(payload=pages[kwargs["params"]["start"]])
            if url.endswith("/api/v1/studies/E-MTAB-1/info"):
                return response(payload={"httpLink": base})
            if url.endswith("/Files/E-MTAB-1.idf.txt"):
                return response(text=IDF)
            if url.endswith("/Files/E-MTAB-1.sdrf.txt"):
                return response(text=SDRF)
            raise AssertionError(url)

        requester.get.side_effect = get

        resolved = AEWebFetcher(
            requester=requester,
            retrieval_policy=public_policy(),
        ).resolve("E-MTAB-1")

        self.assertEqual("E-MTAB-1.idf.txt", resolved.idf.name)
        self.assertEqual(["E-MTAB-1.sdrf.txt"], [item.name for item in resolved.sdrfs])
        file_calls = [
            call.kwargs["params"]
            for call in requester.get.call_args_list
            if call.args[0].endswith("/api/v1/files/E-MTAB-1")
        ]
        self.assertEqual(
            [{"start": 0, "length": 100}, {"start": 1, "length": 100}],
            file_calls,
        )

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
            AEWebFetcher(
                requester=requester,
                retrieval_policy=public_policy(),
            ).resolve("E-MTAB-1")


if __name__ == "__main__":
    unittest.main()
