# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
import io
import os
import sys
import tarfile
import unittest
from unittest.mock import Mock


ROOT = os.path.dirname(os.path.dirname(__file__))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from meta_standards_converter.geo_handlers.geo_webfetcher import GEOWebFetcher  # noqa: E402


def miniml_archive(gse: str, content: str) -> bytes:
    buffer = io.BytesIO()
    encoded = content.encode("utf-8")
    with tarfile.open(fileobj=buffer, mode="w:gz") as tar:
        info = tarfile.TarInfo(f"{gse}_family.xml")
        info.size = len(encoded)
        tar.addfile(info, io.BytesIO(encoded))
    return buffer.getvalue()


class TestGEOWebFetcher(unittest.TestCase):
    def test_fetch_gse_miniml_uses_requester_and_extracts_xml(self):
        response = Mock(content=miniml_archive("GSE1", "<MINiML />"))
        response.raise_for_status = Mock()
        requester = Mock()
        requester.get.return_value = response

        result = GEOWebFetcher(requester=requester).fetch_gse_miniml(gse="GSE1")

        requester.get.assert_called_once_with(
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSEnnn/GSE1/miniml/GSE1_family.xml.tgz"
        )
        response.raise_for_status.assert_called_once()
        self.assertEqual("<MINiML />", result)

    def test_url_gse_miniml_rejects_non_gse_accession(self):
        with self.assertRaises(ValueError):
            GEOWebFetcher().url_gse_miniml("SRX1")


if __name__ == "__main__":
    unittest.main()
