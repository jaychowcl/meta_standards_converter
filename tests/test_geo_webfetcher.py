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
from unittest.mock import Mock, patch


ROOT = os.path.dirname(os.path.dirname(__file__))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from meta_standards_converter.geo_handlers.geo_webfetcher import GEOWebFetcher  # noqa: E402


def miniml_archive(
    gse: str,
    content: str,
    *,
    extra_members: tuple[tuple[str, bytes, bytes | None], ...] = (),
) -> bytes:
    buffer = io.BytesIO()
    encoded = content.encode("utf-8")
    with tarfile.open(fileobj=buffer, mode="w:gz") as tar:
        info = tarfile.TarInfo(f"{gse}_family.xml")
        info.size = len(encoded)
        tar.addfile(info, io.BytesIO(encoded))
        for name, member_type, payload in extra_members:
            info = tarfile.TarInfo(name)
            info.type = member_type
            if member_type in {tarfile.SYMTYPE, tarfile.LNKTYPE}:
                info.linkname = str((payload or b"target").decode("utf-8"))
                tar.addfile(info)
                continue
            encoded_member = payload or b""
            info.size = len(encoded_member)
            tar.addfile(info, io.BytesIO(encoded_member))
    return buffer.getvalue()


class TestGEOWebFetcher(unittest.TestCase):
    def test_fetch_gse_miniml_uses_requester_and_extracts_xml(self):
        content = miniml_archive("GSE1", "<MINiML />")
        response = Mock(
            headers={"Content-Length": str(len(content))},
            iter_content=Mock(return_value=iter([content])),
        )
        response.raise_for_status = Mock()
        requester = Mock()
        requester.get.return_value = response

        result = GEOWebFetcher(requester=requester).fetch_gse_miniml(gse="GSE1")

        requester.get.assert_called_once_with(
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSEnnn/GSE1/miniml/GSE1_family.xml.tgz",
            stream=True,
        )
        response.raise_for_status.assert_called_once()
        self.assertEqual("<MINiML />", result)

    def test_fetch_gse_miniml_allows_safe_auxiliary_archive_members(self):
        content = miniml_archive(
            "GSE1",
            "<MINiML />",
            extra_members=(
                ("tables", tarfile.DIRTYPE, None),
                ("tables/GSM1-tbl-1.txt", tarfile.REGTYPE, b"ID_REF\tVALUE\n"),
            ),
        )
        response = Mock(
            headers={"Content-Length": str(len(content))},
            iter_content=Mock(return_value=iter([content])),
        )
        response.raise_for_status = Mock()
        requester = Mock()
        requester.get.return_value = response

        self.assertEqual(
            "<MINiML />",
            GEOWebFetcher(requester=requester).fetch_gse_miniml(gse="GSE1"),
        )

    def test_fetch_gse_miniml_rejects_unexpected_xml_and_unsafe_members(self):
        cases = (
            (("other.xml", tarfile.REGTYPE, b"<other />"), "unexpected XML"),
            (("../escape.txt", tarfile.REGTYPE, b"escape"), "unsafe path"),
            (("link", tarfile.SYMTYPE, b"target"), "unsafe member type"),
        )
        for extra_member, message in cases:
            with self.subTest(extra_member=extra_member[0]):
                content = miniml_archive(
                    "GSE1",
                    "<MINiML />",
                    extra_members=(extra_member,),
                )
                response = Mock(
                    headers={"Content-Length": str(len(content))},
                    iter_content=Mock(return_value=iter([content])),
                )
                response.raise_for_status = Mock()
                requester = Mock()
                requester.get.return_value = response

                with self.assertRaisesRegex(ValueError, message):
                    GEOWebFetcher(requester=requester).fetch_gse_miniml(gse="GSE1")

    def test_fetch_gse_miniml_rejects_duplicate_names_and_member_overflow(self):
        duplicate = miniml_archive(
            "GSE1",
            "<MINiML />",
            extra_members=(("GSE1_family.xml", tarfile.REGTYPE, b"<MINiML />"),),
        )
        response = Mock(
            headers={"Content-Length": str(len(duplicate))},
            iter_content=Mock(return_value=iter([duplicate])),
        )
        response.raise_for_status = Mock()
        requester = Mock()
        requester.get.return_value = response
        with self.assertRaisesRegex(ValueError, "duplicate member"):
            GEOWebFetcher(requester=requester).fetch_gse_miniml(gse="GSE1")

        overflow = miniml_archive(
            "GSE1",
            "<MINiML />",
            extra_members=(("one.txt", tarfile.REGTYPE, b"1"),),
        )
        response = Mock(
            headers={"Content-Length": str(len(overflow))},
            iter_content=Mock(return_value=iter([overflow])),
        )
        response.raise_for_status = Mock()
        requester.get.return_value = response
        with patch(
            "meta_standards_converter.geo_handlers.geo_webfetcher.MAX_GEO_ARCHIVE_MEMBERS",
            1,
        ):
            with self.assertRaisesRegex(ValueError, "member-count limit"):
                GEOWebFetcher(requester=requester).fetch_gse_miniml(gse="GSE1")

    def test_url_gse_miniml_rejects_non_gse_accession(self):
        with self.assertRaises(ValueError):
            GEOWebFetcher().url_gse_miniml("SRX1")


if __name__ == "__main__":
    unittest.main()
