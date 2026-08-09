# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================

from __future__ import annotations

import pytest

from meta_standards_converter.xml_safety import (
    UnsafeXMLDocumentError,
    XMLSizeLimitError,
    parse_xml,
    read_limited_response,
)


class _Response:
    headers: dict[str, str]

    def __init__(self, chunks: tuple[bytes, ...], length: int | None = None) -> None:
        self.chunks = chunks
        self.headers = {} if length is None else {"Content-Length": str(length)}

    def iter_content(self, chunk_size: int):
        yield from self.chunks


def test_parse_xml_rejects_dtd_and_entity_declarations() -> None:
    payload = b'<!DOCTYPE x [<!ENTITY secret SYSTEM "file:///etc/passwd">]><x>&secret;</x>'

    with pytest.raises(UnsafeXMLDocumentError, match="DTD/entity"):
        parse_xml(payload, max_bytes=1024)


def test_parse_xml_enforces_utf8_byte_limit_before_parsing() -> None:
    with pytest.raises(XMLSizeLimitError, match="exceeds the 8 byte XML limit"):
        parse_xml("<root>é</root>", max_bytes=8)


def test_limited_response_rejects_declared_and_chunked_oversize() -> None:
    with pytest.raises(XMLSizeLimitError, match="declares 9 bytes"):
        read_limited_response(_Response((), length=9), max_bytes=8)

    with pytest.raises(XMLSizeLimitError, match="exceeds the 8 byte response limit"):
        read_limited_response(
            _Response((b"1234", b"5678", b"9")),
            max_bytes=8,
        )
