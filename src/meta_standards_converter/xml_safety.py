# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================

"""Bounded XML response reading and fail-closed external-DTD parsing."""

from __future__ import annotations

import re
from typing import Any
import xml.etree.ElementTree as ET


_DOCTYPE_START = re.compile(br"<!\s*DOCTYPE\b", re.IGNORECASE)
_ENTITY_START = re.compile(br"<!\s*ENTITY\b", re.IGNORECASE)
_XML_NAME = br"[A-Za-z_][A-Za-z0-9_.:-]*"
_QUOTED_LITERAL = br'(?:"[^"]*"|\'[^\']*\')'
_EXTERNAL_DOCTYPE = re.compile(
    br"<!\s*DOCTYPE\s+"
    + _XML_NAME
    + br"\s+(?:SYSTEM\s+"
    + _QUOTED_LITERAL
    + br"|PUBLIC\s+"
    + _QUOTED_LITERAL
    + br"\s+"
    + _QUOTED_LITERAL
    + br")\s*>",
    re.IGNORECASE | re.DOTALL,
)
_XML_PROLOG_PREFIX = re.compile(
    br"(?:\s|<\?[^>]*\?>|<!--(?:(?!-->).)*-->)*",
    re.DOTALL,
)


class XMLSafetyError(ValueError):
    """Base class for rejected external XML content."""


class XMLSizeLimitError(XMLSafetyError):
    """External XML or its transport body exceeded a configured byte ceiling."""


class UnsafeXMLDocumentError(XMLSafetyError):
    """External XML declared an unsafe or malformed DTD/entity construct."""


def _strip_safe_external_doctype(encoded: bytes) -> bytes:
    """Strip one ordinary external DTD without resolving or parsing it.

    NCBI XML commonly declares a SYSTEM or PUBLIC DTD.  ElementTree does not
    need that declaration for the subset of XML consumed here, so accepting
    and removing its identifier is safer than allowing a parser to resolve it.
    Internal subsets and every entity declaration remain forbidden.
    """

    if _ENTITY_START.search(encoded):
        raise UnsafeXMLDocumentError(
            "XML contains an unsafe DTD/entity declaration."
        )
    starts = list(_DOCTYPE_START.finditer(encoded))
    if not starts:
        return encoded
    if len(starts) != 1:
        raise UnsafeXMLDocumentError(
            "XML contains multiple DTD/entity declarations."
        )

    start = starts[0].start()
    quote: int | None = None
    end: int | None = None
    for index in range(starts[0].end(), len(encoded)):
        byte = encoded[index]
        if quote is not None:
            if byte == quote:
                quote = None
            continue
        if byte in (ord('"'), ord("'")):
            quote = byte
        elif byte == ord("["):
            raise UnsafeXMLDocumentError(
                "XML internal DTD/entity subsets are not allowed."
            )
        elif byte == ord(">"):
            end = index + 1
            break

    if end is None:
        raise UnsafeXMLDocumentError(
            "XML contains a malformed DTD/entity declaration."
        )
    declaration = encoded[start:end]
    if _EXTERNAL_DOCTYPE.fullmatch(declaration) is None:
        raise UnsafeXMLDocumentError(
            "XML contains an unsafe or malformed DTD/entity declaration."
        )
    if _XML_PROLOG_PREFIX.fullmatch(encoded[:start]) is None:
        raise UnsafeXMLDocumentError(
            "XML DTD/entity declaration is outside the document prolog."
        )
    return encoded[:start] + encoded[end:]


def _content_length_describes_decoded_body(headers: Any) -> bool:
    """Return whether requests' decoded body should equal Content-Length."""

    raw_encoding = (
        headers.get("Content-Encoding") if hasattr(headers, "get") else None
    )
    return str(raw_encoding or "").strip().casefold() in {"", "identity"}


def parse_xml(value: str | bytes, *, max_bytes: int) -> ET.Element:
    if max_bytes <= 0:
        raise ValueError("max_bytes must be positive")
    encoded = value.encode("utf-8") if isinstance(value, str) else value
    if len(encoded) > max_bytes:
        raise XMLSizeLimitError(
            f"XML document has {len(encoded)} bytes and exceeds the "
            f"{max_bytes} byte XML limit."
        )
    return ET.fromstring(_strip_safe_external_doctype(encoded))


def read_limited_response(
    response: Any,
    *,
    max_bytes: int,
    chunk_size: int = 1024 * 1024,
) -> bytes:
    if max_bytes <= 0:
        raise ValueError("max_bytes must be positive")
    headers = getattr(response, "headers", {})
    raw_length = headers.get("Content-Length") if hasattr(headers, "get") else None
    declared = None
    if raw_length not in (None, ""):
        try:
            declared = int(raw_length)
        except (TypeError, ValueError) as error:
            raise XMLSizeLimitError("Response Content-Length is invalid.") from error
        if declared < 0:
            raise XMLSizeLimitError("Response Content-Length is invalid.")
        if declared > max_bytes:
            raise XMLSizeLimitError(
                f"Response declares {declared} bytes and exceeds the "
                f"{max_bytes} byte response limit."
            )
    body = bytearray()
    for chunk in response.iter_content(chunk_size=chunk_size):
        if not chunk:
            continue
        body.extend(chunk)
        if len(body) > max_bytes:
            raise XMLSizeLimitError(
                f"Response exceeds the {max_bytes} byte response limit."
            )
    if (
        declared is not None
        and _content_length_describes_decoded_body(headers)
        and len(body) != declared
    ):
        raise XMLSizeLimitError(
            f"Response body has {len(body)} bytes; Content-Length declared {declared}."
        )
    return bytes(body)


def stream_limited_response(
    response: Any,
    destination,
    *,
    max_bytes: int,
    chunk_size: int = 1024 * 1024,
) -> int:
    """Stream a response into an open binary file under a hard byte ceiling."""

    if max_bytes <= 0:
        raise ValueError("max_bytes must be positive")
    headers = getattr(response, "headers", {})
    raw_length = headers.get("Content-Length") if hasattr(headers, "get") else None
    declared = None
    if raw_length not in (None, ""):
        try:
            declared = int(raw_length)
        except (TypeError, ValueError) as error:
            raise XMLSizeLimitError("Response Content-Length is invalid.") from error
        if declared < 0 or declared > max_bytes:
            raise XMLSizeLimitError(
                f"Response declares {declared} bytes and exceeds the "
                f"{max_bytes} byte response limit."
            )
    byte_count = 0
    for chunk in response.iter_content(chunk_size=chunk_size):
        if not chunk:
            continue
        byte_count += len(chunk)
        if byte_count > max_bytes:
            raise XMLSizeLimitError(
                f"Response exceeds the {max_bytes} byte response limit."
            )
        destination.write(chunk)
    if (
        declared is not None
        and _content_length_describes_decoded_body(headers)
        and byte_count != declared
    ):
        raise XMLSizeLimitError(
            f"Response body has {byte_count} bytes; Content-Length declared {declared}."
        )
    return byte_count


__all__ = [
    "UnsafeXMLDocumentError",
    "XMLSafetyError",
    "XMLSizeLimitError",
    "parse_xml",
    "read_limited_response",
    "stream_limited_response",
]
