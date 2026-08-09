# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================

"""Bounded XML response reading and DTD/entity-rejecting parsing."""

from __future__ import annotations

import re
from typing import Any
import xml.etree.ElementTree as ET


_UNSAFE_DECLARATION = re.compile(br"<!\s*(?:DOCTYPE|ENTITY)\b", re.IGNORECASE)


class XMLSafetyError(ValueError):
    """Base class for rejected external XML content."""


class XMLSizeLimitError(XMLSafetyError):
    """External XML or its transport body exceeded a configured byte ceiling."""


class UnsafeXMLDocumentError(XMLSafetyError):
    """External XML declared a DTD or entity."""


def parse_xml(value: str | bytes, *, max_bytes: int) -> ET.Element:
    if max_bytes <= 0:
        raise ValueError("max_bytes must be positive")
    encoded = value.encode("utf-8") if isinstance(value, str) else value
    if len(encoded) > max_bytes:
        raise XMLSizeLimitError(
            f"XML document has {len(encoded)} bytes and exceeds the "
            f"{max_bytes} byte XML limit."
        )
    if _UNSAFE_DECLARATION.search(encoded):
        raise UnsafeXMLDocumentError(
            "XML DTD/entity declarations are not allowed."
        )
    return ET.fromstring(encoded)


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
    if declared is not None and len(body) != declared:
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
    if declared is not None and byte_count != declared:
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
