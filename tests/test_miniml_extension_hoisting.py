from __future__ import annotations

import pytest

from meta_standards_converter.miniml import MINiMLCodec, MINiMLModelError


def _package() -> dict:
    return {
        "miniml_schema_version": "3.0",
        "source": {"format": "test"},
        "sample": [],
        "series": {"iid": "GSE1"},
    }


def test_every_legacy_series_extension_is_hoisted_to_package_extensions() -> None:
    payload = _package()
    payload["series"]["extensions"] = {
        "publication_inheritance": {"source": "GSE0"},
        "vendor_note": {"retained": True},
    }
    encoded = MINiMLCodec.encode(MINiMLCodec().decode(payload, strict=True).package)
    assert "extensions" not in encoded["series"]
    assert encoded["extensions"] == payload["series"]["extensions"]


def test_identical_root_and_series_extensions_are_deduplicated() -> None:
    payload = _package()
    payload["extensions"] = {"vendor_note": {"retained": True}}
    payload["series"]["extensions"] = {"vendor_note": {"retained": True}}
    encoded = MINiMLCodec.encode(MINiMLCodec().decode(payload, strict=True).package)
    assert encoded["extensions"] == {"vendor_note": {"retained": True}}
    assert "extensions" not in encoded["series"]


def test_conflicting_root_and_series_extensions_fail_closed() -> None:
    payload = _package()
    payload["extensions"] = {"vendor_note": 1}
    payload["series"]["extensions"] = {"vendor_note": 2}
    with pytest.raises(MINiMLModelError, match="conflicting package and series extension"):
        MINiMLCodec().decode(payload, strict=True)


def test_reserved_harmonization_extension_cannot_come_from_legacy_series() -> None:
    payload = _package()
    payload["series"]["extensions"] = {
        "msc_harmonization": {"schema_version": "1.0", "patches": []}
    }
    with pytest.raises(MINiMLModelError, match="reserved"):
        MINiMLCodec().decode(payload, strict=True)

