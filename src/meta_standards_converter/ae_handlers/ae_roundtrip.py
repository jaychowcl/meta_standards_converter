# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Lossless sidecar support for MAGE-TAB-origin JSON packages."""

from __future__ import annotations

import copy
import hashlib
import json
import logging


logger = logging.getLogger(__name__)


def semantic_sha256(package: dict) -> str:
    semantic = {key: value for key, value in package.items() if key != "mage_tab"}
    payload = json.dumps(
        semantic, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def model_sha256(package: dict) -> str | None:
    mage_tab = package.get("mage_tab") if isinstance(package, dict) else None
    model = mage_tab.get("model") if isinstance(mage_tab, dict) else None
    if not isinstance(model, dict):
        return None
    payload = json.dumps(
        model, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def build_roundtrip(package: dict, idf_rows: list, sdrfs: list[tuple[str, list]]) -> dict:
    result = {
        "schema_version": 1,
        "semantic_sha256": semantic_sha256(package),
        "idf_rows": copy.deepcopy(idf_rows),
        "sdrfs": [
            {"name": name, "rows": copy.deepcopy(rows)} for name, rows in sdrfs
        ],
    }
    fingerprint = model_sha256(package)
    if fingerprint:
        result["model_sha256"] = fingerprint
    return result


def unchanged_magetab(package: dict) -> list | None:
    roundtrip = _roundtrip(package)
    if not roundtrip or roundtrip.get("semantic_sha256") != semantic_sha256(package):
        return None
    expected_model = roundtrip.get("model_sha256")
    if expected_model is not None and expected_model != model_sha256(package):
        return None
    idf_rows = roundtrip.get("idf_rows")
    sdrfs = roundtrip.get("sdrfs")
    if not _table(idf_rows) or not isinstance(sdrfs, list) or len(sdrfs) != 1:
        return None
    sdrf_rows = sdrfs[0].get("rows") if isinstance(sdrfs[0], dict) else None
    if not _table(sdrf_rows):
        return None
    rows = copy.deepcopy(idf_rows)
    for index, row in enumerate(rows):
        if row and _label(row[0]) == "sdrffile":
            rows[index] = ["SDRF File", copy.deepcopy(sdrf_rows), *row[2:]]
            return rows
    return None


def restore_extensions(package: dict, magetab: list) -> list:
    """Restore unsupported source metadata after mapped JSON has been rendered."""
    mage_tab = package.get("mage_tab") if isinstance(package, dict) else None
    if not isinstance(mage_tab, dict):
        return magetab
    rows = copy.deepcopy(magetab)
    known_labels = {_label(row[0]) for row in rows if row}
    sdrf_index = next(
        (index for index, row in enumerate(rows) if row and _label(row[0]) == "sdrffile"),
        None,
    )
    insert_at = sdrf_index if sdrf_index is not None else len(rows)
    for item in mage_tab.get("unmapped_idf_rows", []) or []:
        if not isinstance(item, dict) or not item.get("label"):
            continue
        if _label(item["label"]) in known_labels:
            continue
        rows.insert(insert_at, [item["label"], *(item.get("values") or [])])
        insert_at += 1

    sdrf_index = next(
        (index for index, row in enumerate(rows) if row and _label(row[0]) == "sdrffile"),
        None,
    )
    if sdrf_index is None:
        return rows
    sdrf = rows[sdrf_index][1] if len(rows[sdrf_index]) > 1 else None
    if not _table(sdrf):
        return rows
    for item in mage_tab.get("unmapped_sdrf_columns", []) or []:
        if not isinstance(item, dict) or not item.get("header"):
            continue
        header = item["header"]
        occurrence = int(item.get("occurrence") or 1)
        if sum(_label(value) == _label(header) for value in sdrf[0]) >= occurrence:
            continue
        values = list(item.get("values") or [])
        restored = _match_column_values(package, sdrf, item, values)
        sdrf[0].append(header)
        for index, row in enumerate(sdrf[1:]):
            row.append(restored[index] if index < len(restored) else "")
    return rows


def _match_column_values(package, generated, item, values):
    if len(values) == len(generated) - 1:
        return values
    roundtrip = _roundtrip(package)
    source = None
    for candidate in (roundtrip or {}).get("sdrfs", []) or []:
        if isinstance(candidate, dict) and candidate.get("name") == item.get("file"):
            source = candidate.get("rows")
            break
    if not _table(source):
        logger.warning("Could not safely restore SDRF column %s after JSON edits.", item.get("header"))
        return [""] * (len(generated) - 1)
    source_identity = _identity_index(source[0])
    generated_identity = _identity_index(generated[0])
    if source_identity is None or generated_identity is None:
        return [""] * (len(generated) - 1)
    by_identity = {}
    for index, row in enumerate(source[1:]):
        if index < len(values) and source_identity < len(row):
            by_identity.setdefault(row[source_identity], []).append(values[index])
    restored = []
    for row in generated[1:]:
        distinct = list(dict.fromkeys(by_identity.get(row[generated_identity], [])))
        if len(distinct) <= 1:
            restored.append(distinct[0] if distinct else "")
        else:
            logger.warning(
                "Could not safely restore varying SDRF column %s for %s after JSON edits.",
                item.get("header"), row[generated_identity],
            )
            restored.append("")
    return restored


def _roundtrip(package):
    mage_tab = package.get("mage_tab") if isinstance(package, dict) else None
    value = mage_tab.get("roundtrip") if isinstance(mage_tab, dict) else None
    if not isinstance(value, dict) or value.get("schema_version") != 1:
        return None
    return value


def _table(value):
    return isinstance(value, list) and bool(value) and all(isinstance(row, list) for row in value)


def _label(value):
    return "".join(str(value).split()).casefold()


def _identity_index(header):
    for label in ("Sample Name", "Source Name", "Assay Name"):
        for index, value in enumerate(header):
            if _label(value) == _label(label):
                return index
    return None
