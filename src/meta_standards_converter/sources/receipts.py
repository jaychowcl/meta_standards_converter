# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Bind a dataset's scientific disposition to an unchanged bare MINiML file."""
from __future__ import annotations

from contextlib import contextmanager
from copy import deepcopy
from hashlib import file_digest, sha256
import json
import os
from pathlib import Path
from tempfile import TemporaryDirectory

from meta_standards_converter.runtime_contracts import OperationStatusV2

SUFFIX = ".msc-receipt"
FORMAT = "msc.dataset-receipt/1"
PARTIAL_CODES = {"source_partial", "source_fallback", "source_failed", "conversion_partial", "sample_skipped"}


def status_axes(history):
    """Return conservative source axes while retaining stage detail in the receipt."""
    states = [OperationStatusV2.from_dict(stage["status"]) for stage in history if stage.get("status")]
    if not states:
        return dict(execution="unknown", completeness="unknown", validation="unknown", publication="unknown")
    priorities = {"execution": ("failed", "degraded", "not_run", "succeeded"),
                  "completeness": ("empty", "partial", "unknown", "complete"),
                  "validation": ("invalid", "not_run", "valid"),
                  "publication": ("blocked", "review_required", "publishable", "not_requested")}
    axes = {axis: next(v for v in order if any(getattr(s, axis).value == v for s in states))
            for axis, order in priorities.items()}
    if any(stage.get("stage") == "msc.prior_status_unknown" for stage in history):
        if axes["publication"] in {"publishable", "not_requested"}:
            axes["publication"] = "unknown"
    return axes


def atlas_history(dataset, publications):
    stages = []
    def add(name, operation, policy=None):
        if operation:
            stages.append({"stage": name, "status": deepcopy(operation)})
        if policy and policy.get("review_required"):
            stages.append({"stage": name + ".policy", "status": _status(
                publication="review_required", reason="source_review_required")})
    add("atlas.metadata", dataset.get("metadata_status"))
    review = dataset.get("review") or {}
    add("atlas.review", review.get("status"), review.get("policy"))
    hz = dataset.get("harmonization") or {}
    add("atlas.harmonization", hz.get("status"))
    if hz.get("degraded_stages") and not hz.get("status"):
        stages.append({"stage": "atlas.harmonization", "status": _status(
            partial=True, publication="review_required", reason="legacy_degraded_harmonization")})
    for pid in dataset.get("publication_ids", ()):
        review = publications.get(pid) or {}
        add("atlas.publication." + pid, review.get("status"), review.get("policy"))
    if dataset.get("diagnostics"):
        stages.append({"stage": "atlas.diagnostics", "diagnostics": deepcopy(dataset["diagnostics"])})
    return tuple(stages)


def _status(*, partial=False, invalid=False, publication, reason):
    return {"contract_version": "2.0", "execution": "succeeded", "completeness": "partial" if partial else "complete",
            "validation": "invalid" if invalid else "valid", "evidence_confidence": "not_assessed",
            "publication": publication, "terminal_reason": reason, "errors": []}


def conversion_records(groups, item):
    result = {}
    for group in groups:
        sid = group.dataset_id
        diagnostics = [d for d in item.diagnostics if
                       (d.dataset_id or getattr(d.curation, "dataset_id", None)) == sid
                       or (len(groups) == 1 and not d.dataset_id and not d.curation)]
        invalid = any(d.code == "projection_invalid" or
                      (d.code.startswith("integrity_") and d.severity == "error") for d in diagnostics)
        partial = any(d.code in PARTIAL_CODES for d in diagnostics)
        publication = "blocked" if invalid else "review_required" if partial else "not_requested"
        record = {"stage": "msc.conversion", "status": _status(partial=partial, invalid=invalid,
                  publication=publication, reason="metadata_converted"),
                  "diagnostics": [d.to_dict() for d in diagnostics],
                  "preparation": [deepcopy(x) for x in item.preparation if x.get("dataset_id") == sid]}
        history = list(group.source_history)
        if not history or record != history[-1]:
            history.append(record)
        result[sid] = {"dataset_id": sid, "package_count": len(group.packages), "history": history}
    return result


def receipt_path(path):
    return Path(str(path) + SUFFIX)


@contextmanager
def _writer_lock(path):
    # Do not unlink the lock file. Waiting publishers must compete on one inode.
    lock = path.with_name("." + path.name + ".msc-receipt.lock")
    with lock.open("a+b") as stream:
        if os.name == "nt":
            import msvcrt
            stream.seek(0)
            msvcrt.locking(stream.fileno(), msvcrt.LK_LOCK, 1)
        else:
            import fcntl
            fcntl.flock(stream, fcntl.LOCK_EX)
        yield


def publish_with_receipt(path, payload, records, *, overwrite):
    from meta_standards_converter.converters.archive_results import publish_json
    from meta_standards_converter.operational_events import _guard_output
    path = Path(path)
    receipt = receipt_path(path)
    _guard_output([path, receipt])
    path.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=path.parent, prefix=".msc-save-") as directory:
        data, marker = Path(directory) / "payload", Path(directory) / "receipt"
        publish_json(data, payload)
        with data.open("rb") as stream:
            digest = file_digest(stream, "sha256").hexdigest()
        publish_json(marker, {"format": FORMAT, "sha256": digest, "datasets": records})
        with _writer_lock(path):
            if not overwrite and (path.exists() or receipt.exists()):
                raise FileExistsError(f"Managed output already exists: {path}")
            # An interrupted overwrite can expose data with unknown history. Remove
            # even an equal-payload old receipt before installing the new generation.
            receipt.unlink(missing_ok=True)
            if overwrite:
                os.replace(data, path)
            else:
                os.link(data, path)
            os.replace(marker, receipt)


def read_records(path, raw, packages):
    receipt = receipt_path(path)
    if not receipt.is_file():
        return None
    try:
        data = json.loads(receipt.read_bytes())
        if data["format"] != FORMAT or data["sha256"] != sha256(raw).hexdigest():
            raise ValueError("payload digest or format mismatch")
        records = data["datasets"]
        if not isinstance(records, list) or not records:
            raise ValueError("missing dataset records")
        total, ids = 0, set()
        for record in records:
            key, count, history = record["dataset_id"], record["package_count"], record["history"]
            if not isinstance(key, str) or not key or key in ids or type(count) is not int or count < 1:
                raise ValueError("invalid dataset binding")
            if not isinstance(history, list) or not history:
                raise ValueError("missing dataset history")
            status_axes(history)
            ids.add(key)
            total += count
        if not isinstance(packages, list) or total != len(packages):
            raise ValueError("package cardinality mismatch")
        return records
    except (ValueError, TypeError, KeyError, AttributeError) as exc:
        raise ValueError(f"Invalid MSC receipt for {Path(path).name}: {exc}") from exc
