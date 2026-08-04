# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
import json
import os
from pathlib import Path

import pytest

import meta_standards_converter.artifact_bundle as artifact_bundle

from meta_standards_converter.artifact_bundle import (
    ArtifactRecoveryError,
    DurableArtifactBundlePublisher,
    resolve_current_bundle,
)


def _staged(tmp_path: Path, values: dict[str, str]) -> dict[str, Path]:
    stage = tmp_path / "stage"
    stage.mkdir(exist_ok=True)
    result = {}
    for role, value in values.items():
        path = stage / f"{role}.txt"
        path.write_text(value, encoding="utf-8")
        result[role] = path
    return result


def test_generation_pointer_commits_one_immutable_bundle_and_legacy_views(tmp_path):
    destinations = {
        "table": tmp_path / "out" / "table.tsv",
        "manifest": tmp_path / "out" / "manifest.json",
    }
    publisher = DurableArtifactBundlePublisher()

    first = publisher.publish(
        _staged(tmp_path, {"table": "old-table", "manifest": "old-manifest"}),
        destinations,
        overwrite=False,
    )
    second = publisher.publish(
        _staged(tmp_path, {"table": "new-table", "manifest": "new-manifest"}),
        destinations,
        overwrite=True,
    )

    assert first.generation_path.is_dir()
    assert second.generation_path.is_dir()
    assert first.generation_path != second.generation_path
    assert destinations["table"].read_text() == "new-table"
    assert destinations["manifest"].read_text() == "new-manifest"
    resolved = resolve_current_bundle(second.pointer_path)
    assert resolved.generation_path == second.generation_path
    assert resolved.artifacts["table"].read_text() == "new-table"
    assert resolved.artifacts["manifest"].read_text() == "new-manifest"
    assert first.artifacts["table"].read_text() == "old-table"


def test_pointer_manifest_detects_corrupt_generation_artifact(tmp_path):
    destinations = {"value": tmp_path / "out" / "value.json"}
    published = DurableArtifactBundlePublisher().publish(
        _staged(tmp_path, {"value": "original"}), destinations, overwrite=False
    )
    published.artifacts["value"].write_text("corrupt", encoding="utf-8")

    with pytest.raises(ValueError, match="digest"):
        resolve_current_bundle(published.pointer_path)


def test_pointer_commit_failure_restores_prior_legacy_bundle(tmp_path):
    destinations = {
        "first": tmp_path / "out" / "first.json",
        "second": tmp_path / "out" / "second.json",
    }
    publisher = DurableArtifactBundlePublisher()
    prior = publisher.publish(
        _staged(tmp_path, {"first": "old-1", "second": "old-2"}),
        destinations,
        overwrite=False,
    )

    def fail_pointer(source, destination):
        if Path(destination).name == "current.json":
            raise OSError("pointer commit failed")
        os.replace(source, destination)

    with pytest.raises(OSError, match="pointer commit failed"):
        DurableArtifactBundlePublisher(replace=fail_pointer).publish(
            _staged(tmp_path, {"first": "new-1", "second": "new-2"}),
            destinations,
            overwrite=True,
        )

    assert destinations["first"].read_text() == "old-1"
    assert destinations["second"].read_text() == "old-2"
    assert resolve_current_bundle(prior.pointer_path).generation_path == prior.generation_path


def test_post_pointer_fsync_failure_restores_prior_pointer_and_legacy_views(
    tmp_path, monkeypatch
):
    destinations = {
        "first": tmp_path / "out" / "first.json",
        "second": tmp_path / "out" / "second.json",
    }
    publisher = DurableArtifactBundlePublisher()
    prior = publisher.publish(
        _staged(tmp_path, {"first": "old-1", "second": "old-2"}),
        destinations,
        overwrite=False,
    )
    real_fsync = artifact_bundle._fsync_directory
    pointer_fsyncs = 0

    def fail_new_pointer_fsync(path):
        nonlocal pointer_fsyncs
        if Path(path) == prior.pointer_path.parent:
            pointer_fsyncs += 1
            if pointer_fsyncs == 1:
                raise OSError("pointer fsync failed")
        return real_fsync(path)

    monkeypatch.setattr(artifact_bundle, "_fsync_directory", fail_new_pointer_fsync)

    with pytest.raises(OSError, match="pointer fsync failed"):
        publisher.publish(
            _staged(tmp_path, {"first": "new-1", "second": "new-2"}),
            destinations,
            overwrite=True,
        )

    assert destinations["first"].read_text() == "old-1"
    assert destinations["second"].read_text() == "old-2"
    assert resolve_current_bundle(prior.pointer_path).generation_path == prior.generation_path


def test_failed_restoration_preserves_recovery_paths(tmp_path):
    destinations = {
        "first": tmp_path / "out" / "first.json",
        "second": tmp_path / "out" / "second.json",
    }
    publisher = DurableArtifactBundlePublisher()
    publisher.publish(
        _staged(tmp_path, {"first": "old-1", "second": "old-2"}),
        destinations,
        overwrite=False,
    )

    def fail_publish_and_restore(source, destination):
        source_path = Path(source)
        destination_path = Path(destination)
        if destination_path == destinations["second"] and ".compat." in source_path.name:
            raise OSError("publication failed")
        if destination_path == destinations["first"] and "recovery" in source_path.parts:
            raise OSError("restoration failed")
        os.replace(source, destination)

    with pytest.raises(ArtifactRecoveryError) as exc_info:
        DurableArtifactBundlePublisher(replace=fail_publish_and_restore).publish(
            _staged(tmp_path, {"first": "new-1", "second": "new-2"}),
            destinations,
            overwrite=True,
        )

    error = exc_info.value
    assert "publication failed" in str(error.publication_error)
    assert any("restoration failed" in str(item) for item in error.recovery_errors)
    assert error.recovery_paths
    assert all(path.exists() for path in error.recovery_paths)
    assert any(path.read_text() == "old-1" for path in error.recovery_paths)


def test_pointer_is_json_and_references_relative_generation(tmp_path):
    published = DurableArtifactBundlePublisher().publish(
        _staged(tmp_path, {"value": "ok"}),
        {"value": tmp_path / "out" / "value.json"},
        overwrite=False,
    )

    pointer = json.loads(published.pointer_path.read_text(encoding="utf-8"))
    assert pointer["schema_version"] == "1.0"
    assert not Path(pointer["manifest"]).is_absolute()
