"""Crash-durable publication for related output artifacts.

The immutable generation and its atomically replaced ``current.json`` pointer
are authoritative. Direct destination files are maintained as compatibility
views and are restored on a failed publication attempt.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import shutil
from typing import Any
from uuid import uuid4


class ArtifactRecoveryError(RuntimeError):
    """Publication failed and one or more compatibility views need recovery."""

    def __init__(
        self,
        publication_error: BaseException,
        recovery_errors: list[BaseException],
        recovery_paths: tuple[Path, ...],
    ) -> None:
        self.publication_error = publication_error
        self.recovery_errors = tuple(recovery_errors)
        self.recovery_paths = recovery_paths
        locations = ", ".join(str(path) for path in recovery_paths) or "none"
        super().__init__(
            "artifact publication and recovery failed; preserved recovery paths: "
            f"{locations}; recovery failed: {recovery_errors[0]}"
        )


@dataclass(frozen=True)
class PublishedArtifactBundle:
    """Resolved immutable generation selected by one current pointer."""

    pointer_path: Path
    generation_path: Path
    artifacts: dict[str, Path]


def _digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def _fsync_directory(path: Path) -> None:
    flags = os.O_RDONLY | getattr(os, "O_DIRECTORY", 0)
    descriptor = os.open(path, flags)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _write_json(path: Path, value: Mapping[str, Any]) -> None:
    with path.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, ensure_ascii=False, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())


class DurableArtifactBundlePublisher:
    """Publish immutable generations and one crash-atomic current pointer."""

    def __init__(
        self,
        *,
        replace: Callable[[os.PathLike[str] | str, os.PathLike[str] | str], None] = os.replace,
    ) -> None:
        self._replace = replace

    def publish(
        self,
        staged: Mapping[str, Path],
        destinations: Mapping[str, Path],
        *,
        overwrite: bool,
    ) -> PublishedArtifactBundle:
        """Publish ``staged`` files and return the committed immutable generation."""

        roles = tuple(destinations)
        if not roles or set(roles) != set(staged):
            raise ValueError("staged artifacts and destinations must have identical roles")
        normalized = {role: Path(destinations[role]) for role in roles}
        roots = {path.parent.resolve() for path in normalized.values()}
        if len(roots) != 1:
            raise ValueError("bundle destinations must share one parent directory")
        if len(set(normalized.values())) != len(normalized):
            raise ValueError("bundle destinations must be unique")
        if len({path.name for path in normalized.values()}) != len(normalized):
            raise ValueError("bundle destination basenames must be unique")
        missing = [role for role, path in staged.items() if not Path(path).is_file()]
        if missing:
            raise FileNotFoundError(f"staged artifact is missing: {missing[0]}")
        existing = [path for path in normalized.values() if path.exists()]
        if existing and not overwrite:
            raise FileExistsError(f"Output already exists: {existing[0]}")

        root = next(iter(roots))
        root.mkdir(parents=True, exist_ok=True)
        bundle_key = hashlib.sha256(
            "\0".join(f"{role}:{normalized[role].name}" for role in sorted(roles)).encode()
        ).hexdigest()[:16]
        bundle_root = root / ".artifact-bundles" / bundle_key
        generations = bundle_root / "generations"
        recovery_root = bundle_root / "recovery"
        generation_id = uuid4().hex
        generation = generations / generation_id
        generation.mkdir(parents=True)

        immutable: dict[str, Path] = {}
        manifest_artifacts: dict[str, dict[str, Any]] = {}
        for role in roles:
            destination = normalized[role]
            target = generation / destination.name
            shutil.copyfile(Path(staged[role]), target)
            with target.open("rb") as stream:
                os.fsync(stream.fileno())
            immutable[role] = target
            manifest_artifacts[role] = {
                "path": target.name,
                "legacy_destination": destination.name,
                "sha256": _digest(target),
                "size_bytes": target.stat().st_size,
            }

        manifest_path = generation / ".generation-manifest.json"
        _write_json(
            manifest_path,
            {
                "schema_version": "1.0",
                "generation_id": generation_id,
                "artifacts": manifest_artifacts,
            },
        )
        _fsync_directory(generation)
        _fsync_directory(generations)

        recovery = recovery_root / generation_id
        recovery.mkdir(parents=True)
        compatibility_stages: dict[str, Path] = {}
        backups: dict[str, Path] = {}
        published: list[str] = []
        pointer = bundle_root / "current.json"
        pointer_stage = bundle_root / f".{generation_id}.current.stage"
        recovery_failed = False
        try:
            for role in roles:
                destination = normalized[role]
                compatibility = root / f".{destination.name}.{generation_id}.compat.stage"
                shutil.copyfile(immutable[role], compatibility)
                with compatibility.open("rb") as stream:
                    os.fsync(stream.fileno())
                compatibility_stages[role] = compatibility

            for role in roles:
                destination = normalized[role]
                if destination.exists():
                    backup = recovery / destination.name
                    self._replace(destination, backup)
                    backups[role] = backup
                self._replace(compatibility_stages[role], destination)
                published.append(role)
            _fsync_directory(root)

            manifest_relative = manifest_path.relative_to(bundle_root).as_posix()
            _write_json(
                pointer_stage,
                {
                    "schema_version": "1.0",
                    "generation_id": generation_id,
                    "manifest": manifest_relative,
                    "manifest_sha256": _digest(manifest_path),
                },
            )
            self._replace(pointer_stage, pointer)
            _fsync_directory(bundle_root)
        except BaseException as publication_error:
            recovery_errors: list[BaseException] = []
            for role in reversed(published):
                try:
                    normalized[role].unlink(missing_ok=True)
                except BaseException as error:
                    recovery_errors.append(error)
            for role in reversed(tuple(backups)):
                backup = backups[role]
                if backup.exists():
                    try:
                        self._replace(backup, normalized[role])
                    except BaseException as error:
                        recovery_errors.append(error)
            try:
                _fsync_directory(root)
            except BaseException as error:
                recovery_errors.append(error)
            if recovery_errors:
                recovery_failed = True
                recovery_paths = tuple(
                    [path for path in backups.values() if path.exists()]
                    + [generation, manifest_path]
                )
                raise ArtifactRecoveryError(
                    publication_error, recovery_errors, recovery_paths
                ) from publication_error
            raise
        finally:
            pointer_stage.unlink(missing_ok=True)
            for path in compatibility_stages.values():
                path.unlink(missing_ok=True)
            if not recovery_failed:
                shutil.rmtree(recovery, ignore_errors=True)

        return PublishedArtifactBundle(pointer, generation, immutable)


def resolve_current_bundle(pointer_path: str | Path) -> PublishedArtifactBundle:
    """Validate and resolve the immutable generation selected by ``pointer_path``."""

    pointer = Path(pointer_path)
    payload = json.loads(pointer.read_text(encoding="utf-8"))
    if payload.get("schema_version") != "1.0":
        raise ValueError("unsupported artifact pointer schema version")
    manifest_relative = Path(str(payload.get("manifest") or ""))
    if manifest_relative.is_absolute() or ".." in manifest_relative.parts:
        raise ValueError("artifact pointer manifest must remain inside the bundle")
    manifest_path = pointer.parent / manifest_relative
    if not manifest_path.is_file():
        raise FileNotFoundError(f"artifact generation manifest is missing: {manifest_path}")
    if _digest(manifest_path) != payload.get("manifest_sha256"):
        raise ValueError("artifact generation manifest digest does not match pointer")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema_version") != "1.0":
        raise ValueError("unsupported artifact generation schema version")
    generation = manifest_path.parent
    artifacts: dict[str, Path] = {}
    for role, record in (manifest.get("artifacts") or {}).items():
        relative = Path(str(record.get("path") or ""))
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError("artifact path must remain inside its generation")
        artifact = generation / relative
        if not artifact.is_file():
            raise FileNotFoundError(f"artifact generation file is missing: {artifact}")
        if _digest(artifact) != record.get("sha256"):
            raise ValueError(f"artifact digest does not match manifest: {role}")
        artifacts[str(role)] = artifact
    return PublishedArtifactBundle(pointer, generation, artifacts)


__all__ = [
    "ArtifactRecoveryError",
    "DurableArtifactBundlePublisher",
    "PublishedArtifactBundle",
    "resolve_current_bundle",
]
