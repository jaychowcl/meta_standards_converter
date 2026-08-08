# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Canonical decoding, encoding, and explicit migration for MSC MINiML 2.0."""

from __future__ import annotations

from dataclasses import dataclass
import json
import os
from pathlib import Path
from typing import Any, Mapping, Sequence
from uuid import uuid4

from .model import MINiMLPackage, MINiMLValidationIssue


class MINiMLCompatibilityError(ValueError):
    """Compatibility diagnostics were promoted to a decoding failure."""

    def __init__(self, diagnostics: Sequence[MINiMLValidationIssue]):
        self.diagnostics = tuple(diagnostics)
        rendered = "; ".join(f"{item.code} at {item.path}" for item in diagnostics)
        super().__init__(rendered)


@dataclass(frozen=True)
class MINiMLDecodeResult:
    package: MINiMLPackage
    diagnostics: tuple[MINiMLValidationIssue, ...] = ()


@dataclass(frozen=True)
class MINiMLBatchDecodeResult:
    packages: tuple[MINiMLPackage, ...]
    diagnostics: tuple[MINiMLValidationIssue, ...] = ()


class MINiMLCodec:
    @staticmethod
    def migrate_v1(value: Mapping[str, Any]):
        """Migrate a legacy package without weakening the strict v2 decoder."""
        from .migration import MINiMLV1Migrator

        return MINiMLV1Migrator().migrate(value)

    def decode(self, value: Mapping[str, Any], *, strict: bool = False) -> MINiMLDecodeResult:
        package = MINiMLPackage.from_mapping(value.to_mapping()) if isinstance(value, MINiMLPackage) else MINiMLPackage.from_mapping(value)
        diagnostics = package.validate()
        if strict and diagnostics:
            raise MINiMLCompatibilityError(diagnostics)
        return MINiMLDecodeResult(package, diagnostics)

    def decode_many(self, value: Mapping[str, Any] | Sequence[Mapping[str, Any]], *, strict: bool = False) -> MINiMLBatchDecodeResult:
        values = list(value) if isinstance(value, (list, tuple)) else [value]
        packages = []
        diagnostics = []
        for item in values:
            result = self.decode(item, strict=strict)
            packages.append(result.package)
            diagnostics.extend(result.diagnostics)
        return MINiMLBatchDecodeResult(tuple(packages), tuple(diagnostics))

    @staticmethod
    def encode(package: MINiMLPackage) -> dict[str, Any]:
        if not isinstance(package, MINiMLPackage):
            raise TypeError("package must be a MINiMLPackage")
        canonical = MINiMLPackage.from_mapping(package.to_mapping())
        canonical.validate()
        return canonical.to_mapping()

    def encode_many(self, packages: Sequence[MINiMLPackage]) -> list[dict[str, Any]]:
        return [self.encode(package) for package in packages]

    def load(self, path: str | Path, *, strict: bool = False) -> MINiMLBatchDecodeResult:
        payload = json.loads(Path(path).read_text(encoding="utf-8"))
        return self.decode_many(payload, strict=strict)

    def dump(self, packages: MINiMLPackage | Sequence[MINiMLPackage], path: str | Path) -> None:
        values = [packages] if isinstance(packages, MINiMLPackage) else list(packages)
        destination = Path(path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        temporary = destination.with_name(f".{destination.name}.{uuid4().hex}.tmp")
        try:
            with temporary.open("x", encoding="utf-8") as handle:
                json.dump(
                    self.encode_many(values),
                    handle,
                    indent=2,
                    sort_keys=True,
                    ensure_ascii=False,
                )
                handle.write("\n")
                handle.flush()
                os.fsync(handle.fileno())
            os.replace(temporary, destination)
        finally:
            temporary.unlink(missing_ok=True)
