# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Public requests and outcomes for the unified conversion boundary."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping


@dataclass(frozen=True)
class InputSpec:
    """One logical input. Sources may include explicitly bound companions."""

    sources: Any
    id: str | None = None
    in_type: str | None = None
    companions: Mapping[str, Any] = field(default_factory=dict)
    metadata: Any = None
    input_options: Mapping[str, Any] = field(default_factory=dict)
    output_options: Mapping[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class Diagnostic:
    code: str
    message: str
    stage: str = "input"
    severity: str = "error"


class InputError(ValueError):
    """A safe, actionable input or capability error."""

    def __init__(self, code, message):
        self.code = code
        super().__init__(message)


@dataclass
class LoadedInput:
    metadata: Any = None
    expression: Any = None
    assets: tuple = ()
    direct_magetab: Any = None
    source_path: str | None = None
    content_id: str | None = None
    provider: str | None = None
    enrichment_applied: bool = False
    origins: tuple[str, ...] = ()
    companions: Mapping[str, Any] = field(default_factory=dict)
    diagnostics: list[Diagnostic] = field(default_factory=list)


@dataclass
class ConversionItemResult:
    id: str
    source: str
    in_type: str | None = None
    forced: bool = False
    provider: str | None = None
    content_id: str | None = None
    origins: tuple[str, ...] = ()
    companions: dict[str, Any] = field(default_factory=dict)
    route: tuple[str, ...] = ()
    dataset_ids: tuple[str, ...] = ()
    status: str = "failed"
    execution: str = "failed"
    completeness: str = "unknown"
    validation: str = "unknown"
    payload: Any = None
    artifacts: dict[str, str] = field(default_factory=dict)
    diagnostics: list[Diagnostic] = field(default_factory=list)

    def to_dict(self):
        """Serialize a summary without embedding large or scientific payloads."""
        return {
            name: getattr(self, name)
            for name in (
                "id",
                "source",
                "in_type",
                "forced",
                "provider",
                "content_id",
                "origins",
                "companions",
                "route",
                "dataset_ids",
                "status",
                "execution",
                "completeness",
                "validation",
                "artifacts",
            )
        } | {"diagnostics": [vars(d).copy() for d in self.diagnostics]}


@dataclass
class ConversionBatchResult:
    items: list[ConversionItemResult] = field(default_factory=list)

    @property
    def status(self):
        if not self.items or all(i.status in {"failed", "skipped"} for i in self.items):
            return "failed"
        return (
            "complete" if all(i.status == "complete" for i in self.items) else "partial"
        )

    @property
    def partial(self):
        return self.status == "partial"

    def to_dict(self):
        return {
            "operation": "convert",
            "status": self.status,
            "items": [i.to_dict() for i in self.items],
        }
