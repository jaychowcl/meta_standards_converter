# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================

"""Shared, versioned runtime status, resource, and safe-error contracts."""

from __future__ import annotations

from dataclasses import dataclass, fields, replace
from enum import Enum
from pathlib import Path
import re
import shutil
from typing import Any, Iterable, Mapping
from urllib.parse import urlsplit, urlunsplit
from uuid import uuid4


MIB = 1024**2
GIB = 1024**3
TIB = 1024**4
STATUS_CONTRACT_VERSION = "2.0"
SAFE_ERROR_CONTRACT_VERSION = "1.0"
_REASON_PATTERN = re.compile(r"^[a-z][a-z0-9_.-]{0,127}$")


class ExecutionStatus(str, Enum):
    NOT_RUN = "not_run"
    SUCCEEDED = "succeeded"
    DEGRADED = "degraded"
    FAILED = "failed"


class CompletenessStatus(str, Enum):
    UNKNOWN = "unknown"
    EMPTY = "empty"
    PARTIAL = "partial"
    COMPLETE = "complete"


class EvidenceConfidence(str, Enum):
    NOT_ASSESSED = "not_assessed"
    INSUFFICIENT = "insufficient"
    LOW = "low"
    MEDIUM = "medium"
    HIGH = "high"


class ValidationStatus(str, Enum):
    NOT_RUN = "not_run"
    VALID = "valid"
    INVALID = "invalid"


class PublicationDisposition(str, Enum):
    NOT_REQUESTED = "not_requested"
    PUBLISHABLE = "publishable"
    REVIEW_REQUIRED = "review_required"
    BLOCKED = "blocked"


class RetryCategory(str, Enum):
    UNKNOWN = "unknown"
    RETRYABLE = "retryable"
    TERMINAL = "terminal"


@dataclass(frozen=True)
class SafeErrorEnvelope:
    """A persistence-safe description that intentionally excludes raw messages."""

    error_type: str
    provider: str | None = None
    location: str | None = None
    http_status: int | None = None
    retry_category: RetryCategory = RetryCategory.UNKNOWN
    stage: str | None = None
    item_id: str | None = None
    correlation_id: str = ""
    contract_version: str = SAFE_ERROR_CONTRACT_VERSION

    def __post_init__(self) -> None:
        if self.contract_version != SAFE_ERROR_CONTRACT_VERSION:
            raise ValueError(
                f"Safe error contract_version must be {SAFE_ERROR_CONTRACT_VERSION}."
            )
        if not self.error_type.strip():
            raise ValueError("error_type must not be blank")
        if self.http_status is not None and not 100 <= self.http_status <= 599:
            raise ValueError("http_status must be a valid HTTP status code")
        if not self.correlation_id:
            object.__setattr__(self, "correlation_id", uuid4().hex)

    @classmethod
    def from_exception(
        cls,
        error: BaseException,
        *,
        provider: str | None = None,
        location: str | Path | None = None,
        http_status: int | None = None,
        retry_category: RetryCategory = RetryCategory.UNKNOWN,
        stage: str | None = None,
        item_id: str | None = None,
        correlation_id: str | None = None,
    ) -> "SafeErrorEnvelope":
        return cls(
            error_type=type(error).__name__,
            provider=_safe_label(provider),
            location=_sanitize_location(location),
            http_status=http_status,
            retry_category=retry_category,
            stage=_safe_label(stage),
            item_id=_safe_label(item_id),
            correlation_id=correlation_id or uuid4().hex,
        )

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "SafeErrorEnvelope":
        if value.get("contract_version") != SAFE_ERROR_CONTRACT_VERSION:
            raise ValueError(
                f"Safe error payload requires contract_version "
                f"{SAFE_ERROR_CONTRACT_VERSION}."
            )
        return cls(
            error_type=str(value.get("error_type", "")),
            provider=_optional_string(value.get("provider")),
            location=_optional_string(value.get("location")),
            http_status=value.get("http_status"),
            retry_category=RetryCategory(value.get("retry_category", "unknown")),
            stage=_optional_string(value.get("stage")),
            item_id=_optional_string(value.get("item_id")),
            correlation_id=str(value.get("correlation_id", "")),
            contract_version=str(value["contract_version"]),
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "contract_version": self.contract_version,
            "error_type": self.error_type,
            "provider": self.provider,
            "location": self.location,
            "http_status": self.http_status,
            "retry_category": self.retry_category.value,
            "stage": self.stage,
            "item_id": self.item_id,
            "correlation_id": self.correlation_id,
        }


@dataclass(frozen=True)
class OperationStatusV2:
    """Independent execution, completeness, evidence, validation, and release axes."""

    execution: ExecutionStatus
    completeness: CompletenessStatus
    evidence_confidence: EvidenceConfidence
    validation: ValidationStatus
    publication: PublicationDisposition
    terminal_reason: str
    errors: tuple[SafeErrorEnvelope, ...] = ()
    contract_version: str = STATUS_CONTRACT_VERSION

    def __post_init__(self) -> None:
        if self.contract_version != STATUS_CONTRACT_VERSION:
            raise ValueError(
                f"Operation status requires contract_version "
                f"{STATUS_CONTRACT_VERSION}."
            )
        if not _REASON_PATTERN.fullmatch(self.terminal_reason):
            raise ValueError(
                "terminal_reason must be a nonblank stable lowercase code"
            )
        if (
            self.execution is ExecutionStatus.FAILED
            and self.publication is PublicationDisposition.PUBLISHABLE
        ):
            raise ValueError("failed execution cannot be publishable")
        if (
            self.validation is ValidationStatus.INVALID
            and self.publication is PublicationDisposition.PUBLISHABLE
        ):
            raise ValueError("invalid output cannot be publishable")
        if (
            self.completeness is not CompletenessStatus.COMPLETE
            and self.publication is PublicationDisposition.PUBLISHABLE
        ):
            raise ValueError("incomplete output cannot be publishable")

    @property
    def safe_to_publish(self) -> bool:
        return (
            self.execution is ExecutionStatus.SUCCEEDED
            and self.completeness is CompletenessStatus.COMPLETE
            and self.validation is ValidationStatus.VALID
            and self.publication is PublicationDisposition.PUBLISHABLE
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "contract_version": self.contract_version,
            "execution": self.execution.value,
            "completeness": self.completeness.value,
            "evidence_confidence": self.evidence_confidence.value,
            "validation": self.validation.value,
            "publication": self.publication.value,
            "terminal_reason": self.terminal_reason,
            "errors": [error.to_dict() for error in self.errors],
        }

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "OperationStatusV2":
        if value.get("contract_version") != STATUS_CONTRACT_VERSION:
            raise ValueError(
                f"Operation status requires contract_version "
                f"{STATUS_CONTRACT_VERSION}."
            )
        raw_errors = value.get("errors")
        if not isinstance(raw_errors, list):
            raise ValueError("Operation status errors must be a list.")
        return cls(
            execution=ExecutionStatus(value["execution"]),
            completeness=CompletenessStatus(value["completeness"]),
            evidence_confidence=EvidenceConfidence(value["evidence_confidence"]),
            validation=ValidationStatus(value["validation"]),
            publication=PublicationDisposition(value["publication"]),
            terminal_reason=str(value["terminal_reason"]),
            errors=tuple(SafeErrorEnvelope.from_dict(item) for item in raw_errors),
            contract_version=str(value["contract_version"]),
        )

    @classmethod
    def aggregate(
        cls,
        values: Iterable["OperationStatusV2"],
        *,
        terminal_reason: str = "aggregate_complete",
    ) -> "OperationStatusV2":
        items = tuple(values)
        if not items:
            return cls(
                execution=ExecutionStatus.NOT_RUN,
                completeness=CompletenessStatus.EMPTY,
                evidence_confidence=EvidenceConfidence.NOT_ASSESSED,
                validation=ValidationStatus.NOT_RUN,
                publication=PublicationDisposition.NOT_REQUESTED,
                terminal_reason="no_items",
            )
        execution = _aggregate_execution(items)
        completeness = _aggregate_completeness(items)
        validation = _aggregate_validation(items)
        evidence = min(
            (item.evidence_confidence for item in items),
            key=_evidence_rank,
        )
        publication = _aggregate_publication(
            items,
            execution=execution,
            completeness=completeness,
            validation=validation,
        )
        return cls(
            execution=execution,
            completeness=completeness,
            evidence_confidence=evidence,
            validation=validation,
            publication=publication,
            terminal_reason=terminal_reason,
            errors=tuple(error for item in items for error in item.errors),
        )


@dataclass(frozen=True)
class ResourceProfile:
    """Disk, network, worker, and in-memory admission ceilings."""

    name: str
    max_redirects: int
    connect_timeout_seconds: int
    read_timeout_seconds: int
    max_xml_bytes: int
    max_compressed_archive_bytes: int
    max_expanded_archive_bytes: int
    max_ontology_file_bytes: int
    max_matrix_bytes: int
    max_in_memory_matrix_bytes: int
    max_aggregate_download_bytes: int
    max_cache_bytes: int
    network_workers: int
    ontology_build_workers: int
    disk_headroom_fraction: float = 0.10
    available_memory_fraction: float = 0.70
    force_memory_fraction: float = 0.90

    def __post_init__(self) -> None:
        for field in fields(self):
            if field.name == "name":
                continue
            value = getattr(self, field.name)
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise ValueError(f"Resource limit {field.name} must be numeric.")
            if value <= 0:
                raise ValueError(f"Resource limit {field.name} must be positive.")
        if not 0 < self.disk_headroom_fraction < 1:
            raise ValueError("disk_headroom_fraction must be between zero and one")
        if not 0 < self.available_memory_fraction < 1:
            raise ValueError("available_memory_fraction must be between zero and one")
        if not 0 < self.force_memory_fraction < 1:
            raise ValueError("force_memory_fraction must be between zero and one")
        if self.force_memory_fraction < self.available_memory_fraction:
            raise ValueError(
                "force_memory_fraction cannot be lower than available_memory_fraction"
            )

    def with_overrides(self, overrides: Mapping[str, int | float]) -> "ResourceProfile":
        allowed = {field.name for field in fields(self)} - {"name"}
        unknown = set(overrides) - allowed
        if unknown:
            raise ValueError(
                f"Unknown resource override: {sorted(unknown)[0]}"
            )
        return replace(self, **dict(overrides))


STANDARD_RESOURCE_PROFILE = ResourceProfile(
    name="standard",
    max_redirects=3,
    connect_timeout_seconds=10,
    read_timeout_seconds=60,
    max_xml_bytes=128 * MIB,
    max_compressed_archive_bytes=256 * MIB,
    max_expanded_archive_bytes=1 * GIB,
    max_ontology_file_bytes=3 * GIB,
    max_matrix_bytes=100 * GIB,
    max_in_memory_matrix_bytes=8 * GIB,
    max_aggregate_download_bytes=200 * GIB,
    max_cache_bytes=250 * GIB,
    network_workers=4,
    ontology_build_workers=1,
)

LARGE_RESOURCE_PROFILE = ResourceProfile(
    name="large",
    max_redirects=5,
    connect_timeout_seconds=30,
    read_timeout_seconds=300,
    max_xml_bytes=512 * MIB,
    max_compressed_archive_bytes=1 * GIB,
    max_expanded_archive_bytes=4 * GIB,
    max_ontology_file_bytes=5 * GIB,
    max_matrix_bytes=500 * GIB,
    max_in_memory_matrix_bytes=32 * GIB,
    max_aggregate_download_bytes=1 * TIB,
    max_cache_bytes=500 * GIB,
    network_workers=8,
    ontology_build_workers=2,
)

RESOURCE_PROFILES = {
    "standard": STANDARD_RESOURCE_PROFILE,
    "large": LARGE_RESOURCE_PROFILE,
}


def get_resource_profile(
    name: str | ResourceProfile = "standard",
    *,
    overrides: Mapping[str, int | float] | None = None,
) -> ResourceProfile:
    if isinstance(name, ResourceProfile):
        return name.with_overrides(overrides) if overrides else name
    try:
        profile = RESOURCE_PROFILES[name]
    except KeyError as error:
        raise ValueError(
            f"Unknown resource profile {name!r}; expected standard or large."
        ) from error
    return profile.with_overrides(overrides or {})


@dataclass(frozen=True)
class DiskBudget:
    path: str
    required_bytes: int
    required_with_headroom_bytes: int
    free_bytes: int
    headroom_fraction: float


class DiskBudgetError(RuntimeError):
    """Raised before work starts when its disk envelope cannot be honored."""


def require_disk_headroom(
    path: str | Path,
    *,
    required_bytes: int,
    headroom_fraction: float = 0.10,
) -> DiskBudget:
    if isinstance(required_bytes, bool) or not isinstance(required_bytes, int):
        raise ValueError("required_bytes must be a positive integer")
    if required_bytes <= 0:
        raise ValueError("required_bytes must be a positive integer")
    if not 0 < headroom_fraction < 1:
        raise ValueError("headroom_fraction must be between zero and one")
    destination = Path(path)
    existing = destination
    while not existing.exists() and existing != existing.parent:
        existing = existing.parent
    free = shutil.disk_usage(existing).free
    required_with_headroom = int(required_bytes * (1 + headroom_fraction))
    budget = DiskBudget(
        path=str(destination),
        required_bytes=required_bytes,
        required_with_headroom_bytes=required_with_headroom,
        free_bytes=free,
        headroom_fraction=headroom_fraction,
    )
    if free < required_with_headroom:
        percentage = f"{headroom_fraction:.0%}"
        raise DiskBudgetError(
            f"Disk preflight failed for {destination}: {free} bytes free, "
            f"{required_with_headroom} required including {percentage} headroom."
        )
    return budget


def _sanitize_location(value: str | Path | None) -> str | None:
    if value is None:
        return None
    text = str(value)
    parsed = urlsplit(text)
    if parsed.scheme and parsed.hostname:
        hostname = parsed.hostname
        if ":" in hostname and not hostname.startswith("["):
            hostname = f"[{hostname}]"
        netloc = hostname
        if parsed.port is not None:
            netloc = f"{hostname}:{parsed.port}"
        return urlunsplit((parsed.scheme.lower(), netloc, parsed.path, "", ""))
    name = Path(text).name
    return name or None


def _safe_label(value: str | None) -> str | None:
    if value is None:
        return None
    normalized = re.sub(r"[^A-Za-z0-9_.:-]", "_", str(value).strip())
    return normalized[:128] or None


def _optional_string(value: Any) -> str | None:
    return None if value is None else str(value)


def _aggregate_execution(items: tuple[OperationStatusV2, ...]) -> ExecutionStatus:
    values = {item.execution for item in items}
    if ExecutionStatus.FAILED in values:
        return ExecutionStatus.FAILED
    if ExecutionStatus.DEGRADED in values or ExecutionStatus.NOT_RUN in values:
        return ExecutionStatus.DEGRADED
    return ExecutionStatus.SUCCEEDED


def _aggregate_completeness(
    items: tuple[OperationStatusV2, ...],
) -> CompletenessStatus:
    values = {item.completeness for item in items}
    if values == {CompletenessStatus.EMPTY}:
        return CompletenessStatus.EMPTY
    if values == {CompletenessStatus.COMPLETE}:
        return CompletenessStatus.COMPLETE
    if CompletenessStatus.UNKNOWN in values:
        return CompletenessStatus.UNKNOWN
    return CompletenessStatus.PARTIAL


def _aggregate_validation(items: tuple[OperationStatusV2, ...]) -> ValidationStatus:
    values = {item.validation for item in items}
    if ValidationStatus.INVALID in values:
        return ValidationStatus.INVALID
    if values == {ValidationStatus.VALID}:
        return ValidationStatus.VALID
    return ValidationStatus.NOT_RUN


def _aggregate_publication(
    items: tuple[OperationStatusV2, ...],
    *,
    execution: ExecutionStatus,
    completeness: CompletenessStatus,
    validation: ValidationStatus,
) -> PublicationDisposition:
    values = {item.publication for item in items}
    if (
        PublicationDisposition.BLOCKED in values
        or execution is ExecutionStatus.FAILED
        or validation is ValidationStatus.INVALID
    ):
        return PublicationDisposition.BLOCKED
    if (
        PublicationDisposition.REVIEW_REQUIRED in values
        or execution is not ExecutionStatus.SUCCEEDED
        or completeness is not CompletenessStatus.COMPLETE
        or validation is not ValidationStatus.VALID
    ):
        return PublicationDisposition.REVIEW_REQUIRED
    if values == {PublicationDisposition.PUBLISHABLE}:
        return PublicationDisposition.PUBLISHABLE
    return PublicationDisposition.NOT_REQUESTED


def _evidence_rank(value: EvidenceConfidence) -> int:
    return {
        EvidenceConfidence.NOT_ASSESSED: 0,
        EvidenceConfidence.INSUFFICIENT: 1,
        EvidenceConfidence.LOW: 2,
        EvidenceConfidence.MEDIUM: 3,
        EvidenceConfidence.HIGH: 4,
    }[value]


__all__ = [
    "CompletenessStatus",
    "DiskBudget",
    "DiskBudgetError",
    "EvidenceConfidence",
    "ExecutionStatus",
    "LARGE_RESOURCE_PROFILE",
    "OperationStatusV2",
    "PublicationDisposition",
    "RESOURCE_PROFILES",
    "ResourceProfile",
    "RetryCategory",
    "SAFE_ERROR_CONTRACT_VERSION",
    "STATUS_CONTRACT_VERSION",
    "STANDARD_RESOURCE_PROFILE",
    "SafeErrorEnvelope",
    "ValidationStatus",
    "get_resource_profile",
    "require_disk_headroom",
]
