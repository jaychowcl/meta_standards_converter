# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================

from __future__ import annotations

import json

import pytest

from meta_standards_converter.runtime_contracts import (
    CompletenessStatus,
    DiskBudgetError,
    EvidenceConfidence,
    ExecutionStatus,
    OperationStatusV2,
    PublicationDisposition,
    RetryCategory,
    SafeErrorEnvelope,
    ValidationStatus,
    get_resource_profile,
    require_disk_headroom,
)
from meta_standards_converter.helpers.request_helper import RequestSettings


MIB = 1024**2
GIB = 1024**3
TIB = 1024**4


def test_standard_resource_profile_has_exact_conservative_limits() -> None:
    profile = get_resource_profile("standard")

    assert profile.name == "standard"
    assert profile.max_redirects == 3
    assert profile.connect_timeout_seconds == 10
    assert profile.read_timeout_seconds == 60
    assert profile.max_xml_bytes == 128 * MIB
    assert profile.max_compressed_archive_bytes == 256 * MIB
    assert profile.max_expanded_archive_bytes == 1 * GIB
    assert profile.max_ontology_file_bytes == 3 * GIB
    assert profile.max_matrix_bytes == 100 * GIB
    assert profile.max_aggregate_download_bytes == 200 * GIB
    assert profile.max_cache_bytes == 250 * GIB
    assert profile.network_workers == 4
    assert profile.ontology_build_workers == 1
    assert profile.disk_headroom_fraction == 0.10


def test_large_resource_profile_has_exact_opt_in_limits() -> None:
    profile = get_resource_profile("large")

    assert profile.name == "large"
    assert profile.max_redirects == 5
    assert profile.connect_timeout_seconds == 30
    assert profile.read_timeout_seconds == 300
    assert profile.max_xml_bytes == 512 * MIB
    assert profile.max_compressed_archive_bytes == 1 * GIB
    assert profile.max_expanded_archive_bytes == 4 * GIB
    assert profile.max_ontology_file_bytes == 5 * GIB
    assert profile.max_matrix_bytes == 500 * GIB
    assert profile.max_aggregate_download_bytes == 1 * TIB
    assert profile.max_cache_bytes == 500 * GIB
    assert profile.network_workers == 8
    assert profile.ontology_build_workers == 2


def test_request_settings_are_derived_from_typed_resource_profiles() -> None:
    standard = RequestSettings.from_resource_profile(
        get_resource_profile("standard"),
        request_delay=0.5,
    )
    large = RequestSettings.from_resource_profile(get_resource_profile("large"))

    assert standard.timeout == (10, 60)
    assert standard.max_in_flight == 4
    assert standard.request_delay == 0.5
    assert large.timeout == (30, 300)
    assert large.max_in_flight == 8


def test_resource_profile_api_overrides_are_explicit_and_validated() -> None:
    profile = get_resource_profile(
        "standard",
        overrides={"max_xml_bytes": 64 * MIB, "network_workers": 2},
    )

    assert profile.max_xml_bytes == 64 * MIB
    assert profile.network_workers == 2
    assert get_resource_profile("standard").max_xml_bytes == 128 * MIB

    with pytest.raises(ValueError, match="Unknown resource override"):
        get_resource_profile("standard", overrides={"mystery_limit": 1})
    with pytest.raises(ValueError, match="must be positive"):
        get_resource_profile("standard", overrides={"max_xml_bytes": 0})


def test_disk_preflight_requires_ten_percent_headroom(monkeypatch, tmp_path) -> None:
    class _Usage:
        total = 1_000
        used = 450
        free = 550

    monkeypatch.setattr("shutil.disk_usage", lambda path: _Usage())
    result = require_disk_headroom(tmp_path, required_bytes=500)

    assert result.required_with_headroom_bytes == 550
    assert result.free_bytes == 550

    _Usage.free = 549
    with pytest.raises(DiskBudgetError, match="10% headroom"):
        require_disk_headroom(tmp_path, required_bytes=500)


def test_status_v2_keeps_execution_and_scientific_evidence_independent() -> None:
    status = OperationStatusV2(
        execution=ExecutionStatus.SUCCEEDED,
        completeness=CompletenessStatus.COMPLETE,
        evidence_confidence=EvidenceConfidence.INSUFFICIENT,
        validation=ValidationStatus.VALID,
        publication=PublicationDisposition.REVIEW_REQUIRED,
        terminal_reason="healthy_no_match",
    )

    assert status.to_dict() == {
        "contract_version": "2.0",
        "execution": "succeeded",
        "completeness": "complete",
        "evidence_confidence": "insufficient",
        "validation": "valid",
        "publication": "review_required",
        "terminal_reason": "healthy_no_match",
        "errors": [],
    }
    assert status.safe_to_publish is False


def test_status_v2_rejects_legacy_or_unknown_contracts() -> None:
    with pytest.raises(ValueError, match="contract_version 2.0"):
        OperationStatusV2.from_dict({"status": "complete"})
    with pytest.raises(ValueError, match="contract_version 2.0"):
        OperationStatusV2.from_dict(
            {
                "contract_version": "1.0",
                "execution": "succeeded",
                "completeness": "complete",
                "evidence_confidence": "high",
                "validation": "valid",
                "publication": "publishable",
                "terminal_reason": "complete",
                "errors": [],
            }
        )


def test_safe_error_envelope_never_serializes_credentials_or_query_values() -> None:
    canary = "canary-super-secret"
    error = RuntimeError(
        "request failed for https://user:pass@example.org/data?token="
        f"{canary}#fragment"
    )
    envelope = SafeErrorEnvelope.from_exception(
        error,
        provider="biostudies",
        location=(
            "https://user:pass@example.org/data?token=" f"{canary}#fragment"
        ),
        retry_category=RetryCategory.RETRYABLE,
        stage="download",
        item_id="E-MTAB-1",
        correlation_id="corr-1",
    )

    payload = json.dumps(envelope.to_dict(), sort_keys=True)
    assert canary not in payload
    assert "user" not in payload
    assert "pass" not in payload
    assert "token" not in payload
    assert envelope.location == "https://example.org/data"
    assert envelope.error_type == "RuntimeError"
    assert envelope.correlation_id == "corr-1"


def test_failed_status_cannot_claim_publishable_disposition() -> None:
    with pytest.raises(ValueError, match="failed execution cannot be publishable"):
        OperationStatusV2(
            execution=ExecutionStatus.FAILED,
            completeness=CompletenessStatus.PARTIAL,
            evidence_confidence=EvidenceConfidence.NOT_ASSESSED,
            validation=ValidationStatus.INVALID,
            publication=PublicationDisposition.PUBLISHABLE,
            terminal_reason="provider_failed",
        )
