# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Immutable provider-response records used by the study converters."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib


@dataclass(frozen=True)
class ProviderDocument:
    """One byte-faithful response consumed while constructing a study."""

    kind: str
    name: str
    uri: str | None
    media_type: str
    content: bytes

    def __post_init__(self) -> None:
        if not self.kind.strip() or not self.name.strip():
            raise ValueError("provider document requires nonblank kind and name")
        if not self.media_type.strip():
            raise ValueError("provider document requires a nonblank media type")
        if not isinstance(self.content, bytes):
            raise TypeError("provider document content must be bytes")

    @property
    def sha256(self) -> str:
        return hashlib.sha256(self.content).hexdigest()


@dataclass(frozen=True)
class StudyFetchResult:
    """All provider documents that belong to one resolved INSDC study."""

    provider: str
    requested_accession: str
    study_accession: str
    documents: tuple[ProviderDocument, ...]
    warnings: tuple[dict, ...] = ()

    def __post_init__(self) -> None:
        provider = self.provider.casefold()
        if provider not in {"sra", "ena"}:
            raise ValueError("provider must be 'sra' or 'ena'")
        object.__setattr__(self, "provider", provider)
        if not self.requested_accession.strip() or not self.study_accession.strip():
            raise ValueError("study fetch result requires requested and study accessions")
        if not isinstance(self.documents, tuple) or not all(
            isinstance(item, ProviderDocument) for item in self.documents
        ):
            raise TypeError("study fetch result documents must be ProviderDocument values")
