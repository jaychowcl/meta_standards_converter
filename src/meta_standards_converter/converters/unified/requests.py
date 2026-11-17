# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Immutable configuration at the conversion boundary; one owner per setting."""
from collections.abc import Mapping
from dataclasses import dataclass, fields
from types import MappingProxyType

from meta_standards_converter.runtime_contracts import ResourceProfile, get_resource_profile
from .options import INPUT, OUTPUT, RUNTIME_DEFAULTS, settings


@dataclass(frozen=True)
class PreparationPolicy:
    """Optional enrichment and independent study expansion.

    ``off`` disables enrichment, not loading or expansion. Disable expansion too
    to convert only supplied metadata without supplementary provider retrieval.
    """
    enrichment: str = "standard"
    expand_studies: bool = True
    linked_metadata: bool = True
    peer_provider: bool = True

    def __post_init__(self):
        if self.enrichment not in {"standard", "curators", "off"}:
            raise ValueError("enrichment must be standard, curators or off")
        for name in ("expand_studies", "linked_metadata", "peer_provider"):
            if not isinstance(getattr(self, name), bool):
                raise TypeError(f"{name} must be a boolean")

    def allows(self, operation):
        if self.enrichment == "off":
            return False
        if operation in {"linked_metadata", "peer_provider"}:
            return self.enrichment == "standard" and getattr(self, operation)
        return True


@dataclass(frozen=True)
class RuntimeOptions:
    """Execution controls, shared by every source and destination in a call."""
    overwrite: bool = False
    fail_fast: bool = False
    recursive: bool = False
    allow_processing: bool = False
    resource_profile: str | ResourceProfile = "standard"
    resource_overrides: Mapping[str, int | float] | None = None
    insdc_default: str = "ena"
    allowed_hosts: tuple[str, ...] = ()
    reserved_paths: tuple[str, ...] = ()
    metadata_timeout_seconds: float | None = None
    metadata_max_requests: int | None = None

    def __post_init__(self):
        settings(self.to_dict(), RUNTIME_DEFAULTS, "runtime options")
        get_resource_profile(self.resource_profile, overrides=self.resource_overrides)
        if not isinstance(self.allowed_hosts, (tuple, list)) or not all(isinstance(v, str) and v for v in self.allowed_hosts):
            raise TypeError("allowed_hosts must contain nonempty host names")
        object.__setattr__(self, "allowed_hosts", tuple(self.allowed_hosts))
        object.__setattr__(self, "reserved_paths", tuple(self.reserved_paths))
        if self.resource_overrides is not None:
            object.__setattr__(self, "resource_overrides", MappingProxyType(dict(self.resource_overrides)))

    def to_dict(self):
        return {f.name: (dict(self.resource_overrides) if self.resource_overrides is not None else None)
                if f.name == "resource_overrides" else getattr(self, f.name) for f in fields(self)}


def input_settings(value):
    return settings(value, set().union(*INPUT.values()), "input options")


def output_settings(value, target, *, partial=False):
    return settings(value, OUTPUT[target], "output options", partial=partial)


def merge_input_options(defaults, overrides):
    return dict(defaults) | input_settings(overrides)
