# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from __future__ import annotations


import logging


from dataclasses import dataclass, field


from typing import Any, Mapping, Protocol, Sequence


from meta_standards_converter.expression.assets import Asset

logger = logging.getLogger(__name__)

@dataclass(frozen=True)
class MetadataProjectionContext:
    """Read-only conversion context supplied to metadata projectors."""

    sample: Mapping[str, Any]
    package: Mapping[str, Any]
    study_accession: str
    sample_accession: str
    asset: Asset
    base_metadata: Mapping[str, Any]


@dataclass(frozen=True)
class AnnDataMetadataProjection:
    """Metadata additions returned by an :class:`AnnDataMetadataProjector`."""

    obs: Mapping[str, Any] = field(default_factory=dict)
    var: Mapping[str, Any] = field(default_factory=dict)
    uns: Mapping[str, Any] = field(default_factory=dict)
    obs_renames: Mapping[str, str] = field(default_factory=dict)
    obs_drops: tuple[str, ...] = ()
    warnings: tuple[str, ...] = ()
    errors: tuple[str, ...] = ()


class AnnDataProjectionError(ValueError):
    """Raised when an AnnData projector reports invalid projected metadata."""

    def __init__(self, errors: Sequence[str]):
        self.errors = tuple(str(error) for error in errors)
        super().__init__("; ".join(self.errors))


class AnnDataMetadataProjector(Protocol):
    """Optional extension that adds organization-neutral metadata to AnnData."""

    def project_sample(
        self,
        *,
        adata: Any,
        context: MetadataProjectionContext,
    ) -> AnnDataMetadataProjection:
        """Return metadata additions for one sample AnnData object."""


class AnnDataProjectorRunner:
    def __init__(self, projectors=()):
        self.metadata_projectors = tuple(projectors)

    def project_sample(
        self,
        adata,
        context: MetadataProjectionContext,
        *,
        warnings: list[str],
        errors: list[str],
        allow_invalid: bool,
    ) -> None:
        for projector in self.metadata_projectors:
            callback = getattr(projector, "project_sample", None)
            if callback is None:
                continue
            projection = callback(adata=adata, context=context)
            self._apply_metadata_projection(adata, projection)
            self._extend_warnings(warnings, projection.warnings)
            self._extend_warnings(errors, projection.errors)
            if projection.errors and not allow_invalid:
                raise AnnDataProjectionError(projection.errors)

    def _apply_metadata_projection(
        self,
        adata,
        projection: AnnDataMetadataProjection,
    ) -> None:
        if not isinstance(projection, AnnDataMetadataProjection):
            raise TypeError(
                "metadata projector must return AnnDataMetadataProjection"
            )
        self._apply_obs_transforms(
            adata.obs,
            renames=projection.obs_renames,
            drops=projection.obs_drops,
        )
        self._apply_axis_projection(
            adata.obs,
            projection.obs,
            axis_name="obs",
            axis_length=adata.n_obs,
        )
        self._apply_axis_projection(
            adata.var,
            projection.var,
            axis_name="var",
            axis_length=adata.n_vars,
        )
        for key, value in projection.uns.items():
            if key in adata.uns:
                raise ValueError(f"uns metadata key {key!r} already exists")
            adata.uns[key] = value

    @staticmethod
    def _apply_obs_transforms(
        frame,
        *,
        renames: Mapping[str, str],
        drops: Sequence[str],
    ) -> None:
        renames = dict(renames)
        drops = tuple(drops)
        drop_set = set(drops)
        if len(drop_set) != len(drops):
            raise ValueError("obs drop columns must be unique")
        overlap = sorted(set(renames) & drop_set)
        if overlap:
            raise ValueError(
                f"obs metadata key {overlap[0]!r} cannot be renamed and dropped"
            )
        for source, target in renames.items():
            if not isinstance(source, str) or not isinstance(target, str):
                raise TypeError("obs rename sources and targets must be strings")
            if source not in frame:
                raise ValueError(f"obs rename source {source!r} does not exist")
            if not target:
                raise ValueError(f"obs rename target for {source!r} must be nonblank")
        missing_drops = [column for column in drops if column not in frame]
        if missing_drops:
            raise ValueError(f"obs drop source {missing_drops[0]!r} does not exist")
        targets = list(renames.values())
        if len(set(targets)) != len(targets):
            raise ValueError("obs rename targets must be unique")
        survivors = set(frame.columns) - set(renames) - drop_set
        collisions = sorted(set(targets) & survivors)
        if collisions:
            raise ValueError(
                f"obs rename target {collisions[0]!r} already exists"
            )
        if drops:
            frame.drop(columns=list(drops), inplace=True)
        if renames:
            frame.rename(columns=renames, inplace=True)

    @staticmethod
    def _apply_axis_projection(
        frame,
        values: Mapping[str, Any],
        *,
        axis_name: str,
        axis_length: int,
    ) -> None:
        for key, value in values.items():
            if key in frame:
                raise ValueError(f"{axis_name} metadata key {key!r} already exists")
            if isinstance(value, Sequence) and not isinstance(
                value, (str, bytes, bytearray)
            ):
                if len(value) != axis_length:
                    raise ValueError(
                        f"{axis_name} metadata key {key!r} must contain "
                        f"{axis_length} values"
                    )
                frame[key] = list(value)
            else:
                frame[key] = "" if value is None else value

    @staticmethod
    def _extend_warnings(target: list[str], values: Sequence[str]) -> None:
        for value in values:
            rendered = str(value)
            if rendered and rendered not in target:
                target.append(rendered)
