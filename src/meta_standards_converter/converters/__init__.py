# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================

from .json2h5ad import (
    AnnDataMetadataProjection,
    AnnDataMetadataProjector,
    MetadataProjectionContext,
)

__all__ = [
    "AnnDataMetadataProjection",
    "AnnDataMetadataProjector",
    "MetadataProjectionContext",
]
