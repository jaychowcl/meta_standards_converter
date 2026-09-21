# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""One strict MINiML policy at every metadata acceptance boundary."""
from meta_standards_converter.miniml import MINiMLCodec


def validate_metadata(metadata):
    if metadata is not None:
        codec = MINiMLCodec()
        for group in metadata.groups:
            for package in group.packages:
                codec.decode(package, strict=True)
