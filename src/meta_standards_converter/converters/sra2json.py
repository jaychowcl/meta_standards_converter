# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""NCBI SRA study to MINiML JSON converter."""

from meta_standards_converter.insdc_handlers.study_fetchers import SRAStudyFetcher
from meta_standards_converter.insdc_handlers.study_parser import INSDCStudyParser

from .insdc2json import INSDC2JSONConverter


class sra2json(INSDC2JSONConverter):
    provider = "sra"
    base_suffix = ".sra.json"

    def __init__(self, fetcher=None, parser=None, **kwargs):
        profile = kwargs.get("resource_profile", "standard")
        overrides = kwargs.get("resource_overrides")
        super().__init__(
            fetcher=fetcher or SRAStudyFetcher(
                resource_profile=profile, resource_overrides=overrides
            ),
            parser=parser or INSDCStudyParser(),
            **kwargs,
        )


__all__ = ["sra2json"]
