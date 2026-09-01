# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""ENA study to MINiML JSON converter."""

from meta_standards_converter.insdc_handlers.study_fetchers import ENAStudyFetcher
from meta_standards_converter.insdc_handlers.study_parser import INSDCStudyParser

from .insdc2json import INSDC2JSONConverter


class ena2json(INSDC2JSONConverter):
    provider = "ena"
    base_suffix = ".ena.json"

    def __init__(self, fetcher=None, parser=None, **kwargs):
        profile = kwargs.get("resource_profile", "standard")
        overrides = kwargs.get("resource_overrides")
        super().__init__(
            fetcher=fetcher or ENAStudyFetcher(
                resource_profile=profile, resource_overrides=overrides
            ),
            parser=parser or INSDCStudyParser(),
            **kwargs,
        )


__all__ = ["ena2json"]
