# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Render prepared package groups, including explicit replacement profiles."""
import logging
from meta_standards_converter.conversion_logging import _retained_warning
from meta_standards_converter.miniml import MINiMLCodec, MINiMLPackage
from meta_standards_converter.sources.json import DatasetPackageGroup
from .constructor import AEConstructor
from .technology import series_identity
from ._performance import measured_phase
logger = logging.getLogger(__name__)


class MAGETabExporter:
    def __init__(self, constructor=None):
        self.constructor = constructor if constructor is not None else AEConstructor()

    def render(self, loaded, *, platform_handler=None, replacement_profile=None):
        packages = self._packages_from_loaded(loaded, replacement_profile=replacement_profile)
        options = {"platform_handler": platform_handler} if platform_handler is not None else {}
        return [self.constructor.miniml2magetab(data=package, **options)
                for package in packages]

    @measured_phase("magetab_package_resolution")
    def _packages_from_loaded(self, loaded, *, replacement_profile=None):
        for warning in loaded.warnings:
            _retained_warning(logger, "%s", warning)
        if not loaded.groups:
            raise ValueError("JSON source contains no convertible package groups.")
        packages = []
        for original_group in loaded.groups:
            group = original_group.resolved(replacement_profile=replacement_profile)
            for warning in getattr(group.harmonization_resolution, "warnings", ()):
                _retained_warning(logger, "%s", warning)
            packages.extend(group.packages)

        for index, package in enumerate(packages, start=1):
            if not isinstance(package, MINiMLPackage):
                raise ValueError(f"Parsed MINiML package {index} must decode as MINiMLPackage.")
            if not self._usable_study_accession(package):
                raise ValueError(f"Parsed MINiML package {index} has no usable study accession.")
        return packages

    def _usable_study_accession(self, package: MINiMLPackage) -> str | None:
        return series_identity({"series": {
            "iid": package.series.iid,
            "accession": [{"value": accession.value} for accession in package.series.accessions],
        }})
