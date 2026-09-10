# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from meta_standards_converter.magetab.writer import MAGETabWriter
"""Converter for parsed MINiML JSON to ArrayExpress MAGE-TAB format."""

import logging

from meta_standards_converter.magetab.technology import series_identity
from meta_standards_converter.magetab.constructor import AEConstructor
from meta_standards_converter.sources.json import JSONPackageSource
from meta_standards_converter.metadata.enrichment import MINiMLEnricher
from meta_standards_converter.helpers.json_helper import JSONHandler
from meta_standards_converter.miniml import MINiMLCodec, MINiMLPackage


logger = logging.getLogger(__name__)


class JSON2AEConverter(JSONHandler):
    """Convert parsed MINiML JSON packages into MAGE-TAB payloads."""

    def __init__(self, enricher=None, ae_constructor=None, package_source=None):
        self.enricher = enricher or MINiMLEnricher()
        self.ae_constructor = ae_constructor or AEConstructor()
        self.package_source = package_source or JSONPackageSource()

    def convert(
        self,
        json_path: str,
        out: str = None,
        enrich: bool = True,
        platform_handler: str | None = None,
        use_harmonization_overrides: bool = False,
    ) -> list[list]:
        """Load parsed MINiML JSON and optionally write IDF/SDRF files."""
        packages = self._load_packages(
            json_path=json_path,
            use_harmonization_overrides=use_harmonization_overrides,
        )
        logger.debug("%s: loaded %d parsed package(s)", json_path, len(packages))

        magetabs = []
        codec = MINiMLCodec()
        for index, package in enumerate(packages, start=1):
            converted_package = package
            if enrich:
                logger.info("%s: enriching parsed package %d", json_path, index)
                converted_package = codec.decode(
                    self.enricher.enrich(data=package)
                ).package
            else:
                logger.info("%s: skipping enrichment for parsed package %d", json_path, index)
            logger.info("%s: building MAGE-TAB package %d", json_path, index)
            if platform_handler is None:
                magetab = self.ae_constructor.miniml2magetab(data=converted_package)
            else:
                magetab = self.ae_constructor.miniml2magetab(
                    data=converted_package,
                    platform_handler=platform_handler,
                )
            magetabs.append(magetab)

        if out:
            for index, magetab in enumerate(magetabs, start=1):
                logger.info("%s: writing MAGE-TAB package %d to %s", json_path, index, out)
                MAGETabWriter().write(magetab=magetab, out=out)

        logger.info("%s: conversion produced %d MAGE-TAB package(s)", json_path, len(magetabs))
        return magetabs

    def _load_packages(
        self, json_path: str, *, use_harmonization_overrides: bool = False
    ) -> list[MINiMLPackage]:
        try:
            loaded = self.package_source.load(json_path)
        except FileNotFoundError as error:
            raise FileNotFoundError(
                f"MINiML JSON file not found: {json_path}"
            ) from error

        for warning in loaded.warnings:
            logger.warning("%s", warning)
        if not loaded.groups:
            raise ValueError("JSON source contains no convertible package groups.")
        packages = []
        for original_group in loaded.groups:
            group = original_group.resolved(enabled=use_harmonization_overrides)
            for warning in getattr(group.harmonization_resolution, "warnings", ()):
                logger.warning("%s", warning)
            packages.extend(group.packages)

        for index, package in enumerate(packages, start=1):
            if not isinstance(package, MINiMLPackage):
                raise ValueError(f"Parsed MINiML package {index} must decode as MINiMLPackage.")
            if not self._usable_study_accession(package):
                raise ValueError(f"Parsed MINiML package {index} has no usable study accession.")
        return packages

    def _usable_study_accession(self, package: MINiMLPackage) -> str | None:
        return series_identity(package.to_mapping())
