# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""
Parser for GEO MINiML XML data.
"""

from __future__ import annotations

import re
import logging
import time
from collections import deque
from xml.etree import ElementTree as ET

from meta_standards_converter.miniml import MINiMLPackage, MINiMLV1Migrator
from meta_standards_converter.runtime_contracts import (
    CompletenessStatus,
    EvidenceConfidence,
    ExecutionStatus,
    OperationStatusV2,
    PublicationDisposition,
    SafeErrorEnvelope,
    ValidationStatus,
    get_resource_profile,
)
from meta_standards_converter.xml_safety import parse_xml


logger = logging.getLogger(__name__)




class GEOParser:
    def __init__(
        self,
        resource_profile: str = "standard",
        resource_overrides=None,
    ):
        self.resource_profile = get_resource_profile(
            resource_profile,
            overrides=resource_overrides,
        )
        self.repeated_children = {
            "MINiML": {
                "Organization",
                "Contributor",
                "Database",
                "Platform",
                "Sample",
                "Series",
            },
            "Platform": {
                "Status",
                "Accession",
                "Organism",
                "Web-Link",
                "Pubmed-ID",
                "Citation",
                "Contributor-Ref",
                "Contributor",
                "Contact-Ref",
                "Contact",
                "Supplementary-Data",
                "Relation",
            },
            "Sample": {
                "Status",
                "Accession",
                "Channel",
                "Contact-Ref",
                "Contact",
                "Supplementary-Data",
                "Raw-Data",
                "Relation",
            },
            "Series": {
                "Status",
                "Accession",
                "Pubmed-ID",
                "Citation",
                "Web-Link",
                "Type",
                "Contributor-Ref",
                "Contributor",
                "Contact-Ref",
                "Contact",
                "Sample-Ref",
                "Variable",
                "Repeats",
                "Supplementary-Data",
                "Relation",
                "Data-Table",
            },
            "Status": {"Comment"},
            "Address": {"Line"},
            "Channel": {
                "Organism",
                "Characteristics",
                "Biomaterial-Provider",
            },
            "Data-Table": {"Column"},
            "Variable": {"Sample-Ref"},
            "Repeats": {"Sample-Ref"},
        }

    def parse(
        self,
        miniml: str,
        remove_empty: bool = False,
    ) -> list[MINiMLPackage]:
        started = time.monotonic()
        parsed = self.parse_mapping(miniml=miniml)

        if remove_empty:
            parsed = [self.remove_empty_fields(series_package) for series_package in parsed]

        migrator = MINiMLV1Migrator()
        packages = [migrator.migrate(package).package for package in parsed]

        logger.info(
            "MINiML parse stats packages=%s series=%s samples=%s platforms=%s related_series=%s remove_empty=%s elapsed_seconds=%.3f",
            len(packages),
            len(packages),
            sum(len(package.samples) for package in packages),
            sum(len(package.platforms) for package in packages),
            False,
            remove_empty,
            time.monotonic() - started,
        )

        return packages

    def parse_mapping(self, miniml: str) -> list[dict]:
        root = parse_xml(
            miniml,
            max_bytes=self.resource_profile.max_xml_bytes,
        )
        top_level = self._top_level_nodes(root=root)
        parsed_top_level = {
            name: [self._parse_element(node) for node in nodes]
            for name, nodes in top_level.items()
        }
        indexes = self._build_indexes(parsed_top_level=parsed_top_level)

        return [
            self._series_package(
                root=root,
                series=series,
                indexes=indexes,
            )
            for series in parsed_top_level["series"]
        ]

    def remove_empty_fields(self, data):
        return self._remove_empty_fields(data)



    def _top_level_nodes(self, root: ET.Element) -> dict[str, list[ET.Element]]:
        top_level = {
            "organization": [],
            "contributor": [],
            "database": [],
            "platform": [],
            "sample": [],
            "series": [],
        }
        for child in root:
            key = self._to_snake_case(self._local_name(child.tag))
            if key in top_level:
                top_level[key].append(child)
        return top_level

    def _build_indexes(self, parsed_top_level: dict[str, list[dict]]) -> dict[str, dict]:
        return {
            top_level_key: {
                item["iid"]: item
                for item in items
                if isinstance(item, dict) and item.get("iid")
            }
            for top_level_key, items in parsed_top_level.items()
        }

    def _series_package(
        self,
        root: ET.Element,
        series: dict,
        indexes: dict[str, dict],
    ) -> dict:
        samples = self._resolve_samples(series=series, samples_by_iid=indexes["sample"])
        platforms = self._resolve_platforms(
            samples=samples,
            platforms_by_iid=indexes["platform"],
        )
        contributors = self._resolve_contributors(
            series=series,
            samples=samples,
            platforms=platforms,
            contributors_by_iid=indexes["contributor"],
        )
        databases = self._resolve_databases(
            elements=[series, *samples, *platforms, *contributors],
            databases_by_iid=indexes["database"],
        )
        organizations = self._resolve_organizations(
            elements=[*contributors, *databases],
            organizations_by_iid=indexes["organization"],
        )

        package = {
            "version": root.attrib.get("version"),
            "database": databases,
            "organization": organizations,
            "contributor": contributors,
            "platform": platforms,
            "sample": samples,
            "series": series,
        }
        self._attach_namespaced_root_attributes(root=root, package=package)
        return package

    def _resolve_samples(self, series: dict, samples_by_iid: dict[str, dict]) -> list[dict]:
        values = self._as_list(series.get("sample_ref"))
        values = sorted(values, key=lambda item: self._position_key(item, values.index(item)))
        refs = [
            ref.get("ref")
            for ref in values
            if isinstance(ref, dict)
        ]
        return self._items_for_refs(refs=refs, index=samples_by_iid)

    @staticmethod
    def _position_key(item, fallback):
        try:
            return (0, int(item["position"])) if isinstance(item, dict) and item.get("position") not in (None, "") else (1, fallback)
        except (TypeError, ValueError):
            return (1, fallback)

    def _resolve_platforms(
        self,
        samples: list[dict],
        platforms_by_iid: dict[str, dict],
    ) -> list[dict]:
        refs = []
        for sample in samples:
            platform_ref = sample.get("platform_ref")
            if isinstance(platform_ref, dict) and platform_ref.get("ref"):
                refs.append(platform_ref["ref"])
        return self._items_for_refs(refs=refs, index=platforms_by_iid)

    def _resolve_contributors(
        self,
        series: dict,
        samples: list[dict],
        platforms: list[dict],
        contributors_by_iid: dict[str, dict],
    ) -> list[dict]:
        refs = []
        for element in [series, *samples, *platforms]:
            refs.extend(self._reference_values(element=element, keys={"contributor_ref", "contact_ref"}))
        return self._items_for_refs(refs=refs, index=contributors_by_iid)

    def _resolve_databases(
        self,
        elements: list[dict],
        databases_by_iid: dict[str, dict],
    ) -> list[dict]:
        refs = []
        for element in elements:
            for child in self._walk_dicts(element):
                for accession in self._as_list(child.get("accession")):
                    if isinstance(accession, dict) and accession.get("database"):
                        refs.append(accession["database"])
                for status in self._as_list(child.get("status")):
                    if isinstance(status, dict) and status.get("database"):
                        refs.append(status["database"])
        return self._items_for_refs(refs=refs, index=databases_by_iid)

    def _resolve_organizations(
        self,
        elements: list[dict],
        organizations_by_iid: dict[str, dict],
    ) -> list[dict]:
        refs = []
        for element in elements:
            refs.extend(self._reference_values(element=element, keys={"organization_ref"}))
        return self._items_for_refs(refs=refs, index=organizations_by_iid)

    def _items_for_refs(self, refs: list[str], index: dict[str, dict]) -> list[dict]:
        items = []
        seen = set()
        for ref in refs:
            if not ref or ref in seen or ref not in index:
                continue
            seen.add(ref)
            items.append(index[ref])
        return items

    def _reference_values(self, element: dict, keys: set[str]) -> list[str]:
        refs = []
        for child in self._walk_dicts(element):
            for key, values in child.items():
                if key not in keys:
                    continue
                for value in self._as_list(values):
                    if isinstance(value, dict) and value.get("ref"):
                        refs.append(value["ref"])
        return refs

    def _parse_element(self, node: ET.Element):
        attrs = {
            self._to_snake_case(self._local_name(name)): value
            for name, value in node.attrib.items()
        }
        children = list(node)
        direct_text = self._normalized_text(node.text)

        if not attrs and not children:
            return direct_text

        parsed = dict(attrs)
        if children:
            grouped = {}
            repeatable_keys = set()
            parent_name = self._local_name(node.tag)
            for child in children:
                child_name = self._local_name(child.tag)
                key = self._child_key(parent_name=parent_name, child_name=child_name)
                if child_name in self.repeated_children.get(parent_name, set()):
                    repeatable_keys.add(key)
                grouped.setdefault(key, []).append(self._parse_element(child))

            for key, values in grouped.items():
                parsed[key] = values if key in repeatable_keys or len(values) > 1 else values[0]

        if direct_text is not None:
            parsed["value"] = direct_text

        return parsed

    def _child_key(self, parent_name: str, child_name: str) -> str:
        return self._to_snake_case(child_name)

    def series_accessions(self, series_packages: list[dict]) -> list[str]:
        accessions = []
        for package in series_packages:
            series = package.get("series", {})
            for accession in self._as_list(series.get("accession")):
                if not isinstance(accession, dict):
                    continue
                value = accession.get("value")
                if isinstance(value, str) and re.fullmatch(r"GSE\d+", value, re.IGNORECASE):
                    normalized = value.upper()
                    if normalized not in accessions:
                        accessions.append(normalized)
        return accessions

    def related_accessions(self, series_packages: list[dict]) -> list[str]:
        related_gses = []
        for package in series_packages:
            series = package.get("series", {})
            for relation in self._as_list(series.get("relation")):
                if not isinstance(relation, dict) or not self._is_related_series_relation(relation):
                    continue
                for value in (relation.get("type"), relation.get("target"), relation.get("comment")):
                    if not isinstance(value, str):
                        continue
                    for match in re.findall(r"GSE\d+", value, re.IGNORECASE):
                        normalized = match.upper()
                        if normalized not in related_gses:
                            related_gses.append(normalized)
        return related_gses

    def _is_related_series_relation(self, relation: dict) -> bool:
        relation_text = " ".join(
            value
            for value in (relation.get("type"), relation.get("target"), relation.get("comment"))
            if isinstance(value, str)
        ).lower()
        return "superseries" in relation_text or "subseries" in relation_text

    def _remove_empty_fields(self, value):
        if isinstance(value, dict):
            cleaned = {}
            for key, child_value in value.items():
                cleaned_child = self._remove_empty_fields(child_value)
                if not self._is_empty_value(cleaned_child):
                    cleaned[key] = cleaned_child
            return cleaned

        if isinstance(value, list):
            cleaned = []
            for item in value:
                cleaned_item = self._remove_empty_fields(item)
                if not self._is_empty_value(cleaned_item):
                    cleaned.append(cleaned_item)
            return cleaned

        return value

    def _is_empty_value(self, value) -> bool:
        return value is None or value == "" or value == [] or value == {}

    def _walk_dicts(self, value):
        if isinstance(value, dict):
            yield value
            for child_value in value.values():
                yield from self._walk_dicts(child_value)
        elif isinstance(value, list):
            for item in value:
                yield from self._walk_dicts(item)

    def _as_list(self, value) -> list:
        if value is None:
            return []
        if isinstance(value, list):
            return value
        return [value]

    def _attach_namespaced_root_attributes(self, root: ET.Element, package: dict) -> None:
        for name, value in root.attrib.items():
            key = self._to_snake_case(self._local_name(name))
            if key != "version":
                package[key] = value

    def _normalized_text(self, text: str | None) -> str | None:
        if text is None:
            return None
        normalized = re.sub(r"\s+", " ", text).strip()
        return normalized or None

    def _local_name(self, tag: str) -> str:
        if "}" in tag:
            return tag.split("}", 1)[1]
        return tag

    def _to_snake_case(self, value: str) -> str:
        value = value.replace("-", "_").replace(" ", "_").replace("/", "_")
        value = re.sub(r"([a-z0-9])([A-Z])", r"\1_\2", value)
        return value.lower()
