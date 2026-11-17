# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Verified direct-parent publication inheritance, independent of output expansion."""
from collections.abc import Mapping, Sequence
import logging
from meta_standards_converter.conversion_logging import _retained_warning
import re
_GSE_ACCESSION = re.compile(r"\bGSE\d+\b", re.IGNORECASE)
from meta_standards_converter.miniml import MINiMLCodec, MINiMLPackage
logger = logging.getLogger(__name__)

def _as_list(value):
    if value is None:
        return []
    if isinstance(value, (list, tuple)):
        return list(value)
    return [value]

def _gse_values(value) -> list[str]:
    values = []
    for item in _as_list(value):
        if isinstance(item, Mapping):
            candidates = (item.get("value"), item.get("target"), item.get("comment"))
        else:
            candidates = (item,)
        for candidate in candidates:
            if not isinstance(candidate, str):
                continue
            for match in _GSE_ACCESSION.findall(candidate):
                normalized = match.upper()
                if normalized not in values:
                    values.append(normalized)
    return values

def _series_accessions(series: Mapping) -> list[str]:
    return _gse_values([series.get("iid"), *_as_list(series.get("accession"))])

def _relation_gses(series: Mapping, relation_type: str) -> list[tuple[str, str]]:
    matches = []
    expected = relation_type.strip().casefold()
    for relation in _as_list(series.get("relation")):
        if not isinstance(relation, Mapping):
            continue
        rendered_type = str(relation.get("type") or "").strip()
        if rendered_type.casefold() != expected:
            continue
        for accession in _gse_values(
            [relation.get("target"), relation.get("comment")]
        ):
            item = (accession, rendered_type)
            if item not in matches:
                matches.append(item)
    return matches


class GEOParentPublications:
    def __init__(self, source, session, available=None):
        self.source, self.session = source, session
        self.available = available or {}

    def inherit(
        self,
        *,
        packages: Sequence[MINiMLPackage],
        remove_empty: bool,
    ) -> tuple[list[MINiMLPackage], list[str]]:
        """Hydrate one unambiguous direct parent publication without adding packages."""

        codec = MINiMLCodec()
        canonical = [codec.encode(codec.decode(package).package) for package in packages]
        available = dict(self.available)
        for package in canonical:
            series = package.get("series")
            if not isinstance(series, Mapping):
                continue
            for accession in _series_accessions(series):
                available.setdefault(accession, package)

        hydrated = []
        issues = []
        for package in canonical:
            series = package.get("series")
            if not isinstance(series, dict):
                hydrated.append(codec.decode(package).package)
                continue
            if _as_list(series.get("pubmed_id")) or _as_list(
                series.get("pubmed_publication")
            ):
                hydrated.append(codec.decode(package).package)
                continue
            child_accessions = _series_accessions(series)
            parent_relations = _relation_gses(series, "SubSeries of")
            parent_accessions = list(dict.fromkeys(item[0] for item in parent_relations))
            if len(child_accessions) != 1 or len(parent_accessions) != 1:
                hydrated.append(codec.decode(package).package)
                continue

            child_accession = child_accessions[0]
            parent_accession = parent_accessions[0]
            parent_package = available.get(parent_accession)
            if parent_package is None:
                try:
                    parsed_parents = self.session.load(self.source, parent_accession, remove_empty=remove_empty)
                    parent_candidates = []
                    for parsed_parent in parsed_parents:
                        candidate = codec.encode(codec.decode(parsed_parent, strict=True).package)
                        candidate_series = candidate.get("series")
                        if (
                            isinstance(candidate_series, Mapping)
                            and parent_accession in _series_accessions(candidate_series)
                        ):
                            parent_candidates.append(candidate)
                    if len(parent_candidates) == 1:
                        parent_package = parent_candidates[0]
                except Exception as error:
                    issues.append(f"{child_accession}: parent {parent_accession} publication retrieval failed ({type(error).__name__})")
                    _retained_warning(logger,
                        "%s: parent-series publication lookup degraded "
                        "parent=%s error_type=%s",
                        child_accession,
                        parent_accession,
                        type(error).__name__,
                    )

            parent_series = (
                parent_package.get("series")
                if isinstance(parent_package, Mapping)
                else None
            )
            if not isinstance(parent_series, Mapping):
                hydrated.append(codec.decode(package).package)
                continue
            reciprocal_children = {
                accession
                for accession, _relation_type in _relation_gses(
                    parent_series, "SuperSeries of"
                )
            }
            if child_accession not in reciprocal_children:
                hydrated.append(codec.decode(package).package)
                continue
            pubmed_ids = list(
                dict.fromkeys(
                    str(value).strip()
                    for value in _as_list(parent_series.get("pubmed_id"))
                    if str(value).strip()
                )
            )
            if len(pubmed_ids) != 1:
                hydrated.append(codec.decode(package).package)
                continue

            relation_type = next(
                relation_type
                for accession, relation_type in parent_relations
                if accession == parent_accession
            )
            extensions = package.get("extensions")
            if extensions is None:
                extensions = {}
            if not isinstance(extensions, Mapping):
                hydrated.append(codec.decode(package).package)
                continue
            extensions = dict(extensions)
            extensions["publication_inheritance"] = {
                "relation_type": relation_type,
                "source_series": parent_accession,
                "pubmed_ids": pubmed_ids,
            }
            series["pubmed_id"] = pubmed_ids
            package["extensions"] = extensions
            hydrated.append(codec.decode(package).package)
        return hydrated, issues
