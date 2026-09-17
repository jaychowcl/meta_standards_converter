# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Explicit loaded-data routes; no graph search or provider fallback."""

from dataclasses import dataclass

from .contracts import InputError


@dataclass(frozen=True)
class Route:
    content: str
    destination: str
    service: str


ROUTES = {
    (content, destination): Route(content, destination, service)
    for content, destination, service in (
        ("metadata", "json", "json"),
        ("metadata", "magetab", "magetab"),
        ("metadata", "tsv", "table"),
        ("metadata", "csv", "table"),
        ("metadata", "h5ad", "expression"),
        ("metadata", "obs", "expression"),
        ("expression", "h5ad", "expression"),
        ("expression", "obs", "expression"),
        ("direct_magetab", "magetab", "direct_magetab"),
    )
}


def select_route(loaded, source_kind, destination):
    if loaded.direct_magetab is not None:
        content = "direct_magetab"
    elif loaded.expression is not None and destination in {"h5ad", "obs"}:
        content = "expression"
    elif loaded.metadata is not None:
        content = "metadata"
    else:
        raise InputError(
            "metadata_required",
            "Destination requires supplied or valid embedded MSC metadata",
        )
    route = ROUTES.get((content, destination))
    if route is None:
        raise InputError(
            "unsupported_route", "No registered route for this content and destination"
        )
    steps = (
        (source_kind, "miniml", destination)
        if content == "metadata"
        else (source_kind, destination)
    )
    return route, steps
