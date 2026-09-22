# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Curated canonical aliases; raw characteristic names are not rewritten."""
BIOLOGICAL_ALIASES = {
    "organism_part": {"organism_part", "organismpart", "tissue"},
    "developmental_stage": {"developmental_stage", "developmentalstage"},
}


def biological_characteristics(channels, field, *, slug):
    """Yield explicit source occurrences, including their attached qualifiers."""
    from collections.abc import Mapping
    names = BIOLOGICAL_ALIASES[field]
    for channel in channels:
        items = channel.get("characteristics") or ()
        if isinstance(items, Mapping):
            items = (items,)
        for item in items:
            if isinstance(item, Mapping) and slug(item.get("name", item.get("tag", ""))) in names:
                yield item
