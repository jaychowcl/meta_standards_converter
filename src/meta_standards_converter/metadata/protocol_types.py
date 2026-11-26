# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Explicit protocol families; unknown authored types never imply extraction."""
import re
from collections.abc import Mapping

_TYPES = {
    'extraction': ('extract', 'extraction', 'extract protocol', 'nucleic acid extraction', 'nucleic acid extraction protocol'),
    'library': ('library construction', 'library preparation', 'library construction protocol', 'nucleic acid library construction protocol'),
    'treatment': ('treatment', 'treatment protocol', 'sample treatment protocol'),
    'growth': ('growth', 'growth protocol'),
    'label': ('label', 'labelling', 'labeling', 'label protocol', 'labelling protocol', 'labeling protocol'),
    'hybridization': ('hybridization', 'hybridisation', 'hybridization protocol', 'hybridisation protocol'),
    'scan': ('scan', 'scanning', 'scan protocol', 'array scanning protocol', 'array scanning and feature extraction protocol'),
    'processing': ('feature extraction', 'data processing', 'data processing protocol', 'normalization', 'normalisation', 'normalization data transformation protocol', 'data transformation', 'data transformation protocol'),
}
_FAMILIES = {name: family for family, names in _TYPES.items() for name in names}


def protocol_type_family(value) -> str | None:
    """Normalize separators/case and resolve an exact authored type alias."""
    if isinstance(value, Mapping):
        value = value.get('value', '')
    return _FAMILIES.get(re.sub(r'[\s_-]+', ' ', str(value or '')).strip().casefold())
