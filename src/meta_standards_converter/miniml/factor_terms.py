# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Existing MINiML factor categories shared by archive and MAGE-TAB imports."""
MINIML_VARIABLE_FACTORS = {
    "dose", "time", "tissue", "strain", "gender", "cell line",
    "development stage", "age", "agent", "cell type", "infection",
    "isolate", "metabolism", "shock", "stress", "temperature",
    "speciman", "disease state", "protocol", "growth protocol", "other",
    "genotype/variation", "species", "individual",
}
FACTOR_NORMALIZATION = {
    "compound": "agent",
    "treatment": "agent",
    "drug": "agent",
    "disease": "disease state",
    "organism": "species",
    "sex": "gender",
    "genotype": "genotype/variation",
}


def normalized_factor_category(value):
    normalized = " ".join(value.strip().casefold().split())
    mapped = FACTOR_NORMALIZATION.get(normalized, normalized)
    return mapped if mapped in MINIML_VARIABLE_FACTORS else "other"
