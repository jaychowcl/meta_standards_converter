# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Permitted protocol-description equivalences; original source text is unchanged."""
import re


def comparable_protocol_text(value):
    value = " ".join(value.split()).translate(str.maketrans({"‘": "'", "’": "'", "“": '"', "”": '"'}))
    return re.sub(r"(?<=\d)\s*[uµμ]g\b", " µg", value)

