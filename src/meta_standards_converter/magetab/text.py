# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Authored field text and the final, minimal TSV transport adaptation."""
import re


def literal(value):
    if isinstance(value, dict) and 'value' in value:
        value = value['value']
    return None if value is None else str(value)


def tsv_cell(value):
    return '' if value is None else re.sub(r'\r\n|[\t\r\n]', ' ', str(value))


def geo_origin(data):
    from ..metadata.reconciliation import provider
    return provider(data) == 'GEO'
