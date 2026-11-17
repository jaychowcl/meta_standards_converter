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
