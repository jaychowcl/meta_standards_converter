"""Presence of table identities without rewriting meaningful source text."""


def has_cell_value(value):
    return value is not None and bool(str(value).strip())
