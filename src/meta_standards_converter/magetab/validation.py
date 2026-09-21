# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Validate archive MAGE-TAB references and occurrence-local qualifiers.

These checks do not impose an Atlas processing recipe or infer missing biology.
"""
import re


class MAGETabValidationError(ValueError):
    """A constructed document cannot safely be published as MAGE-TAB."""


def normalized(value):
    return "".join(str(value).split()).casefold()


def text(value):
    return "" if value is None else str(value).strip()


def used_term_sources(rows):
    sources = set()
    for row in rows:
        if not row:
            continue
        label = normalized(row[0])
        if label == "sdrffile" and len(row) > 1 and isinstance(row[1], list):
            table = row[1]
            if table:
                for i, header in enumerate(table[0]):
                    if normalized(header) == "termsourceref":
                        sources.update(text(r[i]) for r in table[1:] if i < len(r) and text(r[i]))
        elif label.endswith("termsourceref"):
            sources.update(text(v) for v in row[1:] if text(v))
    return sources


def _qualifiable(header):
    label = normalized(header)
    if re.fullmatch(r"comment\[hz_.*\]", label):
        from meta_standards_converter.miniml.harmonization import parse_harmonized_key
        return parse_harmonized_key(label[8:-1])[1] == "value"
    return label in {"provider", "materialtype", "arraydesignref", "arraydesignfile",
                     "technologytype", "label", "protocolref", "unit", "performer"} or bool(
        re.fullmatch(r"(?:characteristics|factorvalue|parametervalue|unit)\[.*\](?:\(.*\))?", label))


def validate_magetab(rows):
    """Raise before publication; return the original tables without modifying them."""
    declarations = {normalized(r[0]): r[1:] for r in rows if r}
    known_sources = {text(v) for v in declarations.get("termsourcename", []) if text(v)}
    missing = used_term_sources(rows) - known_sources
    if missing:
        raise MAGETabValidationError(f"Undeclared term source references: {', '.join(sorted(missing))}")
    protocols = {text(v) for v in declarations.get("protocolname", []) if text(v)}
    factors = {normalized(v) for v in declarations.get("experimentalfactorname", []) if text(v)}
    tables = [r[1] for r in rows if r and normalized(r[0]) == "sdrffile" and len(r) > 1]
    if not tables:
        raise MAGETabValidationError("MAGE-TAB requires an SDRF table")
    for table in tables:
        if not isinstance(table, list) or not table or not isinstance(table[0], (tuple, list)):
            raise MAGETabValidationError("SDRF File must contain a nonempty table before publication")
        headers = table[0]
        for i, header in enumerate(headers):
            label = normalized(header)
            if label in {"termsourceref", "termaccessionnumber"}:
                previous = headers[i - 1] if i else ""
                if not (_qualifiable(previous) or (label == "termaccessionnumber" and normalized(previous) == "termsourceref")):
                    raise MAGETabValidationError(f"SDRF column {i + 1}: unbound {header}")
            match = re.fullmatch(r"factorvalue\[(.*?)\](?:\(.*\))?", label)
            if match and match[1] not in factors:
                raise MAGETabValidationError(f"SDRF column {i + 1}: undeclared experimental factor {match[1]!r}")
        for number, row in enumerate(table[1:], 2):
            if len(row) != len(headers):
                raise MAGETabValidationError(f"SDRF row {number}: {len(row)} cells, expected {len(headers)}")
            for i, header in enumerate(headers):
                if normalized(header) == "protocolref" and text(row[i]) and text(row[i]) not in protocols:
                    raise MAGETabValidationError(f"SDRF row {number}, column {i + 1}: undeclared protocol {text(row[i])!r}")
    return rows
