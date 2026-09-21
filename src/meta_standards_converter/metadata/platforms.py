"""Repository namespaces of platform identifiers, without accession rewriting."""
import re


def platform_namespace(identifier, declared=None):
    if declared:
        return str(declared)
    if re.fullmatch(r"GPL\d+", str(identifier), re.I):
        return "GEO"
    if re.fullmatch(r"A-[A-Z]+-\d+", str(identifier), re.I):
        return "ArrayExpress"
    return None
