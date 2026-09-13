# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================

from importlib import import_module

_EXPORTS = {'ENA2JSONConverter': 'ena2json', 'SRA2JSONConverter': 'sra2json', 'GEO2JSONConverter': 'geo2json', 'GEO2AEConverter': 'geo2ae', 'AE2JSONConverter': 'ae2json', 'JSON2AEConverter': 'json2ae', 'JSON2TSVConverter': 'json2tsv', 'JSON2H5ADConverter': 'json2h5ad', 'JSON2OBSConverter': 'json2obs'}
__all__ = list(_EXPORTS)

def __getattr__(name):
    """Load an owning-package export only when requested."""
    if name not in _EXPORTS:
        raise AttributeError(name)
    value = getattr(import_module("." + _EXPORTS[name], __name__), name)
    globals()[name] = value
    return value
