# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================

from importlib import import_module

_EXPORTS = {'ENASource': 'ena', 'SRASource': 'sra', 'GEOXMLParser': 'contracts', 'MAGETabSourceResolver': 'contracts', 'PackageLoader': 'contracts', 'INSDCClient': 'contracts', 'PubMedClient': 'contracts', 'MetricsProvider': 'contracts', 'RequestMetrics': 'contracts', 'JSONPackageSource': 'json', 'DatasetPackageGroup': 'json', 'SourceLoadResult': 'json', 'GEOSource': 'geo', 'GEOWebFetcher': 'geo', 'AEWebFetcher': 'magetab', 'MAGETabInput': 'magetab', 'TextResource': 'magetab'}
__all__ = list(_EXPORTS)

def __getattr__(name):
    """Load an owning-package export only when requested."""
    if name not in _EXPORTS:
        raise AttributeError(name)
    value = getattr(import_module("." + _EXPORTS[name], __name__), name)
    globals()[name] = value
    return value
