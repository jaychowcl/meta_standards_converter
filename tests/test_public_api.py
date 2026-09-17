# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Small public API smoke contract; scientific correctness lives in tests/e2e."""
import importlib


def test_supported_owning_package_exports_import():
    exports = {
        "converters": ("Converter", "InputSpec", "ConversionBatchResult", "GEO2JSONConverter", "GEO2AEConverter", "AE2JSONConverter", "JSON2AEConverter", "JSON2TSVConverter", "JSON2H5ADConverter", "JSON2OBSConverter"),
        "sources": ("GEOXMLParser", "MAGETabSourceResolver", "PackageLoader"),
        "sources.json": ("JSONPackageSource",),
        "sources.insdc": ("INSDCWebfetcher",),
        "sources.magetab": ("AEWebFetcher",),
        "metadata": ("MetadataEnrichment",),
        "metadata.enrichment": ("MINiMLEnricher",),
        "metadata.projection": ("AnnDataMetadataProjection", "TabularMetadataProjection"),
        "expression": ("AssetDiscovery", "AssetReader", "SourcePlanner"),
        "retrieval": ("AssetDownloader",),
    }
    for package, symbols in exports.items():
        module = importlib.import_module("meta_standards_converter." + package)
        for symbol in symbols:
            assert getattr(module, symbol) is not None, (package, symbol)
