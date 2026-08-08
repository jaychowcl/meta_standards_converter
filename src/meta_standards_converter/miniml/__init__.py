# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Versioned data model for the MINiML-compatible JSON package."""

from .model import (
    MINIML_SCHEMA_VERSION,
    Accession,
    Address,
    Channel,
    Characteristics,
    Contributor,
    DataColumn,
    Database,
    DataTable,
    FASTQFile,
    InstrumentModel,
    MINiMLModelError,
    MINiMLPackage,
    MINiMLValidationIssue,
    Organization,
    Organism,
    Person,
    Platform,
    PubMedPublication,
    Reference,
    Relation,
    Repeat,
    Sample,
    Series,
    SRARun,
    Status,
    SupplementLink,
    TableData,
    Variable,
    miniml_schema_path,
)
from .codec import (
    MINiMLBatchDecodeResult,
    MINiMLCodec,
    MINiMLCompatibilityError,
    MINiMLDecodeResult,
)

__all__ = [
    "MINIML_SCHEMA_VERSION",
    "Accession",
    "Address",
    "Channel",
    "Characteristics",
    "Contributor",
    "DataColumn",
    "Database",
    "DataTable",
    "FASTQFile",
    "InstrumentModel",
    "MINiMLModelError",
    "MINiMLPackage",
    "MINiMLValidationIssue",
    "Organization",
    "Organism",
    "Person",
    "Platform",
    "PubMedPublication",
    "Reference",
    "Relation",
    "Repeat",
    "Sample",
    "Series",
    "SRARun",
    "Status",
    "SupplementLink",
    "TableData",
    "Variable",
    "miniml_schema_path",
    "MINiMLBatchDecodeResult",
    "MINiMLCodec",
    "MINiMLCompatibilityError",
    "MINiMLDecodeResult",
]
