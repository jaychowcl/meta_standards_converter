# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Unified conversion interfaces."""

from .contracts import InputSpec, ConversionBatchResult
from .engine import Converter

__all__ = ["Converter", "InputSpec", "ConversionBatchResult"]
