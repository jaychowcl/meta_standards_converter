# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""
Placeholder converter for parsed MINiML JSON to AnnData H5AD output.
"""

import logging
import os


logger = logging.getLogger(__name__)


class json2h5ad:
    def convert(self, json_path: str, out: str | None = None) -> str:
        """
        Validate JSON input path, then fail loudly until H5AD conversion is implemented.
        """
        if not os.path.exists(json_path):
            raise FileNotFoundError(f"MINiML JSON file not found: {json_path}")

        logger.info("%s: validated parsed MINiML JSON input", json_path)
        raise NotImplementedError(
            "json2h5ad is a placeholder: matrix fetching/parsing and AnnData writing "
            "are not implemented yet."
        )
