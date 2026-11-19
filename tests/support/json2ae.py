# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Shared builders for json2ae."""


def package(accession="GSE1"):
    return {
        "miniml_schema_version": "3.1",
        "source": {"format": "test"},
        "series": {
            "accession": [{"value": accession, "database": "GEO"}],
            "title": f"Study {accession}",
        },
        "sample": [{"iid": f"GSM-{accession}"}],
        "platform": [],
    }
