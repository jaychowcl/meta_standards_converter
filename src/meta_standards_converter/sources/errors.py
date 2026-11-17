# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Safe source failures that can cross the conversion/report boundary."""
import re


class MetadataSourceError(ValueError):
    """A classified provider condition, with controlled recovery guidance.

    Never construct messages from response bodies or arbitrary exception text.
    Availability failures are terminal for automatic retries; a user can always
    retry explicitly when repository access or supplied metadata changes.
    """

    MESSAGES = {
        "source_inaccessible": "BioStudies reports Study not accessible. Check repository access or supply the study's IDF and SDRF files.",
        "source_not_found": "BioStudies did not find the study. Check the accession or supply its IDF and SDRF files.",
        "missing_idf": "The accessible study has no IDF file. Supply an IDF and its SDRFs.",
        "missing_sdrf": "The accessible study has no SDRF file. Supply the missing SDRF with its IDF.",
        "ambiguous_idf": "The study exposes multiple IDF files. Select an IDF and its SDRFs explicitly.",
        "ambiguous_sdrf": "The study exposes duplicate SDRF paths. Select an unambiguous local MAGE-TAB bundle.",
        "provider_response_invalid": "The provider returned an invalid or incomplete metadata response. Retry later; this does not establish that files are missing.",
        "fallback_unverified": "No uniquely verified public ENA study links to this ArrayExpress accession. Supply the study's IDF and SDRF files.",
        "fallback_ambiguous": "Multiple distinct public ENA studies link to this ArrayExpress accession. Select the correct source explicitly.",
    }
    FALLBACK_CODES = frozenset({"source_inaccessible", "source_not_found", "missing_idf", "missing_sdrf"})

    def __init__(self, code, accession, *, provider="biostudies", http_status=None):
        self.code = code
        self.accession = accession if re.fullmatch(r"[A-Z0-9-]+", accession) else "study"
        self.provider = provider
        self.http_status = http_status
        self.retryable = code == "provider_response_invalid"
        super().__init__(f"{self.accession}: {self.MESSAGES[code]}")
