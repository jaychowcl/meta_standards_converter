# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Opt-in contracts for MSC public HTTP providers."""

import pytest

from meta_standards_converter.ae_handlers.ae_webfetcher import AEWebFetcher
from meta_standards_converter.helpers.request_helper import RequestSettings
from meta_standards_converter.insdc_handlers.insdc_webfetcher import INSDCWebfetcher
from meta_standards_converter.pubmed_handlers.pubmed_webfetcher import PubmedWebFetcher


pytestmark = pytest.mark.live_api
SETTINGS = RequestSettings(timeout=10, request_delay=0, max_retries=0)


def test_pubmed_esummary_contract():
    doi, authors, title, *_ = PubmedWebFetcher(request_settings=SETTINGS).pubmed_summary("39747812")
    assert doi and authors and title


def test_ncbi_sra_and_ena_file_report_contracts():
    fetcher = INSDCWebfetcher(ncbi_request_settings=SETTINGS, ena_request_settings=SETTINGS)
    runs = fetcher.fetch_sra_runs("SRX32831930")
    report = fetcher.fetch_ena_file_report("SRX32831930")
    assert runs
    assert all(run.get("run") for run in runs)
    assert isinstance(report, list)


def test_biostudies_idf_and_sdrf_contracts():
    resource = AEWebFetcher(request_settings=SETTINGS).resolve("E-MTAB-1")
    assert resource.source_kind == "accession"
    assert resource.idf.text.strip()
    assert resource.sdrfs
    assert all(sdrf.text.strip() for sdrf in resource.sdrfs)
