# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""BioStudies availability must not masquerade as a missing IDF or retry loop."""
import json

import pytest
import requests

from meta_standards_converter.converters.unified.contracts import ConversionItemResult
from meta_standards_converter.provider_issues import collect_provider_issues

from tests.support.ae import ACCESSION, INFO, fetcher, response


IDF = {"Name": ACCESSION + ".idf.txt", "Type": "IDF File"}
SDRF = {"Name": ACCESSION + ".sdrf.txt", "Type": "SDRF File"}


def test_explicit_denial_checks_access_first_and_never_retries(tmp_path):
    source, calls = fetcher(tmp_path, lambda *a: response({"errorMessage": "Study not accessible", "private": "DO_NOT_DISPLAY"}, 403))
    item = ConversionItemResult("one", ACCESSION)
    with collect_provider_issues(item), pytest.raises(ValueError) as caught:
        source.resolve(ACCESSION)
    assert caught.value.code == "source_inaccessible"
    assert len(calls) == 1
    assert calls[0][0].endswith("/studies/" + ACCESSION + "/info")
    issue = item.to_dict()["diagnostics"][0]["provider_issue"]
    assert issue["error"]["provider"] == "biostudies"
    assert issue["error"]["http_status"] == 403
    assert issue["error"]["retry_category"] == "terminal"
    assert "DO_NOT_DISPLAY" not in str(caught.value) + json.dumps(item.to_dict())
    assert "IDF" in str(caught.value)


@pytest.mark.parametrize("status", [403, 429, 500])
def test_other_denials_and_temporary_failures_keep_bounded_retries(tmp_path, status):
    source, calls = fetcher(tmp_path, lambda *a: response({"errorMessage": "temporary gateway problem"}, status))
    with pytest.raises(requests.HTTPError):
        source.resolve(ACCESSION)
    assert len(calls) == 2


def test_not_found_has_specific_terminal_diagnostic(tmp_path):
    source, calls = fetcher(tmp_path, lambda *a: response({}, 404))
    with pytest.raises(ValueError) as caught:
        source.resolve(ACCESSION)
    assert caught.value.code == "source_not_found"
    assert caught.value.http_status == 404
    assert len(calls) == 1


@pytest.mark.parametrize("listing", [{}, {"data": None}, {"data": [None]}, {"data": [], "recordsFiltered": "invalid"}])
def test_invalid_listing_is_not_evidence_of_missing_files(tmp_path, listing):
    source, _ = fetcher(tmp_path, lambda url, _: response(INFO if url.endswith("/info") else listing))
    item = ConversionItemResult("one", ACCESSION)
    with collect_provider_issues(item), pytest.raises(ValueError) as caught:
        source.resolve(ACCESSION)
    assert caught.value.code == "provider_response_invalid"
    issue = item.to_dict()["diagnostics"][0]["provider_issue"]
    assert issue["reason"] == "provider_response_invalid"
    assert issue["error"]["retry_category"] == "retryable"


@pytest.mark.parametrize("rows,code", [([], "missing_idf"), ([IDF], "missing_sdrf"), ([IDF, dict(IDF, Name="second.idf.txt"), SDRF], "ambiguous_idf")])
def test_valid_listings_report_missing_or_ambiguous_files(tmp_path, rows, code):
    source, calls = fetcher(tmp_path, lambda url, _: response(INFO if url.endswith("/info") else {"data": rows}))
    with pytest.raises(ValueError) as caught:
        source.resolve(ACCESSION)
    assert caught.value.code == code
    assert len(calls) == 2
    assert calls[0][0].endswith("/info")


def test_invalid_info_never_proceeds_to_file_discovery(tmp_path):
    source, calls = fetcher(tmp_path, lambda *a: response({}))
    with pytest.raises(ValueError) as caught:
        source.resolve(ACCESSION)
    assert caught.value.code == "provider_response_invalid"
    assert len(calls) == 1


def test_truncated_pagination_is_a_provider_failure(tmp_path):
    pages = iter([{"data": [IDF], "recordsFiltered": 2}, {"data": [], "recordsFiltered": 2}])
    source, _ = fetcher(tmp_path, lambda url, _: response(INFO if url.endswith("/info") else next(pages)))
    with pytest.raises(ValueError) as caught:
        source.resolve(ACCESSION)
    assert caught.value.code == "provider_response_invalid"


def test_malformed_json_has_actionable_safe_diagnostic(tmp_path):
    source, _ = fetcher(tmp_path, lambda *args: response(b'{PRIVATE_RESPONSE'))
    with pytest.raises(ValueError) as caught:
        source.resolve(ACCESSION)
    assert caught.value.code == 'provider_response_invalid'
    assert 'PRIVATE_RESPONSE' not in str(caught.value)


def test_duplicate_sdrf_is_ambiguous(tmp_path):
    source, _ = fetcher(tmp_path, lambda url, _: response(INFO if url.endswith('/info') else {'data': [IDF, SDRF, SDRF]}))
    with pytest.raises(ValueError) as caught:
        source.resolve(ACCESSION)
    assert caught.value.code == 'ambiguous_sdrf'
