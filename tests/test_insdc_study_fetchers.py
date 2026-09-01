# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from __future__ import annotations

import json

import pytest

from meta_standards_converter.insdc_handlers.study_fetchers import (
    MAX_PAGES_HARD,
    MAX_PROVIDER_REQUESTS_HARD,
    MAX_UNIQUE_RUNS_HARD,
    ENAStudyFetcher,
    PaginationSafetyError,
    ProviderEmergencyLimitError,
    SRAStudyFetcher,
    StudyNotFoundError,
    collect_accession_pages,
)


class Response:
    def __init__(self, content: bytes, status_code: int = 200):
        self.content = content
        self.status_code = status_code
        self.headers = {"Content-Length": str(len(content))}

    def raise_for_status(self):
        if self.status_code >= 400:
            raise RuntimeError(f"HTTP {self.status_code}")

    def iter_content(self, chunk_size):
        yield self.content


class Requester:
    def __init__(self, rows):
        self.rows = rows
        self.calls = []

    def get(self, url, **kwargs):
        self.calls.append(("get", url, kwargs))
        return Response(json.dumps(self.rows).encode())


@pytest.mark.parametrize(
    ("accession", "field"),
    [
        ("PRJNA609050", "GPRJ"),
        ("PRJEB24901", "GPRJ"),
        ("SAMN14218700", "BSPL"),
        ("SRP250911", "ACCN"),
        ("SRS6225446", "ACCN"),
        ("SRX7812918", "ACCN"),
        ("SRR11192680", "ACCN"),
        ("DRP000030", "ACCN"),
    ],
)
def test_sra_accessions_select_exact_esearch_field(accession: str, field: str) -> None:
    assert SRAStudyFetcher.search_field(accession) == field


@pytest.mark.parametrize("accession", ["", "GSE18729", "E-MTAB-6486", "SRZ1", "SRP"])
def test_fetchers_reject_unsupported_accession_classes(accession: str) -> None:
    with pytest.raises(ValueError, match="supported INSDC accession"):
        SRAStudyFetcher.validate_accession(accession)
    with pytest.raises(ValueError, match="supported INSDC accession"):
        ENAStudyFetcher.validate_accession(accession)


def test_page_collection_deduplicates_and_stops_on_short_page() -> None:
    pages = {0: ["A", "B"], 2: ["B", "C"], 4: ["D"]}
    assert collect_accession_pages(lambda offset, limit: pages[offset], page_size=2) == ["A", "B", "C", "D"]


def test_page_collection_rejects_repeated_signature() -> None:
    with pytest.raises(PaginationSafetyError, match="repeated page"):
        collect_accession_pages(lambda offset, limit: ["A", "B"], page_size=2)


def test_page_collection_rejects_no_progress_page() -> None:
    pages = {0: ["A", "B"], 2: ["A", "B"]}
    with pytest.raises(PaginationSafetyError, match="repeated page|no new accessions"):
        collect_accession_pages(lambda offset, limit: pages[offset], page_size=2)


def test_emergency_limits_are_internal_high_water_marks() -> None:
    assert MAX_UNIQUE_RUNS_HARD == 1_000_000
    assert MAX_PAGES_HARD == 100_000
    assert MAX_PROVIDER_REQUESTS_HARD == 100_000


def test_not_found_error_names_provider_and_accession() -> None:
    error = StudyNotFoundError("sra", "DRP000158")
    assert str(error) == "sra provider could not resolve accession DRP000158"


def test_ena_project_resolution_preserves_row_level_project_study_pairs() -> None:
    requester = Requester(
        [
            {"study_accession": "PRJDA43743", "secondary_study_accession": "DRP000087"},
            {"study_accession": "PRJDA43743", "secondary_study_accession": "DRP000086"},
        ]
    )
    fetcher = ENAStudyFetcher(requester=requester)

    resolved = fetcher._resolve("PRJDA43743")

    assert resolved == [
        ("DRP000087", "PRJDA43743"),
        ("DRP000086", "PRJDA43743"),
    ]
    params = requester.calls[0][2]["params"]
    assert params["query"] == 'study_accession="PRJDA43743"'
    assert "project_accession" not in params["fields"]
    assert params["limit"] == 0


def test_ena_experiment_resolution_uses_declared_portal_field() -> None:
    requester = Requester(
        [{"study_accession": "PRJNA609050", "secondary_study_accession": "SRP250911"}]
    )
    fetcher = ENAStudyFetcher(requester=requester)

    assert fetcher._resolve("SRX7812918") == [("SRP250911", "PRJNA609050")]

    assert requester.calls[0][2]["params"]["query"] == 'experiment_accession="SRX7812918"'


def test_explicit_secondary_study_does_not_expand_to_project_siblings() -> None:
    requester = Requester(
        [
            {
                "study_accession": "PRJDA39855",
                "secondary_study_accession": "DRP000158;DRP000031",
            }
        ]
    )
    fetcher = ENAStudyFetcher(requester=requester)

    assert fetcher._resolve("DRP000158") == [("DRP000158", "PRJDA39855")]


def test_ena_browser_batch_recursively_splits_oversize_responses() -> None:
    class Fetcher(ENAStudyFetcher):
        def __init__(self):
            super().__init__(requester=Requester([]))

        def _request(self, method, url, **kwargs):
            accessions = kwargs["json"]["accessions"]
            if len(accessions) > 1:
                from meta_standards_converter.xml_safety import XMLSizeLimitError

                raise XMLSizeLimitError("too large")
            accession = accessions[0]
            return f'<STUDY_SET><STUDY accession="{accession}"/></STUDY_SET>'.encode(), url

    documents = Fetcher()._browser_batch(["DRP000086", "DRP000087"], kind="ena_study")

    assert len(documents) == 2
    assert [document.kind for document in documents] == ["ena_study", "ena_study"]


@pytest.mark.parametrize(
    ("result", "field", "accession"),
    [("analysis", "analysis_accession", "ERZ1"), ("assembly", "assembly_accession", "GCA_1")],
)
def test_ena_study_links_can_preserve_non_first_class_objects(result, field, accession) -> None:
    fetcher = ENAStudyFetcher(requester=Requester([{field: accession}]))

    assert fetcher._links("PRJEB1", result) == [accession]


def test_page_collection_has_a_high_emergency_page_guard() -> None:
    with pytest.raises(ProviderEmergencyLimitError, match="page safety guard"):
        collect_accession_pages(
            lambda offset, limit: [f"A{offset}", f"B{offset}"],
            page_size=2,
            max_pages=1,
        )


def test_ncbi_efetch_history_rejects_a_repeated_package_batch() -> None:
    package = b'<EXPERIMENT_PACKAGE_SET><EXPERIMENT_PACKAGE><EXPERIMENT accession="SRX1"/><RUN_SET><RUN accession="SRR1"/></RUN_SET></EXPERIMENT_PACKAGE></EXPERIMENT_PACKAGE_SET>'

    class Fetcher(SRAStudyFetcher):
        def _get_bytes(self, endpoint, params):
            return package, "https://example.test/efetch"

    fetcher = Fetcher(requester=Requester([]))

    with pytest.raises(PaginationSafetyError, match="repeated NCBI EFetch batch"):
        fetcher._history_documents(
            {"count": "200", "querykey": "1", "webenv": "test"},
            prefix="study",
        )
