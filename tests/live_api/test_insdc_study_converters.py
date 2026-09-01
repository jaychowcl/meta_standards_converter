# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Opt-in end-to-end contracts for study-scoped INSDC retrieval."""

from __future__ import annotations

import pytest

from meta_standards_converter.insdc_handlers.study_fetchers import (
    ENAStudyFetcher,
    SRAStudyFetcher,
    StudyNotFoundError,
)
from meta_standards_converter.insdc_handlers.study_parser import INSDCStudyParser
from meta_standards_converter.miniml import MINiMLCodec


CORE_STUDIES = (
    "SRP002056",
    "ERP106767",
    "DRP000030",
    "SRP250911",
    "ERP000265",
)


@pytest.mark.live_api(max_requests=100)
@pytest.mark.parametrize("fetcher_type", [SRAStudyFetcher, ENAStudyFetcher])
@pytest.mark.parametrize("study", CORE_STUDIES)
def test_core_studies_retain_a_complete_provider_hierarchy(fetcher_type, study) -> None:
    results = fetcher_type().fetch(study)

    assert [result.study_accession for result in results] == [study]
    kinds = {document.kind for document in results[0].documents}
    if fetcher_type is SRAStudyFetcher:
        assert "sra_efetch" in kinds
    else:
        assert {"ena_study", "ena_experiment", "ena_run"} <= kinds


@pytest.mark.live_api(max_requests=20)
@pytest.mark.parametrize("fetcher_type", [SRAStudyFetcher, ENAStudyFetcher])
def test_project_probe_resolves_two_ddbj_studies_in_stable_order(fetcher_type) -> None:
    fetcher = fetcher_type()
    if isinstance(fetcher, ENAStudyFetcher):
        studies = sorted(study for study, _project in fetcher._resolve("PRJDA43743"))
    else:
        studies = [result.study_accession for result in fetcher.fetch("PRJDA43743")]

    assert studies == ["DRP000086", "DRP000087"]


@pytest.mark.live_api(max_requests=20)
def test_zero_run_probe_degrades_in_ena_but_is_not_found_in_sra() -> None:
    ena_results = ENAStudyFetcher().fetch("DRP000158")

    assert [result.study_accession for result in ena_results] == ["DRP000158"]
    assert "ena_study" in {document.kind for document in ena_results[0].documents}
    package = INSDCStudyParser().parse(ena_results[0])
    data = MINiMLCodec().encode(package)
    assert not package.samples
    assert "metadata_only_study" in {
        warning["code"] for warning in data["extensions"]["insdc"]["warnings"]
    }
    with pytest.raises(StudyNotFoundError):
        SRAStudyFetcher().fetch("DRP000158")
