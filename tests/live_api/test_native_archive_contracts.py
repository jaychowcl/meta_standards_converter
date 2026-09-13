# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Opt-in provider checks; deterministic semantics are tested in test_native_archive_*.

RUN_LIVE_API_TESTS=1 .venv/bin/python -m pytest tests/live_api/test_native_archive_contracts.py
"""
import pytest
from meta_standards_converter.sources.sra import SRASource
from meta_standards_converter.sources.ena import ENASource
from meta_standards_converter.miniml.sra_parser import SRAParser

pytestmark = pytest.mark.live_api


def test_sra_entrez_run_expansion_and_full_small_study():
    source = SRASource()
    resolution = source.resolve('SRR037073')
    assert not resolution.issues
    assert [s.study for s in resolution.studies] == ['SRP002056']
    records = source.fetch(resolution.studies[0])
    assert not records.issues
    package = SRAParser().parse(records)
    assert package.series.iid == 'SRP002056'
    assert len(package.samples) >= 4
    assert all(s.sra_runs for s in package.samples)


def test_ena_run_identity_count_xml_and_analysis_contracts():
    source = ENASource()
    result = source.resolve('SRR037073')
    assert not result.issues
    assert [(s.study, s.primary) for s in result.studies] == [('SRP002056', 'PRJNA123835')]
    rows = source.search('read_run', 'secondary_study_accession="SRP002056"')
    count = source.http.get(source.portal + 'count', {'result': 'read_run', 'query': 'secondary_study_accession="SRP002056"', 'format': 'json'}, 'json')
    assert len(rows) == int(count['count'])
    root = source.xml(['SRR037073'])
    assert root.find('RUN').get('accession') == 'SRR037073'
    root = source.xml(['ERZ10000609'])
    assert root.find('ANALYSIS') is not None
    root = source.xml(['GCA_000001405.29'])
    assert root.find('ASSEMBLY').get('accession') == 'GCA_000001405.29'
