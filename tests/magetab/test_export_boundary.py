# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Construction must never supplement supplied metadata through a provider."""
from unittest.mock import Mock

from meta_standards_converter.magetab.constructor import AEConstructor
from meta_standards_converter.sources.pubmed import PubmedWebFetcher


def test_constructor_does_not_retrieve_missing_publication(monkeypatch):
    fetch = Mock(return_value=(None,) * 6)
    monkeypatch.setattr(PubmedWebFetcher, "pubmed_summary", fetch)
    package = {
        "miniml_schema_version": "3.1",
        "source": {"format": "test"},
        "series": {"iid": "GSE1", "title": "Supplied metadata", "pubmed_id": ["123"]},
        "sample": [{"iid": "GSM1"}],
        "platform": [],
    }
    result = AEConstructor().miniml2magetab(package)
    assert result
    fetch.assert_not_called()


def test_preparation_and_rendering_do_not_construct_provider_clients(monkeypatch):
    from meta_standards_converter.sources.insdc import INSDCWebfetcher
    from tests.support.archive import legacy

    def forbidden(*args, **kwargs):
        raise AssertionError("provider client created inside export boundary")

    monkeypatch.setattr(PubmedWebFetcher, "__init__", forbidden)
    monkeypatch.setattr(INSDCWebfetcher, "__init__", forbidden)
    constructor = AEConstructor()
    plan = constructor.prepare(legacy())
    assert constructor.render(plan)
