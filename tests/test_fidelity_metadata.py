"""Conservative biological projection and repository-scoped date contracts."""
from copy import deepcopy

import pytest

from meta_standards_converter.metadata.interpretation import MINiMLMetadataService
from meta_standards_converter.magetab.idf import IDFConstructor
from meta_standards_converter.magetab.parser import AEParser
from tests.converters.test_json2tsv import package


@pytest.mark.parametrize("names", [("OrganismPart", "DevelopmentalStage"),
                                    ("organism_part", "developmental_stage"),
                                    ("organism part", "developmental stage")])
def test_biological_aliases_preserve_raw_names_and_values(names):
    value = package()
    sample = value["sample"][0]
    sample["channel"][0]["characteristics"] = [
        {"name": names[0], "value": "whole_organism"},
        {"name": names[1], "value": "embryo"}]
    before = deepcopy(value)
    service = MINiMLMetadataService()
    result = service.sample_metadata_values(sample, value)
    assert result["organism_part"] == ("whole_organism",)
    assert result["developmental_stage"] == ("embryo",)
    assert service.metadata_slug(names[0]) in result["characteristics"]
    assert value == before


def test_tissue_uses_all_explicit_values_and_never_source_description():
    value = package()
    sample = value["sample"][0]
    channel = sample["channel"][0]
    channel.update(source="total RNA", characteristics=[])
    service = MINiMLMetadataService()
    result = service.sample_metadata_values(sample, value)
    assert result["source"] == ("total RNA",)
    assert result["organism_part"] == ()
    channel["characteristics"] = [{"name": "tissue", "value": "lung"},
                                  {"name": "OrganismPart", "value": "heart"},
                                  {"name": "organism part", "value": "lung"}]
    assert service.sample_metadata_values(sample, value)["organism_part"] == ("lung", "heart")


@pytest.mark.parametrize("fmt,accession,provider", [
    ("GEO MINiML", "GSE1", "GEO"), ("MAGE-TAB", "E-MTAB-1", "ArrayExpress"),
    ("ArrayExpress", "custom", "ArrayExpress"), ("MAGE-TAB", "custom", None),
    ("ENA", "ERP1", "ENA"), ("SRA", "SRP1", "SRA")])
def test_date_labels_follow_repository_provenance(fmt, accession, provider, monkeypatch):
    constructor = IDFConstructor()
    monkeypatch.setattr(constructor, "_current_idf_date", lambda: "2026-09-21")
    value = package(accession)
    value["source"]["format"] = fmt
    value["series"]["status"] = [{"submission_date": "2025-01-01", "release_date": "2025-02-01"}]
    rows = dict((r[0], r[1:]) for r in constructor._idf_dates(value))
    assert rows["Comment[ArrayExpressSubmissionDate]"] == ["2026-09-21"]
    assert rows["Public Release Date"] == ["2025-02-01"]
    if fmt == "GEO MINiML": assert rows["Date of Experiment"] == ["2025-01-01"]
    if provider:
        assert rows[f"Comment[{provider}ReleaseDate]"] == ["2025-02-01"]
    if provider != "GEO": assert "Comment[GEOReleaseDate]" not in rows


def test_explicit_date_scopes_and_conflicting_values_survive_import_and_export():
    parser = AEParser()
    statuses = parser._statuses(parser._idf_index([
        ["Public Release Date", "2025-01-01"],
        ["Comment[GEOReleaseDate]", "2025-02-01"],
        ["Comment[ArrayExpressReleaseDate]", "2025-03-01"],
        ["Comment[ArrayExpressLastUpdateDate]", "2025-04-01"],
    ]))
    assert any(s.get("database") == "ArrayExpress" and s.get("release_date") == "2025-03-01" for s in statuses)
    assert any(s.get("database") == "GEO" and s.get("release_date") == "2025-02-01" for s in statuses)
    value = package("E-MTAB-1")
    value["source"]["format"] = "MAGE-TAB"
    value["series"]["status"] = statuses
    rows = dict((r[0], r[1:]) for r in IDFConstructor()._idf_dates(value))
    assert "2025-02-01" in rows["Comment[GEOReleaseDate]"]
    assert "2025-03-01" in rows["Comment[ArrayExpressReleaseDate]"]
    assert rows["Public Release Date"] == ["2025-01-01"]


def test_native_dates_accept_database_declarations_without_optional_iid():
    from meta_standards_converter.miniml.archive_dates import normalize_archive_dates
    value = package()
    value['source']['format'] = 'ENA'
    value['database'] = [{'name': 'ENA'}]
    normalize_archive_dates(value)
    assert value['database'][0] == {'name': 'ENA'}
