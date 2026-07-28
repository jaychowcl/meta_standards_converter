# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
import json

from meta_standards_converter.converters.json_source import JSONPackageSource


def package(study: str, sample: str) -> dict:
    return {
        "series": {"accession": [{"value": study}]},
        "sample": [{"iid": sample}],
    }


def test_miniml_packages_are_grouped_by_study(tmp_path):
    source = tmp_path / "miniml.json"
    source.write_text(
        json.dumps([package("GSE2", "GSM2"), package("GSE1", "GSM1")]),
        encoding="utf-8",
    )

    result = JSONPackageSource().load(source)

    assert [group.dataset_id for group in result.groups] == ["GSE2", "GSE1"]
    assert result.groups[0].packages[0]["sample"][0]["iid"] == "GSM2"
    assert result.warnings == ()


def test_atlas_envelope_yields_only_completed_harmonized_metadata(tmp_path):
    source = tmp_path / "atlas.json"
    source.write_text(
        json.dumps(
            {
                "accessions": [
                    {
                        "datalink_id": "GSE1",
                        "ontology_harmonization_run_status": "completed",
                        "accession_metadata": [package("GSE1", "GSM1")],
                    },
                    {
                        "datalink_id": "GSE2",
                        "ontology_harmonization_run_status": "error",
                        "ontology_harmonization_error": "lookup failed",
                        "accession_metadata": [package("GSE2", "GSM2")],
                    },
                    {
                        "datalink_id": "GSE3",
                        "ontology_harmonization_run_status": "not_run",
                        "accession_metadata": None,
                    },
                ]
            }
        ),
        encoding="utf-8",
    )

    result = JSONPackageSource().load(source)

    assert [group.dataset_id for group in result.groups] == ["GSE1"]
    assert result.groups[0].source_accession == "GSE1"
    assert result.warnings == (
        "GSE2: harmonization error: lookup failed",
        "GSE3: harmonization status not_run; no convertible metadata",
    )


def test_conflicting_duplicate_samples_are_rejected(tmp_path):
    first = package("GSE1", "GSM1")
    second = package("GSE1", "GSM1")
    second["sample"][0]["title"] = "different"
    source = tmp_path / "atlas.json"
    source.write_text(
        json.dumps(
            {
                "accessions": [
                    {
                        "datalink_id": "GSE1",
                        "ontology_harmonization_run_status": "completed",
                        "accession_metadata": [first],
                    },
                    {
                        "datalink_id": "GSE1-copy",
                        "ontology_harmonization_run_status": "completed",
                        "accession_metadata": [second],
                    },
                ]
            }
        ),
        encoding="utf-8",
    )

    try:
        JSONPackageSource().load(source)
    except ValueError as error:
        assert "GSM1" in str(error)
        assert "conflicting" in str(error)
    else:
        raise AssertionError("expected conflicting duplicate sample failure")
