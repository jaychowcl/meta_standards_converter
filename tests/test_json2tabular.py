# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
import csv
import json

import pytest

from meta_standards_converter.converters import (
    TabularMetadataProjection,
)
from meta_standards_converter.converters.json2tabular import (
    JSON2TSVConverter,
    TabularProjectionError,
)


def package(study: str = "GSE1", sample: str = "GSM1") -> dict:
    return {
        "database": [{"public_id": "GEO", "name": "Gene Expression Omnibus"}],
        "series": {"accession": [{"value": study}]},
        "sample": [
            {
                "iid": sample,
                "title": "sample title",
                "channel": [
                    {
                        "organism": [{"value": "Mus musculus", "taxid": "10090"}],
                        "characteristics": [
                            {"tag": "tissue", "value": "brain"},
                        ],
                    }
                ],
            }
        ],
    }


def atlas_v2_dataset(dataset_id: str, metadata: dict) -> dict:
    return {
        "dataset_id": dataset_id,
        "source_repository": "geo",
        "source_ordinal": 0,
        "status": "harmonized",
        "metadata": metadata,
        "publication_ids": [],
        "review": None,
        "harmonization": None,
        "diagnostics": [],
    }


def atlas_v2_payload(datasets: list[dict]) -> dict:
    return {
        "schema_version": "2.0",
        "atlas": {"atlas_id": "atlas-test", "title": "Test", "theme": "test"},
        "run": {
            "run_id": "run-test",
            "created_at": "2026-07-31T00:00:00Z",
            "config": {"queries": [], "metadata_repositories": [], "options": {}},
            "status": "complete",
        },
        "datasets": datasets,
        "publications": [],
        "summary": {
            "dataset_count": len(datasets),
            "publication_count": 0,
            "completed_dataset_count": len(datasets),
            "failed_dataset_count": 0,
        },
    }


def read_rows(path, delimiter):
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter=delimiter)
        return tuple(reader.fieldnames or ()), list(reader)


def test_default_tsv_emits_stable_msc_metadata_and_characteristics(tmp_path):
    source = tmp_path / "input.json"
    output = tmp_path / "output.tsv"
    source.write_text(json.dumps([package()]), encoding="utf-8")

    result = JSON2TSVConverter().convert_source(source, output)
    columns, rows = read_rows(output, "\t")

    assert result.row_count == 1
    assert columns[:4] == (
        "msc.sample.accession",
        "msc.series.accession",
        "msc.sample.title",
        "msc.sample.description",
    )
    assert rows[0]["msc.sample.accession"] == "GSM1"
    assert rows[0]["msc.sample.channel.organism.value"] == "Mus musculus"
    assert rows[0]["msc.sample.channel.organism.taxid"] == "10090"
    assert rows[0]["msc.characteristics.tissue"] == "brain"


def test_atlas_json_is_aggregated_into_one_csv(tmp_path):
    source = tmp_path / "atlas.json"
    output = tmp_path / "output.csv"
    source.write_text(
        json.dumps(
            atlas_v2_payload(
                [
                    atlas_v2_dataset("GSE1", package("GSE1", "GSM1")),
                    atlas_v2_dataset("GSE2", package("GSE2", "GSM2")),
                ]
            )
        ),
        encoding="utf-8",
    )

    result = JSON2TSVConverter(output_format="csv").convert_source(source, output)
    _columns, rows = read_rows(output, ",")

    assert result.row_count == 2
    assert [row["msc.sample.accession"] for row in rows] == ["GSM1", "GSM2"]
    assert result.dataset_ids == ("GSE1", "GSE2")


def test_manifest_converter_rejects_unknown_output_format():
    with pytest.raises(ValueError, match="output_format must be 'tsv' or 'csv'"):
        JSON2TSVConverter(output_format="json")


def test_explicit_projector_replaces_default_contract(tmp_path):
    class Projector:
        def project_sample(self, *, context):
            return TabularMetadataProjection(
                values={
                    "private.sample": context.sample_accession,
                    "private.extra": "value",
                },
                columns=("private.sample",),
            )

    source = tmp_path / "input.json"
    output = tmp_path / "output.tsv"
    source.write_text(json.dumps([package()]), encoding="utf-8")

    JSON2TSVConverter(metadata_projectors=[Projector()]).convert_source(
        source, output
    )
    columns, rows = read_rows(output, "\t")

    assert columns == ("private.sample", "private.extra")
    assert rows == [{"private.sample": "GSM1", "private.extra": "value"}]


def test_projector_column_collisions_fail(tmp_path):
    class Projector:
        def __init__(self, value):
            self.value = value

        def project_sample(self, *, context):
            return TabularMetadataProjection(values={"duplicate": self.value})

    source = tmp_path / "input.json"
    source.write_text(json.dumps([package()]), encoding="utf-8")

    with pytest.raises(ValueError, match="duplicate"):
        JSON2TSVConverter(
            metadata_projectors=[Projector("one"), Projector("two")]
        ).convert_source(source, tmp_path / "output.tsv")


def test_projection_errors_fail_closed_unless_allowed(tmp_path):
    class Projector:
        def project_sample(self, *, context):
            return TabularMetadataProjection(
                values={"sample": context.sample_accession},
                errors=("required value missing",),
            )

    source = tmp_path / "input.json"
    source.write_text(json.dumps([package()]), encoding="utf-8")
    converter = JSON2TSVConverter(metadata_projectors=[Projector()])

    with pytest.raises(TabularProjectionError, match="required value missing"):
        converter.convert_source(source, tmp_path / "strict.tsv")

    result = converter.convert_source(
        source, tmp_path / "inspection.tsv", allow_invalid=True
    )
    assert result.errors == ("GSM1: required value missing",)
