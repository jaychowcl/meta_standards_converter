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

from meta_standards_converter.metadata.projection.tabular import TabularMetadataProjection
from meta_standards_converter.converters.json2tsv import JSON2TSVConverter
from meta_standards_converter.metadata.projection.tabular import TabularProjectionError


def package(study: str = "GSE1", sample: str = "GSM1") -> dict:
    return {
        "miniml_schema_version": "2.0",
        "source": {"format": "test"},
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
                            {"name": "tissue", "value": "brain"},
                        ],
                    }
                ],
            }
        ],
    }


def atlas_v1_dataset(dataset_id: str, metadata: dict) -> dict:
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


def atlas_v1_payload(datasets: list[dict]) -> dict:
    return {
        "schema_version": "1.0",
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


def test_tsv_projects_sample_bound_typed_protocol_and_material_semantics(tmp_path):
    payload = package()
    payload["series"].update({
        "protocols": [
            {
                "name": "P-collect",
                "type": {
                    "value": "bespoke collection protocol",
                    "term_source_ref": "TEST",
                    "term_accession_number": "TEST:001",
                },
            }
        ],
        "assay_paths": [
            {
                "steps": [
                    {
                        "kind": "sample",
                        "name": "GSM1",
                        "sample_ref": "GSM1",
                        "material_type": {"value": "fresh tissue specimen"},
                    },
                    {
                        "kind": "protocol_application",
                        "protocol_ref": "P-collect",
                    },
                ]
            }
        ],
    })
    source = tmp_path / "input.json"
    output = tmp_path / "output.tsv"
    source.write_text(json.dumps([payload]), encoding="utf-8")

    JSON2TSVConverter().convert_source(source, output)
    _columns, rows = read_rows(output, "\t")

    assert rows[0]["msc.sample.channel.material_type"] == "fresh tissue specimen"
    assert rows[0]["msc.protocol.types"] == "bespoke collection protocol"
    assert rows[0]["msc.protocol.term_source_refs"] == "TEST"
    assert rows[0]["msc.protocol.term_accession_numbers"] == "TEST:001"


def test_tabular_converter_uses_injected_neutral_metadata_service(tmp_path):
    class MetadataService:
        def study_accession(self, packages):
            assert len(packages) == 1
            return "GSE-SERVICE"

        def samples(self, package):
            return tuple(package.get("sample", ()))

        def sample_accession(self, sample):
            return "GSM-SERVICE"

        def sample_metadata(self, sample, package):
            return {"title": "from-service"}

        def sample_modality(self, sample):
            return "single_cell"

    source = tmp_path / "input.json"
    output = tmp_path / "output.tsv"
    source.write_text(json.dumps([package()]), encoding="utf-8")

    JSON2TSVConverter(metadata_service=MetadataService()).convert_source(
        source, output
    )
    _columns, rows = read_rows(output, "\t")

    assert rows[0]["msc.sample.accession"] == "GSM-SERVICE"
    assert rows[0]["msc.series.accession"] == "GSE-SERVICE"
    assert rows[0]["msc.sample.title"] == "from-service"
    assert rows[0]["msc.expression.modality"] == "single_cell"


def test_atlas_json_is_aggregated_into_one_csv(tmp_path):
    source = tmp_path / "atlas.json"
    output = tmp_path / "output.csv"
    source.write_text(
        json.dumps(
            atlas_v1_payload(
                [
                    atlas_v1_dataset("GSE1", package("GSE1", "GSM1")),
                    atlas_v1_dataset("GSE2", package("GSE2", "GSM2")),
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


@pytest.mark.parametrize('output_format,delimiter', [('tsv', '\t'), ('csv', ',')])
@pytest.mark.parametrize('virtual', [True, False])
def test_manifest_preserves_injected_source_and_metadata(tmp_path, output_format, delimiter, virtual):
    from meta_standards_converter.sources.json import JSONPackageSource
    from meta_standards_converter.metadata.interpretation import MINiMLMetadataService

    real = tmp_path / 'real.json'
    real.write_text(json.dumps(package()))
    source = tmp_path / 'virtual.json' if virtual else real
    calls = []

    class Source:
        def load(self, path):
            calls.append(path)
            return JSONPackageSource().load(real)

    class Metadata(MINiMLMetadataService):
        def sample_metadata(self, sample, package):
            return {**super().sample_metadata(sample, package), 'title': 'injected\t"title",\nsecond line'}

    converter = JSON2TSVConverter(package_source=Source(), metadata_service=Metadata(), output_format=output_format)
    direct = tmp_path / ('direct.' + output_format)
    converter.convert_source(source, direct)
    result = converter.export_manifest(source, outdir=tmp_path / 'bundle', output_format=output_format)
    from pathlib import Path
    table = Path(result.output_path)
    assert table.read_bytes() == direct.read_bytes()
    assert read_rows(table, delimiter)[1][0]['msc.sample.title'] == 'injected\t"title",\nsecond line'
    assert calls == [source, source]
    manifest = json.loads(Path(result.manifest_path).read_text())
    assert manifest['row_count'] == 1
    assert tuple(manifest['columns']) == read_rows(table, delimiter)[0]
    with pytest.raises(FileExistsError):
        converter.export_manifest(source, outdir=tmp_path / 'bundle', output_format=output_format)
    converter.export_manifest(source, outdir=tmp_path / 'bundle', output_format=output_format, overwrite=True)
    assert table.read_bytes() == direct.read_bytes()


def test_manifest_injected_projector_validation_is_atomic(tmp_path):
    class Projector:
        def project_sample(self, *, context):
            return TabularMetadataProjection(values={'private.sample': context.sample_accession},
                                             columns=('private.sample',), errors=('review required',))
    source = tmp_path / 'input.json'
    source.write_text(json.dumps(package()))
    converter = JSON2TSVConverter(metadata_projectors=[Projector()])
    out = tmp_path / 'bundle'
    with pytest.raises(TabularProjectionError, match='review required'):
        converter.export_manifest(source, outdir=out)
    assert not out.exists()
    result = converter.export_manifest(source, outdir=out, allow_invalid=True)
    from pathlib import Path
    assert read_rows(Path(result.output_path), '\t') == (('private.sample',), [{'private.sample': 'GSM1'}])
    assert json.loads(Path(result.manifest_path).read_text())['errors'] == list(result.errors)
    original = Path(result.output_path).read_bytes()
    with pytest.raises(TabularProjectionError):
        converter.export_manifest(source, outdir=out, overwrite=True)
    assert Path(result.output_path).read_bytes() == original
