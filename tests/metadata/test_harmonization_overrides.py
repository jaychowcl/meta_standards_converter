# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
import json
import csv
import inspect
import anndata
import pandas
from scipy import sparse

from meta_standards_converter.metadata.harmonization_overrides import (
    resolve_harmonization_overrides,
)
from meta_standards_converter.converters.json2h5ad import JSON2H5ADConverter
from meta_standards_converter.converters.json2ae import JSON2AEConverter
from meta_standards_converter.converters.json2tsv import JSON2TSVConverter
from meta_standards_converter.converters.json2obs import JSON2OBSConverter
from meta_standards_converter.sources.json import JSONPackageSource
from meta_standards_converter.cli import json2ae as json2ae_cli
from meta_standards_converter.cli import json2h5ad as json2h5ad_cli
from meta_standards_converter.cli import json2obs as json2obs_cli
from meta_standards_converter.cli import json2tsv as json2tsv_cli


PROFILE = {
    "schema_version": "1.0",
    "replacements": {
        "organism": ["species_name"],
        "disease": ["sample_disease_name", "disease_category"],
        "organism_part": ["tissue_name", "high_level_tissue"],
    },
}


def package():
    return {
        "miniml_schema_version": "2.0",
        "source": {"format": "test"},
        "database": [],
        "organization": [],
        "contributor": [],
        "platform": [],
        "series": {"accession": [{"value": "GSE1"}]},
        "sample": [{
            "iid": "GSM1",
            "channel": [{
                "organism": [{"value": "human", "taxid": "9606"}],
                "characteristics": [
                    {"name": "organism", "value": "human", "annotations": [{
                        "field": "species_name", "value": "Homo sapiens",
                        "term_source_ref": "ncbitaxon",
                        "term_accession_number": "NCBITaxon:9606",
                        "hierarchy_depth": 0,
                    }]},
                    {"name": "disease", "value": "raw case", "annotations": [{
                        "field": "disease_category", "value": "fallback disease",
                        "term_source_ref": "mondo",
                        "term_accession_number": "MONDO:9",
                    }]},
                    {"name": "tissue", "value": "raw lung", "annotations": [
                        {"field": "tissue_name", "value": "lung", "term_source_ref": "uberon", "term_accession_number": "UBERON:0002048", "hierarchy_depth": 0},
                        {"field": "high_level_tissue", "value": "respiratory system", "term_source_ref": "uberon", "term_accession_number": "UBERON:0001004", "hierarchy_depth": 1},
                    ]},
                ],
            }],
        }],
    }


def package_with_magetab_parameter():
    value = package()
    value["series"]["protocols"] = [{"name": "treatment"}]
    value["series"]["assay_paths"] = [{
        "document": "study.sdrf.txt",
        "steps": [
            {"kind": "sample", "name": "GSM1", "sample_ref": "GSM1"},
            {"kind": "protocol_application", "protocol_ref": "treatment", "parameter_values": [{
                "name": "duration", "value": "30",
                "unit": {"value": "minutes", "annotations": [{
                    "field": "unit", "value": "minute",
                    "term_source_ref": "uo",
                    "term_accession_number": "UO:0000031",
                }]},
            }]},
        ],
    }]
    return value


def package_with_fibrosis_ontology_annotations():
    value = package()
    value["sample"][0]["channel"][0]["characteristics"].extend(
        [
            {"name": "exposure", "value": "bleomycin injection", "annotations": [{"field": "exposure_name", "value": "exposure to bleomycin via injection", "term_source_ref": "ecto", "term_accession_number": "ECTO:0900222"}]},
            {"name": "cell state", "value": "Fbl_24", "annotations": [{"field": "cell_state_name", "value": "Fbl_24", "term_source_ref": "pcl", "term_accession_number": "PCL:0015251"}]},
        ]
    )
    return value


def test_agentic_envelope_is_loaded_with_profile_per_group(tmp_path):
    source = tmp_path / "agentic.json"
    source.write_text(json.dumps({
        "miniml_json": package(),
        "harmonization_overrides": PROFILE,
        "harmonization_targets": [],
    }), encoding="utf-8")

    loaded = JSONPackageSource().load(source)

    assert loaded.groups[0].harmonization_overrides == PROFILE
    assert loaded.groups[0].packages[0]["sample"][0]["iid"] == "GSM1"


def test_resolver_replaces_destinations_from_typed_annotations():
    original = package()
    result = resolve_harmonization_overrides([original], PROFILE, enabled=True)
    channel = result.packages[0]["sample"][0]["channel"][0]
    characteristics = {row["name"]: row for row in channel["characteristics"]}

    assert channel["organism"] == [{
        "value": "Homo sapiens", "taxid": "9606",
        "term_source_ref": "ncbitaxon", "term_accession_number": "NCBITaxon:9606",
    }]
    assert characteristics["disease"]["value"] == "fallback disease"
    assert characteristics["organism part"]["value"] == "lung"
    assert characteristics["organism part"]["term_accession_number"] == "UBERON:0002048"
    assert any(
        row.get("name") == "hz_high_level_tissue"
        for row in channel["characteristics"]
    )
    assert original["sample"][0]["channel"][0]["organism"][0]["value"] == "human"
    assert result.selections[0].source_field == "species_name"


def test_invalid_profile_warns_and_uses_unchanged_raw_packages():
    original = package()
    result = resolve_harmonization_overrides(
        [original], {"schema_version": "9", "replacements": {}}, enabled=True
    )

    assert result.packages == (original,)
    assert result.warnings and result.applied is False


def test_resolved_view_drives_canonical_metadata_and_preserves_typed_annotations():
    result = resolve_harmonization_overrides([package()], PROFILE, enabled=True)
    package_view = result.packages[0]
    sample = package_view["sample"][0]

    metadata = JSON2H5ADConverter().normalizer._sample_metadata_values(sample, package_view)

    assert metadata["organism"] == ("Homo sapiens",)
    assert metadata["disease"] == ("fallback disease",)
    assert metadata["organism_part"] == ("lung",)
    assert metadata["characteristics"]["harmonized_tissue_name"] == ("lung",)
    assert metadata["characteristics"]["harmonized_species_name"] == ("Homo sapiens",)


def test_tabular_opt_in_uses_resolved_destinations_and_typed_annotation_columns(tmp_path):
    source = tmp_path / "agentic.json"
    destination = tmp_path / "manifest.tsv"
    source.write_text(json.dumps({
        "miniml_json": package(), "harmonization_overrides": PROFILE
    }), encoding="utf-8")

    JSON2TSVConverter().convert_source(
        source, destination, use_harmonization_overrides=True
    )
    with destination.open(encoding="utf-8", newline="") as stream:
        row = next(csv.DictReader(stream, delimiter="\t"))

    assert row["msc.sample.channel.organism.value"] == "Homo sapiens"
    assert row["msc.sample.channel.disease"] == "fallback disease"
    assert row["msc.characteristics.harmonized_tissue_name"] == "lung"
    assert row["msc.harmonization.organism.source_field"] == "species_name"


def test_tabular_projects_native_assay_parameter_columns(tmp_path):
    source = tmp_path / "agentic.json"
    destination = tmp_path / "manifest.tsv"
    source.write_text(json.dumps(package_with_magetab_parameter()), encoding="utf-8")

    JSON2TSVConverter().convert_source(source, destination)
    with destination.open(encoding="utf-8", newline="") as stream:
        row = next(csv.DictReader(stream, delimiter="\t"))

    assert row["msc.assay.parameter.duration.value"] == "30"
    assert row["msc.assay.parameter.duration.unit"] == "minutes"
    assert row["msc.assay.parameter.duration.harmonized_unit"] == "minute"
    assert row["msc.assay.parameter.duration.harmonized_unit_id"] == "UO:0000031"


def test_tabular_retains_ecto_and_pcl_harmonization_columns(tmp_path):
    source = tmp_path / "fibrosis.json"
    destination = tmp_path / "manifest.tsv"
    source.write_text(
        json.dumps(package_with_fibrosis_ontology_annotations()), encoding="utf-8"
    )

    JSON2TSVConverter().convert_source(source, destination)
    with destination.open(encoding="utf-8", newline="") as stream:
        row = next(csv.DictReader(stream, delimiter="\t"))

    assert row["msc.characteristics.harmonized_exposure_name_id"] == "ECTO:0900222"
    assert row["msc.characteristics.harmonized_cell_state_name_id"] == "PCL:0015251"


def test_magetab_opt_in_replaces_destinations_and_keeps_typed_annotations_internal(tmp_path):
    source = tmp_path / "agentic.json"
    source.write_text(json.dumps({
        "miniml_json": package(), "harmonization_overrides": PROFILE
    }), encoding="utf-8")

    magetab = JSON2AEConverter().convert(
        str(source), enrich=False, use_harmonization_overrides=True
    )[0]
    sdrf = next(row[1] for row in magetab if row and row[0] == "SDRF File")
    header, values = sdrf[0], sdrf[1]

    disease = header.index("Characteristics[disease]")
    assert values[disease] == "fallback disease"
    assert header[disease + 1:disease + 3] == ["Term Source REF", "Term Accession Number"]
    assert values[disease + 1:disease + 3] == ["mondo", "MONDO:9"]
    assert not any("hz_" in str(label) for label in header)


def test_all_json_consumers_expose_explicit_opt_in():
    assert "use_harmonization_overrides" in inspect.signature(JSON2AEConverter.convert).parameters
    assert "use_harmonization_overrides" in inspect.signature(JSON2H5ADConverter.convert).parameters
    assert "use_harmonization_overrides" in inspect.signature(JSON2TSVConverter.convert_source).parameters
    for cli in (json2ae_cli, json2h5ad_cli, json2obs_cli, json2tsv_cli):
        assert any(
            action.dest == "use_harmonization_overrides"
            for action in cli._parser()._actions
        )


def test_h5ad_publishes_resolved_columns_and_raw_and_harmonization_ledgers(tmp_path):
    source = tmp_path / "agentic.json"
    source.write_text(json.dumps({
        "miniml_json": package(), "harmonization_overrides": PROFILE
    }), encoding="utf-8")
    expression = tmp_path / "GSM1.h5ad"
    anndata.AnnData(
        X=sparse.csr_matrix([[1]]),
        obs=pandas.DataFrame(index=["cell-1"]),
        var=pandas.DataFrame(index=["gene-1"]),
    ).write_h5ad(expression)

    result = JSON2H5ADConverter().convert(
        str(source),
        out=str(tmp_path / "out"),
        asset_specs=[f"GSM1={expression}"],
        use_harmonization_overrides=True,
    )
    converted = anndata.read_h5ad(result.sample_h5ads["GSM1"])

    assert converted.obs["msc.sample.channel.organism.value"].iat[0] == "Homo sapiens"
    assert converted.obs["msc.harmonization.organism.source_field"].iat[0] == "species_name"
    assert converted.uns["msc_harmonization"]["schema_version"] == "1.0"
    fields = converted.uns["msc_miniml"]["fields"]
    raw_organism = fields.loc[fields["path"] == "channel[0].organism[0].value", "value"]
    assert raw_organism.tolist() == ["human"]


def test_h5ad_publishes_native_assay_parameter_obs_and_occurrence_ledger(tmp_path):
    source = tmp_path / "agentic.json"
    source.write_text(json.dumps(package_with_magetab_parameter()), encoding="utf-8")
    expression = tmp_path / "GSM1.h5ad"
    anndata.AnnData(
        X=sparse.csr_matrix([[1]]),
        obs=pandas.DataFrame(index=["cell-1"]),
        var=pandas.DataFrame(index=["gene-1"]),
    ).write_h5ad(expression)

    result = JSON2H5ADConverter().convert(
        str(source), out=str(tmp_path / "out"), asset_specs=[f"GSM1={expression}"]
    )
    converted = anndata.read_h5ad(result.sample_h5ads["GSM1"])

    assert converted.obs["msc.assay.parameter.duration.value"].iat[0] == "30"
    assert converted.obs["msc.assay.parameter.duration.harmonized_unit"].iat[0] == "minute"
    ledger = converted.uns["msc_assay"]["parameters"]
    assert ledger["document"].tolist() == ["study.sdrf.txt"]
    assert ledger["harmonized_unit_id"].tolist() == ["UO:0000031"]


def test_h5ad_and_json2obs_retain_ecto_and_pcl_columns(tmp_path):
    source = tmp_path / "fibrosis.json"
    source.write_text(
        json.dumps(package_with_fibrosis_ontology_annotations()), encoding="utf-8"
    )
    expression = tmp_path / "GSM1.h5ad"
    anndata.AnnData(
        X=sparse.csr_matrix([[1]]),
        obs=pandas.DataFrame(index=["cell-1"]),
        var=pandas.DataFrame(index=["gene-1"]),
    ).write_h5ad(expression)

    result = JSON2OBSConverter().convert(
        source,
        outdir=tmp_path / "obs",
        asset_specs=[f"GSM1={expression}"],
    )

    assert result.obs["msc.characteristics.harmonized_exposure_name_id"].iat[0] == (
        "ECTO:0900222"
    )
    assert result.obs["msc.characteristics.harmonized_cell_state_name_id"].iat[0] == (
        "PCL:0015251"
    )
