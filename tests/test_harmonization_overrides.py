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

from meta_standards_converter.converters.harmonization_overrides import (
    resolve_harmonization_overrides,
)
from meta_standards_converter.converters.json2h5ad import JSON2H5ADConverter
from meta_standards_converter.converters.json2ae import json2ae
from meta_standards_converter.converters.json2tabular import JSON2TSVConverter
from meta_standards_converter.converters.json_source import JSONPackageSource
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
        "database": [],
        "organization": [],
        "contributor": [],
        "platform": [],
        "series": {"accession": [{"value": "GSE1"}]},
        "sample": [{
            "iid": "GSM1",
            "channel": [{
                "organism": [{"value": "human", "taxid": "9606"}],
                "hz_species_name": [{
                    "value": "Homo sapiens", "id": "NCBITaxon:9606",
                    "onto": "ncbitaxon", "hierarchy_depth": 0,
                }],
                "characteristics": [
                    {"tag": "disease", "value": "raw case"},
                    {"tag": "tissue", "value": "raw lung"},
                    {"tag": "hz_disease_category", "value": "fallback disease"},
                    {"tag": "hz_disease_category_id", "value": "MONDO:9"},
                    {"tag": "hz_disease_category_onto", "value": "mondo"},
                    {"tag": "hz_tissue_name", "value": "lung"},
                    {"tag": "hz_tissue_name_id", "value": "UBERON:0002048"},
                    {"tag": "hz_tissue_name_onto", "value": "uberon"},
                    {"tag": "hz_tissue_name_hierarchy_depth", "value": 0},
                    {"tag": "hz_high_level_tissue", "value": "respiratory system"},
                    {"tag": "hz_high_level_tissue_id", "value": "UBERON:0001004"},
                    {"tag": "hz_high_level_tissue_onto", "value": "uberon"},
                    {"tag": "hz_high_level_tissue_hierarchy_depth", "value": 1},
                ],
            }],
        }],
    }


def package_with_magetab_parameter():
    value = package()
    value["mage_tab"] = {"model": {
        "schema_version": 1,
        "assay_paths": [{
            "id": "study.sdrf.txt:row:1",
            "sdrf": "study.sdrf.txt",
            "row_index": 1,
            "binding": {"source_name": "GSM1", "sample_name": "GSM1"},
            "steps": [{
                "kind": "attribute",
                "attribute_type": "parameter value",
                "name": "duration",
                "column_index": 3,
                "value": "30",
                "unit": "minutes",
                "hz_unit": "minute",
                "hz_unit_id": "UO:0000031",
                "hz_unit_onto": "uo",
            }],
        }],
    }}
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


def test_resolver_replaces_destinations_and_retains_every_hz_field():
    original = package()
    result = resolve_harmonization_overrides([original], PROFILE, enabled=True)
    channel = result.packages[0]["sample"][0]["channel"][0]
    characteristics = {row["tag"]: row for row in channel["characteristics"]}

    assert channel["organism"] == [{
        "value": "Homo sapiens", "taxid": "9606",
        "term_source_ref": "ncbitaxon", "term_accession_number": "NCBITaxon:9606",
    }]
    assert characteristics["disease"]["value"] == "fallback disease"
    assert characteristics["organism part"]["value"] == "lung"
    assert characteristics["organism part"]["term_accession_number"] == "UBERON:0002048"
    assert characteristics["hz_high_level_tissue"]["value"] == "respiratory system"
    assert original["sample"][0]["channel"][0]["organism"][0]["value"] == "human"
    assert result.selections[0].source_field == "species_name"


def test_invalid_profile_warns_and_uses_unchanged_raw_packages():
    original = package()
    result = resolve_harmonization_overrides(
        [original], {"schema_version": "9", "replacements": {}}, enabled=True
    )

    assert result.packages == (original,)
    assert result.warnings and result.applied is False


def test_resolved_view_drives_canonical_metadata_and_preserves_hz_characteristics():
    result = resolve_harmonization_overrides([package()], PROFILE, enabled=True)
    package_view = result.packages[0]
    sample = package_view["sample"][0]

    metadata = JSON2H5ADConverter()._sample_metadata_values(sample, package_view)

    assert metadata["organism"] == ("Homo sapiens",)
    assert metadata["disease"] == ("fallback disease",)
    assert metadata["organism_part"] == ("lung",)
    assert metadata["characteristics"]["hz_tissue_name"] == ("lung",)
    assert metadata["characteristics"]["hz_species_name"] == ("Homo sapiens",)


def test_tabular_opt_in_uses_resolved_destinations_and_retains_hz_columns(tmp_path):
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
    assert row["msc.characteristics.hz_tissue_name"] == "lung"
    assert row["msc.harmonization.organism.source_field"] == "species_name"


def test_tabular_projects_additive_magetab_parameter_columns(tmp_path):
    source = tmp_path / "agentic.json"
    destination = tmp_path / "manifest.tsv"
    source.write_text(json.dumps(package_with_magetab_parameter()), encoding="utf-8")

    JSON2TSVConverter().convert_source(source, destination)
    with destination.open(encoding="utf-8", newline="") as stream:
        row = next(csv.DictReader(stream, delimiter="\t"))

    assert row["msc.mage_tab.parameter.duration.value"] == "30"
    assert row["msc.mage_tab.parameter.duration.unit"] == "minutes"
    assert row["msc.mage_tab.parameter.duration.hz_unit"] == "minute"
    assert row["msc.mage_tab.parameter.duration.hz_unit_id"] == "UO:0000031"


def test_magetab_opt_in_replaces_destinations_with_companions_and_retains_hz(tmp_path):
    source = tmp_path / "agentic.json"
    source.write_text(json.dumps({
        "miniml_json": package(), "harmonization_overrides": PROFILE
    }), encoding="utf-8")

    magetab = json2ae().convert(
        str(source), enrich=False, use_harmonization_overrides=True
    )[0]
    sdrf = next(row[1] for row in magetab if row and row[0] == "SDRF File")
    header, values = sdrf[0], sdrf[1]

    disease = header.index("Characteristics[disease]")
    assert values[disease] == "fallback disease"
    assert header[disease + 1:disease + 3] == ["Term Source REF", "Term Accession Number"]
    assert values[disease + 1:disease + 3] == ["mondo", "MONDO:9"]
    assert "Characteristics[hz_tissue_name]" in header


def test_all_json_consumers_expose_explicit_opt_in():
    assert "use_harmonization_overrides" in inspect.signature(json2ae.convert).parameters
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


def test_h5ad_publishes_magetab_parameter_obs_and_occurrence_ledger(tmp_path):
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

    assert converted.obs["msc.mage_tab.parameter.duration.value"].iat[0] == "30"
    assert converted.obs["msc.mage_tab.parameter.duration.hz_unit"].iat[0] == "minute"
    ledger = converted.uns["msc_mage_tab"]["parameters"]
    assert ledger["assay_path_id"].tolist() == ["study.sdrf.txt:row:1"]
    assert ledger["hz_unit_id"].tolist() == ["UO:0000031"]
