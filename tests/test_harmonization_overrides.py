# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
import json

from meta_standards_converter.converters.harmonization_overrides import (
    resolve_harmonization_overrides,
)
from meta_standards_converter.converters.json2h5ad import JSON2H5ADConverter
from meta_standards_converter.converters.json_source import JSONPackageSource


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
