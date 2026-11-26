# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Source-level ownership, organism enrichment and protocol origin."""
from copy import deepcopy
import json

import pytest
from meta_standards_converter.metadata.integrity import qualifier_ownership
from meta_standards_converter.metadata.interpretation import MINiMLMetadataService
from meta_standards_converter.metadata.projection.assay import _parameter_rows
from meta_standards_converter.metadata.workflows import prepare_workflows
from pathlib import Path
from meta_standards_converter.miniml.geo_parser import GEOParser


def study():
    samples = [{"iid": "A", "channel": [{"organism": [{"value": "human", "taxid": "9606"}],
               "characteristics": [{"name": "organism part", "value": "lung"}]}]},
               {"iid": "B"}]
    paths = []
    for i, sample in enumerate(samples, 1):
        owner = sample["iid"]
        paths.append({"steps": [
            {"kind": "source", "name": owner, "sample_ref": owner,
             "comments": [{"name": "BioSD_SAMPLE", "value": f"SAMEA{i}"}],
             "characteristics": [{"name": "organism", "value": "human"}]},
            {"kind": "protocol_application", "protocol_ref": owner,
             "parameter_values": [{"name": "duration", "value": str(i * 10),
                                   "unit": {"value": "min"}}]},
            {"kind": "assay", "name": "assay-" + owner, "sample_ref": owner}]})
    return {"miniml_schema_version": "3.1", "source": {"format": "test"},
            "sample": samples, "series": {"iid": "S", "assay_paths": paths,
            "protocols": [{"name": s["iid"], "type": {"value": "protocol-" + s["iid"]}}
                          for s in samples]}}


def test_material_and_anatomy_have_independent_evidence():
    data = study()
    row = MINiMLMetadataService().sample_metadata(data["sample"][0], data)
    assert row["organism_part"] == "lung"
    assert not row["material_type"]


def test_biosample_evidence_and_conflict_are_owned():
    data = study()
    sample = data["sample"][0]
    sample["relation"] = [{"type": "BioSample", "target": "https://www.ncbi.nlm.nih.gov/biosample/SAMN99"},
                          {"type": "Reanalyzed by", "target": "SAMD77"}]
    sample["sra_run"] = [{"biosample": "SAMEA1"}]
    row = MINiMLMetadataService().sample_metadata(sample, data)
    assert set(row["biosample_accession"].split("; ")) == {"SAMN99", "SAMEA1"}
    assert row["biosample_conflict"]
    row = MINiMLMetadataService().sample_metadata(data["sample"][1], data)
    assert row["biosample_accession"] == "SAMEA2"
    assert not row["biosample_conflict"]


def test_source_and_assay_refs_scope_methods_and_preserve_units():
    data = study()
    for index, sample in enumerate(data["sample"], 1):
        row = MINiMLMetadataService().sample_metadata(sample, data)
        assert row["protocol_types"] == "protocol-" + sample["iid"]
        params = _parameter_rows(data, sample)
        assert [(p["protocol_ref"], p["value"], p["unit"]) for p in params] == [
            (sample["iid"], str(index * 10), "min")]


def test_shared_or_unbound_procedure_remains_in_study_scope():
    data = study()
    first = data["series"]["assay_paths"][0]
    first["steps"][-1]["sample_ref"] = "B"
    for step in data["series"]["assay_paths"][1]["steps"]:
        step.pop("sample_ref", None)
    assert all(not _parameter_rows(data, s) for s in data["sample"])
    for s in data["sample"]:
        assert not MINiMLMetadataService().sample_metadata(s, data)["protocol_types"]
    assert [r["ownership"] for r in _parameter_rows(data)] == ["pooled", "unassigned"]


def test_structural_protocols_have_occurrence_evidence():
    source = Path(__file__).parents[1] / "fixtures/studies/GSE100/inputs/geo.xml"
    package = prepare_workflows(GEOParser().parse(source.read_text())[0]).to_mapping()
    for sample in package["sample"]:
        row = MINiMLMetadataService().sample_metadata(sample, package)
        entries = json.loads(row["protocol_evidence"])
        assert {(r["type"], r["origin"]) for r in entries} == {
            ("sample collection protocol", "structural_placeholder"),
            ("hybridization protocol", "technology_inferred")}
    supplied = study()
    row = MINiMLMetadataService().sample_metadata(supplied["sample"][0], supplied)
    assert json.loads(row["protocol_evidence"])[0]["origin"] == "supplied_definition"


@pytest.mark.parametrize("mutation", ["supported", "supported_obo_iri", "wrong_taxon", "changed_unit", "lost_duplicate", "moved", "overwritten"])
def test_taxonomy_addition_preserves_original_occurrences(mutation):
    data = study()
    items = data["series"]["assay_paths"][0]["steps"][0]["characteristics"]
    items[0]["unit"] = {"value": "u"}
    if mutation == "overwritten":
        items[0]["term_source_ref"] = "source"
    items.append(deepcopy(items[0]))
    revised = deepcopy(data)
    steps = revised["series"]["assay_paths"][0]["steps"]
    new = steps[0]["characteristics"]
    new[0].update(term_source_ref="NCBITaxon", term_accession_number="9606")
    if mutation == "supported_obo_iri":
        for entry in revised["series"]["assay_paths"][0]["steps"][0]["characteristics"]:
            if entry.get("term_source_ref") == "NCBITaxon":
                entry["term_accession_number"] = "http://purl.obolibrary.org/obo/NCBITaxon_" + entry["term_accession_number"]
    if mutation == "wrong_taxon":
        new[0]["term_accession_number"] = "10090"
    elif mutation == "changed_unit":
        new[0]["unit"]["value"] = "v"
    elif mutation == "lost_duplicate":
        new.pop()
    elif mutation == "moved":
        steps[-1]["characteristics"] = [new.pop(0)]
    assert bool(qualifier_ownership(data, revised)) == (mutation not in {"supported", "supported_obo_iri"})


def test_distinct_supplied_values_retain_their_relative_order():
    package = study()
    first = package["series"]["assay_paths"][0]["steps"][0]
    first["characteristics"].append({"name": "treatment", "value": "oxygen"})
    generated = deepcopy(package)
    generated["series"]["assay_paths"][0]["steps"][0]["characteristics"].reverse()
    assert qualifier_ownership(package, generated)
