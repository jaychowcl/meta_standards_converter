# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""The public request has one explicit owner for every preparation setting."""
from dataclasses import FrozenInstanceError
import json
import pytest

from meta_standards_converter.converters import Converter, InputSpec, PreparationPolicy, RuntimeOptions
from tests.converters.test_json2tsv import package


def test_typed_interface_and_exact_output(tmp_path):
    policy = PreparationPolicy(enrichment="off", expand_studies=False)
    result = Converter().convert(package(), tmp_path, out_type="json", in_type="miniml", preparation=policy)
    assert result.status == "complete", result.to_dict()
    assert (tmp_path / "GSE1.json").exists()
    result = Converter().convert(package(), out_type="csv", preparation=policy,
                                 outfile=tmp_path / "exact.csv", output_options={"aggregate": True})
    assert result.status == "complete", result.to_dict()
    assert (tmp_path / "exact.csv").read_text().startswith("msc.")


def test_configuration_is_immutable_and_owns_override_snapshot():
    overrides = {"max_xml_bytes": 1234}
    runtime = RuntimeOptions(resource_overrides=overrides)
    overrides["max_xml_bytes"] = 9999
    assert runtime.resource_overrides["max_xml_bytes"] == 1234
    with pytest.raises(TypeError):
        runtime.resource_overrides["max_xml_bytes"] = 9999
    with pytest.raises(FrozenInstanceError):
        runtime.overwrite = True
    with pytest.raises(FrozenInstanceError):
        PreparationPolicy().expand_studies = False


@pytest.mark.parametrize("arguments", [
    {"options": {}}, {"enrichment": "off"}, {"force_in_type": "json"},
    {"runtime_options": {}}, {"preparation": {}}, {"runtime": {}},
    {"input_options": {"enrich": False}}, {"output_options": {"enrich": False}},
    {"input_options": {"related_series": False}},
])
def test_obsolete_or_misplaced_settings_are_rejected(arguments):
    with pytest.raises((ValueError, TypeError)):
        Converter().convert(package(), out_type="json", **arguments)


@pytest.mark.parametrize("settings", [{"enrichment": "unknown"}, {"expand_studies": 1}, {"linked_metadata": None}, {"peer_provider": "yes"}])
def test_invalid_policy_rejected(settings):
    with pytest.raises((ValueError, TypeError)):
        PreparationPolicy(**settings)


def test_manifest_v2_policy_replaces_call_policy(tmp_path):
    (tmp_path / "source.json").write_text(json.dumps(package()))
    manifest = {"schema_version": "2.0", "inputs": [
        {"id": "chosen", "sources": [str(tmp_path / "source.json")],
         "preparation": {"enrichment": "off", "expand_studies": False}}
    ]}
    result = Converter().convert(out_type="json", input_manifest=manifest)
    assert result.status == "complete", result.to_dict()
    assert result.items[0].id == "chosen"
    assert any(p.get("reason") == "disabled" for p in result.items[0].preparation)
    manifest["schema_version"] = "1.0"
    with pytest.raises(ValueError, match="schema_version 2.0"):
        Converter().convert(out_type="json", input_manifest=manifest)


def test_policy_defaults_are_resolved_per_input():
    from meta_standards_converter.converters.unified.contracts import LoadedInput
    from meta_standards_converter.sources.json import JSONPackageSource
    policies = []
    class Handler:
        kind = "miniml"
        def probe(self, source, context): return True
        def validate(self, spec, context): pass
        def load(self, spec, context):
            policies.append(context.preparation_policy)
            return LoadedInput(metadata=JSONPackageSource().decode(spec.sources))
    result = Converter(handlers=[Handler()]).convert([
        InputSpec(package(), id="default"),
        InputSpec(package(), id="custom", preparation=PreparationPolicy(enrichment="off", expand_studies=False)),
    ], out_type="json")
    assert result.status == "complete", result.to_dict()
    assert [(p.enrichment, p.expand_studies) for p in policies] == [("standard", True), ("off", False)]


def test_enrichment_and_expansion_are_independent():
    assert PreparationPolicy(enrichment="off").expand_studies is True
    assert not PreparationPolicy(enrichment="off").allows("publications")
    assert PreparationPolicy(expand_studies=False).allows("publications")
