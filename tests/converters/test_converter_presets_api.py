"""Public preset interface and compatibility option contracts."""
import json
from pathlib import Path

import pytest

from meta_standards_converter.converters import Converter, InputSpec
from tests.converters.test_json2tsv import package


def test_simple_interface_and_exact_file_alias(tmp_path):
    result = Converter().convert(package(), tmp_path, out_type="json",
                                 in_type="miniml", enrichment="off")
    assert result.status == "complete", result.to_dict()
    assert (tmp_path / "GSE1.json").exists()
    target = tmp_path / "exact.csv"
    result = Converter().convert(package(), out_type="csv", enrichment="off",
                                 options={"outfile": target, "aggregate": True})
    assert result.status == "complete", result.to_dict()
    assert target.read_text().startswith("msc.")


def test_equivalent_aliases_and_no_mutation(tmp_path):
    options = {"outfile": tmp_path / "one.json", "overwrite": False}
    result = Converter().convert(package(), out_type="json", in_type="miniml",
        force_in_type="msc_miniml", enrichment="off", options=options,
        outfile=str(options["outfile"]), runtime_options={"overwrite": False})
    assert result.status == "complete", result.to_dict()
    assert options == {"outfile": tmp_path / "one.json", "overwrite": False}


@pytest.mark.parametrize("arguments", [
    {"in_type": "miniml", "force_in_type": "atlas"},
    {"options": {"overwrite": True}, "runtime_options": {"overwrite": False}},
    {"enrichment": "unknown"},
    {"options": {"expand_studies": "yes"}},
    {"options": {"execution_profile": "docker"}},
    {"options": {"overwite": True}},
    {"options": {"outfile": "a.json"}, "outfile": "b.json"},
])
def test_invalid_configuration_precedes_loading(arguments):
    with pytest.raises((ValueError, TypeError)):
        Converter().convert(Path("missing.json"), out_type="json", **arguments)


def test_manifest_alias_and_per_input_preparation_options(tmp_path):
    (tmp_path / "source.json").write_text(json.dumps(package()))
    manifest = tmp_path / "inputs.json"
    manifest.write_text(json.dumps({"schema_version": "1.0", "inputs": [
        {"id": "chosen", "sources": ["source.json"],
         "input_options": {"enrichment": "off", "expand_studies": False}}
    ]}))
    result = Converter().convert(out_type="json", options={"input_manifest": manifest})
    assert result.status == "complete", result.to_dict()
    assert result.items[0].id == "chosen"


def test_policy_defaults_are_resolved_per_input():
    from meta_standards_converter.converters.unified.contracts import LoadedInput
    from meta_standards_converter.sources.json import JSONPackageSource

    policies = []
    class Handler:
        kind = "miniml"
        def probe(self, source, context):
            return True
        def validate(self, spec, context):
            pass
        def load(self, spec, context):
            policies.append(context.preparation_policy)
            return LoadedInput(metadata=JSONPackageSource().decode(spec.sources))

    converter = Converter(handlers=[Handler()])
    result = converter.convert([
        InputSpec(package(), id="default"),
        InputSpec(package(), id="custom", input_options={"enrichment": "off", "expand_studies": False}),
    ], out_type="json")
    assert result.status == "complete", result.to_dict()
    assert [(p.enrichment, p.expand_studies) for p in policies] == [("standard", True), ("off", False)]
