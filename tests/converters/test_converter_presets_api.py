# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
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


def test_flat_matrix_orientation_and_execution_profile_normalization():
    from meta_standards_converter.converters.unified.requests import normalize, AUTO, STANDARD
    def request(options, input_options=None, output_options=None):
        return normalize(target='h5ad', in_type=AUTO, enrichment=STANDARD,
                         options=options, force_in_type=None, outfile=None,
                         input_manifest=None, input_options=input_options,
                         output_options=output_options, runtime_options=None)
    resolved = request({'matrix_orientation': 'genes-by-observations', 'execution_profile': 'docker'},
                       output_options={'matrix_orientation':'genes-by-observations', 'profile':'docker'})
    assert resolved[4]['profile'] == 'docker'
    assert resolved[4]['matrix_orientation'] == 'genes-by-observations'
    with pytest.raises(ValueError, match='Conflicting'):
        request({'execution_profile':'docker'}, output_options={'profile':'singularity'})
    # Interdependent options can be supplied through different aliases.
    assert request({'resume': True}, output_options={'force_memory':True})[4]['force_memory']


def test_per_input_expansion_alias_overrides_default():
    from tests.converters.test_converter_expansion import geo, Geo
    reader = Geo({'GSE1':geo('GSE1',['GSE2'])})
    result = Converter(services={'geo2json':reader}).convert(
        InputSpec('GSE1', input_options={'related_series':False}),
        out_type='json', enrichment='off', options={'expand_studies':True})
    assert result.status == 'complete', result.to_dict()
    assert reader.calls == ['GSE1']


@pytest.mark.parametrize('local,call', [
    ({'input_options':{'enrichment':'off'}}, {'output_options':{'enrich':True}}),
    ({'output_options':{'enrich':False}}, {'enrichment':'standard'}),
])
def test_input_local_preparation_overrides_call_level_aliases(local,call):
    from tests.converters.test_converter_preparation import services
    calls=[]
    result=Converter(services=services(calls)).convert(InputSpec(package(),**local),
        out_type='magetab',**call)
    assert result.status=='complete',result.to_dict()
    assert not calls
