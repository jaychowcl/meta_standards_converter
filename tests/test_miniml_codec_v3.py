# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from __future__ import annotations

from dataclasses import FrozenInstanceError

import pytest

from meta_standards_converter.miniml import (
    FASTQFile,
    MINiMLBatchDecodeResult,
    MINiMLCodec,
    MINiMLCompatibilityError,
    MINiMLDecodeResult,
    MINiMLPackage,
    PubMedPublication,
    SRARun,
    Sample,
    Series,
    SourceInfo,
    iter_harmonized_values,
)


def package_payload() -> dict:
    return {
        "miniml_schema_version": "3.0",
        "source": {"format": "test"},
        "series": {
            "iid": "GSE1",
            "pubmed_publication": [
                {
                    "pubmed_id": "123",
                    "doi": "10.1/example",
                    "title": "Study",
                }
            ],
            "sample_ref": [{"ref": "GSM1"}],
        },
        "sample": [
            {
                "iid": "GSM1",
                "sra_accession": ["SRX1"],
                "ena_accession": ["ERP1"],
                "sra_run": [
                    {
                        "run": "SRR1",
                        "study": "ERP1",
                        "fastq_files": [
                            {
                                "filename": "reads.fastq.gz",
                                "uri": "https://example/reads.fastq.gz",
                                "md5": "0123456789abcdef0123456789abcdef",
                            }
                        ],
                    }
                ],
            }
        ],
    }


def test_codec_decodes_typed_enrichment_from_v3_input() -> None:
    result = MINiMLCodec().decode(package_payload())

    assert isinstance(result, MINiMLDecodeResult)
    assert isinstance(result.package, MINiMLPackage)
    assert isinstance(result.package.series.pubmed_publications[0], PubMedPublication)
    assert isinstance(result.package.samples[0].sra_runs[0], SRARun)
    assert isinstance(result.package.samples[0].sra_runs[0].fastq_files[0], FASTQFile)
    assert MINiMLCodec().encode(result.package)["miniml_schema_version"] == "3.0"


def test_codec_returns_warnings_and_strict_mode_promotes_them() -> None:
    payload = package_payload()
    payload["series"]["sample_ref"].append({"ref": "GSM404"})

    compatible = MINiMLCodec().decode(payload)

    assert [issue.code for issue in compatible.diagnostics] == [
        "unresolved_reference",
    ]
    with pytest.raises(MINiMLCompatibilityError, match="unresolved_reference"):
        MINiMLCodec().decode(payload, strict=True)


def test_codec_decodes_one_or_many_packages() -> None:
    single = MINiMLCodec().decode_many(package_payload())
    multiple = MINiMLCodec().decode_many([package_payload(), package_payload()])

    assert isinstance(single, MINiMLBatchDecodeResult)
    assert len(single.packages) == 1
    assert len(multiple.packages) == 2


def test_typed_package_is_immutable_and_preserves_typed_harmonized_values() -> None:
    payload = package_payload()
    payload['sample'][0]['channel'] = [
        {
            'source': 'lung',
            'characteristics': [
                {
                    'name': 'organism',
                    'value': 'Homo sapiens',
                },
                {
                    'name': 'hz_organism',
                    'value': 'Homo sapiens',
                },
                {
                    'name': 'hz_organism_id',
                    'value': 'NCBITaxon:9606',
                },
                {
                    'name': 'hz_organism_onto',
                    'value': 'NCBITaxon',
                },
            ],
        },
    ]
    package = MINiMLCodec().decode(payload).package

    rows = package.to_mapping()["sample"][0]["channel"][0]["characteristics"]
    annotation = iter_harmonized_values(rows)[0]
    assert annotation.value == "Homo sapiens"
    assert annotation.term_accession_number == "NCBITaxon:9606"
    with pytest.raises(FrozenInstanceError):
        package.version = "changed"  # type: ignore[misc]
    with pytest.raises(FrozenInstanceError):
        annotation.value = "changed"  # type: ignore[misc]


def test_v3_characteristic_annotations_align_duplicate_fields_with_parentheses():
    payload = package_payload()
    payload['sample'][0]['channel'] = [
        {
            'characteristics': [
                {
                    'name': 'tissue',
                    'value': 'kidney',
                },
                {
                    'name': 'hz_tissue_name',
                    'value': 'kidney',
                },
                {
                    'name': 'hz_tissue_name_id',
                    'value': 'UBERON:0002113',
                },
                {
                    'name': 'organism part',
                    'value': 'lung',
                },
                {
                    'name': 'hz_tissue_name(1)',
                    'value': 'lung',
                },
                {
                    'name': 'hz_tissue_name_id(1)',
                    'value': 'UBERON:0002048',
                },
            ],
        },
    ]

    migrated = MINiMLCodec().encode(MINiMLCodec().decode(payload).package)
    rows = migrated["sample"][0]["channel"][0]["characteristics"]

    assert [row["name"] for row in rows] == [
        "tissue",
        "hz_tissue_name",
        "hz_tissue_name_id",
        "organism part",
        "hz_tissue_name(1)",
        "hz_tissue_name_id(1)",
    ]
    assert [value.value for value in iter_harmonized_values(rows)] == [
        "kidney",
        "lung",
    ]


def test_direct_model_construction_deep_freezes_collections_and_extensions() -> None:
    root_extensions = {"nested": {"values": []}}
    series_extensions = {"vendor": []}
    samples = [Sample(iid="GSM1")]
    package = MINiMLPackage(
        series=Series(iid="GSE1", extras=series_extensions),
        source=SourceInfo("test"),
        samples=samples,
        extensions=root_extensions,
    )

    root_extensions["nested"]["values"].append("mutated")
    series_extensions["vendor"].append("mutated")
    samples.append(Sample(iid="GSM2"))

    assert not package.series.extras
    assert package.to_mapping()["extensions"] == {
        "nested": {"values": []},
        "vendor": [],
    }
    assert tuple(item.iid for item in package.samples) == ("GSM1",)
    with pytest.raises(TypeError):
        package.extensions["new"] = "value"  # type: ignore[index]


def test_certified_children_preserve_public_mutation_isolation():
    from meta_standards_converter.miniml import model

    authored = {'nested': {'values': ['original']}}
    child = model.Accession('A1', extras=authored)
    parent = model.Series(iid='S1', accessions=[child])
    assert parent.accessions[0] is child
    authored['nested']['values'].append('changed')
    exported = parent.to_mapping()
    exported['accession'][0]['nested']['values'].append('changed')
    assert child.to_mapping() == {'value': 'A1', 'nested': {'values': ['original']}}
    with pytest.raises(FrozenInstanceError):
        child.value = 'changed'
    with pytest.raises(TypeError):
        child.extras['nested']['values'][0] = 'changed'


def test_unknown_values_and_model_subclasses_keep_defensive_copy_fallback():
    from dataclasses import dataclass
    from meta_standards_converter.miniml import model

    class MutableValue:
        def __init__(self):
            self.values = ['original']

    @dataclass(frozen=True)
    class ForeignFrozen:
        values: list

    class AccessionSubclass(model.Accession):
        pass

    class StringSubclass(str):
        pass

    text = StringSubclass('A1')
    text.values = ['original']
    subclass = AccessionSubclass('A1')
    subclass.extra_state = ['original']
    for value in (MutableValue(), ForeignFrozen(['original']), subclass,
                  model.Accession(MutableValue()), model.Accession(text)):
        assert not model._certified_immutable(value)
        frozen = model._freeze_json(value)
        assert frozen is not value
        if isinstance(value, model.Accession) and type(value) is model.Accession:
            value.value.values.append('changed')
            assert frozen.value.values == ['original']
        elif value is subclass:
            value.extra_state.append('changed')
            assert frozen.extra_state == ['original']
        else:
            value.values.append('changed')
            assert frozen.values == ['original']


def test_unknown_extension_values_and_mapping_subclasses_are_not_certified():
    from meta_standards_converter.miniml import model

    class MutableValue:
        def __init__(self):
            self.values = ['original']

    class FrozenMappingSubclass(model._FrozenJSONMapping):
        pass

    class TupleSubclass(tuple):
        pass

    source = MutableValue()
    mapping = model._FrozenJSONMapping({'custom': source})
    assert not model._certified_immutable(mapping)
    assert not model._certified_immutable(model.Accession('A1', extras=mapping))
    reconstructed = model._freeze_json(mapping)
    assert reconstructed is not mapping
    assert reconstructed['custom'] is not mapping['custom']
    source.values.append('changed')
    assert reconstructed['custom'].values == ['original']
    subclass = FrozenMappingSubclass({'value': 'original'})
    assert not model._certified_immutable(subclass)
    assert type(model._freeze_json(subclass)) is model._FrozenJSONMapping
    assert not model._certified_immutable(TupleSubclass(('original',)))
    assert type(model._freeze_json(TupleSubclass(('original',)))) is tuple


def test_mapping_models_keep_mapping_freeze_behavior():
    from meta_standards_converter.miniml import model

    package = MINiMLCodec().decode(package_payload()).package
    assert type(package) not in model._IMMUTABLE_MODEL_FIELDS
    assert not model._certified_immutable(package)
    frozen = model._freeze_json(package)
    assert type(frozen) is model._FrozenJSONMapping
    assert model._plain(frozen) == package.to_mapping()


def test_construction_errors_and_unexpected_private_state_are_not_certified(monkeypatch):
    from meta_standards_converter.miniml import model

    original = model.NamedComment.__post_init__

    def with_private_state(self):
        original(self)
        object.__setattr__(self, '_custom_state', [])

    monkeypatch.setattr(model.NamedComment, '__post_init__', with_private_state)
    child = model.NamedComment('name', 'value')
    assert not model._certified_immutable(child)
    assert model._freeze_json(child) is not child
    incomplete = object.__new__(model.NamedComment)
    with pytest.raises(model.MINiMLModelError, match='nonblank name'):
        model.NamedComment.__init__(incomplete, '', 'value')
    assert not model._certified_immutable(incomplete)


def test_copy_replace_and_pickle_preserve_value_contract_without_global_deepcopy_override():
    from copy import copy, deepcopy
    from dataclasses import fields, replace
    import pickle
    from meta_standards_converter.miniml import model

    child = model.Accession('A1', extras={'nested': ['original']})
    for duplicate in (copy(child), deepcopy(child), pickle.loads(pickle.dumps(child))):
        assert duplicate == child
        assert duplicate is not child
        assert model._certified_immutable(duplicate)
        assert duplicate.to_mapping() == child.to_mapping()
    assert model._IMMUTABLE_CERTIFICATE not in {item.name for item in fields(child)}
    assert model._IMMUTABLE_CERTIFICATE not in child.to_mapping()
    modified = replace(child, extras={'nested': ['replacement']})
    assert modified.to_mapping()['nested'] == ['replacement']
    assert child.to_mapping()['nested'] == ['original']
    assert model._certified_immutable(modified)


def test_certified_model_construction_preserves_frozen_golden_mapping():
    import json
    from pathlib import Path

    golden = json.loads((Path(__file__).parent / 'fixtures/edge_cases/'
        'magetab-normalization/expected/ae2json/E-MTAB-6486.json').read_text())
    codec = MINiMLCodec()
    decoded = codec.decode_many(golden, strict=True)
    assert [package.to_mapping() for package in decoded.packages] == golden


def test_certification_does_not_bypass_strict_codec_reconstruction():
    from meta_standards_converter.miniml import model

    series = Series(iid='S1', sample_ref=(model.Reference('missing'),))
    assert model._certified_immutable(series)
    package = MINiMLPackage(series, SourceInfo('test'))
    with pytest.raises(MINiMLCompatibilityError, match='unresolved_reference'):
        MINiMLCodec().decode(package, strict=True)


def test_codec_dump_is_deterministic_and_replaces_existing_file(tmp_path) -> None:
    codec = MINiMLCodec()
    package = codec.decode(package_payload()).package
    destination = tmp_path / "package.json"
    destination.write_text("old", encoding="utf-8")

    codec.dump(package, destination)
    first = destination.read_bytes()
    codec.dump(package, destination)

    assert destination.read_bytes() == first
    assert not list(tmp_path.glob(".*.tmp"))


def test_channel_hz_groups_cover_harmonized_channel_scalars() -> None:
    payload = package_payload()
    payload['sample'][0]['channel'] = [{'source': 'lung', 'hz_tissue_name': 'lung', 'hz_tissue_name_id': 'UBERON:0002048', 'hz_tissue_name_onto': 'UBERON'}]

    package = MINiMLCodec().decode(payload).package

    channel = package.to_mapping()["sample"][0]["channel"][0]
    assert channel["hz_tissue_name"] == "lung"
    assert channel["hz_tissue_name_id"] == "UBERON:0002048"
    assert iter_harmonized_values(channel)[0].field == "tissue_name"
