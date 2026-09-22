"""Explicit archive relationships survive without expanding related studies."""
from pathlib import Path

from meta_standards_converter.converters import Converter
from meta_standards_converter.magetab.idf import IDFConstructor


def test_explicit_ae_relationships_survive_facade_roundtrip_without_expansion():
    path = Path(__file__).parents[1] / 'fixtures/fidelity/E-MTAB-14416.idf.txt'
    result = Converter().convert(path, out_type='magetab', enrichment='off', options={'expand_studies': False})
    assert result.status == 'complete', result.to_dict()
    idf = result.items[0].payload[0]
    assert next(r[1:] for r in idf if r[0] == 'Comment[RelatedExperiment]') == [
        'E-MTAB-10619', 'E-MTAB-14415', 'E-MTAB-14417', 'E-MTAB-14419']


def test_explicit_targets_keep_identity_order_and_do_not_promote_other_relations():
    relations = [{'type': 'related experiment', 'target': 'E-MTAB-2'},
                 {'type': 'related experiment', 'target': 'GSE9'},
                 {'type': 'related experiment', 'target': 'E-MTAB-2'},
                 {'type': 'BioProject', 'target': 'PRJNA2'},
                 {'type': 'publication', 'target': 'GSE10'}]
    assert IDFConstructor()._related_experiments({'series': {'relation': relations}}) == ['E-MTAB-2', 'GSE9']
