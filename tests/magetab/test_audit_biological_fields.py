# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Regression evidence from the expanded scientific audit."""
from pathlib import Path
from copy import deepcopy

from meta_standards_converter.converters import Converter
from meta_standards_converter.magetab.sdrf.handlers.base import _BaseSDRFHandler


def test_geo_rna_descriptions_are_not_anatomical_parts():
    source = Path(__file__).parents[1] / 'fixtures/studies/GSE100/inputs/geo.xml'
    result = Converter().convert(source, out_type='magetab', enrichment='off', options={'expand_studies': False})
    assert result.status == 'complete', result.to_dict()
    rows = next(r[1] for r in result.items[0].payload[0] if r[0] == 'SDRF File')
    header = rows[0]
    assert {'Total RNA', 'zzMex67 co-IPed RNA'} <= {r[header.index('Comment[Sample_source_name]')] for r in rows[1:]}
    assert not any(r[header.index('Characteristics[organism part]')] for r in rows[1:])


def test_generated_anatomy_retains_alias_values_and_each_qualifier():
    handler = _BaseSDRFHandler(None, {})
    channel = {'source': 'total RNA', 'characteristics': [
        {'name': 'tissue', 'value': 'lung', 'term_source_ref': 'UBERON', 'term_accession_number': 'UBERON:0002048'},
        {'name': 'OrganismPart', 'value': 'heart', 'term_source_ref': 'UBERON', 'term_accession_number': 'UBERON:0000948'},
        {'name': 'organism_part', 'value': 'lung', 'term_source_ref': 'UBERON', 'term_accession_number': 'UBERON:0002048'},
    ]}
    before = deepcopy(channel)
    attrs = handler.characteristic_attrs(channel)
    anatomy = [a for a in attrs if a.label == 'Characteristics[organism part]']
    assert [a.value for a in anatomy] == ['lung', 'heart']
    assert [[x.value for x in a.attrs] for a in anatomy] == [['UBERON', 'UBERON:0002048'], ['UBERON', 'UBERON:0000948']]
    assert {a.label for a in attrs} >= {'Characteristics[tissue]', 'Characteristics[OrganismPart]', 'Characteristics[organism_part]'}
    assert channel == before
