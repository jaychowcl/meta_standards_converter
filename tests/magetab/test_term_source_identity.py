# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Source-faithful term-source identity and ambiguous-reference round trips."""
import csv
import io
from pathlib import Path
from types import SimpleNamespace

import pytest
from meta_standards_converter.magetab.parser import AEParser
from meta_standards_converter.magetab.idf import IDFConstructor
from meta_standards_converter.magetab.constructor import AEConstructor
from meta_standards_converter.magetab.writer import MAGETabWriter
from meta_standards_converter.sources.magetab import MAGETabInput, TextResource
from meta_standards_converter.converters import AE2JSONConverter
from meta_standards_converter.miniml import MINiMLCodec


def source(declarations, *, secondary=False):
    rows = [['MAGE-TAB Version', '1.1'], ['Investigation Accession', 'E-MTAB-TEST'],
            ['Investigation Title', 'Ontology source identities'],
            ['Experimental Design', 'expression profiling'],
            ['Experimental Design Term Source REF', 'Ontology'],
            ['SDRF File', 'test.sdrf.txt']]
    if secondary:
        rows += [['Comment[SecondaryAccession]', 'CUSTOM1'], ['Comment[SecondaryAccessionTermSourceRef]', 'Ontology']]
    rows += [[label, *[d[index] for d in declarations]] for index, label in enumerate(
        ['Term Source Name', 'Term Source File', 'Term Source Version'])]
    stream = io.StringIO(); csv.writer(stream, delimiter='\t').writerows(rows)
    sdrf = 'Source Name\tCharacteristics[organism]\tTerm Source REF\tTerm Accession Number\tAssay Name\nS1\tHomo sapiens\tOntology\tTEST:1\tA1\n'
    return MAGETabInput(TextResource('test.idf.txt', stream.getvalue(), 'memory:idf'),
                       (TextResource('test.sdrf.txt', sdrf, 'memory:sdrf'),), 'E-MTAB-TEST', 'local')


def triples(package):
    return [(d.get('name'), d.get('url', ''), d.get('version', ''))
            for d in package['database'] if d.get('name', '').startswith('Ontology')]


def render(package):
    return AEConstructor(evidence_resolver=SimpleNamespace(sample_runs=lambda *args: {}, publications=lambda data: [])).miniml2magetab(package)


def table_triples(rows):
    values = {row[0]: row[1:] for row in rows if row[0].startswith('Term Source ')}
    return list(zip(values['Term Source Name'], values['Term Source File'], values['Term Source Version']))


def test_exact_duplicates_merge_without_warning_and_keep_source_hashes():
    parser = AEParser()
    package = parser.parse(source([(' Ontology ', ' https://example.org/o ', ''),
                                   ('Ontology', 'https://example.org/o', ' ')]))
    assert triples(package) == [('Ontology', 'https://example.org/o', '')]
    assert not any('Conflicting term source' in w for w in parser.warnings)
    assert len(package.source.documents) == 2
    assert all(d.sha256 for d in package.source.documents)


@pytest.mark.parametrize('declarations', [
    [('Ontology', 'https://example.org/a', '1'), ('Ontology', 'https://example.org/b', '1')],
    [('Ontology', 'https://example.org/a', '1'), ('Ontology', 'https://example.org/a', '2')],
    [('Ontology', '', ''), ('Ontology', 'https://example.org/a', '')],
    [('Ontology', '', ''), ('Ontology', '', '2')],
])
def test_conflicting_sources_preserved_with_unresolved_refs_and_roundtrip(declarations, tmp_path):
    parser = AEParser()
    package = parser.parse(source(declarations, secondary=True))
    assert triples(package) == declarations
    dbs = [d for d in package['database'] if d.get('name') == 'Ontology']
    assert [d['iid'] for d in dbs] == ['Ontology__msc_1', 'Ontology__msc_2']
    assert not any(d['iid'] == 'Ontology' for d in package['database'])
    assert any('Conflicting term source' in w and 'Ontology' in w for w in parser.warnings)
    assert any('unresolved' in w.lower() and 'Ontology' in w for w in parser.warnings)
    assert next(a for a in package['series']['accession'] if a['value'] == 'CUSTOM1')['database'] == 'Ontology'
    encoded = MINiMLCodec().encode(package)
    package = MINiMLCodec().decode(encoded).package
    rows = render(package)
    assert [tuple(v or '' for v in t) for t in table_triples(rows) if t[0] == 'Ontology'] == declarations
    sdrf = next(r[1] for r in rows if r[0] == 'SDRF File')
    assert sdrf[1][sdrf[0].index('Term Source REF')] == 'Ontology'
    path = MAGETabWriter().write(rows, str(tmp_path))
    again = AE2JSONConverter().convert(path)[0]
    assert triples(again) == declarations
    assert [d['iid'] for d in again['database'] if d.get('name') == 'Ontology'] == ['Ontology__msc_1', 'Ontology__msc_2']
    assert len(again['sample']) == 1


def test_generated_ids_reserve_all_original_names_and_deduplicate_variants():
    declarations = [('Ontology', 'a', ''), ('Ontology', 'b', ''), ('Ontology', 'a', ''),
                    ('Ontology__msc_1', 'reserved', '')]
    package = AEParser().parse(source(declarations))
    dbs = [d for d in package['database'] if d.get('name', '').startswith('Ontology')]
    assert [(d['iid'], d.get('url')) for d in dbs] == [
        ('Ontology__msc_2', 'a'), ('Ontology__msc_3', 'b'), ('Ontology__msc_1', 'reserved')]


def test_idf_does_not_overwrite_same_name_variants():
    rows = IDFConstructor()._idf_term_source([['Experimental Design Term Source REF', 'Ontology']],
        {'database': [{'iid': 'Ontology__msc_1', 'name': 'Ontology', 'url': 'a'},
                      {'iid': 'Ontology__msc_2', 'name': 'Ontology', 'url': 'b'}]})
    assert table_triples(rows) == [('Ontology', 'a', None), ('Ontology', 'b', None)]


def test_full_e_mtab_1_converts_and_retains_sample_and_assay_identities(tmp_path):
    source_dir = Path(__file__).parents[1] / 'fixtures/studies/E-MTAB-1/inputs'
    package = AE2JSONConverter().convert(str(source_dir / 'E-MTAB-1.idf.txt'))[0]
    original = list(csv.DictReader((source_dir / 'E-MTAB-1.sdrf.txt').open(), delimiter='\t'))
    assert len(original) == 176
    assert len(package['sample']) == len({row['Source Name'] for row in original}) == 45
    assert {s['iid'] for s in package['sample']} == {row['Source Name'] for row in original}
    assert len(package['series']['assay_paths']) == 176
    assert len([d for d in package['database'] if d.get('name') == 'The MGED Ontology']) == 1
    path = MAGETabWriter().write(render(package), str(tmp_path))
    again = AE2JSONConverter().convert(path)[0]
    assert {s['iid'] for s in again['sample']} == {s['iid'] for s in package['sample']}
    assert len(again['series']['assay_paths']) == 176


def test_casefolded_alias_of_one_database_remains_unambiguous():
    rows = IDFConstructor()._idf_term_source([['Experimental Design Term Source REF', 'ontology']],
        {'database': [{'iid': 'Ontology', 'name': 'ONTOLOGY', 'url': 'https://example.org/o'}]})
    assert table_triples(rows) == [('Ontology', 'https://example.org/o', None),
                                  ('ontology', 'https://example.org/o', None)]


def test_different_names_are_not_merged_even_with_identical_urls():
    package = AEParser().parse(source([('Ontology', 'same', ''), ('Ontology2', 'same', '')]))
    assert triples(package) == [('Ontology', 'same', ''), ('Ontology2', 'same', '')]
