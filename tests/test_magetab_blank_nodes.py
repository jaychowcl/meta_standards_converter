# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from meta_standards_converter.magetab.parser import AEParser
from meta_standards_converter.sources.magetab import MAGETabInput, TextResource


def parse(header, rows, **kwargs):
    return AEParser().parse(MAGETabInput(
        TextResource('x.idf.txt', 'Investigation Accession\tE-GEOD-1\nInvestigation Title\tExample\nSDRF File\tx.sdrf.txt\n', 'memory:idf'),
        (TextResource('x.sdrf.txt', '\n'.join('\t'.join(r) for r in [header, *rows]), 'memory:sdrf'),),
        'E-GEOD-1', 'accession', **kwargs,
    ))


def test_blank_node_clears_attachment_context_and_preserves_unbound_occurrences():
    data = parse(['Source Name', 'Characteristics[organism]', 'Extract Name',
                  'Characteristics[material]', 'Unit', 'Comment[detail]', 'Assay Name',
                  'Array Data File', 'Comment[File URI]'],
                 [['s1', 'mouse', ' ', 'RNA', ' ng ', 'orphan', 'a1', ' ', 'https://example.org/orphan']])
    steps = data['series']['assay_paths'][0]['steps']
    assert [(s['kind'], s.get('name')) for s in steps] == [('source', 's1'), ('assay', 'a1')]
    assert steps[0]['characteristics'] == [{'name': 'organism', 'value': 'mouse'}]
    assert 'comments' not in steps[1] and 'link' not in steps[1]
    retained = data['extensions']['magetab']['unbound_annotations']
    assert [x['column_index'] for x in retained] == [3, 5, 8]
    assert retained[0]['unit'] == ' ng '
    assert all(x['row_index'] == 1 and x['sdrf'] == 'x.sdrf.txt' for x in retained)


def test_blank_reference_does_not_redirect_parameters_to_preceding_material():
    data = parse(['Source Name', 'Protocol REF', 'Parameter Value[duration]', 'Assay Name'],
                 [['s1', ' ', '20', 'a1']])
    steps = data['series']['assay_paths'][0]['steps']
    assert 'characteristics' not in steps[0]
    assert data['extensions']['magetab']['unbound_annotations'][0]['value'] == '20'


def test_unbound_material_and_description_do_not_populate_sample_scalars():
    data = parse(['Source Name','Extract Name','Material Type','Description','Assay Name'],
                 [['s1',' ','total RNA','orphan description','a1']])
    sample = data['sample'][0]
    assert not sample.get('description')
    assert not sample['channel'][0].get('molecule')
    assert len(data['extensions']['magetab']['unbound_annotations']) == 2


def test_catalogue_resolves_registered_path_and_unique_filename_without_guessing():
    catalog = ({'path': 'results/a.txt', 'uri': 'https://example.org/Files/results/a.txt'},
               {'path': 'one/b.txt', 'uri': 'https://example.org/Files/one/b.txt'},
               {'path': 'two/b.txt', 'uri': 'https://example.org/Files/two/b.txt'})
    data = parse(['Source Name', 'Derived Array Data File'],
                 [['s1', 'a.txt'], ['s1', 'b.txt'], ['s1', 'results/a.txt'], ['s1', 'unknown.txt']],
                 file_catalogue=catalog)
    files = [p['steps'][-1] for p in data['series']['assay_paths']]
    assert files[0]['link']['value'] == catalog[0]['uri']
    assert files[2]['link']['value'] == catalog[0]['uri']
    assert not files[1].get('link') and not files[3].get('link')
    assert files[0]['name'] == 'a.txt'
    assert {'value':catalog[0]['uri']} in data['sample'][0]['supplementary_data']
    biological = parse(['Source Name','Characteristics[label]'], [['s1','a.txt']], file_catalogue=catalog)
    assert biological['sample'][0]['channel'][0]['characteristics'][0]['value'] == 'a.txt'


def test_blank_units_preserve_companions_without_assigning_them_to_value():
    for blank in ('', '   '):
        for label, prefix, destination in [
            ('Characteristics[dose]', ['Source Name'], 'characteristics'),
            ('Factor Value[dose]', ['Source Name'], 'factor_values'),
            ('Parameter Value[dose]', ['Source Name','Protocol REF'], 'parameter_values'),
        ]:
            header = [*prefix, label, 'Unit[TimeUnit]', 'Term Source REF', 'Term Accession Number']
            row = ['s1', *(['P-1'] if len(prefix) == 2 else []), '5', blank, 'UO', 'UO:0000027']
            data = parse(header, [row])
            step = data['series']['assay_paths'][0]['steps'][-1]
            assert step[destination] == [{'name':'dose','value':'5'}]
            chars = data['sample'][0]['channel'][0].get('characteristics', [])
            if destination == 'characteristics': assert chars == [{'name':'dose','value':'5'}]
            unbound = data['extensions']['magetab']['unbound_annotations']
            assert len(unbound) == 1
            unit = unbound[0]
            assert unit['kind'] == 'unbound_unit' and unit['value'] == blank
            assert unit['column_index'] == len(prefix) + 1
            assert unit['unit_term_source_ref'] == 'UO' and unit['unit_term_accession_number'] == 'UO:0000027'
            assert unit['parent_attribute']['value'] == '5'
            assert unit['row_index'] == 1 and unit['sample_ref'] == 's1'
            assert parse(header, [row]) == data


def test_meaningful_unit_remains_a_coupled_core_and_path_group():
    data = parse(['Source Name','Characteristics[dose]','Unit[TimeUnit]','Term Source REF','Term Accession Number'],
                 [['s1','5','second','UO','UO:0000010']])
    expected = {'name':'dose','value':'5','unit':{'value':'second','term_source_ref':'UO','term_accession_number':'UO:0000010'},'unit_type':'TimeUnit'}
    assert data['sample'][0]['channel'][0]['characteristics'] == [expected]
    assert data['series']['assay_paths'][0]['steps'][0]['characteristics'] == [expected]
    assert not data.get('extensions',{}).get('magetab',{}).get('unbound_annotations')


def test_optional_blank_unit_header_alone_does_not_create_orphan_diagnostics():
    data = parse(['Source Name','Characteristics[dose]','Unit[TimeUnit]','Term Source REF','Term Accession Number'],
                 [['s1','5',' ','','']])
    assert not data.get('extensions',{}).get('magetab',{}).get('unbound_annotations')
    assert data['sample'][0]['channel'][0]['characteristics'] == [{'name':'dose','value':'5'}]
