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
