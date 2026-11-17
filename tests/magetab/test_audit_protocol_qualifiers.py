# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Protocol qualifiers belong to their individual application, not a comment."""
from meta_standards_converter.converters import PreparationPolicy
from meta_standards_converter.converters import Converter
from meta_standards_converter.miniml import MINiMLCodec
from meta_standards_converter.miniml.model import ProtocolApplication
from meta_standards_converter.magetab.parser import AEParser
from tests.converters.test_ae2json import resolved_input


def test_protocol_qualifiers_survive_typed_and_archive_boundaries():
    header = ['Source Name', 'Protocol REF', 'Term Source REF', 'Term Accession Number',
              'Parameter Value[temperature]', 'Term Source REF', 'Term Accession Number',
              'Extract Name', 'Protocol REF', 'Term Source REF', 'Assay Name']
    row = ['s1', 'P-MTAB-1', 'ArrayExpress', 'P-MTAB-1', 'cold', 'EFO', 'EFO:1', 'e1', 'P-MTAB-2', 'ArrayExpress', 'a1']
    source = resolved_input(idf='Investigation Accession\tE-MTAB-1\nProtocol Name\tP-MTAB-1\tP-MTAB-2\nTerm Source Name\tArrayExpress\tEFO\n',
        sdrfs=['\t'.join(header)+'\n'+'\t'.join(row)+'\n'])
    package = AEParser().parse(source)
    codec = MINiMLCodec()
    package = codec.decode(package.to_mapping(), strict=True).package
    protocols = [s for s in package.series.assay_paths[0].steps if isinstance(s, ProtocolApplication)]
    assert protocols[0].term_source_ref == 'ArrayExpress'
    assert protocols[0].term_accession_number == 'P-MTAB-1'
    assert protocols[0].parameter_values[0].term_source_ref == 'EFO'
    assert protocols[1].term_source_ref == 'ArrayExpress'
    result = Converter().convert(package, out_type='magetab', preparation=PreparationPolicy(enrichment='off', expand_studies=False))
    assert result.status == 'complete', result.to_dict()
    idf = result.items[0].payload[0]
    table = next(r[1] for r in idf if r[0] == 'SDRF File')
    assert table == [header, row]
    assert 'ArrayExpress' in next(r[1:] for r in idf if r[0] == 'Term Source Name')


def test_existing_positional_protocol_constructor_is_unchanged():
    old = ProtocolApplication('P1', (), 'person', '2026-01-01', ())
    assert old.to_mapping() == {'kind': 'protocol_application', 'protocol_ref': 'P1', 'performer': 'person', 'date': '2026-01-01'}
    qualified = ProtocolApplication.from_mapping({**old.to_mapping(), 'term_source_ref': 'ArrayExpress'})
    assert qualified.to_mapping()['term_source_ref'] == 'ArrayExpress'


def test_recorded_four_study_protocol_qualifier_occurrences_survive():
    import csv
    from collections import Counter
    from pathlib import Path
    total = 0
    def occurrences(table):
        headers = table[0]
        columns, owner, occurrence = [], None, 0
        for i, header in enumerate(headers):
            if header == 'Protocol REF':
                occurrence += 1
                owner = ('protocol', occurrence, i)
            elif header not in {'Term Source REF', 'Term Accession Number'}:
                owner = None
            if owner and header in {'Term Source REF', 'Term Accession Number'}:
                columns.append((owner, header, i))
        return Counter((owner[1], header, row[owner[2]], row[i])
                       for row in table[1:] for owner, header, i in columns if row[i])
    for accession, expected in [('E-GEOD-48034',20), ('E-GEOD-48035',844), ('E-GEOD-52564',34), ('E-GEOD-60450',48)]:
        root = Path(__file__).parents[1] / 'fixtures' / 'fidelity'
        with (root / (accession + '.sdrf.txt')).open() as stream:
            source = list(csv.reader(stream, delimiter='\t'))
        result = Converter().convert(str(root / (accession + '.idf.txt')), out_type='magetab', preparation=PreparationPolicy(enrichment='off', expand_studies=False))
        assert result.status == ('partial' if accession == 'E-GEOD-48035' else 'complete'), result.to_dict()
        assert result.items[0].validation == 'valid'
        assert all('source-text encoding' in d.message for d in result.items[0].diagnostics if d.code == 'source_partial')
        output = next(r[1] for r in result.items[0].payload[0] if r[0] == 'SDRF File')
        original = occurrences(source)
        assert sum(original.values()) == expected
        assert occurrences(output) == original
        assert len(output) == len(source)
        total += expected
    assert total == 946
