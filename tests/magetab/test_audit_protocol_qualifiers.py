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
from tests.support.ae import resolved_input


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
