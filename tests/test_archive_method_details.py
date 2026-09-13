# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Source-faithful protocol detail projection and export contracts."""
import json
import xml.etree.ElementTree as ET

import pytest

from meta_standards_converter.miniml.ena_parser import ENAParser
from meta_standards_converter.miniml.sra_parser import SRAParser
from tests.test_native_archive_parsers import fixture_records

METHOD = 'Library was made using a Hi-C - Arima v2 kit with restriction enzyme motif ^GATC,G^ANTC,C^TNAG,T^TAA.'
CONTEXT = ('Illumina sequencing of sample accession SAMEA7520545 for study accession PRJEB43529. '
           'This is part of an Illumina multiplexed sequencing run (34892_5). '
           'This submission includes reads tagged with the sequence TTATTATG. ')


def records_with_method(provider, short='Hi-C - Arima v2', design=CONTEXT + METHOD):
    records = fixture_records(provider)
    experiment = next(n for root in records.xml for n in root.iter('EXPERIMENT'))
    descriptor = experiment.find('DESIGN/LIBRARY_DESCRIPTOR')
    node = descriptor.find('LIBRARY_CONSTRUCTION_PROTOCOL')
    if node is not None:
        descriptor.remove(node)
    if short is not None:
        ET.SubElement(descriptor, 'LIBRARY_CONSTRUCTION_PROTOCOL').text = short
    node = experiment.find('DESIGN/DESIGN_DESCRIPTION')
    if node is None:
        node = ET.SubElement(experiment.find('DESIGN'), 'DESIGN_DESCRIPTION')
    node.text = design
    return records


@pytest.mark.parametrize('provider,parser', [('ena', ENAParser), ('sra', SRAParser)])
def test_native_parsers_combine_explicit_method_and_keep_assay_scope(provider, parser):
    data = parser().parse(records_with_method(provider)).to_mapping()
    protocol = data['series']['protocols'][0]
    assert protocol['description'] == 'Hi-C - Arima v2\n\n' + METHOD
    assert 'hardware' not in protocol
    assert data['sample'][0]['sra_run'][0]['instrument_model']
    assays = [s for p in data['series']['assay_paths'] for s in p['steps'] if s['kind'] == 'assay']
    assert all(s['description'] == CONTEXT + METHOD for s in assays)
    residual = json.dumps(data['extensions']['insdc'])
    assert 'LIBRARY_CONSTRUCTION_PROTOCOL' not in residual
    assert 'DESIGN_DESCRIPTION' not in residual


@pytest.mark.parametrize('short,design,expected', [
    (None, METHOD, METHOD),
    (None, CONTEXT, None),
    ('Original text', CONTEXT, 'Original text'),
    (METHOD, METHOD, METHOD),
    ('Short', METHOD + ' ' + METHOD, 'Short\n\n' + METHOD),
    ('Short', 'Library was prepared with Kit v2.1 at 1.5 ng input. Library was constructed using Kit B.',
     'Short\n\nLibrary was prepared with Kit v2.1 at 1.5 ng input.\n\nLibrary was constructed using Kit B.'),
    ('Short', 'Library was prepared with Kit A for sample SAMEA7520545.', 'Short'),
    ('Short', 'Library was made using Kit A for multiplexed sequencing run 34892_5.', 'Short'),
    ('Short', 'Library was made using Kit A with barcode TTATTATG.', 'Short'),
    ('Short', 'Library was made using Kit A with reads tagged with TTATTATG.', 'Short'),
    ('Short', 'No library was made using Kit A.', 'Short'),
    ('Short', 'Library was prepared with Kit A. Other experimental information.',
     'Short\n\nLibrary was prepared with Kit A.'),
    ('Short', 'Library was prepared with kits, e.g. Kit v2.1, at 1.5 ng input. Unrelated details.',
     'Short\n\nLibrary was prepared with kits, e.g. Kit v2.1, at 1.5 ng input.'),
    ('Original\ntext', 'Library was prepared with Kit v2.1\n at 1.5 ng input.',
     'Original\ntext\n\nLibrary was prepared with Kit v2.1\n at 1.5 ng input.'),
    ('Short', 'Library was made using Kit A for run 34892_5.', 'Short'),
    ('Short', 'Library was made using Kit A with index sequence ACGT.', 'Short'),
])
def test_explicit_sentence_contract(short, design, expected):
    from meta_standards_converter.miniml.archive_protocols import library_description
    assert library_description(short, design) == expected
    assert library_description(expected, design) == expected


def test_residuals_keep_unmapped_library_siblings():
    records = records_with_method('sra')
    descriptor = records.xml[0].find('.//LIBRARY_DESCRIPTOR')
    ET.SubElement(descriptor, 'SUBMITTER_NOTE').text = 'Unmapped annotation'
    data = SRAParser().parse(records).to_mapping()
    residual = json.dumps(data['extensions']['insdc'])
    assert 'Unmapped annotation' in residual
    assert 'LIBRARY_CONSTRUCTION_PROTOCOL' not in residual
