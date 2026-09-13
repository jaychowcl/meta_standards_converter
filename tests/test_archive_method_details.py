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
from copy import deepcopy
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
    ('Short', 'Library was made using Kit A for sample ABC.', 'Short'),
    ('Short', 'Library was made using Kit A with index ACGTTGCA.', 'Short'),
    (None, 'The libraries were constructed with Kit B v2.1.', 'The libraries were constructed with Kit B v2.1.'),
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


@pytest.mark.parametrize('provider,parser', [('ena', ENAParser), ('sra', SRAParser)])
def test_parser_creates_protocol_only_from_explicit_preparation(provider, parser):
    data = parser().parse(records_with_method(provider, short=None)).to_mapping()
    assert data['series']['protocols'][0]['description'] == METHOD
    assert all(any(s.get('protocol_ref') for s in p['steps']) for p in data['series']['assay_paths'])
    empty = parser().parse(records_with_method(provider, short=None, design=CONTEXT)).to_mapping()
    assert not empty['series'].get('protocols')


def legacy_package():
    data = SRAParser().parse(records_with_method('sra')).to_mapping()
    protocol = data['series']['protocols'][0]
    protocol['description'] = 'Hi-C - Arima v2'
    protocol['hardware'] = [data['sample'][0]['sra_run'][0]['instrument_model']]
    # Avoid unrelated publication lookup in this export contract.
    data['series'].pop('pubmed_id', None)
    data['series'].pop('pubmed_publication', None)
    return data


def test_saved_native_json_export_upgrades_only_in_memory(tmp_path, monkeypatch):
    import requests
    from meta_standards_converter.converters.json2ae import JSON2AEConverter
    from meta_standards_converter.magetab.protocol_export import prepare_protocols
    def forbidden(*args, **kwargs):
        raise AssertionError('no remote lookup allowed')
    monkeypatch.setattr(requests.sessions.Session, 'request', forbidden)
    data = legacy_package()
    original = json.dumps(data)
    path = tmp_path / 'native.json'
    path.write_text(original)
    output = JSON2AEConverter().convert(str(path), enrich=False)[0]
    rows = {r[0]: r[1:] for r in output}
    assert rows['Protocol Description'] == ['Hi-C - Arima v2\n\n' + METHOD]
    assert not any(rows['Protocol Hardware'])
    sdrf = rows['SDRF File'][0]
    header = sdrf[0]
    assert all(r[header.index('Protocol REF')] == rows['Protocol Name'][0] for r in sdrf[1:])
    description_column = header.index('Description', header.index('Assay Name'))
    assert all(r[description_column] == CONTEXT + METHOD for r in sdrf[1:])
    assert all(r[header.index('Comment[INSTRUMENT_MODEL]')] for r in sdrf[1:])
    assert path.read_text() == original
    prepare_protocols(data)
    once = deepcopy(data)
    prepare_protocols(data)
    assert data == once


@pytest.mark.parametrize('change', ['ambiguous_description', 'wrong_assay', 'unreferenced', 'registered', 'geo'])
def test_saved_metadata_skips_unverified_or_enriched_protocols(change):
    from meta_standards_converter.magetab.protocol_export import prepare_protocols
    data = legacy_package()
    protocol = data['series']['protocols'][0]
    if change == 'ambiguous_description':
        other = deepcopy(data['series']['assay_paths'][0])
        next(s for s in other['steps'] if s['kind'] == 'assay')['description'] = 'Library was prepared with conflicting kit.'
        data['series']['assay_paths'].append(other)
    elif change == 'wrong_assay':
        for p in data['series']['assay_paths']:
            next(s for s in p['steps'] if s['kind'] == 'assay')['name'] = 'SRX999999'
    elif change == 'unreferenced':
        for p in data['series']['assay_paths']:
            p['steps'] = [s for s in p['steps'] if s['kind'] != 'protocol_application']
    elif change == 'registered':
        old = protocol['name']; protocol['name'] = 'E-MTAB-1:P-MTAB-123'
        for p in data['series']['assay_paths']:
            for s in p['steps']:
                if s.get('protocol_ref') == old: s['protocol_ref'] = protocol['name']
    else:
        data['source']['format'] = 'GEO'
    prepare_protocols(data)
    assert data['series']['protocols'][0]['description'] == 'Hi-C - Arima v2'


def test_shared_protocol_reference_with_different_assays_is_not_augmented():
    from meta_standards_converter.magetab.protocol_export import prepare_protocols
    data = legacy_package()
    other = deepcopy(data['series']['assay_paths'][0])
    next(s for s in other['steps'] if s['kind'] == 'assay')['name'] = 'SRX999999'
    data['series']['assay_paths'].append(other)
    prepare_protocols(data)
    assert data['series']['protocols'][0]['description'] == 'Hi-C - Arima v2'


def test_preserve_explicit_non_sequencing_hardware():
    from meta_standards_converter.magetab.protocol_export import prepare_protocols
    data = legacy_package()
    data['series']['protocols'][0]['hardware'] = ['Liquid handling robot']
    prepare_protocols(data)
    assert data['series']['protocols'][0]['hardware'] == ['Liquid handling robot']
    assert data['series']['protocols'][0]['description'].endswith(METHOD)


def test_distinct_methods_keep_distinct_export_names_and_references():
    from meta_standards_converter.magetab.protocol_export import prepare_protocols
    data = legacy_package()
    original = data['series']['protocols'][0]
    old_experiment = original['name'].split(':')[0]
    second = deepcopy(original); second['name'] = 'SRX999999:library'
    data['series']['protocols'].append(second)
    # Replicate an explicitly bound experiment with different preparation evidence.
    second_run = deepcopy(data['sample'][0]['sra_run'][0])
    second_run.update(experiment='SRX999999', run='SRR999999')
    data['sample'][0]['sra_run'].append(second_run)
    paths = deepcopy(data['series']['assay_paths'])
    for path in paths:
        for step in path['steps']:
            if step.get('protocol_ref') == original['name']: step['protocol_ref'] = second['name']
            if step.get('kind') == 'assay':
                step['name'] = 'SRX999999'
                step['description'] = 'Library was made using a different kit.'
    data['series']['assay_paths'].extend(paths)
    prepare_protocols(data)
    protocols = {p['name']: p for p in data['series']['protocols']}
    assert len(protocols) == 2
    for path in data['series']['assay_paths']:
        assay = next(s for s in path['steps'] if s['kind'] == 'assay')
        ref = next(s['protocol_ref'] for s in path['steps'] if s['kind'] == 'protocol_application')
        expected = METHOD if assay['name'] == old_experiment else 'Library was made using a different kit.'
        assert protocols[ref]['description'].endswith(expected)


def test_export_deduplicates_equivalent_methods_across_different_instruments():
    from meta_standards_converter.magetab.protocol_export import prepare_protocols
    data = legacy_package()
    first = data['series']['protocols'][0]
    second = deepcopy(first); second.update(name='SRX999999:library', hardware=['Another sequencer'])
    data['series']['protocols'].append(second)
    run = deepcopy(data['sample'][0]['sra_run'][0])
    run.update(experiment='SRX999999', run='SRR999999', instrument_model='Another sequencer')
    data['sample'][0]['sra_run'].append(run)
    path = deepcopy(data['series']['assay_paths'][0])
    for step in path['steps']:
        if step.get('protocol_ref') == first['name']: step['protocol_ref'] = second['name']
        if step['kind'] == 'assay': step['name'] = 'SRX999999'
    data['series']['assay_paths'].append(path)
    prepare_protocols(data)
    assert len(data['series']['protocols']) == 1
    assert data['series']['protocols'][0]['description'].endswith(METHOD)
    assert 'hardware' not in data['series']['protocols'][0]
