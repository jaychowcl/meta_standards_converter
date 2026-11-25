# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Regression coverage for archival findings from the new live metadata-integrity follow-up."""
import json
import csv
from pathlib import Path
import xml.etree.ElementTree as ET

import pytest

from meta_standards_converter.miniml.geo_parser import GEOParser
from meta_standards_converter.metadata.workflows import prepare_workflows
from meta_standards_converter.magetab.constructor import AEConstructor
from meta_standards_converter.magetab.parser import AEParser
from meta_standards_converter.sources.magetab import MAGETabInput, TextResource
from tests.support.live_integrity_semantics import organism_taxids

SOURCE = Path(__file__).parents[1] / 'fixtures/edge_cases/live-integrity-followup/GSE37703_family.xml'


def deposited_package():
    resources = [TextResource('E-GEOD-37703.' + suffix,
        SOURCE.with_name('E-GEOD-37703.' + suffix).read_text(), 'frozen-original')
        for suffix in ('idf.txt', 'sdrf.txt')]
    return AEParser().parse(MAGETabInput(resources[0], (resources[1],), 'E-GEOD-37703', 'local'))


def experiment_relations():
    root = ET.fromstring(SOURCE.read_bytes())
    for node in root.iter():
        node.tag = node.tag.rsplit('}', 1)[-1]
    return [(sample.get('iid'), relation.get('target')) for sample in root.findall('Sample')
            for relation in sample.findall('Relation') if relation.get('type') == 'SRA']


def test_live_geo_json_preserves_all_supplied_experiment_relations():
    package = GEOParser().parse(SOURCE.read_text())[0].to_mapping()
    expected = experiment_relations()
    assert len(expected) == 12
    actual = [(s['iid'], r['target']) for s in package['sample']
              for r in s.get('relation', []) if r['type'] == 'SRA']
    assert actual == expected


@pytest.mark.parametrize('route', ['direct', 'prepared', 'saved-json'])
def test_supplied_experiment_identifiers_reach_their_own_sdrf_assays(route):
    package = GEOParser().parse(SOURCE.read_text())[0]
    if route != 'direct':
        package = prepare_workflows(package)
    if route == 'saved-json':
        package = json.loads(json.dumps(package.to_mapping()))
    document = AEConstructor().miniml2magetab(package)
    sdrf = next(row[1] for row in document if row[0] == 'SDRF File')
    source_index = sdrf[0].index('Source Name')
    missing = []
    for sample, target in experiment_relations():
        accession = target.rsplit('=', 1)[-1]
        owned_rows = [row for row in sdrf[1:] if row[source_index] == sample]
        if not any(target in row or accession in row for row in owned_rows):
            missing.append((sample, accession))
    assert not missing, f'Source-supplied experiment identifiers absent from their own assays: {missing}'


def test_supplied_processing_definition_and_applications_are_retained():
    package = deposited_package().to_mapping()
    supplied = dict((r[0], r[1:]) for r in csv.reader(
        SOURCE.with_name('E-GEOD-37703.idf.txt').open(), delimiter='\t'))
    index = supplied['Protocol Type'].index('feature_extraction')
    name = supplied['Protocol Name'][index]
    protocol = next(p for p in package['series']['protocols'] if p['name'] == name)
    assert protocol['description'] == supplied['Protocol Description'][index]
    assert protocol['type']['value'] == 'feature_extraction'
    assert all(any(s.get('protocol_ref') == name for s in path['steps'])
               for path in package['series']['assay_paths'])


def test_feature_extraction_is_not_nucleic_acid_extraction():
    package = deposited_package().to_mapping()
    supplied = dict((r[0], r[1:]) for r in csv.reader(
        SOURCE.with_name('E-GEOD-37703.idf.txt').open(), delimiter='\t'))
    index = supplied['Protocol Type'].index('feature_extraction')
    computational = supplied['Protocol Description'][index]
    assert 'TopHat' in computational and 'Cuffdiff' in computational
    wrong = [s['iid'] for s in package['sample'] for c in s.get('channel', [])
             if c.get('extract_protocol') == computational]
    assert not wrong, f'Computational processing misclassified as wet-lab extraction: {wrong}'


@pytest.mark.parametrize('route', ['direct', 'prepared', 'saved-json'])
def test_geo_taxid_reaches_the_organism_characteristic(route):
    package = GEOParser().parse(SOURCE.read_text())[0]
    expected = {s['iid']: s['channel'][0]['organism'][0]['taxid'] for s in package.to_mapping()['sample']}
    assert len(expected) == 12 and set(expected.values()) == {'9606'}
    if route != 'direct': package = prepare_workflows(package)
    if route == 'saved-json': package = json.loads(json.dumps(package.to_mapping()))
    document = AEConstructor().miniml2magetab(package)
    sdrf = next(row[1] for row in document if row[0] == 'SDRF File')
    missing = [(row[sdrf[0].index('Source Name')], organism_taxids(sdrf[0],row))
               for row in sdrf[1:] if expected[row[sdrf[0].index('Source Name')]] not in organism_taxids(sdrf[0],row)]
    assert not missing, f'Taxonomy identifiers lost from their organism occurrences: {missing}'


@pytest.mark.parametrize('route', ['direct', 'prepared', 'saved-json'])
def test_multiple_experiments_and_conflicting_taxonomy_preserve_owned_evidence(route):
    from copy import deepcopy
    from meta_standards_converter.miniml import MINiMLCodec
    source = {'miniml_schema_version': '3.1', 'source': {'format': 'GEO'},
        'series': {'iid': 'GSE1', 'assay_paths': [{'steps': [
            {'kind': 'source', 'name': 'GSM1', 'sample_ref': 'GSM1'},
            {'kind': 'assay', 'name': 'SRX99'},
            {'kind': 'scan', 'name': 'SRR99'},
            {'kind': 'array_data_file', 'name': 'reads.fastq', 'link': {'value': 'https://example.org/reads.fastq', 'type': 'FASTQ'}}]}]},
        'sample': [{'iid': 'GSM1', 'library_strategy': 'RNA-Seq',
            'relation': [{'type': 'SRA', 'target': 'https://example.org/?term=SRX1'},
                         {'type': 'SRA', 'target': 'SRX2'}],
            'channel': [{'organism': [{'value': 'Homo sapiens', 'taxid': '9606'},
                {'value': 'Mus musculus', 'taxid': '10090', 'term_source_ref': 'NCBITaxon', 'term_accession_number': '999'}]}]}]}
    original = deepcopy(source)
    prepared = prepare_workflows(source)
    value = prepared.to_mapping()
    MINiMLCodec().decode(value, strict=True)
    assert source == original
    nodes = value['series']['assay_paths'][0]['steps']
    assert [n['name'] for n in nodes if n['kind'] == 'array_data_file'] == ['reads.fastq']
    node = next(n for n in nodes if n['kind'] == 'source')
    assert [c['value'] for c in node['comments'] if c['name'] == 'SRA_EXPERIMENT_RELATION'] == ['SRX1', 'SRX2']
    organisms = [c for c in node['characteristics'] if c['name'] == 'organism']
    assert [(c['value'], c['term_accession_number']) for c in organisms] == [
        ('Homo sapiens', '9606'), ('Mus musculus', '999'), ('Mus musculus', '10090')]
    assert any(d['iid'] == 'NCBITaxon' for d in value['database'])
    assert len(value['extensions']['msc_workflow_preparation']['warnings']) == 2
    if route == 'direct': value = source
    elif route == 'saved-json': value = json.loads(json.dumps(value))
    table = next(r[1] for r in AEConstructor().miniml2magetab(value) if r[0] == 'SDRF File')
    assert len(table) == 2
    assert {'SRX1', 'SRX2', 'SRX99'} <= set(table[1])
    assert {'9606', '999', '10090'} <= organism_taxids(table[0], table[1])


@pytest.mark.parametrize('kind,expected', [('feature_extraction', 'processing'),
    ('Feature-Extraction', 'processing'), ('nucleic acid extraction protocol', 'extraction'),
    ('unknown extraction analysis', None), ('novel library algorithm', None)])
def test_protocol_classifier_is_explicit_and_preparation_excludes_computation(kind, expected):
    from meta_standards_converter.metadata.protocol_types import protocol_type_family
    from meta_standards_converter.metadata.preparation import bound_protocols
    assert protocol_type_family(kind) == expected
    series = {'protocols': [{'name': 'P', 'type': {'value': kind}, 'description': 'authored'}],
        'assay_paths': [{'steps': [{'kind': 'source', 'name': 'S', 'sample_ref': 'S'},
            {'kind': 'protocol_application', 'protocol_ref': 'P'},
            {'kind': 'assay', 'name': 'SRX1'}, {'kind': 'scan', 'name': 'SRR1'}]}]}
    result = bound_protocols({'iid': 'S'}, {'run': 'SRR1', 'experiment': 'SRX1'}, series)
    assert bool(result) == (expected == 'extraction')


@pytest.mark.parametrize('accession', ['10090', 'NCBITaxon:10090', 'NCBITaxon_10090',
    'http://purl.obolibrary.org/obo/NCBITaxon_10090', 'https://purl.obolibrary.org/obo/NCBITaxon_10090'])
def test_equivalent_authored_taxonomy_identifier_is_not_a_conflict(accession):
    from meta_standards_converter.metadata.source_fidelity import organism_annotations
    organism = {'value': 'Mus musculus', 'taxid': '10090',
        'term_source_ref': 'NCBITaxon', 'term_accession_number': accession}
    warnings = []
    annotations = organism_annotations(organism, warnings)
    assert annotations == [{'name': 'organism', **{k:v for k,v in organism.items() if k != 'taxid'}}]
    assert not warnings
