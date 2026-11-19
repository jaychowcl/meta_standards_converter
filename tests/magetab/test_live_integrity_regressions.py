# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Unfixed archival findings from the new live metadata-integrity follow-up."""
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


@pytest.mark.xfail(strict=True, reason='MSC-LIVE-001: supplied sample experiment relations need network-free SDRF projection')
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


@pytest.mark.xfail(strict=True, reason='MSC-LIVE-002: feature_extraction must not project into wet-lab extract_protocol')
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
