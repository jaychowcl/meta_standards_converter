from copy import deepcopy
import json

import pytest

from meta_standards_converter.miniml.archive_entities import local_platforms
from meta_standards_converter.miniml.ena_parser import ENAParser
from meta_standards_converter.miniml.sra_parser import SRAParser
from tests.test_native_archive_parsers import fixture_records


@pytest.mark.parametrize('provider,parser', [('ena', ENAParser), ('sra', SRAParser)])
@pytest.mark.parametrize('family,model', [('ILLUMINA', 'HiSeq X Ten'), ('PACBIO_SMRT', 'Sequel II')])
def test_xml_platform_family_is_mapped_with_model_and_references(provider, parser, family, model):
    records = fixture_records(provider)
    experiment = next(e for root in records.xml for e in root.iter('EXPERIMENT'))
    platform = experiment.find('PLATFORM')[0]
    platform.tag = family
    platform.find('INSTRUMENT_MODEL').text = model
    data = parser().parse(records).to_mapping()
    run = data['sample'][0]['sra_run'][0]
    assert run['instrument_platform'] == family
    assert run['instrument_model'] == model
    p = next(p for p in data['platform'] if p['iid'] == run['platform_ref']['ref'])
    assert p['instrument_platform'] == family
    assert p['instrument_model'] == model
    assert data['sample'][0]['platform_ref'] == run['platform_ref']


def test_ena_indexed_family_fills_missing_xml_evidence_without_indexed_wrapper():
    records = fixture_records('ena')
    records.xml = [root for root in records.xml if root.tag != 'EXPERIMENT_SET']
    records.indexed['read_experiment'] = [{'experiment_accession': 'SRX7812918',
        'study_accession': 'SRP250911', 'sample_accession': 'SAMN14218700',
        'instrument_model': 'Illumina MiSeq', 'instrument_platform': 'ILLUMINA'}]
    data = ENAParser().parse(records).to_mapping()
    assert data['sample'][0]['sra_run'][0]['instrument_platform'] == 'ILLUMINA'
    assert all(p.get('instrument_platform') != 'INDEXED' for p in data['platform'])
    residual = [r for r in data['extensions']['insdc']['records'] if r['kind'] == 'read_experiment']
    assert all('instrument_platform' not in r['metadata'] for r in residual)


def test_conflicting_indexed_family_remains_residual_not_overwriting_xml():
    records = fixture_records('ena')
    records.indexed['read_run'][0]['instrument_platform'] = 'PACBIO_SMRT'
    data = ENAParser().parse(records).to_mapping()
    assert data['sample'][0]['sra_run'][0]['instrument_platform'] == 'ILLUMINA'
    assert any(r['metadata'].get('instrument_platform') == 'PACBIO_SMRT'
               for r in data['extensions']['insdc']['records'] if r['kind'] == 'read_run')


def test_compatible_existing_local_platform_is_completed_without_new_id():
    data = {'source': {'format': 'SRA'}, 'platform': [{'iid': 'sra:platform:existing',
        'instrument_model': 'Model', 'title': 'Model', 'technology': 'high-throughput sequencing'}],
        'sample': [{'iid': 'SRS1', 'sra_run': [{'run': 'SRR1', 'experiment': 'SRX1',
            'instrument_model': 'Model', 'instrument_platform': 'ILLUMINA',
            'platform_ref': {'ref': 'sra:platform:existing'}}]}]}
    local_platforms(data)
    assert len(data['platform']) == 1
    assert data['platform'][0]['instrument_platform'] == 'ILLUMINA'
    assert data['sample'][0]['platform_ref']['ref'] == 'sra:platform:existing'
    before = deepcopy(data)
    local_platforms(data)
    assert data == before


def test_family_only_and_conflicting_platforms_keep_mixed_scope_and_gpl():
    data = {'source': {'format': 'ENA'}, 'platform': [{'iid': 'GPL1', 'title': 'Registered'}],
        'sample': [{'iid': 'SAM1', 'platform_ref': {'ref': 'GPL1'}, 'sra_run': [
            {'run': 'ERR1', 'experiment': 'ERX1', 'instrument_platform': 'ILLUMINA'},
            {'run': 'ERR2', 'experiment': 'ERX2', 'instrument_platform': 'PACBIO_SMRT'}]}]}
    local_platforms(data)
    assert len(data['platform']) == 3
    assert not any(p.get('instrument_model') for p in data['platform'])
    assert data['sample'][0]['platform_ref']['ref'] == 'GPL1'
    assert len({r['platform_ref']['ref'] for r in data['sample'][0]['sra_run']}) == 2
    data['sample'][0].pop('platform_ref')
    local_platforms(data)
    assert not data['sample'][0].get('platform_ref')


def test_unmapped_xml_family_survives_when_model_is_mapped():
    records = fixture_records('sra')
    package = SRAParser().parse(records)
    data = package.to_mapping()
    from meta_standards_converter.miniml.archive_residuals import finalize, source_records
    for sample in data['sample']:
        for run in sample['sra_run']:
            run['instrument_platform'] = 'OTHER'
    result = finalize(data, source_records(package)).to_mapping()
    records = [r for r in result['extensions']['insdc']['records'] if r['kind'] == 'EXPERIMENT']
    assert any('ILLUMINA' in json.dumps(r['metadata']) for r in records)


def test_saved_residual_family_is_recovered_only_for_explicit_experiment():
    data = {'source': {'format': 'ENA'}, 'sample': [{'iid': 'SAM1', 'sra_run': [
        {'run': 'ERR1', 'experiment': 'ERX1', 'instrument_model': 'Model'},
        {'run': 'ERR2', 'experiment': 'ERX2', 'instrument_model': 'Model'}]}],
        'extensions': {'insdc': {'version': '2.0', 'records': [
            {'provider': 'ena', 'kind': 'read_experiment', 'accession': 'ERX1',
             'metadata': {'instrument_platform': 'ILLUMINA'}}]}}}
    local_platforms(data)
    one, two = data['sample'][0]['sra_run']
    assert one['instrument_platform'] == 'ILLUMINA'
    assert not two.get('instrument_platform')
    assert one['platform_ref'] != two['platform_ref']


def test_two_explicit_families_cannot_mutate_one_existing_model_declaration():
    data = {'source': {'format': 'SRA'}, 'platform': [{'iid': 'sra:platform:old',
        'title': 'Model', 'instrument_model': 'Model'}], 'sample': [{'iid': 'SRS1', 'sra_run': [
            {'run': 'SRR1', 'instrument_model': 'Model', 'instrument_platform': 'ILLUMINA', 'platform_ref': {'ref': 'sra:platform:old'}},
            {'run': 'SRR2', 'instrument_model': 'Model', 'instrument_platform': 'OTHER', 'platform_ref': {'ref': 'sra:platform:old'}}]}]}
    local_platforms(data)
    assert not data['platform'][0].get('instrument_platform')
    assert len({r['platform_ref']['ref'] for r in data['sample'][0]['sra_run']}) == 2
    before = deepcopy(data)
    local_platforms(data)
    assert data == before
