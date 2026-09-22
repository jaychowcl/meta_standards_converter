# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Identity-based platform reconciliation must survive strict peer boundaries."""
from copy import deepcopy

import pytest

from meta_standards_converter.metadata.archive_workflows import merge_declarations
from meta_standards_converter.metadata.archive_enrichment import merge_archive_metadata
from meta_standards_converter.miniml import MINiMLCodec
from tests.test_native_archive_enrichment import native


def platform(iid):
    return {'iid': iid, 'title': 'NextSeq 500', 'instrument_model': 'NextSeq 500',
            'instrument_platform': 'ILLUMINA', 'technology': 'high-throughput sequencing'}


def test_equivalent_native_platforms_merge_and_remap_samples():
    data = native().to_mapping()
    data['platform'][0]['iid'] = 'sra:platform:one'
    for sample in data['sample']: sample['platform_ref'] = {'ref': 'sra:platform:one'}
    extra = deepcopy(data)
    extra['source']['format'] = 'ENA'
    extra['platform'][0]['iid'] = 'ena:platform:two'
    extra['platform'][0]['description'] = 'Sequencing instrument'
    for sample in extra['sample']: sample['platform_ref'] = {'ref': 'ena:platform:two'}
    codec = MINiMLCodec()
    result, issues = merge_archive_metadata(codec.decode(data, strict=True).package,
                                           codec.decode(extra, strict=True).package)
    merged = codec.decode(result, strict=True).package.to_mapping()
    assert not issues
    assert len(merged['platform']) == 1
    assert merged['platform'][0]['description'] == 'Sequencing instrument'
    assert merged['sample'][0]['platform_ref'] == {'ref': 'sra:platform:one'}


@pytest.mark.parametrize('different', ['accession', 'model', 'missing_model'])
def test_same_title_is_not_platform_identity(different):
    a, b = platform('a'), platform('b')
    if different == 'accession':
        a['accession'] = [{'database': 'GEO', 'value': 'GPL1'}]
        b['accession'] = [{'database': 'GEO', 'value': 'GPL2'}]
    elif different == 'model': b['instrument_model'] = 'MiSeq'
    else: b.pop('instrument_model')
    data, extra = {'platform': [a]}, {'platform': [b]}
    merge_declarations(data, extra, 'peer:')
    assert len(data['platform']) == 2


def test_verified_platform_accession_reuses_id_with_complementary_fields():
    a = {'iid': 'one', 'accession': [{'database': 'GEO', 'value': 'GPL221'}], 'technology': 'other'}
    b = {**deepcopy(a), 'iid': 'two', 'description': 'Deposited design'}
    data, extra = {'platform': [a]}, {'platform': [b], 'sample': [{'platform_ref': {'ref': 'two'}}]}
    mappings = merge_declarations(data, extra, 'peer:')
    assert mappings['two'] == 'one'
    assert len(data['platform']) == 1
    assert data['platform'][0]['description'] == 'Deposited design'
    assert extra['sample'][0]['platform_ref'] == {'ref': 'one'}
