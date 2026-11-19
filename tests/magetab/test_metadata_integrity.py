# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Integrity regressions derived from independently captured repository records."""
from copy import deepcopy
import json
from pathlib import Path
import xml.etree.ElementTree as ET

import pytest

from meta_standards_converter.miniml.geo_parser import GEOParser
from meta_standards_converter.magetab.constructor import AEConstructor
from meta_standards_converter.sources.insdc import INSDCWebfetcher

SOURCE = json.loads((Path(__file__).parents[1] / 'fixtures/edge_cases/metadata-integrity/source-cases.json').read_text())


def geo(study):
    return GEOParser().parse(SOURCE['geo'][study])[0]


def tables(package):
    idf = AEConstructor().miniml2magetab(package)
    return idf, next(r[1] for r in idf if r[0] == 'SDRF File')


@pytest.mark.xfail(strict=True, raises=AssertionError, reason='MSC integrity: channel material identity collapses')
@pytest.mark.parametrize('column', ['Source Name', 'Extract Name'])
def test_two_channel_materials_have_distinct_identities(column):
    package = geo('GSE13398')
    assert len(package.samples[0].channels) == 2
    _, sdrf = tables(package)
    assert len(sdrf) == 3
    assert len({r[sdrf[0].index(column)] for r in sdrf[1:]}) == 2
    assert len({r[sdrf[0].index('Assay Name')] for r in sdrf[1:]}) == 1


def test_whitespace_bearing_cel_link_is_exported_as_raw_array_file():
    package = geo('GSE31356')
    source = ET.fromstring(SOURCE['geo']['GSE31356']).find('Sample/Supplementary-Data')
    assert source.get('type') == 'CEL'
    assert source.text != source.text.strip()
    _, sdrf = tables(package)
    assert 'Array Data File' in sdrf[0]
    assert sdrf[1][sdrf[0].index('Array Data File')] == source.text.strip()
    assert package.samples[0].supplementary_data[0].value == source.text


@pytest.mark.xfail(strict=True, raises=AssertionError, reason='MSC integrity: GEO treatment precedes growth')
def test_geo_growth_precedes_drug_exposure_without_reordering_authored_paths():
    package = geo('GSE52778')
    source = ET.fromstring(SOURCE['geo']['GSE52778']).find('Sample/Channel')
    growth, treatment = source.findtext('Growth-Protocol'), source.findtext('Treatment-Protocol')
    assert 'passage' in growth.lower() and '18' in treatment
    idf, sdrf = tables(package)
    refs = next(r[1:] for r in idf if r[0] == 'Protocol Name')
    descriptions = dict(zip(refs, next(r[1:] for r in idf if r[0] == 'Protocol Description')))
    ordered = [' '.join(descriptions.get(sdrf[1][i], '').split()) for i, h in enumerate(sdrf[0]) if h == 'Protocol REF']
    assert ordered.index(' '.join(growth.split())) < ordered.index(' '.join(treatment.split()))


def test_generated_fastq_is_not_presented_as_submitted_file(monkeypatch):
    fetcher = INSDCWebfetcher()
    source = ET.fromstring(SOURCE['sra'])
    originals = {n.get('filename') for n in source.findall('.//SRAFile') if n.get('supertype') == 'Original'}
    monkeypatch.setattr(fetcher, 'fetch_sra_xml', lambda **_: source)
    monkeypatch.setattr(fetcher, 'fetch_ena_file_report', lambda **_: SOURCE['ena_report'])
    record = fetcher.fetch_sra_runs('SRP033351')[0]
    assert originals == {'LL01_R1.fastq', 'LL01_R2.fastq'}
    assert not record.get('submitted_file_name') or record['submitted_file_name'] in originals
    assert originals <= {f['filename'] for f in record.get('files', [])}
    assert record['fastq_files'][0]['filename'] not in originals


def test_equivalent_file_formats_complete_primary_checksum_and_size():
    from meta_standards_converter.magetab.native_files import _combine
    source = SOURCE['fastq']
    primary = [{k: source[k] for k in ('URI', 'NAME')} | {'FORMAT': 'FASTQ'}]
    _combine(primary, source)
    assert len(primary) == 1
    assert primary[0]['MD5'] == source['MD5']
    assert primary[0]['BYTES'] == source['BYTES']


def test_conflicting_file_checksums_do_not_merge():
    from meta_standards_converter.magetab.native_files import _combine
    source = SOURCE['fastq']
    primary = [dict(source, MD5='0' * 32)]
    _combine(primary, source)
    assert len(primary) == 2


@pytest.mark.xfail(strict=True, raises=AssertionError, reason='MSC integrity: peer paths retain unselected taxonomy')
def test_peer_generated_taxonomy_uses_selected_sample_facts():
    from meta_standards_converter.metadata.archive_workflows import merge_workflows, _biological_node
    def sample(taxid, name):
        return {'iid': 'S1', 'channel': [{'organism': [{'value': name, 'taxid': taxid}]}],
                'sra_run': [{'run': 'ERR1', 'experiment': 'ERX1'}]}
    selected = sample('83332', 'Mycobacterium tuberculosis H37Rv')
    peer = sample('1773', 'Mycobacterium tuberculosis')
    def path(record):
        return {'steps': [_biological_node(record), {'kind': 'assay', 'name': 'ERX1'}, {'kind': 'scan', 'name': 'ERR1'}]}
    data = {'source': {'format': 'SRA'}, 'sample': [selected], 'series': {'iid': 'ERP185509', 'assay_paths': [path(selected)]}}
    extra = {'source': {'format': 'ENA'}, 'sample': [peer], 'series': {'iid': 'ERP185509', 'assay_paths': [path(peer)]}}
    original = deepcopy(extra)
    issues = []
    merge_workflows(data, extra, {'S1': 'S1'}, {}, False, issues, original_channels={'S1': selected['channel']})
    for row in data['series']['assay_paths']:
        assert row['steps'][0]['characteristics'] == _biological_node(selected)['characteristics']
    assert extra == original


@pytest.mark.xfail(strict=True, raises=AssertionError, reason='MSC integrity: channel identity depends on label')
@pytest.mark.parametrize('labels', [('Cy3', 'Cy3'), ('', '')])
def test_channel_occurrence_identity_does_not_depend_on_label(labels):
    from meta_standards_converter.miniml import MINiMLCodec
    value = geo('GSE13398').to_mapping()
    for channel, label in zip(value['sample'][0]['channel'], labels):
        channel['label'] = label
    _, sdrf = tables(MINiMLCodec().decode(value).package)
    assert len({r[sdrf[0].index('Source Name')] for r in sdrf[1:]}) == 2
