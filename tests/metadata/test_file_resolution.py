# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""File interpretation does not conflate format, provenance and read roles."""
from copy import deepcopy
import pytest

from meta_standards_converter.metadata.file_resolution import (
    classify, combine, download_assets, ena_read_role, structural_uri,
)


@pytest.mark.parametrize('value, expected', [
    ('\nhttps://example.org/a.CEL.gz\n  ', 'array_raw'),
    ({'value': 'https://example.org/opaque.bin', 'type': 'CEL'}, 'array_raw'),
    ({'value': 'https://example.org/a.fastq.gz', 'role': 'processed'}, 'matrix_or_derived'),
    ('https://example.org/a.txt', 'supplementary'),
    ('https://example.org/a.unknown', 'supplementary'),
    ('https://example.org/a.fq.gz', 'sequencing_raw'),
])
def test_explicit_file_evidence_precedes_extension(value, expected):
    assert classify(value) == expected


def test_uri_framing_does_not_rewrite_literal_path():
    assert structural_uri(' \nftp://host/path%23one/read 1.CEL.gz\t') == 'ftp://host/path%23one/read 1.CEL.gz'


@pytest.mark.parametrize('name, expected', [
    ('SRR1_1.fastq.gz', {'read_role': 'paired', 'read_index': '1'}),
    ('SRR1_2.fastq.gz', {'read_role': 'paired', 'read_index': '2'}),
    ('SRR1.fastq.gz', {'read_role': 'unpaired'}),
    ('AA_ATAC_S1_R2.fastq.gz', {}),
])
def test_read_roles_require_archive_origin_and_exact_run(name, expected):
    file = {'filename': name, 'role': 'GENERATED_FILE'}
    assert ena_read_role(file, 'SRR1', 'PAIRED') == expected
    assert ena_read_role(dict(file, role='Original'), 'SRR1', 'PAIRED') == {}
    assert ena_read_role(file, 'SRR2', 'PAIRED') == {}
    assert ena_read_role(file, 'SRR1', None) == {}


def test_download_view_preserves_occurrences_and_conflicts():
    record = {'uri': 'ftp://host/a.fastq', 'md5': 'a' * 32, 'bytes': '100'}
    records = [record, deepcopy(record), dict(record, bytes='200')]
    before = deepcopy(records)
    assets = download_assets(records)
    assert [a.occurrences for a in assets] == [(0, 1), (2,)]
    assert records == before


def test_filename_alone_does_not_join_distinct_locations():
    records = [{'NAME': 'same.fastq', 'URI': 'ftp://one/a', 'FORMAT': 'fastq'}]
    combine(records, {'NAME': 'same.fastq', 'URI': 'ftp://two/a', 'FORMAT': 'fastq'})
    assert len(records) == 2


def test_sparse_descriptor_cannot_bridge_incompatible_files():
    from meta_standards_converter.magetab.native_file_selection import _coalesced
    a = {'NAME': 'same.fastq', 'URI': 'ftp://one/a', 'FORMAT': 'FASTQ', 'MD5': 'a' * 32}
    b = {'NAME': 'same.fastq', 'FORMAT': 'fastq'}
    c = dict(a, MD5='b' * 32)
    assert len(_coalesced([a, b, c], combine)) == 3


def test_unknown_sample_link_retains_owner_in_export():
    from tests.magetab.test_metadata_integrity import geo, tables
    from meta_standards_converter.miniml import MINiMLCodec
    data = geo('GSE31356').to_mapping()
    uri = 'https://example.org/unclassified.txt'
    data['sample'][0]['supplementary_data'] = [{'value': uri, 'type': 'TXT'}]
    _, sdrf = tables(MINiMLCodec().decode(data).package)
    assert sdrf[1][sdrf[0].index('Comment[Unclassified File]')] == uri
    assert sdrf[1][sdrf[0].index('Source Name')] == data['sample'][0]['iid']


def test_file_occurrence_is_not_a_read_role():
    from meta_standards_converter.magetab.sdrf.constructor import SDRFConstructor
    from tests.magetab.test_metadata_integrity import SOURCE, geo
    data = geo('GSE52778').to_mapping()
    handler = SDRFConstructor().create_handler(data, technology_type='sequencing')
    supplied = SOURCE['fastq']
    attrs = handler.sequencing_file_attrs(data['sample'][0], {'fastq_files': [
        {'filename': supplied['NAME'], 'uri': supplied['URI'], 'md5': supplied['MD5']}]})
    assert not any(a.label == 'Comment[read1 file]' for a in attrs)
    assert any(a.value == supplied['NAME'] for a in attrs)
