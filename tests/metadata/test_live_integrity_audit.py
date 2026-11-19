# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Timing labels must not turn reused responses or failures into successful live runs."""
from scripts.audit_live_integrity import request_counts, summary, COHORT
from tests.support.live_integrity_semantics import compare_rows, organism_taxids
from scripts.audit_live_integrity_review import exported_files


def test_request_reuse_and_repeated_attempts_are_distinct():
    event = {'identity': {'url': 'https://example.org/x'}, 'physical': True, 'attempt': {'status': 503}}
    result = request_counts([event, dict(event, attempt={'status': 200}), dict(event, physical=False)])
    assert result == {'physical_requests': 2, 'reused_responses': 1,
        'repeated_physical_identities': 1, 'response_statuses': {'503': 1, '200': 1}}


def test_timing_statistics_keep_failed_outcomes_and_detect_variation():
    rows = [{'wall_seconds': wall, 'cpu_seconds': 1, 'peak_rss_bytes': 3,
        'payload_sha256': 'same', 'items': [{'execution': 'failed'}]} for wall in (3, 1, 2)]
    result = summary(rows)
    assert (result['median_seconds'], result['min_seconds'], result['max_seconds']) == (2, 1, 3)
    assert result['outcomes_identical']
    rows[-1]['items'] = [{'execution': 'succeeded'}]
    assert not summary(rows)['outcomes_identical']


def test_live_cohort_has_five_new_provider_inputs():
    assert len(COHORT['inputs']) == 5
    assert {s['evidence'] for s in COHORT['inputs']} == {'live'}
    assert len({s['reference'] for s in COHORT['inputs']}) == 3


def test_sparse_protocol_slots_and_header_aliases_are_representation():
    source = [['Source Name', 'Protocol REF', 'Term Source REF', 'Protocol REF', 'Term Source REF',
        'Extract Name', 'Characteristics [treatment]', 'FactorValue [TREATMENT]'],
        ['S1', 'P1', 'EFO', '  ', '  ', 'E1', 'drug', 'drug'],
        ['S2', '  ', '  ', 'P2', 'EFO', 'E2', 'control', 'control']]
    output = [['Source Name', 'Protocol REF', 'Term Source REF', 'Extract Name',
        'Characteristics[treatment]', 'Factor Value[TREATMENT]'],
        ['S2', 'P2', 'EFO', 'E2', 'control', 'control'],
        ['S1', 'P1', 'EFO', 'E1', 'drug', 'drug']]
    assert compare_rows(source, output)['equivalent']


def test_repeated_applications_and_blank_owners_cannot_disappear():
    source = [['Source Name', 'Protocol REF', 'Protocol REF', 'Extract Name',
        'Characteristics[exposure]', 'Comment[evidence]'], ['S1', 'P1', 'P1', 'E1', '', 'known']]
    lost_protocol = [[*source[0][:2], *source[0][3:]], [*source[1][:2], *source[1][3:]]]
    lost_owner = [source[0][:-2], source[1][:-2]]
    assert not compare_rows(source, lost_protocol)['equivalent']
    assert not compare_rows(source, lost_owner)['equivalent']


def test_qualifier_swap_changes_semantics_and_duplicate_rows_count():
    source = [['Source Name', 'Characteristics[a]', 'Term Source REF', 'Characteristics[b]', 'Term Source REF'],
        ['S1', 'x', 'EFO', 'y', 'NCBITaxon']]
    swapped = [source[0], ['S1', 'x', 'NCBITaxon', 'y', 'EFO']]
    assert not compare_rows(source, swapped)['equivalent']
    assert not compare_rows([*source, source[1]], source)['equivalent']


def test_repeated_qualifiers_retain_their_ordinal_binding():
    headers = ['Source Name', 'Characteristics[a]', 'Term Source REF', 'Term Accession Number',
        'Term Source REF', 'Term Accession Number']
    source = [headers, ['S1', 'x', 'EFO', 'one', 'OBI', 'two']]
    swapped = [headers, ['S1', 'x', 'EFO', 'two', 'OBI', 'one']]
    assert not compare_rows(source, swapped)['equivalent']


def test_fastq_uri_retains_basename_and_its_owned_checksum_only():
    table = [['Assay Name','Comment[MD5]','Scan Name','Comment[ENA_RUN]',
        'Comment[FASTQ_URI]','Comment[MD5]','Comment[FASTQ_BYTES]'],
        ['A','unrelated','SRR1','SRR1','ftp://example.org/read1.fastq.gz','right','100']]
    result = exported_files(table)['SRR1','ftp://example.org/read1.fastq.gz'][0]
    assert result == {'URI':'ftp://example.org/read1.fastq.gz','MD5':'right','BYTES':'100','NAME':'read1.fastq.gz'}
    table[0].pop(5); table[1].pop(5)
    assert 'MD5' not in exported_files(table)['SRR1','ftp://example.org/read1.fastq.gz'][0]


def test_taxonomy_identifier_cannot_be_borrowed_from_another_characteristic():
    headers=['Source Name','Characteristics[organism]','Characteristics[cell line]',
             'Term Source REF','Term Accession Number']
    row=['S','Homo sapiens','cells','NCBITaxon','9606']
    assert not organism_taxids(headers,row)
    headers.pop(2);row.pop(2)
    assert organism_taxids(headers,row)=={'9606'}
