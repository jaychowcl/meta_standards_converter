"""Negative controls for full-field and source/reference distinction."""
from copy import deepcopy
import pytest
from tests.support.authorship_fidelity import (
    EQUIVALENT, LOSS, UNSUPPORTED, CURATION, ENRICHMENT, UNRESOLVED,
    assert_idf_preserved, idf_cells, classify, align_sdrf, pair_sdrf,
)


def idf():
    return [['PubMed ID', '1', '', '3'], ['Publication DOI', '10.1/a', '10.1/b', '10.1/c'],
            ['Publication Author List', 'Smith AB, Börsch A', 'Example Consortium', 'Jones A Jr'],
            ['Person First Name', 'Anna', 'Anna'], ['Person Last Name', 'Smith', 'Smith'],
            ['Person Email', 'a@example.org', 'b@example.org'], ['Person Affiliation', 'A', 'B'],
            ['Person Roles', 'submitter', 'investigator'],
            ['Protocol Name', 'P1', 'P2'], ['Protocol Description', 'extract', 'sequence'],
            ['Protocol Term Source REF', 'EFO', 'OBI'],
            ['Experimental Factor Name', 'dose'], ['Experimental Factor Type', 'compound dose'],
            ['Experimental Factor Term Source REF', 'EFO'], ['Comment[author supplied]', 'verbatim']]


@pytest.mark.parametrize('field', ['Publication DOI', 'Publication Author List', 'Person Email',
    'Person Affiliation', 'Person Roles', 'Protocol Description', 'Protocol Term Source REF',
    'Experimental Factor Type', 'Experimental Factor Term Source REF', 'Comment[author supplied]'])
def test_every_idf_field_is_checked_and_identity_swaps_fail(field):
    source, actual = idf(), idf()
    row = next(row for row in actual if row[0] == field)
    if len(row) > 2:
        row[1], row[2] = row[2], row[1]
    else:
        row[1] = 'altered'
    with pytest.raises(AssertionError): assert_idf_preserved(source, actual)


def test_whole_publication_groups_can_reorder_but_internal_blank_is_not_shifted():
    source, actual = idf(), idf()
    for row in actual:
        if row[0].startswith('Publication ') or row[0] == 'PubMed ID':
            row[1:] = [row[3], row[1], row[2]]
    actual.reverse()
    assert_idf_preserved(source, actual)
    next(row for row in actual if row[0] == 'PubMed ID')[2:] = ['', '1']
    with pytest.raises(AssertionError): assert_idf_preserved(source, actual)


def test_all_blank_publication_and_trailing_padding_are_equivalent():
    source = [['PubMed ID', '', ''], ['Investigation Title', 'Study', '']]
    actual = [['Investigation Title', 'Study'], ['Publication Author List', None]]
    assert_idf_preserved(source, actual)


def test_repeated_protocol_definitions_and_contact_occurrences_are_not_sets():
    source = [['Protocol Name', 'P1', 'P1'], ['Protocol Description', 'one', 'two'],
              ['Person First Name', 'Anna', 'Anna'], ['Person Last Name', 'Smith', 'Smith']]
    actual = [row[:2] for row in source]
    with pytest.raises(AssertionError): assert_idf_preserved(source, actual)


def test_unknown_idf_annotation_is_not_silently_ignored():
    with pytest.raises(AssertionError):
        assert_idf_preserved([['Comment[unusual]', 'curator annotation']], [])


def test_transformation_requires_exact_cell_identity_and_both_values():
    source = [['Publication DOI', ' 10.1/a ', '10.1/b']]
    actual = [['Publication DOI', '10.1/a', '10.1/b']]
    with pytest.raises(AssertionError): assert_idf_preserved(source, actual)
    # An unrelated waiver must not exempt a DOI mutation.
    with pytest.raises(AssertionError):
        assert_idf_preserved(source, actual, transformations={('study', 'x', 1, 'x'): ('a', 'b')})


@pytest.mark.parametrize('source,generated,reference,options,expected', [
    ('source', 'source', 'source', {}, EQUIVALENT),
    ('source', None, None, {}, LOSS),
    (None, 'invented', None, {'source_complete': True}, UNSUPPORTED),
    (None, 'maybe supported', None, {}, UNRESOLVED),
    (None, None, 'curated annotation', {}, CURATION),
    (None, 'supported', 'supported', {'evidence': 'PubMed response bound to this PMID'}, ENRICHMENT),
    ('old', 'new', 'new', {'transformation': 'Exact declared encoding correction'}, EQUIVALENT),
    ('source', 'source', 'other', {}, UNRESOLVED),
])
def test_evidence_classification_does_not_count_curation_only_as_failure(source, generated, reference, options, expected):
    assert classify(source, generated, reference, **options)[0] == expected


def two_channel():
    return [['Source Name', 'Labeled Extract Name', 'Label', 'Assay Name', 'Array Data File'],
            ['GSM1', 'extract1', 'Cy3', 'SRX1', 'a.txt'],
            ['GSM1', 'extract2', 'Cy5', 'SRX1', 'b.txt']]


def test_pairing_uses_channel_and_file_identity_without_physical_row_order():
    source = two_channel()
    actual = [source[0], source[2], source[1]]
    report = align_sdrf(source, actual)
    assert [(r['reference_row'], r['generated_row']) for r in report['matches']] == [(1, 2), (2, 1)]
    assert not report['unresolved']


def test_sample_run_subset_and_unmatched_reference_are_explicit():
    source = two_channel()
    report = align_sdrf(source, source[:2])
    assert len(report['matches']) == 1
    assert len(report['unresolved']) == 1
    assert report['reference_rows'] == 2 and report['generated_rows'] == 1


def test_ambiguous_sample_identity_does_not_get_a_positional_join():
    source = [['Source Name', 'Comment[tissue]'], ['GSM1', 'heart'], ['GSM1', 'lung']]
    result = align_sdrf(source, source)
    assert not result['matches']
    assert len(result['unresolved']) == 2


def test_pairing_retains_exact_duplicate_occurrences():
    source = [['Source Name', 'Comment[x]'], ['GSM1', 'same'], ['GSM1', 'same']]
    assert len(align_sdrf(source, source)['matches']) == 2
    assert len(align_sdrf(source, source[:2])['unresolved']) == 1


def test_qualifier_swap_cannot_count_as_reference_agreement():
    source = [['Source Name', 'Characteristics[dose]', 'Unit', 'Term Source REF'], ['GSM1', '1', 'hour', 'UO']]
    actual = deepcopy(source)
    actual[1][2], actual[1][3] = actual[1][3], actual[1][2]
    result = pair_sdrf(source, actual, 'synthetic')
    assert sum(item['classification'] == UNRESOLVED for item in result['comparisons']) == 2


def test_related_experiment_alias_and_known_6486_blank_qualifiers_are_narrow():
    from tests.support.deposited_fidelity import assert_deposited_sdrf
    assert_idf_preserved([['Comment [RelatedExperiment]', 'E-MTAB-1', 'E-MTAB-2']],
                         [['Comment[RelatedExperiment]', 'E-MTAB-1', 'E-MTAB-2']])
    source = [['Source Name', 'Factor Value[compound]', 'Factor Value[dose]', 'Unit[dose unit]'],
              ['s1', 'none', '  ', '        ']]
    actual = [source[0], ['s1', 'none', '', '']]
    assert_deposited_sdrf(source, actual, accession='E-MTAB-6486')
    with pytest.raises(AssertionError): assert_deposited_sdrf(source, actual, accession='synthetic')
    actual[1][1] = ''
    with pytest.raises(AssertionError): assert_deposited_sdrf(source, actual, accession='E-MTAB-6486')


def test_related_studies_reorder_without_losing_duplicate_occurrences():
    source = [['Comment[RelatedExperiment]', 'E-MTAB-1'], ['Comment[RelatedExperiment]', 'E-MTAB-2'],
              ['Comment[RelatedExperiment]', 'E-MTAB-1']]
    assert_idf_preserved(source, list(reversed(source)))
    assert_idf_preserved(source, [['Comment[RelatedExperiment]', 'E-MTAB-2', 'E-MTAB-1', 'E-MTAB-1']])
    with pytest.raises(AssertionError): assert_idf_preserved(source, source[:2])


def test_accession_mentioned_in_prose_does_not_create_a_sample_join():
    reference = [['Source Name', 'Comment[description]'], ['GSM1', 'study sample']]
    actual = [['Source Name', 'Comment[description]'], ['GSM2', 'Compared to GSM1']]
    result = align_sdrf(reference, actual)
    assert not result['matches'] and len(result['unresolved']) == 1
