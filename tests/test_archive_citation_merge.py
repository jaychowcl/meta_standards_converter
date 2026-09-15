# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from copy import deepcopy
import pytest
from meta_standards_converter.metadata.archive_enrichment import _merge_entity, merge_archive_metadata
from meta_standards_converter.miniml import MINiMLCodec
from tests.test_native_archive_enrichment import native, linked


@pytest.mark.parametrize('prefer', [True, False])
def test_citation_merge_preserves_distinct_papers_and_completes_identified_citations(prefer):
    original = {'pubmed_id': ['19213877', '12345'], 'pubmed_publication': [
        {'pubmed_id': '19213877', 'pmcid': 'PMC2746483', 'doi': '10.1234/source', 'title': 'Native title'},
        {'pubmed_id': '12345', 'title': 'Other explicitly linked paper'}]}
    incoming = {'pubmed_id': ['19213877', '98765'], 'pubmed_publication': [
        {'pubmed_id': '19213877', 'title': 'Preferred title', 'author_list': 'Supplied author'},
        {'pubmed_id': '98765', 'title': 'Additional linked paper'}]}
    before = deepcopy(incoming)
    _merge_entity(original, incoming, prefer)
    assert original['pubmed_id'] == ['19213877', '12345', '98765']
    assert [p['pubmed_id'] for p in original['pubmed_publication']] == ['19213877', '12345', '98765']
    merged = original['pubmed_publication'][0]
    assert merged['title'] == ('Preferred title' if prefer else 'Native title')
    assert merged['pmcid'] == 'PMC2746483' and merged['doi'] == '10.1234/source'
    assert merged['author_list'] == 'Supplied author'
    assert incoming == before
    snapshot = deepcopy(original)
    _merge_entity(original, incoming, prefer)
    assert original == snapshot


@pytest.mark.parametrize('identifier', ['', 'None', 'invalid'])
def test_titles_and_invalid_identifiers_never_establish_publication_matches(identifier):
    data = {'pubmed_publication': [{'pubmed_id': identifier, 'title': 'Same title', 'author_list': 'One person'}]}
    _merge_entity(data, {'pubmed_publication': [{'pubmed_id': identifier, 'title': 'Same title', 'doi': '10.1234/unmatched'}]}, True)
    assert len(data['pubmed_publication']) == 2
    assert 'doi' not in data['pubmed_publication'][0]
    assert 'author_list' not in data['pubmed_publication'][1]


@pytest.mark.parametrize('identifier,value', [('doi', '10.1234/Paper'), ('pmcid', 'PMC2746483')])
def test_exact_doi_or_pmcid_joins_complete_missing_pmid(identifier, value):
    data = {'pubmed_publication': [{'pubmed_id': '19213877', identifier: value, 'title': 'Native'}]}
    _merge_entity(data, {'pubmed_publication': [{'pubmed_id': '', identifier: value.lower() if identifier == 'doi' else value, 'title': 'Preferred'}]}, True)
    assert len(data['pubmed_publication']) == 1
    assert data['pubmed_publication'][0]['pubmed_id'] == '19213877'
    assert data['pubmed_publication'][0]['title'] == 'Preferred'


@pytest.mark.parametrize('extra', [
    {'pubmed_id': '1', 'doi': '10.1234/conflict'},
    {'pubmed_id': '1', 'pmcid': 'PMC999'},
    {'pubmed_id': '1', 'doi': '10.1234/second'},
])
def test_conflicting_identifier_groups_do_not_bridge_papers(extra, caplog):
    values = [{'pubmed_id': '1', 'doi': '10.1234/first', 'pmcid': 'PMC1', 'title': 'First'},
              {'pubmed_id': '2', 'doi': '10.1234/second', 'title': 'Second'}]
    data = {'pubmed_publication': deepcopy(values)}
    _merge_entity(data, {'pubmed_publication': [extra]}, True)
    assert data['pubmed_publication'][:2] == values
    assert data['pubmed_publication'][2] == extra
    assert 'publication' in caplog.text.lower() and 'conflict' in caplog.text.lower()


def test_repeated_native_occurrences_survive_and_status_annotations_stay_coupled():
    citation = {'pubmed_id': '1', 'title': 'Native', 'status': 'published',
                'status_term_source_ref': 'EFO', 'status_term_accession_number': 'EFO:1'}
    data = {'pubmed_publication': [deepcopy(citation), deepcopy(citation)]}
    _merge_entity(data, {'pubmed_publication': [{'pubmed_id': '1', 'title': 'Preferred', 'status': 'published'}]}, True)
    assert len(data['pubmed_publication']) == 2
    assert all(p['status_term_accession_number'] == 'EFO:1' and p['title'] == 'Preferred' for p in data['pubmed_publication'])
    _merge_entity(data, {'pubmed_publication': [{'pubmed_id': '1', 'status': 'in press'}]}, True)
    assert all(p['status'] == 'in press' and 'status_term_accession_number' not in p and 'status_term_source_ref' not in p
               for p in data['pubmed_publication'])


def test_publication_enrichment_keeps_entity_scope_measurements_and_residual_conflicts():
    data = native().to_mapping()
    data['series']['pubmed_publication'] = [{'pubmed_id': '1', 'title': 'Native study title', 'pmcid': 'PMC1'}]
    data['sample'][0]['pubmed_publication'] = [{'pubmed_id': '2', 'title': 'Native sample title'}]
    run = data['sample'][0]['sra_run'][0]
    run['pubmed_publication'] = [{'pubmed_id': '3', 'title': 'Native run title'}]
    run['total_spots'] = 123
    package = MINiMLCodec().decode(data).package
    extra = linked(package, 'GSE1', 'Example').to_mapping()
    extra['series']['pubmed_publication'] = [{'pubmed_id': '1', 'title': 'Preferred study title'}]
    extra['sample'][0]['pubmed_publication'] = [{'pubmed_id': '2', 'title': 'Preferred sample title'}]
    extra['sample'][0]['sra_run'][0]['pubmed_publication'] = [{'pubmed_id': '3', 'title': 'Preferred run title'}]
    extra['sample'][0]['sra_run'][0]['total_spots'] = 999
    extra['series'].setdefault('relation', []).append({'type': 'cited_by', 'target': '4', 'publication': {'pubmed_id': '4', 'title': 'Citing paper'}})
    result, issues = merge_archive_metadata(package, MINiMLCodec().decode(extra).package, prefer=True)
    merged = result.to_mapping()
    assert len(merged['series']['pubmed_publication']) == 1
    assert {k: merged['series']['pubmed_publication'][0][k] for k in ('pubmed_id', 'title', 'pmcid')} == {
        'pubmed_id': '1', 'title': 'Preferred study title', 'pmcid': 'PMC1'}
    assert merged['sample'][0]['pubmed_publication'][0]['title'] == 'Preferred sample title'
    assert merged['sample'][0]['sra_run'][0]['pubmed_publication'][0]['title'] == 'Preferred run title'
    assert merged['sample'][0]['sra_run'][0]['total_spots'] == 123
    assert any(r['type'] == 'cited_by' and r['publication']['pubmed_id'] == '4' for r in merged['series']['relation'])
    assert 'Native study title' in str(merged['extensions']['insdc'])
    assert not issues


def test_conflicting_incoming_occurrences_are_retained_without_order_winners(caplog):
    original = {'pubmed_id': '1', 'title': 'Native'}
    incoming = [{'pubmed_id': '1', 'title': 'First supplied title', 'status': 'published'},
                {'pubmed_id': '1', 'title': 'Second supplied title', 'status': 'in press'}]
    for values in (incoming, incoming[::-1]):
        data = {'pubmed_publication': [deepcopy(original)]}
        _merge_entity(data, {'pubmed_publication': values}, True)
        assert data['pubmed_publication'][0] == original
        assert {p['title'] for p in data['pubmed_publication']} == {'Native', 'First supplied title', 'Second supplied title'}
        before = deepcopy(data)
        _merge_entity(data, {'pubmed_publication': values}, True)
        assert data == before
    assert 'conflict' in caplog.text.lower()


@pytest.mark.parametrize('preferred', [
    {'status': 'published', 'status_term_source_ref': 'OTHER'},
    {'status': 'published', 'status_term_accession_number': 'OTHER:1'},
])
def test_status_namespace_conflict_and_disjoint_annotations_are_not_combined(preferred):
    data = {'pubmed_publication': [{'pubmed_id': '1', 'status': 'published', 'status_term_source_ref': 'EFO'}]}
    _merge_entity(data, {'pubmed_publication': [{'pubmed_id': '1', **preferred}]}, True)
    assert {k: v for k, v in data['pubmed_publication'][0].items() if k.startswith('status')} == preferred


def test_identifierless_record_retention_is_idempotent_without_title_donation():
    data = {'pubmed_publication': [{'pubmed_id': '', 'title': 'Shared title', 'author_list': 'Native author'}]}
    extra = {'pubmed_publication': [{'pubmed_id': '', 'title': 'Shared title', 'status': 'published'}]}
    _merge_entity(data, extra, True)
    assert len(data['pubmed_publication']) == 2
    before = deepcopy(data)
    _merge_entity(data, extra, True)
    assert data == before


def test_publication_residual_prunes_mapped_identifiers_and_keeps_displaced_text():
    from meta_standards_converter.miniml.archive_residuals import diff
    original = [{'pubmed_id': '1', 'pmcid': 'PMC1', 'doi': '10.1234/one', 'title': 'Native title'}]
    merged = [{'pubmed_id': '1', 'pmcid': 'PMC1', 'doi': '10.1234/one', 'title': 'Preferred title'}]
    assert diff(original, merged, field='pubmed_publication') == [{'pubmed_id': '1', 'title': 'Native title'}]


def test_incoming_disjoint_status_identifiers_do_not_select_by_order(caplog):
    values = [{'pubmed_id': '1', 'status': 'published', 'status_term_source_ref': 'OTHER'},
              {'pubmed_id': '1', 'status': 'published', 'status_term_accession_number': 'EFO:1'}]
    for incoming in (values, values[::-1]):
        data = {'pubmed_publication': [{'pubmed_id': '1', 'title': 'Native'}]}
        _merge_entity(data, {'pubmed_publication': incoming}, True)
        assert len(data['pubmed_publication']) == 3
        assert data['pubmed_publication'][1:] == incoming
    assert 'conflict' in caplog.text.lower()


def test_repeated_conflicting_citation_does_not_multiply_and_reports_scoped_outcome():
    data = native().to_mapping()
    data['series']['pubmed_publication'] = [{'pubmed_id': '1', 'doi': '10.1234/native'}]
    package = MINiMLCodec().decode(data).package
    incoming = linked(package, 'GSE1', 'Example').to_mapping()
    incoming['series']['pubmed_publication'] = [{'pubmed_id': '1', 'doi': '10.1234/conflict'}]
    extra = MINiMLCodec().decode(incoming).package
    once, issues = merge_archive_metadata(package, extra, prefer=True)
    twice, _ = merge_archive_metadata(once, extra, prefer=True)
    assert len(once.to_mapping()['series']['pubmed_publication']) == 2
    assert twice.to_mapping()['series']['pubmed_publication'] == once.to_mapping()['series']['pubmed_publication']
    assert any('SRP250911' in issue and 'conflicting publication' in issue for issue in issues)


@pytest.mark.parametrize('reverse', [False, True])
def test_complementary_incoming_identifiers_cannot_bridge_conflicting_native_papers(reverse):
    native = [{'pubmed_id': '1', 'doi': '10.1234/one'}, {'pubmed_id': '2', 'doi': '10.1234/two'}]
    incoming = [{'pubmed_id': '1', 'pmcid': 'PMC99'}, {'pmcid': 'PMC99', 'doi': '10.1234/two'}]
    if reverse:
        incoming.reverse()
    data = {'pubmed_publication': deepcopy(native)}
    _merge_entity(data, {'pubmed_publication': incoming}, True)
    assert data['pubmed_publication'][:2] == native
    assert data['pubmed_publication'][2:] == incoming


def test_same_priority_citation_conflicts_joined_through_native_identity_have_no_order_winner():
    original = {'pubmed_id': '1', 'doi': '10.1234/one', 'title': 'Native'}
    incoming = [{'pubmed_id': '1', 'title': 'First'}, {'doi': '10.1234/one', 'title': 'Second'}]
    for values in (incoming, incoming[::-1]):
        data = {'pubmed_publication': [deepcopy(original)]}
        _merge_entity(data, {'pubmed_publication': values}, True)
        assert data['pubmed_publication'] == [original, *values]


@pytest.mark.parametrize('native_count,incoming_count', [(0, 2), (1, 2), (3, 2), (1, 3)])
def test_incoming_citation_occurrences_preserve_multiplicity_idempotently(native_count, incoming_count):
    original = {'pubmed_id': '1', 'title': 'Native', 'pmcid': 'PMC1'}
    incoming = {'pubmed_id': '1', 'title': 'Preferred'}
    data = {'pubmed_publication': [deepcopy(original) for _ in range(native_count)]}
    extra = {'pubmed_publication': [deepcopy(incoming) for _ in range(incoming_count)]}
    _merge_entity(data, extra, True)
    assert len(data['pubmed_publication']) == max(native_count, incoming_count)
    assert all(v['title'] == 'Preferred' for v in data['pubmed_publication'])
    if native_count:
        assert all(v['pmcid'] == 'PMC1' for v in data['pubmed_publication'])
    before = deepcopy(data)
    _merge_entity(data, extra, True)
    assert data == before


def test_added_citation_occurrence_does_not_intersect_conflicting_status_groups():
    data = {'pubmed_publication': [
        {'pubmed_id': '1', 'status': 'published', 'status_term_source_ref': 'EFO', 'status_term_accession_number': 'EFO:1'},
        {'pubmed_id': '1', 'status': 'in press', 'status_term_source_ref': 'EFO', 'status_term_accession_number': 'EFO:2'}]}
    value = {'pubmed_id': '1', 'title': 'Preferred'}
    _merge_entity(data, {'pubmed_publication': [deepcopy(value) for _ in range(3)]}, True)
    assert data['pubmed_publication'][2] == value


def test_partial_identifier_occurrences_count_the_whole_verified_paper_for_idempotence():
    data = {'pubmed_publication': [
        {'pubmed_id': '1', 'doi': '10.1234/one', 'title': 'Native'}, {'pubmed_id': '1'}]}
    extra = {'pubmed_publication': [
        {'pubmed_id': '1', 'author_list': 'Supplied'}, {'pubmed_id': '1', 'author_list': 'Supplied'},
        {'doi': '10.1234/one', 'author_list': 'Supplied'}]}
    _merge_entity(data, extra, True)
    assert len(data['pubmed_publication']) == 3
    before = deepcopy(data)
    for _ in range(3):
        _merge_entity(data, extra, True)
        assert data == before


def test_verified_paper_component_updates_compatible_partial_identifier_occurrences():
    data = {'pubmed_publication': [
        {'pubmed_id': '1', 'doi': '10.1234/one', 'title': 'Native'},
        {'pubmed_id': '1', 'title': 'Native', 'source_note': 'Separate occurrence'}]}
    _merge_entity(data, {'pubmed_publication': [{'doi': '10.1234/one', 'title': 'Preferred'}]}, True)
    assert [p['title'] for p in data['pubmed_publication']] == ['Preferred', 'Preferred']
    assert data['pubmed_publication'][1]['source_note'] == 'Separate occurrence'
