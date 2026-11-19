"""Curator findings describe final metadata without changing scientific outcomes."""
from copy import deepcopy
from unittest.mock import Mock

import pytest
import requests

from meta_standards_converter.converters import Converter, PreparationPolicy
from meta_standards_converter.metadata.enrichment import MINiMLEnricher
from meta_standards_converter.miniml import MINiMLCodec
from tests.support.conversion import package


def findings(data, **kwargs):
    from meta_standards_converter.metadata.curation import assess_curation
    return assess_curation(data, dataset_id='GSE1', **kwargs)


def fields(values):
    return {d.curation.field for d in values}


def test_assessment_is_pure_and_reports_study_sample_and_single_missing_citation():
    data = package()
    before = deepcopy(data)
    result = findings(data)
    assert data == before
    assert {'study.title', 'study.description', 'study.overall_design', 'study.experimental_type',
            'study.contributors', 'study.contacts', 'publication.citation',
            'sample.source', 'sample.developmental_stage', 'sample.disease',
            'sample.genotype', 'sample.molecule'} <= fields(result)
    assert 'sample.organism_part' not in fields(result)  # explicit tissue alias
    pubs = [d for d in result if d.curation.field.startswith('publication.')]
    assert len(pubs) == 1 and pubs[0].location == '/series/pubmed_publication'
    assert all(d.severity == 'warning' and d.action and d.location for d in result)


@pytest.mark.parametrize('value,reason', [(None, 'missing'), ('  ', 'missing'),
                                         ('unknown', 'reported_unknown'), ('not provided', 'reported_unknown'),
                                         ('none', None), ('not applicable', None), (0, None), (False, None)])
def test_missing_tokens_do_not_erase_legitimate_biology(value, reason):
    data = package()
    data['sample'][0]['channel'][0]['characteristics'].append({'name': 'disease', 'value': value})
    selected = [d for d in findings(data) if d.curation.field == 'sample.disease']
    assert [d.curation.reason for d in selected] == ([reason] if reason else [])


def test_publication_fields_identifiers_and_ontology_companions():
    data = package()
    data['series']['pubmed_publication'] = [{'doi': '10.1234/paper', 'title': 'Paper',
        'status': 'published', 'status_term_source_ref': 'EFO'}]
    selected = [d for d in findings(data) if d.curation.field.startswith('publication.')]
    assert fields(selected) == {'publication.authors', 'publication.status_term_accession_number'}
    data['series']['pubmed_publication'][0] = {'title': 'Paper'}
    assert 'publication.identifier' in fields(findings(data))


def test_declared_values_factors_and_partial_unit_groups():
    data = package()
    data['series']['variable'] = [{'name': 'dose', 'type': 'dose'}]
    data['sample'][0]['channel'][0]['characteristics'] += [
        {'name': 'age', 'value': '', 'unit': {'value': 'day', 'term_source_ref': 'UO'}},
        {'name': 'sex', 'value': 'female'}]
    result = findings(data)
    assert {'characteristic.age', 'factor.dose', 'ontology.term_accession_number'} <= fields(result)
    assert not any(d.location.endswith('sex/term_source_ref') for d in result)


def test_mixed_assay_checks_are_scoped_to_the_relevant_sample():
    data = package()
    array = deepcopy(data['sample'][0]); array['iid'] = 'GSM2'
    array['platform_ref'] = {'ref': 'GPL1'}
    data['platform'] = [{'iid': 'GPL1', 'technology': 'in situ oligonucleotide'}]
    data['sample'].append(array)
    data['sample'][0]['library_strategy'] = 'RNA-Seq'
    data['sample'][0]['sra_run'] = [{'run': 'SRR1'}]
    result = findings(data)
    seq = fields([d for d in result if d.curation.entity_id.startswith('GSM1')])
    arr = fields([d for d in result if d.curation.entity_id.startswith('GSM2')])
    assert {'assay.library_layout', 'assay.instrument_model', 'protocol.sequencing', 'data.raw_file'} <= seq
    assert 'assay.library_layout' not in arr
    assert 'protocol.hybridization' in arr


def test_bound_protocols_and_factors_satisfy_checks_without_cross_sample_leakage():
    data = package()
    data['sample'][0]['library_strategy'] = 'RNA-Seq'
    data['series']['protocols'] = [{'name': 'P1', 'type': 'nucleic acid extraction', 'description': 'Extract RNA.'}]
    data['series']['variable'] = [{'name': 'dose'}]
    data['series']['assay_paths'] = [{'steps': [
        {'kind': 'source', 'name': 'source1', 'sample_ref': 'GSM1',
         'factor_values': [{'name': 'dose', 'value': '0'}]},
        {'kind': 'protocol_application', 'protocol_ref': 'P1'},
        {'kind': 'extract', 'name': 'extract1'}]}]
    second = deepcopy(data['sample'][0]); second['iid'] = 'GSM2'; data['sample'].append(second)
    result = findings(data)
    first = fields([d for d in result if d.curation.entity_id.startswith('GSM1')])
    second = fields([d for d in result if d.curation.entity_id.startswith('GSM2')])
    assert 'protocol.extraction' not in first and 'factor.dose' not in first
    assert {'protocol.extraction', 'factor.dose'} <= second


@pytest.mark.parametrize('target', ['json', 'magetab'])
def test_converter_retains_advisories_without_changing_success_or_input(target):
    data = package(); before = deepcopy(data)
    result = Converter().convert(data, out_type=target,
        preparation=PreparationPolicy(enrichment='off', expand_studies=False))
    item = result.items[0]
    assert (result.status, item.execution, item.validation) == ('complete', 'succeeded', 'valid')
    assert data == before
    advisories = [d for d in item.diagnostics if getattr(d, 'curation', None)]
    assert 'publication.citation' in fields(advisories)
    assert all(d.to_dict()['curation']['dataset_id'] == 'GSE1' for d in advisories)
    if target == 'json':
        assert item.payload[0].to_mapping() == MINiMLCodec().decode(data).package.to_mapping()


@pytest.mark.parametrize('response,status,lookup', [
    (requests.Timeout('private-secret'), 'partial', 'failed'),
    ((None,) * 6, 'partial', 'empty'),
    (('10.1234/paper', None, 'Title', None, None, None), 'complete', 'incomplete'),
    (('10.1234/paper', 'Authors', 'Title', 'published', None, None), 'complete', 'complete'),
    (('10.1234/paper', 'Authors', 'Title', 'published', 'EFO', 'EFO_0001796'), 'complete', 'complete')])
def test_lookup_outcomes_and_remaining_fields(response, status, lookup):
    data = package(); data['series']['pubmed_id'] = ['123']
    fetch = Mock()
    if isinstance(response, Exception): fetch.pubmed_summary.side_effect = response
    else: fetch.pubmed_summary.return_value = response
    service = MINiMLEnricher(pubmed_fetcher=fetch, insdc_fetcher=Mock())
    result = Converter(services={'enricher': service}).convert(data, out_type='json',
        preparation=PreparationPolicy(enrichment='curators', expand_studies=False))
    item = result.items[0]
    assert result.status == status
    assert fetch.pubmed_summary.call_count == 1
    records = [r for r in item.preparation if r['operation'] == 'publication_lookup']
    assert len(records) == 1 and records[0]['status'] == lookup and records[0]['identifier'] == '123'
    pubs = [d for d in item.diagnostics if getattr(d, 'curation', None) and d.curation.field.startswith('publication.')]
    if lookup == 'complete': assert not pubs
    else: assert 'publication.authors' in fields(pubs)
    assert 'private-secret' not in str(result.to_dict())


def test_no_match_and_disabled_enrichment_are_distinct():
    data = package(); data['series']['pubmed_publication'] = [{'doi': '10.1234/paper', 'title': 'Title'}]
    service = MINiMLEnricher(pubmed_fetcher=Mock(), insdc_fetcher=Mock(),
                            publication_identifier_resolver=lambda *_: None)
    result = Converter(services={'enricher': service}).convert(data, out_type='json',
        preparation=PreparationPolicy(enrichment='curators', expand_studies=False))
    assert result.status == 'complete'
    assert any(d.curation and d.curation.reason == 'lookup_no_match' for d in result.items[0].diagnostics)
    service.pubmed_fetcher.pubmed_summary.assert_not_called()
    off = Converter(services={'enricher': service}).convert(data, out_type='json',
        preparation=PreparationPolicy(enrichment='off', expand_studies=False))
    assert any(d.curation and d.curation.lookup_status == 'disabled' for d in off.items[0].diagnostics)
    assert not any(d.curation and d.curation.reason == 'lookup_no_match' for d in off.items[0].diagnostics)


def test_native_values_keep_exact_paths_and_empty_declared_factors_are_checked():
    data = package()
    data['series']['assay_paths'] = [{'steps': [
        {'kind': 'source', 'name': 'source1', 'sample_ref': 'GSM1',
         'characteristics': [{'name': 'disease', 'value': 'unknown'}],
         'factor_values': [{'name': 'time', 'value': ''}]},
        {'kind': 'array_data_file', 'name': 'reads.fastq'}]}]
    data['sample'][0]['library_strategy'] = 'RNA-Seq'
    result = findings(data)
    assert any(d.curation.field == 'sample.disease' and
               d.location == '/series/assay_paths/0/steps/0/characteristics/0/value' for d in result)
    assert any(d.curation.field == 'factor.time' and
               d.location == '/series/assay_paths/0/steps/0/factor_values/0/value' for d in result)
    assert 'data.raw_file' not in fields(result)


def test_every_empty_declared_characteristic_even_with_populated_alias():
    data = package()
    data['sample'][0]['channel'][0]['characteristics'] += [
        {'name': 'organism_part', 'value': ''}, {'name': 'tissue', 'value': 'unknown'}]
    selected = [d for d in findings(data) if d.curation.field == 'sample.organism_part']
    assert {d.curation.reason for d in selected} == {'missing', 'reported_unknown'}
    assert len(selected) == 2


def test_scoped_factor_declaration_does_not_require_other_samples():
    data = package()
    data['series']['variable'] = [{'factor': 'dose', 'sample_ref': [{'ref': 'GSM1'}]}]
    data['sample'].append(dict(deepcopy(data['sample'][0]), iid='GSM2'))
    selected = [d for d in findings(data) if d.curation.field == 'factor.dose']
    assert [d.curation.entity_id for d in selected] == ['GSM1']


def test_empty_biological_factor_keeps_its_declared_name_and_canonical_biology():
    data = package()
    data['series']['variable'] = [{'factor': 'tissue'}]
    data['sample'][0]['channel'][0]['characteristics'] = [{'name': 'tissue', 'value': ''}]
    selected = fields(findings(data))
    assert {'factor.tissue', 'sample.organism_part'} <= selected
    assert 'factor.organism_part' not in selected


def test_sequencing_without_runs_and_extraction_only_need_run_and_library_protocol():
    data = package()
    sample = data['sample'][0]
    sample['library_strategy'] = 'RNA-Seq'
    sample['channel'][0]['extract_protocol'] = 'RNA was extracted with TRIzol.'
    result = fields(findings(data))
    assert {'run.identifier', 'protocol.library_preparation'} <= result
    assert 'protocol.extraction' not in result


def test_run_files_are_not_borrowed_from_a_different_run():
    data = package()
    sample = data['sample'][0]
    sample['library_strategy'] = 'RNA-Seq'
    sample['sra_run'] = [{'run': 'SRR1', 'fastq_files': [{'uri': 'reads1.fastq'}]}, {'run': 'SRR2'}]
    sample['raw_data'] = [{'value': 'reads1.fastq'}]
    selected = [d for d in findings(data) if d.curation.field == 'data.raw_file']
    assert [d.curation.entity_id for d in selected] == ['GSM1/run[SRR2]']


@pytest.mark.parametrize('kind', ['derived_array_data_file', 'derived_array_data_matrix_file'])
def test_processed_paths_require_processing_description(kind):
    data = package()
    data['series']['assay_paths'] = [{'steps': [
        {'kind': 'sample', 'name': 'GSM1', 'sample_ref': 'GSM1'}, {'kind': kind, 'name': 'matrix.tsv'}]}]
    assert 'protocol.processing' in fields(findings(data))


def test_relation_publication_context_and_status_companions():
    data = package()
    data['series']['relation'] = [{'type': 'literature cross-reference', 'target': '123',
        'experiment_ref': 'SRX1', 'publication': {'pubmed_id': '123', 'title': 'Paper',
        'status': 'published', 'status_term_source_ref': 'EFO'}}]
    result = findings(data)
    assert 'publication.citation' in fields(result)  # experiment citation is not study attribution
    selected = [d for d in result if d.curation.field == 'publication.status_term_accession_number']
    assert len(selected) == 1 and selected[0].curation.entity_type == 'publication'
    assert selected[0].location == '/series/relation/0/publication/status_term_accession_number'
    assert 'SRX1' in selected[0].curation.entity_id


def test_free_text_citation_is_not_absent_and_reports_unstructured_fields():
    data = package(); data['series']['citation'] = ['A source citation without identifiers']
    result = findings(data)
    assert 'publication.citation' not in fields(result)
    assert 'publication.identifier' in fields(result)
    assert any(d.location == '/series/citation/0' for d in result)


def test_empty_people_references_do_not_count_as_contact_information():
    data = package()
    data['series']['contributor_ref'] = [{'ref': 'person1'}]
    data['series']['contact'] = [{'iid': 'person2'}]
    data['contributor'] = [{'iid': 'person1'}]
    assert {'study.contributors', 'study.contacts'} <= fields(findings(data))


def test_successful_empty_run_lookup_is_advisory_with_requested_identifier():
    data = package()
    data['sample'][0]['relation'] = [{'type': 'SRA', 'target': 'SRX1'}]
    fetch = Mock()
    fetch.extract_sra_accessions.return_value = ['SRX1']
    fetch.fetch_sra_runs.return_value = []
    service = MINiMLEnricher(pubmed_fetcher=Mock(), insdc_fetcher=fetch)
    result = Converter(services={'enricher': service}).convert(data, out_type='json',
        preparation=PreparationPolicy(expand_studies=False))
    assert result.status == 'complete'
    selected = [d for d in result.items[0].diagnostics if d.curation and d.curation.reason == 'lookup_no_match']
    assert len(selected) == 1 and selected[0].provider == 'sra'
    assert selected[0].curation.lookup_identifier == 'SRX1'
    fetch.fetch_sra_runs.assert_called_once()


def test_native_array_and_sequencing_assays_on_one_sample_use_their_own_expectations():
    data = package()
    data['series']['assay_paths'] = [{'steps': [
        {'kind': 'source', 'name': 'source1', 'sample_ref': 'GSM1'},
        {'kind': 'assay', 'name': name, 'technology_type': technology}]}
        for name, technology in [('array1', 'array assay'), ('seq1', 'sequencing assay')]]
    result = findings(data)
    arrays = [d for d in result if '/assay[array1]' in d.curation.entity_id]
    sequencing = [d for d in result if '/assay[seq1]' in d.curation.entity_id]
    assert {'assay.array_design', 'protocol.hybridization'} <= fields(arrays)
    assert 'assay.library_layout' not in fields(arrays)
    assert {'assay.library_layout', 'run.identifier'} <= fields(sequencing)
    assert 'assay.array_design' not in fields(sequencing)


def test_empty_publication_placeholder_is_one_absent_citation():
    data = package(); data['series']['pubmed_publication'] = [{}]
    selected = [d for d in findings(data) if d.curation.field.startswith('publication.')]
    assert fields(selected) == {'publication.citation'} and len(selected) == 1


def test_shared_failed_lookup_has_one_finding_per_final_occurrence_and_no_false_association():
    sources = [package('GSE1', 'GSM1'), package('GSE2', 'GSM2')]
    for data in sources:
        data['series']['pubmed_id'] = ['123', '123']
    fetch = Mock(); fetch.pubmed_summary.side_effect = requests.Timeout('secret')
    service = MINiMLEnricher(pubmed_fetcher=fetch, insdc_fetcher=Mock())
    result = Converter(services={'enricher': service}).convert(sources, out_type='json',
        preparation=PreparationPolicy(enrichment='curators', expand_studies=False))
    assert result.status == 'partial'
    fetch.pubmed_summary.assert_called_once()
    findings_ = [d for item in result.items for d in item.diagnostics if d.curation]
    locations = [(d.curation.dataset_id, d.curation.package_index, d.curation.field, d.location) for d in findings_]
    assert len(locations) == len(set(locations))
    authors = [d for d in findings_ if d.curation.field == 'publication.authors']
    assert len(authors) == 2 and all(d.curation.lookup_status == 'failed' for d in authors)
    assert all(d.curation.lookup_identifier == '123' and d.curation.lookup_provider == 'pubmed' for d in authors)
    assert all(d.curation.lookup_status is None for d in findings_ if d.curation.field.startswith('sample.'))


def test_lookup_context_never_crosses_packages_in_one_study():
    data = package(); data['series']['pubmed_publication'] = [{'pubmed_id': '123'}]
    preparation = [dict(operation='publication_lookup', dataset='GSE1', dataset_id='GSE1',
                        package_index=1, identifier_kind='pubmed_id', identifier='123', status='failed', provider='pubmed')]
    assert all(d.curation.lookup_status is None for d in findings(data, preparation=preparation, package_index=0))
    assert any(d.curation.lookup_status == 'failed' for d in findings(data, preparation=preparation, package_index=1))


def test_magetab_import_is_assessed_after_validation_without_enrichment(tmp_path):
    idf = tmp_path/'study.idf.txt'
    idf.write_text('Investigation Title\tStudy\nComment[ArrayExpressAccession]\tE-MTAB-1\nSDRF File\tstudy.sdrf.txt\n')
    (tmp_path/'study.sdrf.txt').write_text('Source Name\tCharacteristics[organism part]\tAssay Name\tTechnology Type\n'
                                        'sample1\tunknown\tarray1\tarray assay\n')
    result = Converter().convert(idf, out_type='json', preparation=PreparationPolicy(enrichment='off', expand_studies=False))
    assert result.status == 'complete', result.to_dict()
    selected = [d for d in result.items[0].diagnostics if d.curation]
    assert 'assay.array_design' in fields(selected)
    assert any(d.curation.field == 'sample.organism_part' and d.curation.reason == 'reported_unknown' for d in selected)


def test_exhausted_http_retries_produce_one_set_of_final_citation_findings(tmp_path):
    from meta_standards_converter.helpers.request_helper import HostRequestGate, RateLimitedRequester, RequestSettings
    from meta_standards_converter.sources.pubmed import PubmedWebFetcher
    get = Mock(side_effect=requests.Timeout('private retry failure'))
    requester = RateLimitedRequester('ncbi_eutils', get=get, host_gate=HostRequestGate(tmp_path),
        settings=RequestSettings(max_retries=2, request_delay=0, backoff_base=0, backoff_max=0))
    service = MINiMLEnricher(pubmed_fetcher=PubmedWebFetcher(requester=requester), insdc_fetcher=Mock())
    data = package(); data['series']['pubmed_id'] = ['123']
    result = Converter(services={'enricher': service}).convert(data, out_type='json',
        preparation=PreparationPolicy(enrichment='curators', expand_studies=False))
    assert result.status == 'partial'
    assert get.call_count == 3 and requester.retry_count == 2
    authors = [d for d in result.items[0].diagnostics if d.curation and d.curation.field == 'publication.authors']
    assert len(authors) == 1 and authors[0].curation.lookup_status == 'failed'
    provider_failures = [d for d in result.items[0].diagnostics if d.provider_issue]
    assert len(provider_failures) == 1
    record = next(r for r in result.items[0].preparation if r['operation'] == 'publication_lookup')
    assert record['correlation_id'] == provider_failures[0].correlation_id


def test_recorded_replay_retains_curation_without_changing_frozen_scientific_goldens():
    from pathlib import Path
    from scripts.audit_recorded_provider import audit
    report = audit(Path(__file__).parents[1]/'fixtures/recorded_provider/success.manifest.json')
    assert report['classification'] == 'complete-captured-coverage'
    for call in report['scenarios'][0]['calls']:
        assert any(d.get('curation') for item in call['result']['items'] for d in item['diagnostics'])
