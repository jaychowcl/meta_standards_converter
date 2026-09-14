# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""A refresh must not invalidate already resolved experimental evidence."""
from copy import deepcopy
from unittest.mock import Mock

import pytest
import requests

from meta_standards_converter.metadata.enrichment import MINiMLEnricher


def fixture(response):
    fetcher = Mock()
    fetcher.extract_sra_accessions.return_value = ['SRX1']
    if isinstance(response, Exception):
        fetcher.fetch_sra_runs.side_effect = response
    else:
        fetcher.fetch_sra_runs.return_value = response
    data = {'sample': [{'iid': 'GSM1', 'relation': [{'type': 'SRA', 'target': 'SRX1'}],
                       'sra_accession': ['SRR1'], 'ena_accession': ['SRP1'],
                       'sra_run': [{'run': 'SRR1', 'experiment': 'SRX1', 'sample': 'SRS1',
                                    'study': 'SRP1', 'fastq_files': [{'uri': 'https://reads/1', 'md5': 'a'}]}]}]}
    return MINiMLEnricher(insdc_fetcher=fetcher, pubmed_fetcher=Mock()), data


@pytest.mark.parametrize('response', [requests.ConnectionError('offline'), []])
def test_failed_or_empty_refresh_preserves_saved_runs_and_accessions(response):
    enricher, data = fixture(response)
    original = deepcopy(data['sample'][0])
    enricher.enrich_sra(data)
    sample = data['sample'][0]
    assert sample['sra_run'] == original['sra_run']
    assert sample['sra_accession'] == ['SRR1', 'SRX1']
    assert sample['ena_accession'] == ['SRP1']


def test_partial_refresh_completes_exact_files_without_losing_runs_or_conflicts():
    new = {'run': 'SRR1', 'experiment': 'SRX1', 'sample': 'SRS1', 'instrument_model': 'machine',
           'fastq_files': [{'uri': 'https://reads/1', 'md5': 'a', 'bytes': '10'},
                           {'uri': 'https://reads/2', 'md5': 'b'},
                           {'uri': 'https://reads/1', 'md5': 'conflict'}]}
    enricher, data = fixture([new, {'run': 'SRR3', 'experiment': 'SRX1'}])
    data['sample'][0]['sra_run'].append({'run': 'SRR2', 'experiment': 'SRX1'})
    enricher.enrich_sra(data)
    runs = data['sample'][0]['sra_run']
    assert [r['run'] for r in runs] == ['SRR1', 'SRR2', 'SRR3']
    assert runs[0]['study'] == 'SRP1' and runs[0]['instrument_model'] == 'machine'
    assert runs[0]['fastq_files'] == new['fastq_files']
    once = deepcopy(data); enricher.enrich_sra(data); assert data == once


@pytest.mark.parametrize('key,value', [('sample', 'SRS2'), ('experiment', 'SRX2'), ('biosample', 'SAMN2')])
def test_same_run_with_conflicting_owner_cannot_supply_metadata(key, value):
    enricher, data = fixture([{'run': 'SRR1', key: value, 'instrument_model': 'wrong'}])
    data['sample'][0]['sra_run'][0].setdefault('biosample', 'SAMN1')
    before = deepcopy(data['sample'][0]['sra_run'])
    enricher.enrich_sra(data)
    assert data['sample'][0]['sra_run'] == before


def test_run_reconciliation_is_sample_local_and_does_not_join_missing_file_identity():
    enricher, data = fixture([{'run': 'SRR1', 'fastq_files': [{'filename': 'same.fastq', 'md5': 'x'}]}])
    data['sample'].append({'iid': 'GSM2', 'sra_run': [{'run': 'SRR1', 'instrument_model': 'other'}]})
    before = deepcopy(data['sample'][1])
    enricher.enrich_sra(data)
    assert data['sample'][1] == before
    assert data['sample'][0]['sra_run'][0]['fastq_files'][0]['uri'] == 'https://reads/1'
    assert len(data['sample'][0]['sra_run'][0]['fastq_files']) == 2


@pytest.mark.parametrize('response', [requests.ConnectionError('offline'), (None,) * 6,
                                    ('new-doi', 'new authors', 'new title', 'published', 'EFO', 'EFO:1')])
def test_saved_citation_facts_and_status_groups_survive_refresh(response):
    fetcher = Mock()
    if isinstance(response, Exception): fetcher.pubmed_summary.side_effect = response
    else: fetcher.pubmed_summary.return_value = response
    data = {'series': {'pubmed_id': ['1'], 'pubmed_publication': [
        {'pubmed_id': '1', 'title': 'source title', 'doi': 'source-doi', 'author_list': 'source authors',
         'status': 'in press', 'status_term_source_ref': 'OTHER', 'unknown_detail': 'retain'}]}}
    before = deepcopy(data)
    MINiMLEnricher(pubmed_fetcher=fetcher, insdc_fetcher=Mock()).enrich_pubmed(data)
    assert data == before


def test_saved_citation_fills_only_missing_compatible_details():
    fetcher = Mock(); fetcher.pubmed_summary.return_value = ('doi', 'authors', 'remote title', 'published', 'EFO', 'EFO:1')
    data = {'series': {'pubmed_id': ['1'], 'pubmed_publication': [{'pubmed_id': '1', 'title': 'source title', 'status': 'published'}]}}
    MINiMLEnricher(pubmed_fetcher=fetcher, insdc_fetcher=Mock()).enrich_pubmed(data)
    assert data['series']['pubmed_publication'] == [{'pubmed_id': '1', 'title': 'source title', 'doi': 'doi',
        'author_list': 'authors', 'status': 'published', 'status_term_source_ref': 'EFO', 'status_term_accession_number': 'EFO:1'}]


@pytest.mark.parametrize('reverse', [False, True])
def test_sparse_file_cannot_choose_between_conflicting_refresh_versions(reverse):
    files=[{'uri':'https://reads/1','md5':'a'}, {'uri':'https://reads/1','md5':'b'}]
    if reverse:files.reverse()
    enricher,data=fixture([{'run':'SRR1','experiment':'SRX1','fastq_files':files}])
    sparse={'uri':'https://reads/1','read_index':'1'}
    data['sample'][0]['sra_run'][0]['fastq_files']=[deepcopy(sparse)]
    enricher.enrich_sra(data)
    assert data['sample'][0]['sra_run'][0]['fastq_files']==[sparse,*files]


def test_sparse_run_cannot_choose_between_conflicting_refresh_owners():
    enricher,data=fixture([{'run':'SRR1','sample':'SRS1'}, {'run':'SRR1','sample':'SRS2'}])
    data['sample'][0]['sra_run']=[{'run':'SRR1'}]
    enricher.enrich_sra(data)
    assert data['sample'][0]['sra_run']==[{'run':'SRR1'}]


def test_full_enrich_preserves_input_package_on_failed_refresh():
    from meta_standards_converter.miniml import MINiMLCodec
    enricher,data=fixture(requests.ConnectionError('offline'))
    data.update(miniml_schema_version='3.0',source={'format':'GEO'},series={'iid':'GSE1'})
    package=MINiMLCodec().decode(data).package;before=package.to_mapping()
    out=enricher.enrich(package).to_mapping()
    assert package.to_mapping()==before
    assert out['sample'][0]['sra_run']==before['sample'][0]['sra_run']


def test_explicit_empty_saved_file_list_can_be_completed():
    enricher,data=fixture([{'run':'SRR1','fastq_files':[{'uri':'https://reads/1'}]}])
    data['sample'][0]['sra_run'][0]['fastq_files']=None
    enricher.enrich_sra(data)
    assert data['sample'][0]['sra_run'][0]['fastq_files']==[{'uri':'https://reads/1'}]
