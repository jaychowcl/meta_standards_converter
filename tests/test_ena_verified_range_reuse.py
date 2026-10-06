# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""ENA range evidence reuse never changes native membership or link literals."""
from copy import deepcopy
import xml.etree.ElementTree as ET

import pytest

from meta_standards_converter.sources.archive_support import StudyRecords, StudySeed
from meta_standards_converter.sources.ena import ENASource
from tests.test_native_archive_sources import HTTP


def range_records(literals=('ERR001-ERR003',), database='ENA-RUN', evidence=''):
    links = ''.join(f'<XREF_LINK><DB>{database}</DB><ID>{literal}</ID></XREF_LINK>'
                    for literal in literals)
    records = StudyRecords(StudySeed('ERP1', 'PRJEB1'))
    records.xml = [ET.fromstring(f'<SAMPLE_SET><SAMPLE accession="ERS1">{links}</SAMPLE></SAMPLE_SET>')]
    if evidence:
        records.xml.append(ET.fromstring(evidence))
    return records


def range_http(tag='RUN'):
    def handler(url, params, fmt):
        ids = url.rsplit('/', 1)[-1].split(',')
        return ET.fromstring(f'<{tag}_SET>' + ''.join(f'<{tag} accession="{acc}"/>' for acc in ids) + f'</{tag}_SET>')
    return HTTP(handler)


def requested(http):
    return [url.rsplit('/', 1)[-1].split(',') for url, _, _ in http.calls]


def test_reuses_core_and_overlapping_range_records_without_membership_mutation():
    literals = ('ERR001-ERR003', 'ERR002-ERR004,ERR001')
    records = range_records(literals, evidence='<RUN_SET><RUN accession="ERR001"/></RUN_SET>')
    original = [ET.tostring(root) for root in records.xml]
    http = range_http()
    ENASource(http=http).fetch_reference_ranges(records)
    assert requested(http) == [['ERR002', 'ERR003'], ['ERR004']]
    assert records.linked == [
        {'provider': 'ena', 'kind': 'accession_range', 'accession': 'ERS1',
         'metadata': {'database': 'ENA-RUN', 'literal': literal, 'accessions': ids}}
        for literal, ids in zip(literals, [
            ['ERR001', 'ERR002', 'ERR003'], ['ERR002', 'ERR003', 'ERR004', 'ERR001']])]
    assert [ET.tostring(root) for root in records.xml] == original
    assert records.issues == []


@pytest.mark.parametrize(('database', 'tag', 'prefix'), [
    ('ENA-RUN', 'RUN', 'ERR'), ('ENA-EXPERIMENT', 'EXPERIMENT', 'ERX'), ('ENA-STUDY', 'STUDY', 'ERP')])
def test_exact_primary_and_provider_secondary_aliases_are_valid_evidence(database, tag, prefix):
    records = range_records((f'{prefix}001-{prefix}003',), database, f'''<{tag}_SET>
        <{tag} accession="{prefix}001"/>
        <{tag}><IDENTIFIERS><PRIMARY_ID>{prefix}002</PRIMARY_ID></IDENTIFIERS></{tag}>
        <{tag} accession="{prefix}900"><IDENTIFIERS><SECONDARY_ID>{prefix}003</SECONDARY_ID></IDENTIFIERS></{tag}>
        </{tag}_SET>''')
    http = range_http(tag)
    ENASource(http=http).fetch_reference_ranges(records)
    assert http.calls == []
    assert len(records.linked) == 1


@pytest.mark.parametrize('evidence', [
    '<RUN_SET><RUN_REF accession="ERR001"/></RUN_SET>',
    '<EXPERIMENT_SET><EXPERIMENT accession="ERR001"/></EXPERIMENT_SET>',
    '<RUN_SET><RUN accession="ERR001"><ERROR>unavailable</ERROR></RUN></RUN_SET>',
    '<ERROR><RUN accession="ERR001"/></ERROR>',
    '<RUN_SET><RUN accession="ERR001"><e:ERROR xmlns:e="error">unavailable</e:ERROR></RUN></RUN_SET>',
    '<RUN_SET><RUN accession="ERR001.1"/></RUN_SET>',
    '<RUN_SET><RUN accession="SRR001"/></RUN_SET>',
    '<RUN_SET xmlns="unexpected"><RUN accession="ERR001"/></RUN_SET>',
    '<RUN_SET><RUN accession="ERR900"><IDENTIFIERS><SUBMITTER_ID>ERR001</SUBMITTER_ID></IDENTIFIERS></RUN></RUN_SET>',
    '<RUN_SET><RUN accession="ERR900"><IDENTIFIERS><EXTERNAL_ID namespace="guessed">ERR001</EXTERNAL_ID></IDENTIFIERS></RUN></RUN_SET>',
])
def test_non_record_or_unverified_identity_cannot_suppress_fetch(evidence):
    records = range_records(('ERR001-ERR001',), evidence=evidence)
    records.indexed['read_run'] = [{'run_accession': 'ERR001'}]
    http = range_http()
    ENASource(http=http).fetch_reference_ranges(records)
    assert requested(http) == [['ERR001']]
    assert len(records.linked) == 1


@pytest.mark.parametrize('response', [
    '<RUN_SET><RUN_REF accession="ERR001"/></RUN_SET>',
    '<EXPERIMENT_SET><EXPERIMENT accession="ERR001"/></EXPERIMENT_SET>',
    '<RUN_SET><RUN accession="ERR001"><ERROR>unavailable</ERROR></RUN></RUN_SET>',
    '<ERROR><RUN accession="ERR001"/></ERROR>',
    '<RUN_SET><RUN accession="ERR001"><e:ERROR xmlns:e="error">unavailable</e:ERROR></RUN></RUN_SET>',
    '<RUN_SET><RUN accession="ERR900"><IDENTIFIERS><SUBMITTER_ID>ERR001</SUBMITTER_ID></IDENTIFIERS></RUN></RUN_SET>',
])
def test_invalid_fetched_nodes_do_not_validate_range(response):
    records = range_records(('ERR001-ERR001',))
    http = HTTP(lambda *_: ET.fromstring(response))
    ENASource(http=http).fetch_reference_ranges(records)
    assert records.linked == []
    assert any('ERR001' in issue for issue in records.issues)


def test_new_provider_aliases_validate_missing_range_union():
    records = range_records(('ERR001-ERR002',), evidence='<RUN_SET><RUN accession="ERR001"/></RUN_SET>')
    http = HTTP(lambda *_: ET.fromstring('<RUN_SET><RUN accession="ERR900"><IDENTIFIERS><SECONDARY_ID>ERR002</SECONDARY_ID></IDENTIFIERS></RUN></RUN_SET>'))
    ENASource(http=http).fetch_reference_ranges(records)
    assert requested(http) == [['ERR002']]
    assert len(records.linked) == 1
    assert records.issues == []


def test_range_evidence_is_invocation_local_and_missing_members_stay_unresolved():
    records = range_records(('ERR001-ERR002', 'ERR001-ERR003'))
    http = HTTP(lambda *_: ET.fromstring('<RUN_SET><RUN accession="ERR001"/></RUN_SET>'))
    source = ENASource(http=http)
    source.fetch_reference_ranges(records)
    # Complete XML that omits an identity receives one bounded missing-only
    # recovery request. Previously received identities are never downloaded again.
    assert requested(http) == [['ERR001', 'ERR002'], ['ERR002'],
                               ['ERR002', 'ERR003'], ['ERR002', 'ERR003']]
    assert records.linked == []
    fresh = range_records(('ERR001-ERR001',))
    source.fetch_reference_ranges(fresh)
    assert requested(http)[-1] == ['ERR001']


def test_missing_records_keep_existing_hundred_record_chunks_and_literal_order():
    records = range_records(('ERR001-ERR250',), evidence='<RUN_SET><RUN accession="ERR001"/></RUN_SET>')
    http = range_http()
    ENASource(http=http).fetch_reference_ranges(records)
    assert list(map(len, requested(http))) == [100, 100, 49]
    assert requested(http)[0][0] == 'ERR002'
    assert records.linked[0]['metadata']['accessions'] == [f'ERR{i:03d}' for i in range(1, 251)]


def test_inventory_and_779_core_records_eliminate_fourteen_redundant_requests():
    """GeoMx-shaped inventory: all 785 Xref targets remain in unchanged order."""
    run_ids = [f'ERR{i:06d}' for i in range(1, 780)]
    inventories = {
        'study': [], 'sample': [{'sample_accession': f'ERS{i}'} for i in (1, 2)],
        'read_experiment': [{'experiment_accession': f'ERX{i}', 'sample_accession': f'ERS{i}'} for i in (1, 2)],
        'read_run': [{'run_accession': acc, 'sample_accession': 'ERS1'} for acc in run_ids],
        'analysis': [], 'assembly': [],
    }
    def handler(url, params, fmt):
        if 'returnFields' in url:
            return 'columnId\nrun_accession\n'
        if '/xml/' in url:
            nodes = []
            for acc in url.rsplit('/', 1)[-1].split(','):
                tag = {'ERP': 'STUDY', 'PRJ': 'PROJECT', 'ERS': 'SAMPLE', 'ERX': 'EXPERIMENT', 'ERR': 'RUN'}[acc[:3]]
                link = ('<XREF_LINK><DB>ENA-RUN</DB><ID>ERR000001-ERR000779</ID></XREF_LINK>' if acc == 'ERS1' else '')
                nodes.append(f'<{tag} accession="{acc}">{link}</{tag}>')
            return ET.fromstring('<ROOT>' + ''.join(nodes) + '</ROOT>')
        if url.endswith('/search') and '/portal/' in url:
            return deepcopy(inventories[params['result']])
        if url.endswith('/filereport'):
            return deepcopy(inventories['read_run'])
        if url.endswith('/count'):
            return {'count': '779'}
        if '/xref/' in url:
            return ([{'Target Primary Accession': 'ERP1', 'Target Secondary Accession': '',
                      'Source': 'ArrayExpress', 'Source Primary Accession': 'E-MTAB-1'}]
                    if params['accession'] == 'ERP1' else [])
        raise AssertionError((url, params, fmt))
    http = HTTP(handler)
    records = ENASource(http=http).fetch(StudySeed('ERP1', 'PRJEB1'))
    assert not any('returnFields' in url for url, _, _ in http.calls)
    fetched_run_batches = [url for url, _, _ in http.calls if '/xml/ERR' in url]
    assert len(fetched_run_batches) == 8
    targets = sorted(['ERP1', 'PRJEB1', 'ERS1', 'ERS2', 'ERX1', 'ERX2', *run_ids])
    assert [params['accession'] for url, params, _ in http.calls if '/xref/' in url] == targets
    assert len(targets) == 785
    assert records.indexed['read_run'] == inventories['read_run']
    assert records.linked[0]['metadata']['accessions'] == run_ids
    assert records.linked[1]['metadata'][0]['Source Primary Accession'] == 'E-MTAB-1'
    assert records.issues == []
