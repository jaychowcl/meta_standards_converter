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

from tests.test_native_file_layout import package, file, projected, comments
from tests.test_protocol_export import render


def with_origin(uri, repository, **fields):
    result = file(uri, **fields)
    result['link'].update(repository=repository, source_accession='study')
    return result


def test_select_complete_linked_set_preserve_index_and_archive_alternatives():
    data = package(); base = data['series']['assay_paths'][0]['steps'][:-1]
    files = [with_origin('https://ae/'+n, 'ArrayExpress', READ_INDEX=n) for n in ('R1.fastq', 'R2.fastq', 'I1.fastq')]
    files += [with_origin('s3://native/'+n, 'SRA', ROLE='Original') for n in ('one.bz2', 'two.bz2')]
    data['series']['assay_paths'] = [{'steps':deepcopy(base)+[f]} for f in files]
    before = deepcopy(data)
    rows = render(data); table = next(r[1] for r in rows if r[0]=='SDRF File')
    assert len(table) == 4
    index = table[0].index('Comment[FASTQ_URI]')
    assert {r[index] for r in table[1:]} == {'https://ae/'+n for n in ('R1.fastq','R2.fastq','I1.fastq')}
    assert all('s3://native/one.bz2' in r and 's3://native/two.bz2' in r for r in table[1:])
    assert data == before


def test_incomplete_preferred_set_cannot_displace_available_reads():
    data = package(); base = data['series']['assay_paths'][0]['steps'][:-1]
    files = [with_origin('https://ae/R1.fastq', 'ArrayExpress'), with_origin('', 'ArrayExpress')]
    files[-1]['name'] = 'R2.fastq'
    files += [with_origin('https://native/'+n, 'ENA', ROLE='GENERATED_FILE') for n in ('R1.fastq','R2.fastq')]
    data['series']['assay_paths'] = [{'steps':deepcopy(base)+[f]} for f in files]
    paths = projected(data)
    assert {v for p in paths for k,v in comments(p) if k=='FASTQ_URI'} == {'https://native/R1.fastq','https://native/R2.fastq'}


@pytest.mark.parametrize('reverse',[False,True])
def test_preferred_exact_subset_cannot_hide_native_set_membership(reverse):
    data=package();base=data['series']['assay_paths'][0]['steps'][:-1]
    files=[with_origin('https://ena/run_1.fastq.gz','ENA',ROLE='GENERATED_FILE'),
           with_origin('https://ena/run_2.fastq.gz','ENA',ROLE='GENERATED_FILE'),
           with_origin('https://ena/run_1.fastq.gz','ArrayExpress')]
    if reverse:files.reverse()
    data['series']['assay_paths']=[{'steps':deepcopy(base)+[f]} for f in files]
    before=deepcopy(data);paths=projected(data)
    assert {v for p in paths for k,v in comments(p) if k=='FASTQ_URI'}=={'https://ena/run_1.fastq.gz','https://ena/run_2.fastq.gz'}
    assert not any(k=='ARCHIVE_FILE_URI' for p in paths for k,v in comments(p))
    assert before['series']['assay_paths']!=paths


def test_distinct_one_file_representation_is_not_a_proven_partial_pair():
    data=package();base=data['series']['assay_paths'][0]['steps'][:-1]
    files=[with_origin('https://ena/run_1.fastq.gz','ENA',ROLE='GENERATED_FILE'),
           with_origin('https://ena/run_2.fastq.gz','ENA',ROLE='GENERATED_FILE'),
           with_origin('https://ae/interleaved.fastq.gz','ArrayExpress')]
    data['series']['assay_paths']=[{'steps':deepcopy(base)+[f]} for f in files]
    paths=projected(data)
    assert {v for p in paths for k,v in comments(p) if k=='FASTQ_URI'}=={'https://ae/interleaved.fastq.gz'}


def test_subset_inventory_survives_raw_prefix_absorption_into_explicit_result():
    data=package();base=data['series']['assay_paths'][0]['steps'][:-1]
    files=[with_origin('https://ena/run_1.fastq.gz','ENA',ROLE='GENERATED_FILE'),
           with_origin('https://ena/run_2.fastq.gz','ENA',ROLE='GENERATED_FILE'),
           with_origin('https://ena/run_1.fastq.gz','ArrayExpress')]
    data['series']['assay_paths']=[{'steps':deepcopy(base)+[f]} for f in files]
    data['series']['assay_paths'][-1]['steps'] += [
        {'kind':'protocol_application','protocol_ref':'analysis'},
        {'kind':'derived_array_data_file','name':'result.tsv','link':{'value':'https://results/result.tsv'}}]
    paths=projected(data)
    assert {v for p in paths for k,v in comments(p) if k=='FASTQ_URI'}=={'https://ena/run_1.fastq.gz','https://ena/run_2.fastq.gz'}
    results=[p for p in paths if any(s['kind']=='derived_array_data_file' for s in p['steps'])]
    assert len(results)==1 and ('FASTQ_URI','https://ena/run_1.fastq.gz') in comments(results[0])


@pytest.mark.parametrize('case',['complete','checksum_conflict','ambiguous'])
def test_subset_requires_unique_compatible_file_correspondence(case):
    data=package();base=data['series']['assay_paths'][0]['steps'][:-1]
    files=[with_origin('https://ena/run_1.fastq.gz','ENA',ROLE='GENERATED_FILE',MD5='a'*32),
           with_origin('https://ena/run_2.fastq.gz','ENA',ROLE='GENERATED_FILE',MD5='b'*32),
           with_origin('https://ena/run_1.fastq.gz','ArrayExpress')]
    if case=='complete':files.append(with_origin('https://ena/run_2.fastq.gz','ArrayExpress'))
    if case=='checksum_conflict':files[-1]['comments'].append({'name':'MD5','value':'c'*32})
    if case=='ambiguous':files.insert(0,with_origin('https://ena/run_1.fastq.gz','ENA',ROLE='GENERATED_FILE',MD5='c'*32))
    data['series']['assay_paths']=[{'steps':deepcopy(base)+[f]} for f in files]
    before=deepcopy(data);rows=render(data);table=next(r[1] for r in rows if r[0]=='SDRF File')
    uri=table[0].index('Comment[FASTQ_URI]')
    assert {r[uri] for r in table[1:]}==({'https://ena/run_1.fastq.gz','https://ena/run_2.fastq.gz'} if case=='complete' else {'https://ena/run_1.fastq.gz'})
    if case=='ambiguous':
        assert not any(r[i] for r in table[1:] for i,h in enumerate(table[0]) if h=='Comment[FASTQ_MD5]')
        assert {'a'*32,'c'*32}<={r[i] for r in table[1:] for i,h in enumerate(table[0]) if h=='Comment[ARCHIVE_FILE_MD5]'}
    assert data==before


def test_sparse_occurrence_cannot_bridge_conflicting_versions_within_one_source_set():
    data=package();base=data['series']['assay_paths'][0]['steps'][:-1]
    files=[with_origin('https://ena/run_1.fastq.gz','ENA',ROLE='GENERATED_FILE',MD5='a'*32),
           with_origin('https://ena/run_1.fastq.gz','ENA',ROLE='GENERATED_FILE',MD5='b'*32),
           with_origin('https://ena/run_1.fastq.gz','ENA',ROLE='GENERATED_FILE',READ_INDEX='1')]
    data['series']['assay_paths']=[{'steps':deepcopy(base)+[f]} for f in files]
    paths=projected(data)
    facts={(dict(comments(p)).get('FASTQ_MD5'),dict(comments(p)).get('FASTQ_READ_INDEX')) for p in paths}
    assert facts=={('a'*32,None),('b'*32,None),(None,'1')}


@pytest.mark.parametrize('reverse',[False,True])
def test_excluded_archive_versions_do_not_borrow_sparse_lane_by_order(reverse):
    data=package();base=data['series']['assay_paths'][0]['steps'][:-1]
    files=[with_origin('https://ena/run_1.fastq.gz','ENA',ROLE='GENERATED_FILE',MD5='a'*32),
           with_origin('https://ena/run_1.fastq.gz','ENA',ROLE='GENERATED_FILE',MD5='b'*32),
           with_origin('https://ena/run_1.fastq.gz','ENA',ROLE='GENERATED_FILE',LANE='1')]
    if reverse:files.reverse()
    files.append(with_origin('https://ae/other.fastq.gz','ArrayExpress'))
    data['series']['assay_paths']=[{'steps':deepcopy(base)+[f]} for f in files]
    paths=projected(data);records=[];current={}
    for key,value in comments(paths[0]):
        if key.startswith('ARCHIVE_FILE_'):
            key=key.removeprefix('ARCHIVE_FILE_')
            if key=='NAME' and current:records.append(current);current={}
            current[key]=value
    if current:records.append(current)
    assert {(r.get('MD5') or None,r.get('LANE') or None) for r in records}=={('a'*32,None),('b'*32,None),(None,'1')}


def test_sparse_prefix_cannot_supply_facts_to_conflicting_result_file_versions():
    data=package();base=data['series']['assay_paths'][0]['steps'][:-1]
    data['series']['assay_paths']=[{'steps':deepcopy(base)+[with_origin('https://ena/run.fastq.gz','ENA',LANE='1')]}]
    for checksum in ('a'*32,'b'*32):
        data['series']['assay_paths'].append({'steps':deepcopy(base)+[
            with_origin('https://ena/run.fastq.gz','ENA',MD5=checksum),
            {'kind':'derived_array_data_file','name':checksum+'.tsv','link':{'value':'https://results/'+checksum+'.tsv'}}]})
    paths=projected(data)
    results=[p for p in paths if any(s['kind']=='derived_array_data_file' for s in p['steps'])]
    assert len(results)==2
    assert not any(k=='FASTQ_LANE' for p in results for k,v in comments(p))
    assert any(('FASTQ_LANE','1') in comments(p) for p in paths if p not in results)


def test_sample_results_are_source_annotations_without_arbitrary_run_binding():
    data = package(); source = deepcopy(data['series']['assay_paths'][0]['steps'][0])
    result = {'kind':'derived_array_data_file','name':'counts.tsv','link':{'value':'https://results/counts.tsv','type':'tsv'},
              'comments':[{'name':'CHECKSUM','value':'x'},{'name':'CHECKSUM_METHOD','value':'SHA256'}]}
    data['series']['assay_paths'].append({'steps':[source,result]})
    before = deepcopy(data)
    paths = projected(data)
    assert len(paths) == 3
    assert all(any(c['name']=='SAMPLE_FILE_URI' and c['value']=='https://results/counts.tsv'
                   for c in p['steps'][0].get('comments',[])) for p in paths)
    assert not any(s['kind']=='derived_array_data_file' for p in paths for s in p['steps'])
    once = deepcopy(data); projected(data); assert data == once
    assert len(before['series']['assay_paths']) > len(paths)


def test_only_globally_empty_optional_archive_columns_removed():
    data = package(); rows = render(data); table = next(r[1] for r in rows if r[0]=='SDRF File')
    assert 'Comment[SUBMITTED_FILE_CHECKSUM_METHOD]' not in table[0]
    assert 'Comment[SUBMITTED_FILE_URI]' in table[0]
    assert all(len(r)==len(table[0]) for r in table)


def test_enrichment_marks_record_origin_only_after_valid_identity_join():
    from meta_standards_converter.metadata.archive_enrichment import merge_archive_metadata
    from meta_standards_converter.miniml import MINiMLCodec
    from tests.test_archive_fidelity_enrichment import enriched_native, workflow
    extra = workflow().to_mapping()
    scan = next(s for s in extra['series']['assay_paths'][0]['steps'] if s['kind']=='scan')
    scan.setdefault('comments',[]).append({'name':'FASTQ_URI','value':'https://ae/reads.fastq'})
    result, issues = merge_archive_metadata(enriched_native(), MINiMLCodec().decode(extra).package, prefer=True)
    raw = [s for p in result.to_mapping()['series']['assay_paths'] for s in p['steps'] if s['kind']=='array_data_file']
    assert any(s.get('link',{}).get('repository')=='ArrayExpress' for s in raw)


def test_explicit_missing_mate_or_index_prevents_primary_set_selection():
    data=package();base=data['series']['assay_paths'][0]['steps'][:-1]
    files=[with_origin('https://ae/R1.fastq','ArrayExpress',READ_INDEX='1')]
    files += [with_origin('https://ena/'+n+'.fastq','ENA',ROLE='GENERATED_FILE',READ_INDEX=n) for n in ('1','2','index1')]
    data['series']['assay_paths']=[{'steps':deepcopy(base)+[f]} for f in files]
    paths=projected(data)
    assert {v for p in paths for k,v in comments(p) if k=='FASTQ_URI'}=={'https://ena/'+n+'.fastq' for n in ('1','2','index1')}


def test_sample_annotations_preserve_source_facts_and_sequential_processing():
    data=package();source=deepcopy(data['series']['assay_paths'][0]['steps'][0])
    source['characteristics'].append({'name':'aliquot','value':'processed portion'})
    first={'kind':'derived_array_data_file','name':'normalized.tsv','link':{'value':'https://results/normalized.tsv'}}
    second={'kind':'derived_array_data_file','name':'filtered.tsv','link':{'value':'https://results/filtered.tsv'}}
    data['series']['assay_paths'].append({'document':'explicit-results.sdrf','steps':[source,
        {'kind':'protocol_application','protocol_ref':'normalize'},first,
        {'kind':'protocol_application','protocol_ref':'filter'},second]})
    paths=projected(data)
    scoped=[p for p in paths if 'processed portion' in str(p)]
    assert len(scoped)==1 and len(scoped[0]['steps'])==1
    assert all('SAMPLE_FILE_URI' not in str(p) for p in paths if p not in scoped)
    c=comments(scoped[0]);processing=[v for k,v in c if k=='SAMPLE_FILE_PROCESSING']
    assert 'filter' not in processing[0] and 'filter' in processing[1]
    assert 'explicit-results.sdrf' in str(c)


def test_verified_mirror_promotion_keeps_repository_with_its_location():
    from meta_standards_converter.magetab.native_files import _combine, _record
    ena=_record(with_origin('https://ena/R1.fastq','ENA',MD5='a'*32))
    ae=_record(with_origin('https://ae/R1.fastq','ArrayExpress',MD5='a'*32))
    records=[ena];_combine(records,ae)
    assert len(records)==1
    assert records[0]['URI']=='https://ae/R1.fastq' and records[0]['_repository']=='ArrayExpress'
    assert records[0]['_alternatives'][0]['URI']=='https://ena/R1.fastq'


def test_excluded_set_keeps_nested_verified_mirrors_in_run_annotations():
    data=package();base=data['series']['assay_paths'][0]['steps'][:-1]
    files=[with_origin('https://ena/reads.fastq','ENA',MD5='a'*32,ROLE='GENERATED_FILE'),
           with_origin('https://sra/reads.fastq','SRA',MD5='a'*32,ROLE='GENERATED_FILE'),
           with_origin('https://ae/reads.fastq','ArrayExpress',MD5='b'*32)]
    data['series']['assay_paths']=[{'steps':deepcopy(base)+[f]} for f in files]
    paths=projected(data)
    assert len(paths)==1
    c=comments(paths[0])
    assert ('FASTQ_URI','https://ae/reads.fastq') in c
    assert {v for k,v in c if k=='ARCHIVE_FILE_URI'}=={'https://ena/reads.fastq','https://sra/reads.fastq'}


@pytest.mark.parametrize('key,value',[('repository',{'name':'legacy repository'}),('repository',['legacy']),('source_accession',{'accession':'legacy'}),('source_accession',['legacy'])])
def test_legacy_record_origin_extras_are_opaque_not_ranked(key,value):
    data=package();data['series']['assay_paths'][0]['steps'][-1]['link'][key]=value
    before=deepcopy(data);rows=render(data)
    assert str(value.get('name','legacy') if isinstance(value,dict) else 'legacy') in str(rows)
    assert data==before
