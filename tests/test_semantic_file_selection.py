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
