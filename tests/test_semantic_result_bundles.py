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
from tests.test_native_file_layout import package
from tests.test_protocol_export import render


def mex():
    data=package()
    from meta_standards_converter.miniml.archive_dates import normalize_archive_dates
    normalize_archive_dates(data)
    prefix=data['series']['assay_paths'][0]['steps']
    data['series']['protocols']=[{'name':'align','type':{'value':'alignment protocol'},
                                'description':'Cell Ranger v6.1.1 generates the gene expression matrix.'}]
    steps=deepcopy(prefix)
    for name in ('barcodes.tsv.gz','features.tsv.gz','matrix.mtx.gz'):
        steps += [{'kind':'protocol_application','protocol_ref':'align','performer':'operator','parameter_values':[{'name':'reference','value':'mm10'}]},
                  {'kind':'derived_array_data_file','name':'sample_'+name,'link':{'value':'https://results/sample_'+name,'type':'txt'},
                   'comments':[{'name':'MD5','value':'a'*32}]}]
    data['series']['assay_paths']=[{'steps':steps}]
    return data


def test_verified_mex_is_one_result_with_companions_and_minimal_residual():
    from meta_standards_converter.miniml.archive_residuals import finalize
    data=mex();original=deepcopy(data)
    records=[{'provider':'MAGE-TAB','kind':'MINiML','accession':data['series']['iid'],'metadata':deepcopy(data)}]
    result=finalize(data,records).to_mapping()
    steps=result['series']['assay_paths'][0]['steps']
    files=[s for s in steps if s['kind'].startswith('derived_')]
    assert len(files)==1 and files[0]['name']=='sample_matrix.mtx.gz'
    assert [c['role'] for c in files[0]['link']['companion_files']]==['barcodes','features']
    assert sum(s.get('protocol_ref')=='align' for s in steps)==1
    residual=result['extensions']['insdc']['records']
    assert any(r['kind']=='result_layout' for r in residual)
    assert all(set(r['metadata']) <= {'database'} for r in residual if r['kind']=='MINiML')
    assert finalize(result).to_mapping()==result
    table=next(r[1] for r in render(result) if r[0]=='SDRF File')
    assert table[0].count('Derived Array Data File')==1
    assert 'Comment[MATRIX_BARCODES_URI]' in table[0] and 'Comment[MATRIX_FEATURES_MD5]' in table[0]
    assert 'sample_barcodes.tsv.gz' in str(table) and 'sample_features.tsv.gz' in str(table)
    assert len(original['series']['assay_paths'][0]['steps'])>len(steps)


@pytest.mark.parametrize('change',['method','negated','negated_clause','negated_tool','existing','parameter','performer','directory','prefix','missing'])
def test_mex_requires_common_result_identity_and_processing_event(change):
    from meta_standards_converter.miniml.archive_results import normalize_result_bundles
    data=mex();steps=data['series']['assay_paths'][0]['steps']
    if change=='method':data['series']['protocols'][0]['description']='Unrelated processing'
    if change=='negated':data['series']['protocols'][0]['description']='Processed with STAR, not Cell Ranger.'
    if change=='negated_clause':data['series']['protocols'][0]['description']='These data were not processed with Cell Ranger.'
    if change=='negated_tool':data['series']['protocols'][0].update(description='Cell Ranger was not used.',software=['Cell Ranger 6.1.1'])
    if change=='existing':steps[-1]['link']['companion_files']=[{'role':'qc','node':{'kind':'derived_array_data_file','name':'qc.json'}}]
    if change=='parameter':steps[-4]['parameter_values'][0]['value']='different-reference'
    if change=='performer':steps[-4]['performer']='different-operator'
    if change=='directory':steps[-3]['link']['value']='https://other/sample_features.tsv.gz'
    if change=='prefix':steps[-3]['name']='another_features.tsv.gz'
    if change=='missing':del steps[-4:-2]
    before=deepcopy(data);assert normalize_result_bundles(data)==[] and data==before


def test_saved_native_export_copy_normalizes_without_changing_input():
    data=mex();before=deepcopy(data);table=next(r[1] for r in render(data) if r[0]=='SDRF File')
    assert table[0].count('Derived Array Data File')==1
    assert data==before


@pytest.mark.parametrize('field,value', [('link','unparsed source value'),('comments',['unparsed']),('name',{'unknown':'value'})])
def test_malformed_companion_extras_are_preserved_without_interpretation(field,value):
    from meta_standards_converter.miniml.archive_results import normalize_result_bundles, comparison_view
    data=mex();normalize_result_bundles(data)
    data['series']['assay_paths'][0]['steps'][-1]['link']['companion_files'][0]['node'][field]=value
    before=deepcopy(data)
    render(data);assert data==before
    assert comparison_view(data)==data


def test_assembly_support_annotations_are_idempotent_and_do_not_borrow_aliquot_facts():
    from meta_standards_converter.magetab.native_files import project_native_files
    data=package();sid=data['sample'][0]['iid']
    data['sample'][0]['supplementary_data']=[{'value':'https://assembly/stats.txt','type':'assembly report','assembly_accession':'GCA_1.2'}]
    for i,path in enumerate(data['series']['assay_paths']):
        source=path['steps'][0];source['name']='aliquot'+str(i%2)
        source['characteristics']=[{'name':'aliquot','value':str(i%2)}]
    reversed_data=deepcopy(data);reversed_data['series']['assay_paths'].reverse()
    project_native_files(data);once=deepcopy(data);project_native_files(data);assert data==once
    project_native_files(reversed_data)
    def supported(d):
        return [s for p in d['series']['assay_paths'] for s in p['steps'] if any(c.get('name')=='SAMPLE_FILE_URI' for c in s.get('comments',[]))]
    assert supported(data)==supported(reversed_data)
    assert len(supported(data))==1
    assert supported(data)[0]['sample_ref']==sid and not supported(data)[0].get('characteristics')


def test_legacy_native_without_paths_keeps_assembly_support_on_source():
    data=package();data['series']['assay_paths']=[]
    data['sample'][0]['supplementary_data']=[{'value':'https://data/stats.txt','type':'assembly report','assembly_accession':'GCA1'}]
    before=deepcopy(data);table=next(r[1] for r in render(data) if r[0]=='SDRF File')
    assert 'Comment[SAMPLE_FILE_URI]' in table[0] and 'https://data/stats.txt' in str(table)
    assert not any(r[i]=='https://data/stats.txt' or r[i]=='stats.txt' for r in table[1:] for i,h in enumerate(table[0]) if h.startswith('Derived'))
    assert data==before
