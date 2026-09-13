# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from copy import deepcopy
from unittest.mock import patch
from meta_standards_converter.magetab.semantics import _miniml_path_columns, _insert_retained_patch_comments
from meta_standards_converter.magetab.parser import AEParser
from meta_standards_converter.miniml.ena_parser import ENAParser
from tests.converters.test_ae2json import resolved_input
from tests.test_native_archive_parsers import fixture_records


def test_generic_file_uri_roundtrip():
    node={'kind':'array_data_file','name':'33308_4#4.cram','link':{'value':'ftp://example.org/33308_4%234.cram','type':'CRAM'}}
    pairs=_miniml_path_columns([node])
    assert ('Comment[File URI]',node['link']['value']) in pairs
    assert not any('FASTQ' in label for label,_ in pairs)
    header=['Source Name']+[p[0] for p in pairs];row=['sample']+[p[1] for p in pairs]
    data=AEParser().parse(resolved_input(sdrfs=['\t'.join(header)+'\n'+'\t'.join(row)+'\n'])).to_mapping()
    file=data['series']['assay_paths'][0]['steps'][-1]
    assert file['name']==node['name'] and file['link']['value']==node['link']['value']
    assert data['sample'][0]['raw_data'][0]['value']==node['link']['value']


def test_no_harmonization_does_not_decode_per_sample():
    data={'sample':[{'iid':f'S{i}'} for i in range(100)]}
    with patch('meta_standards_converter.metadata.provenance.iter_harmonization_operations', side_effect=AssertionError('unnecessary decode')):
        _insert_retained_patch_comments(data, [['SDRF File', [['Source Name'], *[[s['iid']] for s in data['sample']]]]])


def test_ena_indexed_counts_are_run_scoped_not_read_lengths():
    records=fixture_records('ena')
    records.indexed['read_run'][0].update(read_count='500',base_count='10000')
    data=ENAParser().parse(records).to_mapping()
    run=data['sample'][0]['sra_run'][0]
    assert run['indexed_statistics']=={'read_count':'500','base_count':'10000'}
    assert '500' not in run.get('read_lengths',[])


def test_retained_patches_decode_once_and_follow_samples_across_file_rows():
    from tests.converters.test_retained_patch_converters import _applied_package
    from meta_standards_converter.miniml import MINiMLCodec
    data=_applied_package()
    data['sample'].extend({'iid':f'GSM{i}','channel':[]} for i in range(2,52))
    rows=[['SDRF File',[['Source Name'],['GSM1'],['GSM1'],['GSM2']]]]
    original=MINiMLCodec.decode
    with patch.object(MINiMLCodec,'decode',autospec=True,side_effect=original) as decode:
        _insert_retained_patch_comments(data,rows)
        assert decode.call_count<=1
    table=rows[0][1]; index=table[0].index('Comment[msc_harmonization_sample_disease_name_value]')
    assert [r[index] for r in table[1:]]==['disease','disease','']
