# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""File sentinels are not biological missing values or filenames."""
from copy import deepcopy
from meta_standards_converter.miniml.geo_parser import GEOParser
from tests.miniml.test_geo_parser import miniml_body
from tests.test_native_file_layout import package, file, render


def test_geo_none_file_references_do_not_become_files():
    xml=miniml_body('<Sample iid="GSM1"><Supplementary-Data type="unknown">NONE</Supplementary-Data><Supplementary-Data type="TXT">NONE.txt</Supplementary-Data><Supplementary-Data type="TXT">https://example.org/NONE</Supplementary-Data><Channel><Characteristics tag="treatment">NONE</Characteristics></Channel></Sample><Series iid="GSE1"><Sample-Ref ref="GSM1"/><Supplementary-Data> none </Supplementary-Data></Series>')
    data=GEOParser().parse(xml)[0].to_mapping()
    assert [x['value'] for x in data['sample'][0]['supplementary_data']]==['NONE.txt','https://example.org/NONE']
    assert not data['series'].get('supplementary_data')
    assert data['sample'][0]['channel'][0]['characteristics'][0]['value']=='NONE'


def test_saved_native_none_paths_do_not_add_rows_or_change_input():
    data=package();before_rows=render(deepcopy(data));source=deepcopy(data['series']['assay_paths'][0]['steps'][0])
    data['sample'][0]['supplementary_data']=[{'value':'NONE','type':'unknown'}]
    data['series']['assay_paths'].append({'steps':[source,{**file('NONE','unknown'),'kind':'derived_array_data_file'}]})
    before=deepcopy(data);rows=render(data)
    assert rows==before_rows
    assert data==before
