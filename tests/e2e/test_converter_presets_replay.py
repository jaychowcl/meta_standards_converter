# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Exercise the facade with recorded scientific evidence, across loading routes."""
from pathlib import Path
import json

import pytest

from meta_standards_converter.converters import Converter, GEO2AEConverter
from meta_standards_converter.miniml.geo_parser import GEOParser
from tests.support.contracts import GEO, assert_expected


@pytest.mark.parametrize('preset', ['standard','curators','off'])
@pytest.mark.parametrize('source_kind', ['accession','file','object'])
@pytest.mark.parametrize('target', ['json','magetab','tsv','csv'])
def test_geo_facade_routes_have_scientific_evidence_and_one_preparation_pass(preset, source_kind, target, workspace, replay):
    path=GEO/'inputs/geo.xml'
    source={'accession':'GSE328265','file':path,'object':GEOParser().parse(path.read_text())}[source_kind]
    options={'expand_studies':False}
    if target=='magetab':options['platform_handler']='single_cell_sequencing'
    result=Converter().convert(source, workspace/'out', out_type=target, enrichment=preset, options=options)
    assert result.status=='complete', result.to_dict()
    expected=['geo:GSE328265'] if source_kind=='accession' else []
    if preset!='off':expected+=['pubmed:42129775']
    if preset=='standard':expected+=['sra:SRX32831930','ena:SRX32831930']
    assert replay==expected
    assert result.items[0].dataset_ids==(('E-GEOD-328265',) if target=='magetab' else ('GSE328265',))
    assert result.items[0].provider=='geo'
    if target=='json':
        data=result.items[0].payload[0].to_mapping()
        assert data['sample'][0]['iid']=='GSM9651991'
        assert data['sample'][0]['channel'][0]['organism'][0]['value']=='Homo sapiens'
        if preset=='standard':assert data['sample'][0]['sra_run'][0]['run']
        if preset!='off':assert data['series']['pubmed_publication'][0]['title']
    if preset=='standard' and target=='magetab':
        assert_expected(GEO,'geo2ae',workspace/'out',workspace)


def test_direct_geo_magetab_equals_facade_with_equivalent_preparation(workspace,replay):
    direct=GEO2AEConverter().convert('GSE328265',platform_handler='single_cell_sequencing')
    direct_calls=list(replay);replay.clear()
    result=Converter().convert('GSE328265',out_type='magetab',
        options={'expand_studies':False,'platform_handler':'single_cell_sequencing'})
    assert result.status=='complete',result.to_dict()
    assert result.items[0].payload==direct
    assert replay==direct_calls
