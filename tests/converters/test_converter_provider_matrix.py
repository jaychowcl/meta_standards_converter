# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Recorded native metadata preserves scientific facts across presets and readers."""
from copy import deepcopy
from types import SimpleNamespace
import json

import pytest

from meta_standards_converter.converters import Converter, InputSpec
from meta_standards_converter.converters.archive_results import ArchiveImportResult, StudyImportOutcome
from meta_standards_converter.miniml.sra_parser import SRAParser
from meta_standards_converter.miniml.ena_parser import ENAParser
from meta_standards_converter.converters.sra2json import SRA2JSONConverter
from meta_standards_converter.converters.ena2json import ENA2JSONConverter
from meta_standards_converter.sources.archive_support import Resolution
from tests.test_native_archive_parsers import fixture_records


@pytest.mark.parametrize('provider', ['sra','ena'])
@pytest.mark.parametrize('preset', ['standard','curators','off'])
@pytest.mark.parametrize('reader', ['accession','records','file','object'])
@pytest.mark.parametrize('target', ['json','magetab','csv','tsv'])
def test_native_provider_reader_preset_matrix(provider,preset,reader,target,tmp_path):
    records=fixture_records(provider)
    parser=SRAParser() if provider=='sra' else ENAParser()
    original=parser.parse(records)
    calls=[]
    class Source:
        def resolve(self,accession):
            calls.append('resolve')
            return Resolution(studies=[records.seed])
        def fetch(self,seed):
            calls.append('fetch')
            return deepcopy(records)
    class Enricher:
        def enrich_selected(self,data,*,publications,run_metadata):
            assert publications and not run_metadata
            calls.append('publications')
            return data
    class Linked:
        def enrich(self,data):
            calls.append('linked')
            return data,[]
    class Peer:
        def convert(self,accession,**kwargs):
            assert accession==records.seed.study
            calls.append('peer')
            return ArchiveImportResult(accession,[StudyImportOutcome(records.seed.study,records.seed.primary,'complete',original)])
    cls=SRA2JSONConverter if provider=='sra' else ENA2JSONConverter
    services={provider+'2json':cls(source=Source(),publication_enricher=Enricher()),
              ('ena' if provider=='sra' else 'sra')+'2json':Peer(),
              'enricher':Enricher(),'linked_enricher':Linked()}
    path=tmp_path/'source.json';path.write_text(json.dumps(original.to_mapping()))
    source={'accession':InputSpec(records.seed.study,in_type=provider+'_accession'),
            'records':InputSpec(records,in_type=provider+'_records'),'file':path,'object':original}[reader]
    result=Converter(services=services).convert(source,out_type=target,enrichment=preset,options={'expand_studies':False})
    item=result.items[0]
    assert item.execution=='succeeded',result.to_dict()
    assert not any(r['status'] in {'failed','partial'} for r in item.preparation),result.to_dict()
    assert item.provider==provider
    assert calls.count('publications')==(preset!='off')
    assert calls.count('peer')==(preset=='standard')
    assert ('linked' not in calls) if preset!='standard' else True
    assert calls.count('fetch')==(reader=='accession')
    if target=='json':
        sample=item.payload[0].to_mapping()['sample'][0]
        assert sample['channel'][0]['organism'][0]['taxid']=='408170'
        assert any(x['name']=='host' and x['value']=='Homo sapiens' for x in sample['channel'][0]['characteristics'])
        assert sample['sra_run'][0]['run']=='SRR11192680'
        assert sample['sra_run'][0]['library_layout']=='PAIRED'
    assert original==parser.parse(records)


@pytest.mark.parametrize('preset',['standard','curators','off'])
@pytest.mark.parametrize('reader',['accession','file','object'])
@pytest.mark.parametrize('target',['json','magetab'])
def test_biostudies_preset_reader_equivalence(preset,reader,target):
    from pathlib import Path
    from meta_standards_converter.converters.ae2json import AE2JSONConverter
    from meta_standards_converter.sources.magetab import AEWebFetcher
    path=Path(__file__).parents[1]/'fixtures/studies/E-MTAB-6486/inputs/E-MTAB-6486.idf.txt'
    resolved=AEWebFetcher().resolve(str(path))
    calls=[]
    class Fetcher:
        def resolve(self,*args,**kwargs):return resolved
    class Enricher:
        def enrich_selected(self,data,*,publications,run_metadata):
            calls.append((publications,run_metadata));return data
    ae=AE2JSONConverter(fetcher=Fetcher())
    original=ae.convert('E-MTAB-6486')
    source={'accession':'E-MTAB-6486','file':path,'object':original}[reader]
    result=Converter(services={'ae2json':ae,'enricher':Enricher()}).convert(source,
        out_type=target,enrichment=preset,options={'expand_studies':False})
    assert result.status=='complete',result.to_dict()
    assert result.items[0].provider=='biostudies'
    assert calls==([] if preset=='off' else [(True,preset=='standard')])
    if target=='json':
        assert result.items[0].payload==original
        assert len(original[0].samples)==6
