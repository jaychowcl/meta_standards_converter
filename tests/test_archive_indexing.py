# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Relationship lookups must scale with supplied evidence, preserving ambiguity."""
from copy import deepcopy
from collections import Counter
import pytest

from meta_standards_converter.miniml.archive_dates import normalize_archive_dates
from meta_standards_converter.miniml.archive_administration import normalize_administration
from meta_standards_converter.miniml.archive_entities import local_platforms


class CountedList(list):
    visits = 0
    def __iter__(self):
        for value in super().__iter__():
            type(self).visits += 1
            yield value


@pytest.mark.parametrize('normalize,attribute', [
    (normalize_archive_dates, {'name':'ENA first public','value':'2020-01-01'}),
    (normalize_administration, {'name':'INSDC status','value':'public'})])
def test_administrative_path_binding_is_linear_and_leaves_ambiguous_steps(normalize, attribute):
    data={'source':{'format':'ENA'},'sample':[], 'series':{'assay_paths':[]}}
    n=48
    for i in range(n):
        iid=f'SAMEA{i}'
        data['sample'].append({'iid':iid,'channel':[]})
        data['series']['assay_paths'].append({'steps':CountedList([
            {'kind':'source','sample_ref':iid,'characteristics':[deepcopy(attribute)]},
            {'kind':'sample','characteristics':[deepcopy(attribute),deepcopy(attribute)]}])})
    ambiguous={'kind':'sample','characteristics':[deepcopy(attribute)]}
    data['series']['assay_paths'].append({'steps':CountedList([
        {'kind':'source','sample_ref':'SAMEA0'}, {'kind':'source','sample_ref':'SAMEA1'}, ambiguous])})
    CountedList.visits=0
    normalize(data)
    assert CountedList.visits < n*12
    assert ambiguous['characteristics']==[attribute]
    assert all(len(s['status'])==2 for s in data['sample'])
    before=deepcopy(data); normalize(data); assert data==before


def test_platform_association_looks_up_runs_once_and_preserves_order():
    counts=Counter()
    class Run(dict):
        def get(self,key,*args):
            counts[key]+=1
            return super().get(key,*args)
    n=48
    data={'source':{'format':'ENA'}, 'sample':[{'iid':f'SAMEA{i}','sra_run':[
        Run(run=f'ERR{i}',experiment=f'ERX{i}',instrument_model='HiSeq',instrument_platform='ILLUMINA')]} for i in range(n)]}
    local_platforms(data)
    assert counts['platform_ref'] < 6*n
    assert len(data['platform'])==1
    assert [s['sra_run'][0]['run'] for s in data['sample']]==[f'ERR{i}' for i in range(n)]


def test_sample_identity_join_extracts_aliases_once(monkeypatch):
    from meta_standards_converter.metadata import archive_enrichment as enrichment
    from meta_standards_converter.miniml import MINiMLCodec
    from tests.test_native_archive_enrichment import native
    base=native().to_mapping()
    base['series']['assay_paths']=[]
    template=base['sample'][0];base['sample']=[]
    n=24
    for i in range(n):
        sample=deepcopy(template);sample['iid']=f'SAMN{i}'
        sample['accession']=[{'value':sample['iid'],'database':'BioSample'}]
        sample['sra_accession']=[];sample['sra_run']=[];sample['relation']=[]
        base['sample'].append(sample)
    base['series']['sample_ref']=[{'ref':s['iid']} for s in base['sample']]
    package=MINiMLCodec().decode(base).package
    actual=enrichment.entity_ids;calls=[]
    def counted(*args,**kwargs):
        calls.append(1);return actual(*args,**kwargs)
    monkeypatch.setattr(enrichment,'entity_ids',counted)
    result,issues=enrichment.merge_archive_metadata(package,package)
    assert len(calls)<n*8
    assert [s.iid for s in result.samples]==[s.iid for s in package.samples]
    assert not any('ambiguous' in issue or 'unmatched' in issue for issue in issues)


def test_ena_experiment_inventory_is_indexed_once(monkeypatch):
    import xml.etree.ElementTree as ET
    from tests.test_ena_missing_record_recovery import missing_sample_records
    from meta_standards_converter.miniml.ena_parser import ENAParser
    from meta_standards_converter.miniml import insdc_support
    records=missing_sample_records();records.linked=[]
    row=records.indexed['read_run'][0]
    experiment=ET.tostring(records.xml[1][0],encoding='unicode')
    run=ET.tostring(records.xml[2][0],encoding='unicode')
    records.xml[1].clear();records.xml[2].clear()
    rows=[];n=48
    for i in range(1,n+1):
        def replace(value):
            for prefix in ('SAMEA','ERS','ERX','ERR'):
                value=value.replace(prefix+'1',prefix+str(i))
            return value
        records.xml[1].append(ET.fromstring(replace(experiment)))
        records.xml[2].append(ET.fromstring(replace(run)))
        rows.append({key:replace(value) for key,value in row.items()})
    records.indexed={'read_experiment':CountedList(rows),'read_run':rows}
    monkeypatch.setattr(insdc_support,'finish',lambda provider,records,series,samples,*args:samples)
    CountedList.visits=0
    samples=ENAParser().parse(records)
    assert CountedList.visits < n*12
    assert len(samples)==n
    assert all(len(sample['sra_run'])==1 for sample in samples)
