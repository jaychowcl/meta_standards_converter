# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""New library conflict evidence is stable across interpreter hash seeds."""
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest


PROBE = r'''
from copy import deepcopy
import json, socket
def no_network(*args,**kwargs): raise AssertionError('network access forbidden')
socket.socket.connect=socket.create_connection=no_network
from meta_standards_converter.metadata.reconciliation import SourceSelection
from meta_standards_converter.miniml.archive_libraries import resolve_library_facts
native={'library_source':'TRANSCRIPTOMIC','library_strategy':'RNA-Seq','library_selection':'other'}
incoming={'library_source':'NON GENOMIC','library_strategy':'partial cDNA sequencing','library_selection':'3 prime mRNA pull down'}
data={'source':{'format':'SRA'},'series':{'iid':'SRP1'}}
extra={'source':{'format':'MAGE-TAB'},'series':{'iid':'E-MTAB-1'}}
before=deepcopy([native,incoming,data,extra]);issues=[]
selection=SourceSelection(data,extra,issues,mode='standard')
selected=[]
for run in ('SRR2','SRR1'):
 selected.append(resolve_library_facts(native,[incoming],issues,run,selection=selection,path=f'/sample/SRS1/sra_run/{run}'))
assert [native,incoming,data,extra]==before
same_source_issues=[]
unchanged=resolve_library_facts(native,[incoming,native],same_source_issues,'SRR3')
assert unchanged==native
print(json.dumps({'selected':selected,'records':selection.records,'issues':issues,'evidence':selection.evidence(),'same_source_issues':same_source_issues},sort_keys=True))
'''


@pytest.mark.fake_process
@pytest.mark.parametrize('hash_seed', ['0','1','2','42'])
def test_library_conflicts_have_stable_full_evidence_and_run_encounter_order(hash_seed):
    completed=subprocess.run([sys.executable,'-c',PROBE],cwd=Path(__file__).parents[2],
        env={**os.environ,'PYTHONHASHSEED':hash_seed,'PATH':''},capture_output=True,text=True,check=True,timeout=30)
    actual=json.loads(completed.stdout)
    values=[('library_selection','other','3 prime mRNA pull down'),
            ('library_source','TRANSCRIPTOMIC','NON GENOMIC'),
            ('library_strategy','RNA-Seq','partial cDNA sequencing')]
    selected=dict((field,new) for field,_,new in values)
    reason='AE > GEO > native; retain current value for equal source rank'
    records,issues,selections=[],[],{}
    for run in ('SRR2','SRR1'):
        for field,old,new in values:
            path=f'/sample/SRS1/sra_run/{run}/{field}'
            left={'provider':'SRA','accession':'SRP1','value':old}
            right={'provider':'ArrayExpress','accession':'E-MTAB-1','value':new}
            records.append({'provider':'ArrayExpress','kind':'source_conflict','accession':'E-MTAB-1',
                'metadata':{'field':field,'path':path,'statements':[left,right],
                            'selected':right,'reason':reason,'mode':'standard'}})
            issues.append(f"source conflict at {path}: SRA SRP1 {old!r} versus ArrayExpress E-MTAB-1 {new!r}; selected ArrayExpress {new!r} ({reason})")
            selections[path]=right
    assert actual=={'selected':[selected,selected],'records':records,'issues':issues,
        'evidence':records+[{'provider':'SRA','kind':'source_selection','accession':'SRP1','metadata':{'selections':selections}}],
        'same_source_issues':[f'SRR3: conflicting same-source {field} facts' for field,_,_ in values]}


CURATOR_PROBE = r'''
from copy import deepcopy
import json, socket
from meta_standards_converter.metadata.reconciliation import inspect_recorded_conflicts
from meta_standards_converter.miniml import MINiMLCodec
from tests.converters.test_audit_source_precedence import package
def no_network(*args,**kwargs): raise AssertionError('network access forbidden')
socket.socket.connect=socket.create_connection=no_network
data=package();extra=package('MAGE-TAB')
data['sample'][0]['sra_run'][0].update(library_source='TRANSCRIPTOMIC',library_strategy='RNA-Seq',library_selection='other')
extra['sample'][0]['sra_run'][0].update(library_source='NON GENOMIC',library_strategy='partial cDNA sequencing',library_selection='3 prime mRNA pull down')
extra['series']['assay_paths']=[{'steps':[{'kind':kind,'name':name,'sample_ref':'SRS1',
    'comments':[{'name':'LIBRARY_SELECTION','value':value}]}]}
    for kind,name,value in [('assay','SRX1','experiment-bound'),('scan','SRR1','run-first'),('scan','SRR1','run-second')]]
def recorded(field):
 return {'provider':'ArrayExpress','kind':'source_conflict','accession':'E-MTAB-1','metadata':{
  'field':field,'path':'/retained/'+field,'statements':[{'value':'before'},{'value':'after'}],
  'selected':{'value':'before'},'reason':'retained','mode':'curators'}}
prefix=[recorded('z_existing'),recorded('a_existing'),{'provider':'ArrayExpress','kind':'MINiML','accession':'E-MTAB-1','metadata':extra}]
data['extensions']={'insdc':{'version':'1.0','records':prefix}}
p=MINiMLCodec().decode(data,strict=True).package;before=deepcopy(p.to_mapping())
result,issues=inspect_recorded_conflicts(p,'curators');after=result.to_mapping()
assert p.to_mapping()==before
assert after['sample']==before['sample']
assert after['extensions']['insdc']['records'][:3]==before['extensions']['insdc']['records']
print(json.dumps({'records':after['extensions']['insdc']['records'][3:],'issues':issues},sort_keys=True))
'''


@pytest.mark.fake_process
@pytest.mark.parametrize('hash_seed', ['0','1','2','42'])
def test_curators_new_conflicts_keep_field_scope_order_and_imported_prefix(hash_seed):
    completed=subprocess.run([sys.executable,'-c',CURATOR_PROBE],cwd=Path(__file__).parents[2],
        env={**os.environ,'PYTHONHASHSEED':hash_seed,'PATH':''},capture_output=True,text=True,check=True,timeout=30)
    actual=json.loads(completed.stdout)
    reason='primary input value retained'
    values=[('library_selection','other',value) for value in
            ('3 prime mRNA pull down','run-first','run-second','experiment-bound')]
    values += [('library_source','TRANSCRIPTOMIC','NON GENOMIC'),
               ('library_strategy','RNA-Seq','partial cDNA sequencing')]
    records=[]
    issues=[f"recorded source conflict at /retained/{field}: [{{'value': 'before'}}, {{'value': 'after'}}]; selected {{'value': 'before'}}"
            for field in ('z_existing','a_existing')]
    for field,old,new in values:
        path=f'/sample/SRS1/sra_run/SRR1/{field}'
        left={'provider':'SRA','accession':'SRP1','value':old}
        right={'provider':'ArrayExpress','accession':'E-MTAB-1','value':new}
        records.append({'provider':'ArrayExpress','kind':'source_conflict','accession':'E-MTAB-1',
            'metadata':{'field':field,'path':path,'statements':[left,right],
                        'selected':left,'reason':reason,'mode':'curators'}})
        issues.append(f"source conflict at {path}: SRA SRP1 {old!r} versus ArrayExpress E-MTAB-1 {new!r}; selected SRA {old!r} ({reason})")
    assert actual=={'records':records,'issues':issues}
