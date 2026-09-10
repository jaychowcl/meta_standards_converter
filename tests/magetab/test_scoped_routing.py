# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Observable contracts for structured chemistry and sample-specific routing."""
from copy import deepcopy
import pytest
from meta_standards_converter.magetab.chemistry import resolve_chemistry
from meta_standards_converter.magetab.technology import detect_ae_technology
from meta_standards_converter.magetab.sdrf.constructor import SDRFConstructor


def sample(title='', protocol='', code=None, iid='GSM1'):
    channel={'source': 'source', 'extract_protocol': protocol}
    if code is not None:
        channel['characteristics']=[{'name':'singlecell_type','value':code}]
    return {'iid':iid,'title':title,'type':'SRA','library_strategy':'RNA-Seq',
            'accession':[{'database':'GEO','value':iid}], 'channel':[channel]}


@pytest.mark.parametrize('code,family,versions', [
    ('SC3Pv1','3 prime',('1',)), ('SC3Pv2','3 prime',('2',)),
    ('SC3Pv3','3 prime',('3',)), ('SC3Pv4','3 prime',('4',)),
    ('SC3Pv3HT','3 prime',('3.1',)), ('SC5P-PE','5 prime',()),
    ('SC5P-R2','5 prime',()), ('SC5P-PE-v3','5 prime',('3',)),
    ('SC5P-R2-v3','5 prime',('3',)), ('SC5PHT','5 prime',('2',)),
])
def test_documented_identifiers(code,family,versions):
    s=sample(code=code)
    original=deepcopy(s)
    result=resolve_chemistry(s)
    assert (result.manufacturer,result.family,result.versions)==('10x Genomics',family,versions)
    assert not any('read' in name or 'barcode' in name or 'offset' in name for name,_ in result.attributes)
    assert any(e.text==code and 'characteristics[0]' in e.path for e in result.evidence)
    assert s==original


@pytest.mark.parametrize('tag',['singlecell_type','Singlecell Type','SINGLECELL-TYPE','chemistry','Library Chemistry','library-chemistry'])
def test_identifier_keys_and_case(tag):
    s=sample(code='  sc3pV2  ')
    s['channel'][0]['characteristics'][0]['name']=tag
    assert resolve_chemistry(s).versions==('2',)


@pytest.mark.parametrize('code',['auto','SC3Pv2-other','prefix SC3Pv2','SC3Pv2/SC3Pv3','SFRP'])
def test_unknown_identifiers_do_not_guess(code):
    s=sample(code=code)
    result=resolve_chemistry(s)
    assert result.family is None
    assert result.attributes==()
    assert any(e.text==code for e in result.evidence)


def test_identifier_narrows_only_compatible_alternatives():
    s=sample(code='SC3Pv2',protocol="Chromium Single Cell 3' Kit v2 or v3")
    result=resolve_chemistry(s)
    assert (result.family,result.versions)==('3 prime',('2',))
    assert any(e.value=='3' for e in result.evidence)  # source alternatives retained
    assert not any(d.code=='ambiguous_version' for d in result.diagnostics)
    s['channel'][0]['extract_protocol']="Chromium Single Cell 3' Kit v3"
    result=resolve_chemistry(s)
    assert result.versions==('2','3')
    assert any(d.code=='ambiguous_version' for d in result.diagnostics)
    s['channel'][0]['extract_protocol']="Chromium Single Cell 5' Kit v2"
    result=resolve_chemistry(s)
    assert result.family is None
    assert ('end bias','3 prime tag') not in result.attributes


def test_channel_identifiers_do_not_leak():
    s=sample(code='SC3Pv2')
    s['channel'].append(sample(code='SC5PHT')['channel'][0])
    assert resolve_chemistry(s,channel=s['channel'][1]).family=='5 prime'
    assert resolve_chemistry(s).family is None


SHARED="Library preparation for snRNA-seq: Chromium Single Cell 3' Kit v3.1. The spatial libraries were prepared using Visium Spatial Gene Expression Reagent Kits."


def test_sample_identity_over_shared_protocol():
    s=sample('Heart, scRNA-seq, Sham',SHARED)
    assert detect_ae_technology({'sample':[s]})=='droplet_single_cell_sequencing'


def test_conflicting_identity_is_generic_with_diagnostic():
    from meta_standards_converter.magetab.technology import resolve_technology
    decision=resolve_technology(sample('scRNA-seq and Visium spatial library',SHARED))
    assert decision.handler=='sequencing'
    assert any(d.code=='ambiguous_technology' for d in decision.diagnostics)
    assert decision.evidence


@pytest.mark.parametrize('reverse',[False,True])
def test_mixed_sample_rows_keep_their_technology_and_protocols(reverse):
    samples=[sample('scRNA-seq',SHARED,iid='GSM1'),sample('Visium spatial',SHARED,iid='GSM2')]
    if reverse:samples.reverse()
    data={'series':{'iid':'GSE1','sample_ref':[{'ref':s['iid']} for s in samples]},'sample':samples}
    original=deepcopy(data)
    rows=SDRFConstructor()._miniml2sdrf(data)
    values={r[0]:r for r in rows[1:]}
    column=rows[0].index('Comment[library construction]')
    assert values['GSM1'][column]=='10x 3 prime v3.1'
    assert values['GSM2'][column]=='10x Visium'
    assert [r[0] for r in rows[1:]]==[s['iid'] for s in samples]
    assert data==original


def test_run_identity_refines_dispatch():
    s=sample(protocol=SHARED)
    s['sra_run']=[{'accession':'SRR1','library_name':'scRNA-seq'}, {'accession':'SRR2','library_name':'Visium spatial'}]
    rows=SDRFConstructor()._miniml2sdrf({'sample':[s]})
    col=rows[0].index('Comment[library construction]')
    assert [r[col] for r in rows[1:]]==['10x 3 prime v3.1','10x Visium']


def test_shared_visium_cannot_name_other_spatial_library():
    s=sample('Spatial transcriptomics')
    rows=SDRFConstructor()._miniml2sdrf({'sample':[s], 'series':{'summary':'Other libraries used Visium'}},technology_type='spatial_sequencing')
    assert 'Comment[library construction]' not in rows[0]


def test_full_construction_shares_lookup_and_registry_across_mixed_samples():
    from meta_standards_converter.magetab.constructor import AEConstructor
    from meta_standards_converter.miniml import MINiMLCodec
    from meta_standards_converter.sources.insdc import INSDCWebfetcher

    class RecordedINSDC(INSDCWebfetcher):
        def __init__(self):
            self.calls=[]
        def fetch_sra_runs(self,accession):
            self.calls.append(accession)
            return [{'run':'SRR1','library_strategy':'RNA-Seq','fastq_files':[]}]

    samples=[sample('scRNA-seq',SHARED,iid='GSM1'),sample('Visium spatial',SHARED,iid='GSM2')]
    for s in samples:
        s['relation']=[{'type':'SRA','target':'SRX1'}]
        s['data_processing']='shared processing'
    data={'miniml_schema_version':'2.0','source':{'format':'test'},'series':{'iid':'GSE1','title':'mixed','sample_ref':[{'ref':s['iid']} for s in samples]},'sample':samples}
    provider=RecordedINSDC()
    converter=AEConstructor(insdc_client=provider)
    rows=converter.miniml2magetab(MINiMLCodec().decode(data).package)
    assert provider.calls==['SRX1']
    table=next(r[1] for r in rows if r[0]=='SDRF File')
    assert [r[0] for r in table[1:]]==['GSM1','GSM2']
    assert [r[table[0].index('Comment[ENA_RUN]')] for r in table[1:]]==['SRR1','SRR1']
    assert [r[table[0].index('Comment[library construction]')] for r in table[1:]]==['10x 3 prime v3.1','10x Visium']
    names=next(r[1:] for r in rows if r[0]=='Protocol Name')
    descriptions=next(r[1:] for r in rows if r[0]=='Protocol Description')
    assert descriptions.count(' '.join(SHARED.split()))==1
    assert descriptions[0]=='shared processing'
    refs=[value for row in table[1:] for label,value in zip(table[0],row) if label=='Protocol REF' and value]
    assert set(refs)<=set(names)
    assert next(r for r in rows if r[0]=='Comment[AEExperimentType]')==['Comment[AEExperimentType]']
    # Reuse of the constructor is a new operation, with fresh registry/audit/cache.
    again=converter.miniml2magetab(MINiMLCodec().decode(data).package)
    assert provider.calls==['SRX1','SRX1']
    assert rows==again


def test_conflict_audit_and_explicit_override():
    s=sample('scRNA-seq and Visium spatial library',SHARED)
    converter=SDRFConstructor()
    table=converter._miniml2sdrf({'sample':[s]})
    assert 'Comment[library construction]' not in table[0]
    assert any('ambiguous_technology' in w for w in converter.last_sdrf_audit.warnings)
    table=converter._miniml2sdrf({'sample':[s]},technology_type='tenx_v2_droplet_single_cell_sequencing')
    assert table[1][table[0].index('Comment[library construction]')]=='10x 3 prime v3.1'
    assert not any('ambiguous_technology' in w for w in converter.last_sdrf_audit.warnings)


def test_channel_identity_refines_routing():
    s=sample(protocol=SHARED)
    s['channel'][0]['characteristics']=[{'name':'assay_type','value':'scRNA-seq'}]
    other=deepcopy(s['channel'][0])
    other['characteristics'][0]['value']='Visium spatial'
    s['channel'].append(other)
    table=SDRFConstructor()._miniml2sdrf({'sample':[s]})
    col=table[0].index('Comment[library construction]')
    assert [r[col] for r in table[1:]]==['10x 3 prime v3.1','10x Visium']


def test_neighbouring_samples_and_shared_text_do_not_assign_each_other():
    from meta_standards_converter.magetab.technology import resolve_technology
    unknown=sample()
    data={'sample':[unknown,sample('Visium spatial')], 'series':{'summary':'scRNA-seq and Visium spatial libraries'}}
    result=resolve_technology(unknown,data=data)
    assert result.handler=='sequencing'
    assert result.diagnostics


@pytest.mark.parametrize('field,text',[
    ('data_processing','Cell Ranger processed 10x Visium files'),
    ('description','Not a Visium spatial library'),
])
def test_non_preparation_mentions_cannot_reenter_via_fallback(field,text):
    from meta_standards_converter.magetab.technology import resolve_technology
    s=sample()
    s[field]=text
    assert resolve_technology(s).handler=='bulk_sequencing'


def test_study_text_assigned_to_another_sample_is_not_a_fallback():
    from meta_standards_converter.magetab.technology import resolve_technology
    s=sample(iid='GSM1')
    assert resolve_technology(s,data={'series':{'overall_design':'GSM2 was prepared using Visium'}}).handler=='bulk_sequencing'


def test_structured_alternative_narrowing_does_not_cross_channels():
    s=sample(code='SC3Pv2')
    s['channel'].append(sample(protocol="Chromium Single Cell 3' Kit v2 or v3")['channel'][0])
    assert resolve_chemistry(s).versions==('2','3')


def test_homogeneous_operation_preserves_audit_and_rendering():
    s=sample('scRNA-seq',"Chromium Single Cell 3' Kit v2")
    s['channel'].append(deepcopy(s['channel'][0]))
    data={'sample':[s]}
    automatic=SDRFConstructor();forced=SDRFConstructor()
    assert automatic._miniml2sdrf(data)==forced._miniml2sdrf(data,technology_type='droplet_single_cell_sequencing')
    assert automatic.last_sdrf_audit==forced.last_sdrf_audit


@pytest.mark.parametrize('technology,key',[('expression array','array'),('high-throughput sequencing','bulk_sequencing')])
def test_empty_study_keeps_platform_rendering(technology,key):
    data={'platform':[{'iid':'GPL1','technology':technology}]}
    assert SDRFConstructor()._miniml2sdrf(data)==SDRFConstructor()._miniml2sdrf(data,technology_type=key)


def test_single_cell_array_does_not_become_sequencing():
    from meta_standards_converter.magetab.technology import resolve_technology
    s={'iid':'GSM1','title':'single-cell expression','platform_ref':{'ref':'GPL1'},'channel':[{'source':'cell'}]}
    data={'sample':[s], 'platform':[{'iid':'GPL1','technology':'expression array'}]}
    assert resolve_technology(s,data=data).handler=='array'


def test_mixed_array_does_not_suppress_sequencing_retrieval():
    from meta_standards_converter.magetab.constructor import AEConstructor
    from meta_standards_converter.miniml import MINiMLCodec
    from meta_standards_converter.sources.insdc import INSDCWebfetcher

    class RecordedINSDC(INSDCWebfetcher):
        def __init__(self):
            self.calls=[]
        def fetch_sra_runs(self,accession):
            self.calls.append(accession)
            return [{'run':'SRR1','fastq_files':[]}]

    seq=sample('scRNA-seq',SHARED)
    seq['relation']=[{'type':'SRA','target':'SRX1'}]
    array={'iid':'GSM2','platform_ref':{'ref':'GPL1'},'channel':[{'source':'array'}]}
    data={'miniml_schema_version':'2.0','source':{'format':'test'},'series':{'iid':'GSE1','title':'mixed'},
          'platform':[{'iid':'GPL1','technology':'expression array'}], 'sample':[seq,array]}
    provider=RecordedINSDC()
    rows=AEConstructor(insdc_client=provider).miniml2magetab(MINiMLCodec().decode(data).package)
    assert provider.calls==['SRX1']
    table=next(r[1] for r in rows if r[0]=='SDRF File')
    assert table[1][table[0].index('Comment[ENA_RUN]')]=='SRR1'
