# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Conflict provenance survives artifacts, scopes and repeated imports."""
from meta_standards_converter.converters import PreparationPolicy
import csv,io,json
from copy import deepcopy
import pytest
from meta_standards_converter.converters import Converter
from meta_standards_converter.magetab.parser import AEParser
from meta_standards_converter.sources.magetab import MAGETabInput,TextResource
from meta_standards_converter.miniml import MINiMLCodec
from tests.converters.test_audit_source_precedence import package


def record(path='/sample/SRS1/sra_run/SRR1/library_layout'):
    return {'provider':'ArrayExpress','kind':'source_conflict','accession':'E-MTAB-1','metadata':{
        'field':'library_layout','path':path,'statements':[{'provider':'SRA','accession':'SRP1','value':'SINGLE'},
        {'provider':'ArrayExpress','accession':'E-MTAB-1','value':'PAIRED'}],
        'selected':{'provider':'ArrayExpress','accession':'E-MTAB-1','value':'PAIRED'},'reason':'AE > GEO > native','mode':'standard'}}


def conflicted():
    data=package();data['extensions']={'insdc':{'version':'1.0','records':[record(),record()]}}
    return data


def import_tables(tables):
    idf=[];sdrf=next(r[1] for r in tables if r[0]=='SDRF File')
    def text(rows):
        out=io.StringIO();csv.writer(out,delimiter='\t').writerows(rows);return out.getvalue()
    for row in tables:idf.append(['SDRF File','study.sdrf.txt'] if row[0]=='SDRF File' else row)
    source=MAGETabInput(TextResource('study.idf.txt',text(idf),'recorded'),
        (TextResource('study.sdrf.txt',text(sdrf),'recorded'),),'recorded','local')
    parser=AEParser();value=parser.parse(source)
    return value,parser.warnings


def records(data):
    if hasattr(data,'to_mapping'):data=data.to_mapping()
    return [r for r in data.get('extensions',{}).get('insdc',{}).get('records',[]) if r.get('kind')=='source_conflict']


def test_facade_roundtrip_keeps_conflicts_without_refetch_or_reselection():
    data=conflicted()
    for _ in range(2):
        result=Converter().convert(data, out_type='magetab', preparation=PreparationPolicy(enrichment='off', expand_studies=False))
        assert result.status=='partial'
        assert result.items[0].validation=='valid'
        assert any('recorded source conflict' in d.message for d in result.items[0].diagnostics)
        tables=result.items[0].payload[0]
        assert any(r[0]=='Comment[MSCConflict]' for r in tables)
        data,warnings=import_tables(tables)
        assert len(records(data))==1
        assert records(data)[0]==record()
        MINiMLCodec().decode(data,strict=True)


@pytest.mark.parametrize('value',['not JSON',json.dumps({'schema_version':'2.0','record':record()}),json.dumps({'schema_version':'1.0','record':{}})])
def test_unrecognized_conflict_comments_survive_with_diagnostics(value):
    tables=Converter().convert(package(), out_type='magetab', preparation=PreparationPolicy(enrichment='off')).items[0].payload[0]
    tables.append(['Comment[MSCConflict]',value])
    parsed,warnings=import_tables(tables)
    assert not records(parsed)
    assert any('MSCConflict' in w for w in warnings)
    emitted=Converter().convert(parsed, out_type='magetab', preparation=PreparationPolicy(enrichment='off')).items[0].payload[0]
    assert value in [v for row in emitted if row[0]=='Comment[MSCConflict]' for v in row[1:]]


@pytest.mark.parametrize('fmt',['tsv','csv'])
def test_conflicts_project_only_to_affected_sample_or_study(fmt):
    data=conflicted();other=deepcopy(data['sample'][0]);other['iid']='SRS2';other['accession']=[{'value':'SRS2'}];other['sra_run'][0].update(run='SRR2',experiment='SRX2',sample='SRS2')
    data['sample'].append(other);data['series']['sample_ref'].append({'ref':'SRS2'})
    data['extensions']['insdc']['records'].append(record('/series/title'))
    result=Converter().convert(data, out_type=fmt, preparation=PreparationPolicy(enrichment='off'))
    table=next(iter(result.items[0].payload.values()));rows={r['msc.sample.accession']:r for r in table.rows}
    assert len(json.loads(rows['SRS1']['msc.provenance.conflicts']))==2
    assert len(json.loads(rows['SRS2']['msc.provenance.conflicts']))==1


def test_anndata_conflicts_are_portable(tmp_path):
    import anndata,numpy as np
    from meta_standards_converter.converters.json2h5ad import JSON2H5ADConverter
    from meta_standards_converter.expression.assets import Asset
    data=conflicted();source=tmp_path/'matrix.tsv';source.write_text('gene\tcell\na\t1\n')
    adata=anndata.AnnData(np.ones((1,1)))
    JSON2H5ADConverter().normalizer.normalize(adata,data['sample'][0],data,'SRP1',Asset('SRS1',str(source),'matrix'),[],tmp_path)
    assert json.loads(adata.obs['msc.provenance.conflicts'].iloc[0])==[record()]
    assert json.loads(adata.uns['meta_standards_converter']['conflicts'])==[record()]
    adata.write_h5ad(tmp_path/'out.h5ad')
    assert json.loads(anndata.read_h5ad(tmp_path/'out.h5ad').obs['msc.provenance.conflicts'].iloc[0])==[record()]


def test_conflict_values_are_not_whitespace_normalized_by_projection():
    data=conflicted();data['extensions']['insdc']['records'][0]['metadata']['statements'][0]['value']='two  spaces'
    result=Converter().convert(data, out_type='tsv', preparation=PreparationPolicy(enrichment='off'))
    rows=next(iter(result.items[0].payload.values())).rows
    assert json.loads(rows[0]['msc.provenance.conflicts'])[0]['metadata']['statements'][0]['value']=='two  spaces'


def test_unrecognized_conflict_summary_remains_source_text():
    tables=Converter().convert(package(), out_type='magetab', preparation=PreparationPolicy(enrichment='off')).items[0].payload[0]
    tables.append(['Comment[MSCConflict]',json.dumps({'schema_version':'99.0','record':record()})])
    tables.append(['Comment[MSCConflictSummary]','Original summary for the unknown version'])
    for _ in range(2):
        parsed,warnings=import_tables(tables)
        tables=Converter().convert(parsed, out_type='magetab', preparation=PreparationPolicy(enrichment='off')).items[0].payload[0]
        assert [v for r in tables if r[0]=='Comment[MSCConflictSummary]' for v in r[1:]]==['Original summary for the unknown version']
