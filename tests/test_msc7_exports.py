# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Observable MSC 7 export contracts, independent of provider access."""
import copy
import csv
import inspect
import json
import importlib

import pytest

from meta_standards_converter.converters import JSON2AEConverter, JSON2TSVConverter, JSON2H5ADConverter, JSON2OBSConverter
from meta_standards_converter.miniml import MINiMLCodec, iter_harmonized_values
from meta_standards_converter.magetab.constructor import AEConstructor
from meta_standards_converter.sources.json import JSONPackageSource
from meta_standards_converter.expression.checkpoints import ProcessedCheckpointStore

PROFILE = {"schema_version": "1.0", "replacements": {"disease": ["missing", "disease"]}}

def package():
    return {
        "miniml_schema_version": "3.0", "source": {"format": "test"},
        "series": {"iid": "GSE1", "sample_ref": [{"ref": "GSM1"}]},
        "sample": [{"iid": "GSM1", "channel": [{
            "source": {"value": "lung", "hz_tissue": "lung", "hz_tissue_id": "UBERON:0002048"},
            "characteristics": [
                {"name": "disease", "value": "IPF"},
                {"name": "hz_disease", "value": "idiopathic pulmonary fibrosis"},
                {"name": "hz_disease_id", "value": "MONDO:0002771"},
                {"name": "hz_disease_onto", "value": "MONDO"},
                {"name": "hz_disease_hierarchy_depth", "value": 0},
                {"name": "hz_disease(1)", "value": "pulmonary fibrosis"},
                {"name": "hz_disease_id(1)", "value": "MONDO:0003782"},
            ],
        }]}],
    }

def table(rows):
    return next(row[1] for row in rows if row[0] == "SDRF File")

def test_generated_sdrf_exposes_local_groups_and_companions():
    original = package()
    before = copy.deepcopy(original)
    header, row = table(AEConstructor().miniml2magetab(MINiMLCodec().decode(original).package, platform_handler="generic"))
    i = header.index("Characteristics[hz_disease]")
    assert header[i:i+4] == ["Characteristics[hz_disease]", "Term Source REF", "Term Accession Number", "Comment[hz_disease_hierarchy_depth]"]
    assert row[i:i+4] == ["idiopathic pulmonary fibrosis", "MONDO", "MONDO:0002771", 0]
    assert row[header.index("Characteristics[hz_disease(1)]")] == "pulmonary fibrosis"
    assert row[header.index("Characteristics[hz_tissue]")] == "lung"
    assert row[header.index("Characteristics[disease]")] == "IPF"
    assert original == before

@pytest.mark.parametrize("converter,method", [(JSON2AEConverter,"convert"),(JSON2TSVConverter,"convert_source"),(JSON2TSVConverter,"export_manifest"),(JSON2H5ADConverter,"convert"),(JSON2H5ADConverter,"convert_source"),(JSON2OBSConverter,"convert")])
def test_direct_profile_is_keyword_only(converter, method):
    parameters = inspect.signature(getattr(converter, method)).parameters
    assert parameters["replacement_profile"].kind == inspect.Parameter.KEYWORD_ONLY
    assert "use_harmonization_overrides" not in parameters

def test_native_profile_changes_standard_value_retaining_hz(tmp_path):
    source = tmp_path / "input.json"
    source.write_text(json.dumps(package()))
    rows = JSON2AEConverter().convert(str(source), enrich=False, platform_handler="generic", replacement_profile=PROFILE)
    header, row = table(rows[0])
    assert row[header.index("Characteristics[disease]")] == "idiopathic pulmonary fibrosis"
    assert "Characteristics[hz_disease]" in header
    assert json.loads(source.read_text()) == package()
    destination = tmp_path / "out.tsv"
    JSON2TSVConverter().convert_source(source, destination, replacement_profile=PROFILE)
    record = next(csv.DictReader(destination.open(), delimiter="\t"))
    assert "msc.characteristics.hz_disease" in record
    assert not any("harmonized_" in key for key in record)

def test_embedded_profile_rejected(tmp_path):
    source = tmp_path / "old.json"
    source.write_text(json.dumps({"miniml_json": package(), "harmonization_overrides": PROFILE}))
    with pytest.raises(ValueError, match="replacement.profile"):
        JSONPackageSource().load(source)

@pytest.mark.parametrize("command", ["json2ae", "json2tsv", "json2h5ad", "json2obs"])
def test_cli_direct_profile_and_invalid_input(command, tmp_path):
    extra = ["--outdir", str(tmp_path / "output")] if command == "json2obs" else []
    cli = importlib.import_module("meta_standards_converter.cli." + command)
    with pytest.raises(SystemExit) as error:
        cli.main([str(tmp_path / "missing.json"), "--replacement-profile", "[]", *extra])
    assert error.value.code == 2
    args = cli._parser().parse_args(["input.json", "--replacement-profile", json.dumps(PROFILE), *extra])
    assert args.replacement_profile == json.dumps(PROFILE)


def reparse(rows):
    import io
    from meta_standards_converter.magetab.parser import AEParser
    from meta_standards_converter.sources.magetab import MAGETabInput, TextResource
    def tsv(items):
        stream = io.StringIO()
        csv.writer(stream, delimiter='\t', lineterminator='\n').writerows(items)
        return stream.getvalue()
    idf = copy.deepcopy(rows)
    sdrf = table(idf)
    next(row for row in idf if row[0] == 'SDRF File')[1] = 'study.sdrf.txt'
    return AEParser().parse(MAGETabInput(TextResource('study.idf.txt', tsv(idf), 'local'), (TextResource('study.sdrf.txt', tsv(sdrf), 'local'),), 'study', 'local'))


def test_explicit_assay_values_units_and_factors_round_trip():
    p = package()
    p['series']['protocols'] = [{'name': 'P1', 'type': {'value': 'treatment'}}]
    p['series']['assay_paths'] = [{'steps': [
        {'kind': 'sample', 'name': 'GSM1', 'sample_ref': 'GSM1',
         'factor_values': [{'name': 'condition', 'value': 'disease', 'hz_condition': 'fibrosis', 'hz_condition_id': 'MONDO:1'}]},
        {'kind': 'protocol_application', 'protocol_ref': 'P1', 'parameter_values': [
            {'name': 'duration', 'value': '30', 'hz_value': '30', 'hz_value_hierarchy_depth': 0,
             'unit': {'value': 'min', 'hz_unit': 'minute', 'hz_unit_id': 'UO:0000031', 'hz_unit_onto': 'UO'}}]},
        {'kind': 'assay', 'name': 'A1'},
    ]}]
    result = AEConstructor().miniml2magetab(MINiMLCodec().decode(p).package, platform_handler='generic')
    header = table(result)[0]
    assert 'Characteristics[hz_disease]' in header
    assert 'Factor Value[hz_condition]' in header
    assert 'Parameter Value[hz_value]' in header
    assert 'Comment[hz_unit]' in header
    parsed = reparse(result)
    applications = [s for path in parsed.series.assay_paths for s in path.steps if s.kind == 'protocol_application']
    parameter = next(v for a in applications for v in a.parameter_values if v.name == 'duration')
    assert list(iter_harmonized_values(parameter.to_mapping()))[0].value == '30'
    assert list(iter_harmonized_values(parameter.unit.to_mapping()))[0].term_accession_number == 'UO:0000031'
    second = AEConstructor().miniml2magetab(parsed, platform_handler='generic')
    assert table(second)[0].count('Parameter Value[hz_value]') == 1


def test_generated_round_trip_keeps_hz_identifiers_and_depth():
    rows = AEConstructor().miniml2magetab(MINiMLCodec().decode(package()).package, platform_handler='generic')
    parsed = reparse(rows)
    values = list(iter_harmonized_values(parsed['sample'][0]['channel'][0]['characteristics']))
    disease = [v for v in values if v.field == 'disease']
    assert [(v.index, v.term_accession_number, v.hierarchy_depth) for v in disease] == [(0,'MONDO:0002771',0),(1,'MONDO:0003782',None)]


def test_explicit_paths_keep_missing_companions_local():
    from meta_standards_converter.magetab.semantics import render_miniml_assay_documents
    paths = [{'steps': [{'kind': 'sample', 'name': name, 'characteristics': rows}]} for name, rows in [
        ('one', [{'name':'hz_a','value':'A'}, {'name':'hz_b','value':'B'}, {'name':'hz_b_id','value':'B:1'}]),
        ('two', [{'name':'hz_a','value':'AA'}, {'name':'hz_a_id','value':'A:1'}, {'name':'hz_b','value':'BB'}]),
    ]]
    h, a, b = render_miniml_assay_documents(paths)['study.sdrf.txt']
    assert a[h.index('Characteristics[hz_a]')+1] == ''
    assert b[h.index('Characteristics[hz_a]')+1] == 'A:1'
    assert a[h.index('Characteristics[hz_b]')+1] == 'B:1'
    assert b[h.index('Characteristics[hz_b]')+1] == ''


def test_replacement_preserves_inline_source_evidence():
    from meta_standards_converter.metadata.harmonization_overrides import resolve_harmonization_overrides
    p = package()
    profile = {'schema_version':'1.0', 'replacements':{'source':['tissue']}}
    result = resolve_harmonization_overrides([p], profile, enabled=True)
    h, r = table(AEConstructor().miniml2magetab(MINiMLCodec().decode(result.packages[0]).package, platform_handler='generic'))
    assert r[h.index('Characteristics[hz_tissue]')] == 'lung'


def test_explicit_replacements_and_assay_characteristics_stay_local(tmp_path):
    p = package()
    p['series']['assay_paths'] = [{'steps':[
        {'kind':'sample', 'name':'GSM1', 'sample_ref':'GSM1', 'characteristics':[{'name':'disease','value':'IPF'}]},
        {'kind':'assay','name':'A1','characteristics':[{'name':'hz_instrument','value':'machine'}]},
    ]}]
    source=tmp_path/'input.json'; source.write_text(json.dumps(p))
    output=JSON2AEConverter().convert(str(source), enrich=False, platform_handler='generic', replacement_profile=PROFILE)[0]
    h,r=table(output)
    assert r[h.index('Characteristics[disease]')] == 'idiopathic pulmonary fibrosis'
    parsed=reparse(output)
    assert not any(v.field=='instrument' for v in iter_harmonized_values(parsed.samples[0].channels[0].to_mapping()['characteristics']))
    assay=next(s for s in parsed.series.assay_paths[0].steps if s.kind=='assay')
    assert list(iter_harmonized_values(assay.to_mapping()['characteristics']))[0].field=='instrument'


def test_generated_characteristic_units_roundtrip():
    p=package()
    p['sample'][0]['channel'][0]['characteristics'].insert(0, {
        'name': 'duration',
        'value': '30',
        'unit': {
            'value': 'min',
            'hz_unit': 'minute',
            'hz_unit_id': 'UO:1',
            'hz_unit_hierarchy_depth': 0,
        },
    })
    output=AEConstructor().miniml2magetab(MINiMLCodec().decode(p).package,platform_handler='generic')
    assert 'Comment[hz_unit]' in table(output)[0]
    parsed=reparse(output)
    duration=next(v for s in parsed.series.assay_paths[0].steps for v in getattr(s,'characteristics',[]) if v.name=='duration')
    assert list(iter_harmonized_values(duration.unit.to_mapping()))[0].hierarchy_depth==0


@pytest.mark.parametrize('atlas', [False, True])
@pytest.mark.parametrize('enabled', [False, True])
def test_all_python_outputs_preserve_evidence_and_optional_replacements(tmp_path, atlas, enabled):
    import anndata
    import pandas as pd
    from scipy import sparse
    from pathlib import Path
    p=package()
    document=p
    if atlas:
        document=json.loads((Path(__file__).parent/'fixtures/edge_cases/atlas-groups/inputs/atlas.json').read_text())
        document['datasets']=document['datasets'][:1]
        document['datasets'][0].update(dataset_id='GSE1',metadata=p)
        document['summary'].update(dataset_count=1,failed_dataset_count=0)
    source=tmp_path/'input.json';source.write_text(json.dumps(document))
    expression=tmp_path/'input.h5ad'
    anndata.AnnData(sparse.csr_matrix([[1]]),obs=pd.DataFrame(index=['cell1']),var=pd.DataFrame(index=['gene1'])).write_h5ad(expression)
    options={'replacement_profile':PROFILE} if enabled else {}
    expected='idiopathic pulmonary fibrosis' if enabled else 'IPF'
    rows=JSON2AEConverter().convert(str(source),enrich=False,platform_handler='generic',**options)
    h,r=table(rows[0]); assert r[h.index('Characteristics[disease]')]==expected
    assert r[h.index('Characteristics[hz_disease]')]=='idiopathic pulmonary fibrosis'
    for format in ('tsv','csv'):
        out=tmp_path/format
        JSON2TSVConverter().export_manifest(source,outdir=out,output_format=format,**options)
        row=next(csv.DictReader(next(out.glob('*.'+format)).open(),delimiter='\t' if format=='tsv' else ','))
        assert row['msc.sample.channel.disease'].split('; ')[0]==expected
        assert 'msc.characteristics.hz_disease' in row
    for kind in ('h5ad','obs'):
        out=tmp_path/kind
        if kind=='h5ad':
            JSON2H5ADConverter().convert_source(str(source),out=str(out),asset_specs=[f'GSM1={expression}'],**options)
            result=anndata.read_h5ad(next(out.rglob('GSM1.h5ad')))
        else:
            result=JSON2OBSConverter().convert(source,outdir=out,include_uns=True,asset_specs=[f'GSM1={expression}'],**options)
        assert result.obs['msc.sample.channel.disease'].iloc[0].split('; ')[0]==expected
        assert 'msc.characteristics.hz_disease' in result.obs
        assert not any('harmonized_' in k for k in result.obs)
        assert 'msc_miniml' in result.uns
    assert json.loads(source.read_text())==document


def test_checkpoint_identity_includes_profile_and_version(tmp_path):
    from meta_standards_converter.expression.assets import Asset
    args=dict(sample_id='GSM1',source_json_sha256='source',sample=package()['sample'][0],asset=Asset('GSM1','test.h5ad','h5ad'),orientation='auto')
    old=ProcessedCheckpointStore(lambda:'7.0.0').key(tmp_path,**args)
    new=ProcessedCheckpointStore(lambda:'8.0.0')
    raw=new.key(tmp_path,**args)
    profiled=new.key(tmp_path,**args,replacement_profile=PROFILE)
    reordered=new.key(tmp_path,**args,replacement_profile=dict(reversed(list(PROFILE.items()))))
    assert old[0]!=raw[0]
    assert raw[2]!=profiled[2]
    assert reordered==profiled


@pytest.mark.parametrize('command', ['json2ae','json2tsv','json2h5ad','json2obs'])
@pytest.mark.parametrize('bad', ['{','[]','missing-file','both'])
def test_cli_policy_failure_precedes_output_creation(command,bad,tmp_path):
    cli=importlib.import_module('meta_standards_converter.cli.'+command)
    out=tmp_path/'not-created'
    args=['input.json', '--outdir' if command=='json2obs' else '--out',str(out)]
    args+=['--replacement-profile-file',str(tmp_path/'missing')] if bad=='missing-file' else ['--replacement-profile',bad]
    if bad=='both':args+=['--replacement-profile-file','second.json']
    with pytest.raises(SystemExit) as error:cli.main(args)
    assert error.value.code==2
    assert not out.exists()


def test_indexed_parameter_and_unit_groups_keep_original_indexes():
    p=package();p['series']['protocols']=[{'name':'P'}]
    p['series']['assay_paths']=[{'steps':[{'kind':'sample','name':'GSM1'},
        {'kind':'protocol_application','protocol_ref':'P','parameter_values':[{'name':'time','value':'5','hz_value(3)':'five','hz_value_hierarchy_depth(3)':0,'unit':{'value':'min','hz_unit(2)':'minute','hz_unit_hierarchy_depth(2)':0}}]},
        {'kind':'assay','name':'A'}]}]
    rendered=AEConstructor().miniml2magetab(MINiMLCodec().decode(p).package,platform_handler='generic')
    parsed=reparse(rendered)
    v=next(s for s in parsed.series.assay_paths[0].steps if s.kind=='protocol_application').parameter_values[0]
    assert [(a.index,a.hierarchy_depth) for a in iter_harmonized_values(v.to_mapping())]==[(3,0)]
    assert [(a.index,a.hierarchy_depth) for a in iter_harmonized_values(v.unit.to_mapping())]==[(2,0)]


def test_generated_channels_reparse_separately_with_repeated_source_names():
    p=package();p['sample'][0]['channel']=[copy.deepcopy(p['sample'][0]['channel'][0]) for _ in range(2)]
    second=p['sample'][0]['channel'][1]
    second['characteristics'][1]['value']='healthy'
    second['characteristics'][2]['value']='MONDO:healthy'
    result=AEConstructor().miniml2magetab(MINiMLCodec().decode(p).package,platform_handler='generic')
    parsed=reparse(result)
    assert len(parsed.samples[0].channels)==2
    values=[list(iter_harmonized_values(c.to_mapping()['characteristics']))[0].value for c in parsed.samples[0].channels]
    assert values==['idiopathic pulmonary fibrosis','healthy']
    second_render=AEConstructor().miniml2magetab(parsed,platform_handler='generic')
    h,*rows=table(second_render)
    assert [r[h.index('Characteristics[hz_disease]')] for r in rows]==values


def test_explicit_multichannel_binding_uses_labels_and_reports_ambiguity(caplog):
    p=package();first=p['sample'][0]['channel'][0];first['label']={'value':'Cy3'}
    second=copy.deepcopy(first);second['label']={'value':'Cy5'};second['characteristics'][1]['value']='healthy'
    p['sample'][0]['channel'].append(second)
    p['series']['assay_paths'] = [{'steps': [{'kind': 'sample', 'name': 'GSM1', 'sample_ref': 'GSM1'}, {'kind': 'labeled_extract', 'name': label, 'label': {'value': label}}, {'kind': 'assay', 'name': 'A' + label}]} for label in ['Cy5', 'Cy3']]
    result=AEConstructor().miniml2magetab(MINiMLCodec().decode(p).package,platform_handler='generic')
    h,*rows=table(result)
    assert [r[h.index('Characteristics[hz_disease]')] for r in rows]==['healthy','idiopathic pulmonary fibrosis']
    for path in p['series']['assay_paths']:path['steps'][1].pop('label')
    ambiguous=AEConstructor().miniml2magetab(MINiMLCodec().decode(p).package,platform_handler='generic')
    assert 'Characteristics[hz_disease]' not in table(ambiguous)[0]
    assert 'Ambiguous channel association' in caplog.text


def test_explicit_profile_updates_organism_and_source_standard_occurrences(tmp_path):
    p=package();p['sample'][0]['channel'][0]['organism']=[{'value':'human','hz_species':'Homo sapiens','hz_species_id':'NCBITaxon:9606','hz_species_onto':'ncbitaxon'}]
    p['series']['assay_paths'] = [
        {
            'steps': [
                {
                    'kind': 'sample',
                    'name': 'GSM1',
                    'sample_ref': 'GSM1',
                    'characteristics': [
                        {
                            'name': 'organism',
                            'value': 'human',
                        },
                    ],
                    'comments': [
                        {
                            'name': 'Sample_source_name',
                            'value': 'raw source',
                        },
                    ],
                },
            ],
        },
    ]
    profile={'schema_version':'1.0','replacements':{'organism':['species'],'source':['tissue']}}
    source=tmp_path/'input.json';source.write_text(json.dumps(p))
    output=JSON2AEConverter().convert(str(source),enrich=False,platform_handler='generic',replacement_profile=profile)[0]
    h,r=table(output)
    assert r[h.index('Characteristics[organism]')]=='Homo sapiens'
    assert r[h.index('Comment[Sample_source_name]')]=='lung'
    assert 'Characteristics[hz_species]' in h
    assert 'Characteristics[hz_tissue]' in h


@pytest.mark.parametrize('command', ['json2ae','json2tsv','json2h5ad','json2obs'])
@pytest.mark.parametrize('atlas', [False,True])
def test_cli_direct_profile_file_reaches_written_artifact(command, atlas, tmp_path):
    from pathlib import Path
    import anndata
    import pandas as pd
    from scipy import sparse
    p=package();document=p
    if atlas:
        document=json.loads((Path(__file__).parent/'fixtures/edge_cases/atlas-groups/inputs/atlas.json').read_text())
        document['datasets']=document['datasets'][:1];document['datasets'][0].update(dataset_id='GSE1',metadata=p)
        document['summary'].update(dataset_count=1,failed_dataset_count=0)
    source=tmp_path/'input.json';source.write_text(json.dumps(document))
    profile=tmp_path/'profile.json';profile.write_text(json.dumps(PROFILE))
    out=tmp_path/'out'
    args=[str(source),'--outdir' if command=='json2obs' else '--out',str(out),'--replacement-profile-file',str(profile)]
    if command=='json2ae':args+=['--no-enrich','--platform-handler','generic']
    if command in {'json2h5ad','json2obs'}:
        expression=tmp_path/'counts.h5ad'
        anndata.AnnData(sparse.csr_matrix([[1]]),obs=pd.DataFrame(index=['cell']),var=pd.DataFrame(index=['gene'])).write_h5ad(expression)
        args+=['--asset',f'GSM1={expression}']
    assert importlib.import_module('meta_standards_converter.cli.'+command).main(args) in (0,None)
    if command=='json2ae':
        h,r=list(csv.reader(next(out.glob('*.sdrf.txt')).open(),delimiter='\t'))
        assert r[h.index('Characteristics[disease]')]=='idiopathic pulmonary fibrosis'
        assert 'Characteristics[hz_disease]' in h
    else:
        if command=='json2h5ad':frame=anndata.read_h5ad(next(out.rglob('GSM1.h5ad'))).obs
        else:
            path=next(out.glob('*.tsv' if command=='json2tsv' else '*.obs.csv'))
            frame=pd.read_csv(path,sep='\t' if command=='json2tsv' else ',')
        assert frame['msc.sample.channel.disease'].iloc[0].split('; ')[0]=='idiopathic pulmonary fibrosis'
        assert 'msc.characteristics.hz_disease' in frame


def test_atlas_profile_applies_to_every_selected_group(tmp_path):
    from pathlib import Path
    atlas=json.loads((Path(__file__).parent/'fixtures/edge_cases/atlas-groups/inputs/atlas.json').read_text())
    atlas['datasets']=[copy.deepcopy(atlas['datasets'][0]) for _ in range(2)]
    for index,dataset in enumerate(atlas['datasets'],1):
        p=package();p['series']['iid']=f'GSE{index}';p['series']['sample_ref']=[{'ref':f'GSM{index}'}];p['sample'][0]['iid']=f'GSM{index}'
        dataset.update(dataset_id=f'GSE{index}',metadata=p,source_ordinal=index-1)
    atlas['summary'].update(completed_dataset_count=2,failed_dataset_count=0)
    source=tmp_path/'atlas.json';source.write_text(json.dumps(atlas))
    outputs=JSON2AEConverter().convert(str(source),enrich=False,platform_handler='generic',replacement_profile=PROFILE)
    assert len(outputs)==2
    for output in outputs:
        h,r=table(output)
        assert r[h.index('Characteristics[disease]')]=='idiopathic pulmonary fibrosis'


def test_hz_ontology_sources_are_declared_in_idf():
    rows=AEConstructor().miniml2magetab(MINiMLCodec().decode(package()).package,platform_handler='generic')
    names=next(r[1:] for r in rows if r[0]=='Term Source Name')
    assert 'MONDO' in names


def test_reparsed_channel_identity_binds_new_evidence_before_rendering():
    p=package();p['sample'][0]['channel']=[copy.deepcopy(p['sample'][0]['channel'][0]) for _ in range(2)]
    first=AEConstructor().miniml2magetab(MINiMLCodec().decode(p).package,platform_handler='generic')
    parsed=reparse(first).to_mapping()
    parsed['sample'][0]['channel'][1]['characteristics'].extend([
        {
            'name': 'hz_new',
            'value': 'second only',
        },
    ])
    rendered=AEConstructor().miniml2magetab(MINiMLCodec().decode(parsed).package,platform_handler='generic')
    h,*rows=table(rendered)
    assert [r[h.index('Characteristics[hz_new]')] for r in rows]==['','second only']
