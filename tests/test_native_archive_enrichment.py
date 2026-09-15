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
from meta_standards_converter.miniml import MINiMLCodec
from meta_standards_converter.miniml.sra_parser import SRAParser
from meta_standards_converter.metadata.archive_enrichment import merge_archive_metadata, LinkedArchiveEnricher
from tests.test_native_archive_parsers import fixture_records


def native():
    return SRAParser().parse(fixture_records('sra'))


def linked(package, accession, title):
    data = package.to_mapping()
    data['source'] = {'format': 'MAGE-TAB' if accession.startswith('E-') else 'GEO'}
    data['series']['iid'] = accession
    data['series']['accession'].append({'value': accession, 'database': 'ArrayExpress' if accession.startswith('E-') else 'GEO'})
    data['series']['title'] = title
    data['sample'][0]['title'] = title
    old_iid = data['sample'][0]['iid']
    data['sample'][0]['iid'] = 'GSM1'
    for path in data['series'].get('assay_paths', []):
        for step in path['steps']:
            if step.get('sample_ref') == old_iid: step['sample_ref'] = 'GSM1'
    data['sample'][0]['channel'][0]['characteristics'] = [{'name': 'host', 'value': title, 'term_accession_number': 'ONT:1', 'term_source_ref': 'ONT'}]
    # This fixture clones generated native paths: update their same host projection too.
    for path in data['series'].get('assay_paths', []):
        for step in path['steps']:
            if step['kind'] in ('source', 'sample'):
                step['characteristics'] = [c for c in step.get('characteristics', []) if c['name'] != 'host'] + deepcopy(data['sample'][0]['channel'][0]['characteristics'])
    return MINiMLCodec().decode(data).package


def test_coherent_priority_identity_and_membership():
    package = native()
    other = linked(package, 'GSE1', 'enriched')
    merged, issues = merge_archive_metadata(package, other, prefer=True)
    data = merged.to_mapping()
    assert data['series']['title'] == 'enriched'
    assert data['series']['iid'] == package.series.iid
    assert data['sample'][0]['iid'] == package.samples[0].iid
    assert data['sample'][0]['channel'][0]['characteristics'][0]['value'] != 'enriched' or any(a['name'] == 'host' for a in data['sample'][0]['channel'][0]['characteristics'])
    host = next(a for a in data['sample'][0]['channel'][0]['characteristics'] if a['name'] == 'host')
    assert host['value'] == 'enriched'
    assert any(s.get('characteristics') and any(a.get('value') == 'enriched' for a in s['characteristics']) for p in data['series']['assay_paths'] for s in p['steps'])
    assert data['extensions']['insdc']['records'][:1] == package.to_mapping()['extensions']['insdc']['records'][:1]
    assert not issues


def test_missing_values_do_not_displace_informative_and_titles_do_not_join():
    package = native()
    other = linked(package, 'GSE1', 'not provided').to_mapping()
    other['sample'][0]['accession'] = [{'value': 'GSM9', 'database': 'GEO'}]
    other['sample'][0]['sra_accession'] = []
    other['sample'][0]['relation'] = []
    other['sample'][0]['sra_run'] = []
    merged, issues = merge_archive_metadata(package, MINiMLCodec().decode(other).package, prefer=True)
    assert merged.series.title == package.series.title
    assert merged.samples[0].title == package.samples[0].title
    assert len(merged.samples) == 1
    assert any('unmatched' in i for i in issues)


def test_ae_over_geo_and_unavailable_enrichment():
    package = native().to_mapping()
    package['series']['relation'] = [{'type': 'GEO', 'target': 'GSE1'}, {'type': 'ArrayExpress', 'target': 'E-MTAB-1'}]
    package = MINiMLCodec().decode(package).package
    class Converter:
        def convert(self, accession, **kwargs):
            return [linked(package, accession, 'AE' if accession.startswith('E-') else 'GEO')]
    enriched, issues = LinkedArchiveEnricher(geo_converter=Converter(), ae_converter=Converter()).enrich(package)
    assert enriched.series.title == 'AE'
    class Broken:
        def convert(self, *args, **kwargs):
            raise OSError('private detail')
    enriched, issues = LinkedArchiveEnricher(geo_converter=Broken(), ae_converter=Broken()).enrich(package)
    assert enriched.series.iid == package.series.iid
    assert issues and 'private detail' not in str(issues)


def test_ambiguous_sample_join_and_unrelated_peer_are_rejected():
    package = native()
    other = linked(package, 'GSE1', 'enriched').to_mapping()
    other['sample'].append(deepcopy(other['sample'][0]))
    other['sample'][1]['iid'] = 'GSM2'
    merged, issues = merge_archive_metadata(package, MINiMLCodec().decode(other).package, prefer=True)
    assert merged.samples[0].title == package.samples[0].title
    assert any('ambiguous' in i for i in issues)
    other['series']['accession'] = [{'value': 'SRP999', 'database': 'SRA'}]
    other['series']['iid'] = 'SRP999'
    merged, issues = merge_archive_metadata(package, MINiMLCodec().decode(other).package)
    assert merged == package
    assert issues


def test_enrichment_files_and_sample_protocols_reach_native_assay_paths():
    package = native()
    other = linked(package, 'GSE1', 'GEO').to_mapping()
    other['sample'][0]['supplementary_data'] = [{'value': 'https://example.org/counts.tsv', 'type': 'tsv'}]
    other['sample'][0]['channel'][0]['extract_protocol'] = 'Explicit extraction method'
    merged, issues = merge_archive_metadata(package, MINiMLCodec().decode(other).package, prefer=True)
    data = merged.to_mapping()
    assert any(s.get('link', {}).get('value') == 'https://example.org/counts.tsv' for p in data['series']['assay_paths'] for s in p['steps'])
    assert any(p.get('description') == 'Explicit extraction method' for p in data['series']['protocols'])
    assert any(s.get('protocol_ref', '').endswith('extract_protocol') for p in data['series']['assay_paths'] for s in p['steps'])


def test_peer_only_runs_are_added_on_exact_study_and_sample_match():
    package = native()
    other = package.to_mapping()
    other['sample'][0]['sra_run'][0]['run'] = 'SRR99'
    for path in other['series']['assay_paths']:
        for step in path['steps']:
            if step['kind'] == 'scan':
                step['name'] = 'SRR99'
                for comment in step.get('comments', []):
                    if comment['name'] in ('ENA_RUN', 'SRA_RUN'): comment['value'] = 'SRR99'
    merged, issues = merge_archive_metadata(package, MINiMLCodec().decode(other).package)
    assert {r['run'] for r in merged.to_mapping()['sample'][0]['sra_run']} == {'SRR11192680', 'SRR99'}
    assert any(s.get('name') == 'SRR99' for p in merged.to_mapping()['series']['assay_paths'] for s in p['steps'])


def test_processing_protocol_priority_is_coherent_and_scoped():
    package = native()
    geo = linked(package, 'GSE1', 'GEO').to_mapping()
    geo['sample'][0]['data_processing'] = 'GEO processing'
    merged, _ = merge_archive_metadata(package, MINiMLCodec().decode(geo).package, prefer=True)
    ae = linked(package, 'E-MTAB-1', 'AE').to_mapping()
    ae['sample'][0]['data_processing'] = 'AE processing'
    merged, _ = merge_archive_metadata(merged, MINiMLCodec().decode(ae).package, prefer=True)
    data = merged.to_mapping()
    names = {p['name']: p for p in data['series']['protocols']}
    for path in data['series']['assay_paths']:
        applied = [names[s['protocol_ref']].get('description') for s in path['steps'] if s.get('protocol_ref') in names]
        assert 'AE processing' not in applied
        assert 'GEO processing' not in applied
    assert any(p.get('description') == 'AE processing' for p in names.values())


def test_experiment_protocol_enrichment_does_not_leak_to_other_experiments():
    package = native().to_mapping()
    original = deepcopy(package['series']['assay_paths'][0])
    for step in original['steps']:
        if step['kind'] == 'assay': step['name'] = 'SRX99'
        if step['kind'] == 'scan': step['name'] = 'SRR99'; step['comments'] = []
    package['series']['assay_paths'].append(original)
    package = MINiMLCodec().decode(package).package
    extra = linked(package, 'E-MTAB-1', 'AE').to_mapping()
    path = deepcopy(extra['series']['assay_paths'][0])
    for step in path['steps']:
        if step.get('sample_ref'): step['sample_ref'] = 'GSM1'
    path['steps'].insert(2, {'kind': 'protocol_application', 'protocol_ref': 'specific'})
    extra['series']['assay_paths'] = [path]
    extra['series']['protocols'] = [{'name': 'specific', 'description': 'Specific experiment preparation', 'type': {'value': 'library construction protocol'}}]
    merged, issues = merge_archive_metadata(package, MINiMLCodec().decode(extra).package, prefer=True)
    for path in merged.to_mapping()['series']['assay_paths']:
        specific = any(s.get('protocol_ref', '').endswith(':specific') for s in path['steps'])
        if any(s.get('name') == 'SRX99' for s in path['steps']): assert not specific
        else: assert specific


def test_geo_link_discovers_verified_egeod_and_applies_ae_last():
    package=native().to_mapping()
    package['series']['relation']=[{'type':'GEO','target':'GSE1'}]
    package=MINiMLCodec().decode(package).package
    calls=[]
    class Converter:
        def convert(self, accession, **kwargs):
            calls.append(accession)
            return [linked(package,accession,'AE' if accession.startswith('E-') else 'GEO')]
    result,issues=LinkedArchiveEnricher(geo_converter=Converter(),ae_converter=Converter()).enrich(package)
    assert not issues
    assert calls==['GSE1','E-GEOD-1']
    assert result.series.title=='AE'
    assert result.series.iid==package.series.iid


def test_unrelated_egeod_candidate_cannot_redefine_native_study():
    package=native().to_mapping();package['series']['relation']=[{'type':'GEO','target':'GSE1'}]
    package=MINiMLCodec().decode(package).package
    class Geo:
        def convert(self, accession, **kwargs):return [linked(package,accession,'GEO')]
    class AE:
        def convert(self, accession, **kwargs):
            other=linked(package,accession,'wrong').to_mapping()
            other['series']['accession']=[{'value':accession,'database':'ArrayExpress'},{'value':'SRP999999','database':'SRA'}]
            other['series']['relation']=[]
            return [MINiMLCodec().decode(other).package]
    result,issues=LinkedArchiveEnricher(geo_converter=Geo(),ae_converter=AE()).enrich(package)
    assert result.series.title=='GEO'
    assert any('identity' in i for i in issues)


def test_egeod_probe_cannot_join_sibling_read_study_by_project_alone():
    package=native().to_mapping();package['series']['relation']=[{'type':'GEO','target':'GSE1'}]
    package=MINiMLCodec().decode(package).package
    class Geo:
        def convert(self,*a,**k):return []
    class AE:
        def convert(self,accession):
            other=linked(package,accession,'wrong sibling').to_mapping()
            project=[a for a in package.to_mapping()['series']['accession'] if a['database']=='BioProject']
            other['series']['accession']=[*project,{'value':accession,'database':'ArrayExpress'},{'value':'SRP999999','database':'SRA'}]
            other['series']['relation']=[]
            return [MINiMLCodec().decode(other).package]
    result,issues=LinkedArchiveEnricher(geo_converter=Geo(),ae_converter=AE()).enrich(package)
    assert result.series.title==package.series.title
    assert any('identity' in i for i in issues)


def test_shared_organization_ids_do_not_create_sample_identity_matches():
    from meta_standards_converter.metadata.archive_enrichment import entity_ids
    data=native().to_mapping();sample=data['sample'][0]
    sample.setdefault('relation',[]).extend([
        {'type':'archive center','target':'ena:BioSample:SAMN999:organization-2'},
        {'type':'derived from','target':'SAMN998'}])
    assert not {'SAMN999','SAMN998'} & entity_ids(sample,sample=True)


def test_enrichment_preserves_native_and_explicit_linked_contact_scope():
    data=native().to_mapping();sid=data['sample'][0]['iid']
    data['contributor']=[{'iid':'native-person','person':{'last':'Native'}}]
    data['sample'][0]['contact_ref']=[{'ref':'native-person'}]
    data['series']['contributor_ref']=[{'ref':'native-person'}]
    package=MINiMLCodec().decode(data).package
    extra=linked(package,'GSE1','linked').to_mapping()
    extra['contributor']=[{'iid':'linked-person','person':{'last':'Linked'}}]
    extra['sample'][0]['contact_ref']=[{'ref':'linked-person'}]
    extra['series']['contact_ref']=[{'ref':'linked-person'}]
    extra['series']['contributor_ref']=[{'ref':'linked-person'}]
    result,issues=merge_archive_metadata(package,MINiMLCodec().decode(extra).package,prefer=True)
    assert not issues
    mapped=result.to_mapping()
    assert {r['ref'] for r in mapped['sample'][0]['contact_ref']}=={'native-person','GSE1:linked-person'}
    assert {r['ref'] for r in mapped['series']['contributor_ref']}=={'native-person','GSE1:linked-person'}
    assert mapped['series']['contact_ref']==[{'ref':'GSE1:linked-person'}]
    assert mapped['sample'][0]['iid']==sid


def test_sample_only_geo_contact_is_not_promoted_to_study_and_repeated_join_is_idempotent():
    package=native();extra=linked(package,'GSE1','linked').to_mapping()
    extra['contributor']=[{'iid':'sample-person','person':{'last':'Sample only'}}]
    extra['series'].pop('contributor_ref',None);extra['series'].pop('contact_ref',None)
    extra['sample'][0]['contact_ref']=[{'ref':'sample-person'}]
    incoming=MINiMLCodec().decode(extra).package
    for _ in range(2):
        package,issues=merge_archive_metadata(package,incoming,prefer=True);assert not issues
        package=MINiMLCodec().decode(package.to_mapping()).package
    d=package.to_mapping()
    assert [r['ref'] for r in d['sample'][0]['contact_ref']].count('GSE1:sample-person')==1
    assert not any(r['ref']=='GSE1:sample-person' for r in d['series'].get('contributor_ref',[]))
    assert not d['series'].get('contact_ref')


def test_ambiguous_sample_contacts_remain_unassigned_to_native_sample():
    package=native();extra=linked(package,'GSE1','linked').to_mapping()
    extra['contributor']=[{'iid':'linked-person','person':{'last':'Linked'}}]
    extra['sample'][0]['contact_ref']=[{'ref':'linked-person'}]
    extra['sample'].append(deepcopy(extra['sample'][0]));extra['sample'][1]['iid']='GSM2'
    result,issues=merge_archive_metadata(package,MINiMLCodec().decode(extra).package,prefer=True)
    assert any('ambiguous' in x for x in issues)
    d=result.to_mapping()
    assert not any(r['ref']=='GSE1:linked-person' for r in d['sample'][0]['contact_ref'])
    assert any('GSE1:linked-person' in str(r['metadata']) for r in d['extensions']['insdc']['records'])


def test_plain_enrichment_completes_unique_native_characteristic_groups_and_paths():
    from meta_standards_converter.miniml import MINiMLCodec
    original=native().to_mapping();other=linked(native(),'GSE1','linked').to_mapping()
    item={'name':'tissue','value':'brain','term_source_ref':'UBERON','term_accession_number':'UBERON:0000955',
          'unit':{'value':'mg','term_source_ref':'UO','term_accession_number':'UO:0000022'}}
    original['sample'][0]['channel'][0]['characteristics']=[item]
    other['sample'][0]['channel'][0]['characteristics']=[{'name':'tissue','value':'brain'}]
    package,issues=merge_archive_metadata(MINiMLCodec().decode(original).package,MINiMLCodec().decode(other).package,prefer=True,linked_accession='GSE1')
    data=package.to_mapping()
    assert data['sample'][0]['channel'][0]['characteristics']==[item]
    for p in data['series']['assay_paths']:
        assert item in p['steps'][0]['characteristics']


@pytest.mark.parametrize('old,new', [
    ([{'name':'x','value':'brain','term_accession_number':'U:1'}],[{'name':'x','value':'heart'}]),
    ([{'name':'x','value':'brain','term_accession_number':'U:1'}]*2,[{'name':'x','value':'brain'}]),
    ([{'name':'x','value':'brain','term_accession_number':'U:1'}],[{'name':'x','value':'brain'}]*2),
    ([{'name':'x','value':'5','unit':{'value':'mg'},'term_accession_number':'U:1'}],[{'name':'x','value':'5','unit':{'value':'g'}}]),
    ([{'name':'x','value':'brain','term_source_ref':'U','term_accession_number':'U:1'}],[{'name':'x','value':'brain','term_accession_number':'V:2'}]),
])
def test_characteristic_completion_never_bridges_conflicting_or_ambiguous_groups(old,new):
    from meta_standards_converter.metadata.archive_enrichment import _merge_entity
    target={'characteristics':deepcopy(old)}
    _merge_entity(target,{'characteristics':deepcopy(new)},True)
    assert target['characteristics']==new


def test_characteristic_completion_preserves_compatible_unit_ontology():
    from meta_standards_converter.metadata.archive_enrichment import _merge_entity
    target={'characteristics':[{'name':'dose','value':'2','unit':{'value':'mg','term_source_ref':'UO','term_accession_number':'UO:0000022'}}]}
    _merge_entity(target,{'characteristics':[{'name':'dose','value':'2','unit':{'value':'mg'}}]},True)
    assert target['characteristics'][0]['unit']['term_accession_number']=='UO:0000022'


def test_disjoint_partial_ontology_groups_do_not_create_namespace_pairs():
    from meta_standards_converter.metadata.archive_enrichment import _merge_entity
    target={'characteristics':[{'name':'x','value':'same','term_source_ref':'U'}]}
    incoming={'characteristics':[{'name':'x','value':'same','term_accession_number':'V:2'}]}
    _merge_entity(target,incoming,True)
    assert target==incoming


@pytest.mark.parametrize('where', ['native_path','native_path_conflict','preferred_path','preferred_path_conflict','preferred_channel'])
def test_workflow_only_annotations_survive_at_their_explicit_biological_scope(where):
    original=native().to_mapping();other=linked(native(),'E-MTAB-1','linked').to_mapping()
    plain={'name':'tissue','value':'brain'};term={**plain,'term_source_ref':'UBERON','term_accession_number':'UBERON:0000955'}
    original['sample'][0]['channel'][0]['characteristics']=[deepcopy(plain)]
    other['sample'][0]['channel'][0]['characteristics']=[deepcopy(plain)]
    for d in (original,other):
        for path in d['series']['assay_paths']:
            path['steps'][0]['characteristics']=[deepcopy(plain)]
    if where in ('native_path','native_path_conflict'):
        if where=='native_path_conflict':original['sample'][0]['channel'][0]['characteristics']=[{**plain,'term_source_ref':'OTHER','term_accession_number':'OTHER:2'}]
        for path in original['series']['assay_paths']:path['steps'][0]['characteristics']=[deepcopy(term)]
    elif where=='preferred_channel':
        other['sample'][0]['channel'][0]['characteristics']=[deepcopy(term)]
        for path in original['series']['assay_paths']:path['steps'][0]['characteristics']=[{**plain,'term_source_ref':'OTHER','term_accession_number':'OTHER:2'}]
    else:
        for path in other['series']['assay_paths']:path['steps'][0]['characteristics']=[deepcopy(term)]
        if where=='preferred_path_conflict':original['sample'][0]['channel'][0]['characteristics']=[{**plain,'term_source_ref':'OTHER','term_accession_number':'OTHER:2'}]
    package,issues=merge_archive_metadata(MINiMLCodec().decode(original).package,MINiMLCodec().decode(other).package,prefer=True)
    assert not issues
    for path in package.to_mapping()['series']['assay_paths']:
        value=next(v for v in path['steps'][0]['characteristics'] if v['name']=='tissue')
        assert value==term


def test_bound_biological_fallback_preserves_names_but_not_ambiguous_or_extract_annotations():
    from meta_standards_converter.metadata.archive_workflows import _bind_native
    plain={'name':'tissue','value':'brain'}
    term={**plain,'term_source_ref':'UBERON','term_accession_number':'UBERON:0000955'}
    old={'kind':'source','name':'native','characteristics':[term,{'name':'sex','value':'female'}]}
    for incoming in ([],[plain]):
        node={'kind':'source','name':'linked','characteristics':deepcopy(incoming)}
        result=_bind_native([node],'sample',[old])
        assert sorted(result[0]['characteristics'],key=lambda v:v['name'])==[{'name':'sex','value':'female'},term]
    for native_nodes in ([old,deepcopy(old)],[{**old,'kind':'extract'}]):
        node={'kind':'source','name':'linked','characteristics':[deepcopy(plain)]}
        assert _bind_native([node],'sample',native_nodes)[0]['characteristics']==[plain]


def test_missing_biological_workflow_values_do_not_displace_informative_fallback():
    from meta_standards_converter.metadata.archive_workflows import _biological_characteristics
    fallback=[{'name':'tissue','value':'brain','term_source_ref':'UBERON','term_accession_number':'UBERON:0000955'}]
    assert _biological_characteristics([{'name':'tissue','value':'not provided'}],fallback)==fallback


def test_explicit_preferred_workflow_annotation_takes_priority_over_channel_copy():
    original=native();other=linked(original,'E-MTAB-1','linked').to_mapping()
    term={'name':'tissue','value':'brain','term_source_ref':'UBERON','term_accession_number':'UBERON:0000955'}
    other['sample'][0]['channel'][0]['characteristics']=[{**term,'term_source_ref':'OTHER','term_accession_number':'OTHER:2'}]
    for path in other['series']['assay_paths']:
        path['steps'][0]['characteristics']=[deepcopy(term)]
    package,issues=merge_archive_metadata(original,MINiMLCodec().decode(other).package,prefer=True)
    assert not issues
    for path in package.to_mapping()['series']['assay_paths']:
        assert term in path['steps'][0]['characteristics']


@pytest.mark.parametrize('separate', [False,True])
def test_accepted_taxonomy_updates_generated_nodes_and_preserves_separate_occurrences(separate):
    from tests.test_protocol_export import render
    original=native().to_mapping()
    old=deepcopy(original['sample'][0]['channel'][0]['organism'][0])
    literal={'name':'organism','value':old['value']}
    if separate:
        original['sample'][0]['channel'][0]['characteristics'].append(literal)
        for path in original['series']['assay_paths']:
            path['steps'][0]['characteristics'].insert(-1,deepcopy(literal))
    package=MINiMLCodec().decode(original).package
    other=linked(package,'GSE1','linked').to_mapping();other['series']['assay_paths']=[]
    other['sample'][0]['channel'][0]['organism']=[{'value':'Mus musculus','taxid':'10090'}]
    merged,issues=merge_archive_metadata(package,MINiMLCodec().decode(other).package,prefer=True)
    assert not issues
    mapped=merged.to_mapping()
    for path in mapped['series']['assay_paths']:
        groups=[v for v in path['steps'][0]['characteristics'] if v['name']=='organism']
        assert {'name':'organism','value':'Mus musculus','term_source_ref':'NCBITaxon','term_accession_number':'10090'} in groups
        assert not any(v.get('term_accession_number')==old['taxid'] for v in groups)
        assert (literal in groups)==separate
    table=next(row[1] for row in render(mapped) if row[0]=='SDRF File')
    for row in table[1:]:
        values=[row[i] for i,h in enumerate(table[0]) if h=='Characteristics[organism]']
        assert any('Mus musculus' in v for v in values)
    assert package.to_mapping()==original


def test_authored_workflow_taxonomy_is_not_replaced_by_generated_sample_projection():
    package=native();other=linked(package,'E-MTAB-1','linked').to_mapping()
    other['sample'][0]['channel'][0]['organism']=[{'value':'Mus musculus','taxid':'10090'}]
    supplied={'name':'organism','value':'Danio rerio','term_source_ref':'NCBITaxon','term_accession_number':'7955'}
    for path in other['series']['assay_paths']:
        path['steps'][0]['characteristics']=[deepcopy(supplied)]
    merged,issues=merge_archive_metadata(package,MINiMLCodec().decode(other).package,prefer=True)
    assert not issues
    for path in merged.to_mapping()['series']['assay_paths']:
        assert [v for v in path['steps'][0]['characteristics'] if v['name']=='organism']==[supplied]
