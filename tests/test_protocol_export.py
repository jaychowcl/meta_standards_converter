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
from meta_standards_converter.magetab.constructor import AEConstructor
from tests.test_archive_export_cleanup import legacy
from tests.test_archive_publications import PubMed, NoRuns


def render(data):
    return AEConstructor(pubmed_client=PubMed(), insdc_client=NoRuns()).miniml2magetab(MINiMLCodec().decode(data).package)


def definition(name, description="20 ug RNA’s preparation", **kwargs):
    return dict(name=name, description=description, type={"value": "library construction protocol"}, hardware=["machine"], **kwargs)


def test_protocol_export_keeps_unused_deduplicates_and_binds_references():
    data=legacy()
    data['series']['protocols']=[definition('ERX1:library'),definition('ERX2:library', "20 µg RNA's preparation"),
        definition('linked:P-MTAB-9', 'other'),definition('linked:P-MTAB-10','other'),
        definition('P-'+data['series']['iid']+'-1','reserved')]
    data['series']['assay_paths'][0]['steps'].insert(1, {'kind':'protocol_application','protocol_ref':'ERX2:library'})
    before=deepcopy(data);rows=render(data);values={r[0]:r[1:] for r in rows}
    names=values['Protocol Name']
    assert len(names)==4 and len(set(names))==4
    assert {'P-MTAB-9','P-MTAB-10'} <= set(names)
    generated='P-'+data['series']['iid']+'-2'
    assert generated in names
    assert values['Protocol Description'][names.index(generated)]=="20 ug RNA’s preparation"
    table=values['SDRF File'][0]
    assert generated in table[1] and 'ERX2:library' not in str(table)
    indices=[i for i,r in enumerate(rows) if r[0].startswith('Protocol ')]
    assert indices==list(range(min(indices),max(indices)+1))
    assert all(len(rows[i])==len(names)+1 for i in indices)
    assert data==before


@pytest.mark.parametrize('field,value',[('hardware',['other']),('software',['other']),('parameters',['temperature']),('contacts',['person']),('performers',['person']),('type',{'value':'different'})])
def test_conflicting_protocol_metadata_is_not_merged(field,value):
    from meta_standards_converter.magetab.protocol_export import prepare_protocols
    data=legacy();a=definition('first');b=definition('second');a[field]=['native'] if field!='type' else {'value':'native'};b[field]=value
    data['series']['protocols']=[a,b];prepare_protocols(data)
    assert len(data['series']['protocols'])==2


def test_registered_collision_and_no_speculative_repairs(caplog):
    from meta_standards_converter.magetab.protocol_export import prepare_protocols
    data=legacy();data['series']['protocols']=[definition('one:P-MTAB-1','x'),definition('two:P-MTAB-1','y'),definition('third',"manufacturerÕs"),definition('fourth',"manufacturer's")]
    prepare_protocols(data)
    names=[p['name'] for p in data['series']['protocols']]
    assert len(names)==len(set(names))==4 and 'P-MTAB-1' in names
    assert 'conflicting' in caplog.text.lower()
    once=deepcopy(data);prepare_protocols(data);assert data==once
