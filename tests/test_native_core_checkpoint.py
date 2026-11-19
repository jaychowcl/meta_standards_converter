# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Reuse constructed native entities without weakening valid-core recovery."""
from copy import deepcopy
from types import SimpleNamespace
import pytest

from meta_standards_converter.converters.archive_results import _load_native
from meta_standards_converter.metadata.retrieval_budget import _budget_scope, _BudgetExhausted
from meta_standards_converter.miniml import insdc_support as support
from meta_standards_converter.miniml.ena_parser import ENAParser
from meta_standards_converter.miniml.sra_parser import SRAParser
from meta_standards_converter.magetab.constructor import AEConstructor
from tests.support.archive import fixture_records


@pytest.mark.parametrize('provider,parser_type',[('ena',ENAParser),('sra',SRAParser)])
@pytest.mark.parametrize('interrupted',[False,True])
def test_native_core_is_constructed_once_and_supplements_preserve_goldens(monkeypatch,provider,parser_type,interrupted):
    records=fixture_records(provider)
    extra={'provider':'biosamples','kind':'sample','accession':'SAMN14218700',
           'metadata':{'accession':'SAMN14218700','characteristics':{'checkpoint fixture':[{'text':'preserved'}]}}}
    complete=deepcopy(records);complete.linked.append(extra)
    expected=parser_type().parse(complete)
    expected_bundle=AEConstructor().miniml2magetab(expected,platform_handler='generic')
    sample_record=support.sample_record;calls=[]
    def counted(*args):calls.append(1);return sample_record(*args)
    monkeypatch.setattr(support,'sample_record',counted)
    def fetch(seed,*,_on_core=None):
        _on_core(records);records.linked.append(deepcopy(extra))
        if interrupted:raise _BudgetExhausted('supplementary provider interrupted')
        return records
    with _budget_scope(max_requests=100):
        actual,issues=_load_native(SimpleNamespace(fetch=fetch),parser_type(),records.seed)
    assert len(calls)==1
    assert actual.to_mapping()==expected.to_mapping()
    assert AEConstructor().miniml2magetab(actual,platform_handler='generic')==expected_bundle
    assert any('valid native candidate' in issue for issue in issues)==interrupted


def test_changed_core_identity_invalidates_construction_checkpoint(monkeypatch):
    import xml.etree.ElementTree as ET
    records=fixture_records('ena')
    sample_record=support.sample_record;calls=[]
    def counted(*args):calls.append(1);return sample_record(*args)
    monkeypatch.setattr(support,'sample_record',counted)
    def fetch(seed,*,_on_core=None):
        _on_core(records)
        sample=next(n for root in records.xml for n in root if n.tag=='SAMPLE')
        node=sample.find('TITLE')
        if node is None:node=ET.SubElement(sample,'TITLE')
        node.text='updated source title'
        return records
    with _budget_scope(max_requests=100):
        actual,issues=_load_native(SimpleNamespace(fetch=fetch),ENAParser(),records.seed)
    assert len(calls)==2
    assert actual.samples[0].title=='updated source title'


def test_invalid_supplementary_projection_keeps_the_valid_checkpoint(monkeypatch):
    from meta_standards_converter.miniml.model import MINiMLModelError
    records=fixture_records('sra')
    finish=support.finish
    def project(provider,records,*args):
        if getattr(records,'invalid_supplement',False):
            raise MINiMLModelError('invalid supplementary value')
        return finish(provider,records,*args)
    monkeypatch.setattr(support,'finish',project)
    def fetch(seed,*,_on_core=None):
        _on_core(records);records.invalid_supplement=True;return records
    with _budget_scope(max_requests=100):
        actual,issues=_load_native(SimpleNamespace(fetch=fetch),SRAParser(),records.seed)
    assert len(actual.samples)==1 and not actual.validate()
    assert any('projection' in issue and 'retained' in issue for issue in issues)
    coverage=actual.to_mapping()['extensions']['insdc']['records'][-1]
    assert coverage['kind']=='retrieval_coverage' and coverage['metadata']['complete'] is False
