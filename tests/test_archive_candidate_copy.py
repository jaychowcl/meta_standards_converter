# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Share certified immutable branches while keeping mutable candidate state isolated."""
from meta_standards_converter.miniml import model
from meta_standards_converter.metadata.session import PreparationSession
from meta_standards_converter.converters.archive_results import ArchiveImportResult, StudyImportOutcome


def test_construction_reuses_owned_frozen_fields_without_thawing():
    frozen=model._FrozenJSONMapping({'nested':[{'value':'keep'}]})
    accession=model.Accession('A1',extras=frozen)
    package=model.MINiMLPackage(model.Series(iid='S1'),model.SourceInfo('ENA'),extensions=frozen)
    assert accession.extras is frozen
    assert package.extensions is frozen
    editable=package.to_mapping();editable['extensions']['nested'][0]['value']='changed'
    assert frozen['nested'][0]['value']=='keep'


def test_operation_memo_shares_only_certified_model_children():
    from tests.test_native_archive_enrichment import native
    original=native()
    result=ArchiveImportResult('SRP1',[StudyImportOutcome('SRP1','PRJ1','complete',package=original)])
    with PreparationSession() as session:
        first=session.retrieve('source',lambda:result)
        second=session.retrieve('source',lambda:None)
    a,b=first.packages[0],second.packages[0]
    assert a is not b and a is not original
    assert a.samples is b.samples
    first.studies[0].issues.append('changed')
    assert second.studies[0].issues==[] and result.studies[0].issues==[]
    a._archive_source_records[0]['metadata']['private_change']='changed'
    assert 'private_change' not in b._archive_source_records[0]['metadata']
    assert 'private_change' not in original._archive_source_records[0]['metadata']
    assert a.to_mapping()==b.to_mapping()


def test_foreign_mutable_values_do_not_gain_sharing_from_a_frozen_container():
    class Mutable:
        def __init__(self):self.values=['original']
    child=model.Accession(Mutable())
    assert not model._certified_immutable(child)
    with PreparationSession() as session:
        first=session.retrieve('foreign',lambda:{'child':child})
        second=session.retrieve('foreign',lambda:None)
    first['child'].value.values.append('changed')
    assert second['child'].value.values==['original']
