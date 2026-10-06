# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Indexed residual matching keeps greedy source occurrence semantics."""
from copy import deepcopy
import pytest
from meta_standards_converter.miniml import archive_residuals as residuals


@pytest.mark.parametrize('paths',[False,True])
def test_reversed_identity_lists_do_not_compare_all_pairs(monkeypatch,paths):
    n=64
    source=[{'iid':f'S{i}','characteristics':[{'name':'tissue','value':'lung'}]} for i in range(n)]
    if paths:
        source=[{'steps':[{'kind':'source','sample_ref':f'S{i}'},
                         {'kind':'assay','name':f'ERX{i}'}]} for i in range(n)]
    target=list(reversed(deepcopy(source)))
    for node in target:node['extension']='richer'
    original=residuals.contains;calls=[]
    def counted(*args):
        calls.append(1);return original(*args)
    monkeypatch.setattr(residuals,'contains',counted)
    assert residuals.contains(source,target)
    assert len(calls)<n*30


def test_residual_diff_prunes_exact_and_keeps_changed_occurrences(monkeypatch):
    n=64
    source=[{'iid':f'S{i}','title':'original','extra':'retain'} for i in range(n)]
    target=[{'iid':f'S{i}','title':'changed'} for i in reversed(range(n))]
    original=residuals.contains;calls=[]
    def counted(*args):calls.append(1);return original(*args)
    monkeypatch.setattr(residuals,'contains',counted)
    assert residuals.diff(source,target)==source
    assert len(calls)<n*40


def test_duplicate_richer_matches_keep_greedy_order_and_do_not_reuse_occurrences():
    # Matching a less specific first occurrence consumes the first richer target,
    # exactly as the historical occurrence algorithm did.
    source=[{'name':'a'},{'name':'a','value':'x'}]
    target=[{'name':'a','value':'x'},{'name':'a','value':'y'}]
    assert residuals.contains(source,target) is False
    assert residuals.diff(source,target)==[{'name':'a','value':'x'}]
    assert not residuals.contains([{'iid':'S1'}]*2,[{'iid':'S1'}])
    assert residuals.diff([{'iid':'S1','title':'supplied'}],
                          [{'iid':'S1','title':'a'},{'iid':'S1','title':'b'}])==[{'iid':'S1','title':'supplied'}]


def test_unhashable_and_unbound_records_keep_full_comparison_semantics():
    source=[{'name':{'value':'complex'},'members':[2,1]},{'unknown':{'nested':[1,2]}}]
    target=list(reversed(deepcopy(source)))
    before=deepcopy((source,target))
    assert residuals.contains(source,target)
    assert residuals.diff(source,target) is None
    assert (source,target)==before


def test_empty_expected_fields_do_not_narrow_richer_candidates():
    source=[{'steps':[], 'value':'', 'iid':None}]*16
    target=[{'steps':[{'kind':'sample','sample_ref':str(i)}]} for i in range(16)]
    assert residuals.contains(source,target)
    assert residuals.diff(source,target) is None
