# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from copy import deepcopy
from meta_standards_converter.metadata.archive_diagnostics import consistency_issues


def test_explicit_modality_conflict_reports_run_and_fields_without_rewriting():
    data = {'series':{'iid':'ERP1'},'sample':[{'iid':'s','sra_run':[{'run':'ERR1',
        'library_strategy':'ChIP-Seq','library_source':'GENOMIC',
        'library_protocol':'single-cell RNA preparation with reverse transcription of mRNA'}]}]}
    before = deepcopy(data)
    issues = consistency_issues(data)
    assert len(issues) == 1
    assert 'ERR1' in issues[0] and 'library_source' in issues[0] and 'library_protocol' in issues[0]
    assert data == before


def test_explicit_associated_antibody_targets_conflict_but_titles_do_not_guess_targets():
    data = {'series':{},'sample':[{'iid':'s','title':'H3K27ac assay', 'channel':[{'characteristics':[
        {'name':'antibody target','value':'H3K4me3'}, {'name':'chip target','value':'H3K27ac'}]}]}]}
    assert 'target' in consistency_issues(data)[0]
    data['sample'][0]['channel'][0]['characteristics'].pop()
    assert not consistency_issues(data)
