# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Review transformations cannot excuse material or qualifier ownership loss."""
from copy import deepcopy
import pytest
from tests.support.integrity_legacy import reviewed_legacy_tables
from tests.support.deposited_fidelity import assert_deposited_sdrf


def example():
    package = {'source': {'format': 'GEO MINiML'}, 'series': {}, 'sample': [{'iid': 'S1', 'channel': [{}, {}]}]}
    old = [['Source Name', 'Protocol REF', 'Extract Name', 'Assay Name'], ['S1', '', 'S1', 'A'], ['S1', '', 'S1', 'A']]
    new = [['Source Name', 'Comment[msc_sample]', 'Comment[msc_channel]', 'Extract Name', 'Assay Name'],
           ['S1:channel1', 'S1', '0', 'S1:channel1', 'A'], ['S1:channel2', 'S1', '1', 'S1:channel2', 'A']]
    return package, old, new


def test_review_accepts_exact_channel_occurrences_and_empty_unqualified_column():
    package, old, new = example()
    saved = deepcopy((package, old, new))
    a, b, review = reviewed_legacy_tables(package, old, new, {})
    assert_deposited_sdrf(a, b, accession='example')
    assert len(review) == 2
    assert (package, old, new) == saved


def test_review_rejects_wrong_channel_marker():
    package, old, new = example()
    new[1][2] = '1'
    with pytest.raises(AssertionError): reviewed_legacy_tables(package, old, new, {})


def test_review_does_not_remove_qualified_blank_protocol():
    package, old, new = example()
    old[0].insert(2, 'Comment[duration]')
    for row in old[1:]: row.insert(2, '4')
    a, b, _ = reviewed_legacy_tables(package, old, new, {})
    with pytest.raises(AssertionError): assert_deposited_sdrf(a, b, accession='example')


@pytest.mark.parametrize('value', ['lost.txt', 'https://other/file.txt'])
def test_review_cannot_excuse_an_unsupported_file_link(value):
    data = {'series': {}, 'sample': [{'iid': 'S1', 'channel': [{}], 'supplementary_data': [{'value': 'https://source/file.txt'}]}]}
    old = [['Source Name'], ['S1']]
    new = [['Source Name', 'Comment[Unclassified File]'], ['S1', value]]
    with pytest.raises(AssertionError): reviewed_legacy_tables(data, old, new, {})


def test_review_cannot_reorder_an_authored_protocol_path():
    data = {'source': {'format': 'GEO MINiML'}, 'series': {'assay_paths': [{}]},
            'sample': [{'iid': 'S1', 'channel': [{'growth_protocol': 'growth', 'treatment_protocol': 'treat'}]}]}
    old = [['Source Name', 'Protocol REF', 'Protocol REF'], ['S1', 'T', 'G']]
    new = [['Source Name', 'Protocol REF', 'Protocol REF'], ['S1', 'G', 'T']]
    a, b, _ = reviewed_legacy_tables(data, old, new, {'T': 'treat', 'G': 'growth'})
    with pytest.raises(AssertionError): assert_deposited_sdrf(a, b, accession='example')
