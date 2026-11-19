# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Canonical workflows are shared by fresh conversion and legacy JSON export."""
from copy import deepcopy
from dataclasses import replace

import pytest

from meta_standards_converter.metadata.workflows import prepare_workflows
from meta_standards_converter.converters import Converter, PreparationPolicy
from tests.magetab.test_metadata_integrity import geo, tables


@pytest.mark.parametrize('study', ['GSE13398', 'GSE31356', 'GSE52778'])
def test_preparation_is_network_free_idempotent_and_input_immutable(study):
    source = geo(study)
    original = deepcopy(source.to_mapping())
    prepared = prepare_workflows(source)
    assert prepared.series.assay_paths
    assert prepare_workflows(prepared) == prepared
    assert source.to_mapping() == original
    assert tables(source) == tables(prepared) == tables(prepared.to_mapping())


def test_prepared_authored_repeated_applications_and_qualifiers_stay_ordered():
    prepared = prepare_workflows(geo('GSE52778'))
    path = prepared.series.assay_paths[0]
    applications = [s for s in path.steps if s.kind == 'protocol_application']
    # Authored ordering can differ from a provider's generated default.
    revised = replace(path, steps=(path.steps[0], applications[-1], applications[0],
                                  applications[-1], *path.steps[1:]), document='authored.sdrf.txt')
    authored = replace(prepared, series=replace(prepared.series, assay_paths=(revised,)))
    assert prepare_workflows(authored) == authored


def test_unified_json_contains_the_same_prepared_graph_as_direct_export():
    source = geo('GSE13398')
    policy = PreparationPolicy(enrichment='off', expand_studies=False)
    result = Converter().convert(source, out_type='json', preparation=policy)
    assert result.status == 'complete', result.to_dict()
    saved = result.items[0].payload[0]
    assert saved.series.assay_paths == prepare_workflows(source).series.assay_paths
    assert tables(saved) == tables(source)
