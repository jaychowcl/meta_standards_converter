# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Extend the deposited gate to every supplied IDF field and saved JSON route."""
import json
from pathlib import Path
import pytest

from meta_standards_converter.converters import Converter, PreparationPolicy
from tests.support.authorship_fidelity import assert_idf_preserved, table
from tests.support.deposited_fidelity import DEPOSITED_PANEL, load_deposited_reference, without_derived_files, without_added_protocol_placeholders, assert_deposited_sdrf
from tests.support.recorded_replay import ReplayTransport

FIXTURES = Path(__file__).parents[1] / 'fixtures'


@pytest.mark.acceptance
@pytest.mark.parametrize('accession', [*DEPOSITED_PANEL, 'E-MTAB-6486'])
def test_full_idf_and_sdrf_preserve_supplied_fields(accession, tmp_path):
    if accession in DEPOSITED_PANEL:
        reference = load_deposited_reference(FIXTURES / 'fidelity', accession)
        original = FIXTURES / 'fidelity' / (accession + '.idf.txt')
    else:
        original = FIXTURES / 'studies/E-MTAB-6486/inputs' / (accession + '.idf.txt')
        reference = {kind: table(original.with_name(accession + '.' + kind + '.txt')) for kind in ('idf', 'sdrf')}
    isolation = ReplayTransport(FIXTURES / 'recorded_provider/success.manifest.json')
    policy = PreparationPolicy(enrichment='off', expand_studies=False, linked_metadata=False, peer_provider=False)
    with isolation.offline():
        imported = Converter().convert(original, out_type='json', preparation=policy)
        assert imported.items[0].execution == 'succeeded'
        saved = tmp_path / 'supplied.json'
        saved.write_text(json.dumps(imported.items[0].payload[0].to_mapping(), ensure_ascii=False))
        result = Converter().convert(saved, out_type='magetab', preparation=policy)
    isolation.assert_complete()
    item = result.items[0]
    assert item.execution == 'succeeded' and item.validation == 'valid', result.to_dict()
    document = item.payload[0]
    assert_idf_preserved(reference['idf'], document, documented=True)
    sdrf = without_added_protocol_placeholders(reference['idf'], document)
    assert_deposited_sdrf(without_derived_files(reference['sdrf']), sdrf, accession=accession)
