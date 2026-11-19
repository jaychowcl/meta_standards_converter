"""Extend the deposited gate to every supplied IDF field and saved JSON route."""
import json
from pathlib import Path
import pytest

from meta_standards_converter.converters import Converter, PreparationPolicy
from tests.support.authorship_fidelity import assert_idf_preserved, idf_comparison, LOSS, table
from tests.support.deposited_fidelity import DEPOSITED_PANEL, load_deposited_reference, without_derived_files, assert_deposited_sdrf
from tests.support.recorded_replay import ReplayTransport

FIXTURES = Path(__file__).parents[1] / 'fixtures'


@pytest.mark.acceptance
@pytest.mark.parametrize('accession', [*DEPOSITED_PANEL, 'E-MTAB-6486'])
def test_full_idf_and_sdrf_with_explicit_known_failure_budget(accession, tmp_path):
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
    known = {}
    if accession == 'E-MTAB-6486':
        # MSC-AUTH-003 is independently asserted as a strict xfail. Keep the
        # rest of this real study fully gated; do not xfail an entire dataset.
        known = {('study', 'Comment[RelatedExperiment]:' + value, 1, 'Comment[RelatedExperiment]'): (value, None)
                 for value in ['E-MTAB-6483', 'E-MTAB-6432']}
        observed = {(*row['identity'], row['field']): (row['source'], row['generated'])
                    for row in idf_comparison(reference['idf'], document, documented=True)
                    if row['classification'] == LOSS}
        assert observed == known, 'The exact known loss budget changed; review the reproduction, never widen automatically.'
    assert_idf_preserved(reference['idf'], document, documented=True, transformations=known)
    sdrf = next(row[1] for row in document if row[0] == 'SDRF File')
    assert_deposited_sdrf(without_derived_files(reference['sdrf']), sdrf, accession=accession)
