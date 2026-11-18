# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Shared validation of Atlas 1.2 review evidence, independent of its producer."""
from collections.abc import Mapping

from meta_standards_converter.runtime_contracts import OperationStatusV2


REVIEW_DETAIL_FIELDS = {
    'status', 'reasoning', 'criteria_definitions', 'accession_assessments',
    'review_revision', 'evidence_contract_version',
}


def validate_review_details(review: Mapping, name: str = 'review') -> None:
    """Validate transport and consistency; source grounding remains producer-owned."""
    def require(condition, message):
        if not condition:
            raise ValueError(f'{name}: {message}')

    def text(value):
        return isinstance(value, str) and bool(value.strip())

    status = review.get('status')
    if status is not None:
        status = OperationStatusV2.from_dict(status)
        if status.execution.value in {'failed', 'not_run'}:
            require(review.get('judgement') is None and review.get('confidence') is None,
                    'unsuccessful review cannot carry an accepted judgement or confidence')
            require(not review.get('accession_assessments'), 'unsuccessful review cannot carry accepted assessments')
            require(all(c.get('passed') is None and c.get('confidence') is None
                        for c in review.get('criteria', [])), 'unsuccessful review carries stale criteria')
    if review.get('judgement') is None:
        require(status is not None and status.execution.value in {'failed', 'not_run'},
                'null judgement requires explicit failed or not_run status')
    require(review.get('reasoning') is None or isinstance(review['reasoning'], str), 'reasoning must be text or null')
    for key in ('review_revision', 'evidence_contract_version'):
        value = review.get(key)
        require(value is None or text(value) or type(value) is int and value >= 0, f'invalid {key}')
    definitions = review.get('criteria_definitions', [])
    require(isinstance(definitions, list), 'criteria_definitions must be an array')
    identifiers = []
    for definition in definitions:
        require(isinstance(definition, Mapping) and set(definition) == {'id', 'description'}, 'invalid criterion definition')
        require(text(definition['id']) and text(definition['description']), 'empty criterion definition')
        identifiers.append(definition['id'])
    require(len(identifiers) == len(set(identifiers)), 'duplicate criterion definition')
    assessments = review.get('accession_assessments', [])
    require(isinstance(assessments, list), 'accession_assessments must be an array')
    accessions = []
    for assessment in assessments:
        require(isinstance(assessment, Mapping) and set(assessment) == {
            'accession', 'criteria', 'decision', 'exclusion_basis', 'confidence', 'reason'}, 'invalid accession assessment')
        require(text(assessment['accession']) and text(assessment['reason']), 'empty accession assessment')
        accessions.append(assessment['accession'])
        require(assessment['decision'] in {'qualifies', 'exclude', 'uncertain'}, 'invalid accession decision')
        require(assessment['confidence'] in {'low', 'medium', 'high'}, 'invalid accession confidence')
        require(assessment['exclusion_basis'] in {'none', 'explicit_theme_mismatch', 'unrelated_external_dataset'}, 'invalid exclusion basis')
        criteria = assessment['criteria']
        require(isinstance(criteria, list), 'assessment criteria must be an array')
        seen = []
        for criterion in criteria:
            require(isinstance(criterion, Mapping) and set(criterion) == {'criterion_id', 'status', 'evidence'}, 'invalid criterion assessment')
            seen.append(criterion['criterion_id'])
            require(criterion['status'] in {'meets', 'fails', 'uncertain'}, 'invalid criterion status')
            evidence = criterion['evidence']
            if evidence is not None:
                require(isinstance(evidence, Mapping) and set(evidence) == {'source_kind', 'source_id', 'quote'}, 'invalid evidence reference')
                require(evidence['source_kind'] in {'title', 'publication_text', 'metadata'} and
                        text(evidence['source_id']) and text(evidence['quote']), 'invalid evidence reference values')
            require(criterion['status'] == 'uncertain' or evidence is not None, 'decisive criterion requires evidence')
        require(len(seen) == len(set(seen)) and set(seen) == set(identifiers), 'assessment must cover every criterion exactly once')
        if assessment['decision'] == 'exclude':
            require(assessment['confidence'] == 'high' and assessment['exclusion_basis'] != 'none'
                    and any(c['status'] == 'fails' for c in criteria), 'exclusion lacks supported failed criterion')
        else:
            require(assessment['exclusion_basis'] == 'none', 'non-exclusion cannot carry exclusion basis')
    require(len(accessions) == len(set(accessions)), 'duplicate accession assessment')
