# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Timing labels must not turn reused responses or failures into successful live runs."""
from scripts.audit_live_integrity import request_counts, summary, COHORT


def test_request_reuse_and_repeated_attempts_are_distinct():
    event = {'identity': {'url': 'https://example.org/x'}, 'physical': True, 'attempt': {'status': 503}}
    result = request_counts([event, dict(event, attempt={'status': 200}), dict(event, physical=False)])
    assert result == {'physical_requests': 2, 'reused_responses': 1,
        'repeated_physical_identities': 1, 'response_statuses': {'503': 1, '200': 1}}


def test_timing_statistics_keep_failed_outcomes_and_detect_variation():
    rows = [{'wall_seconds': wall, 'cpu_seconds': 1, 'peak_rss_bytes': 3,
        'payload_sha256': 'same', 'items': [{'execution': 'failed'}]} for wall in (3, 1, 2)]
    result = summary(rows)
    assert (result['median_seconds'], result['min_seconds'], result['max_seconds']) == (2, 1, 3)
    assert result['outcomes_identical']
    rows[-1]['items'] = [{'execution': 'succeeded'}]
    assert not summary(rows)['outcomes_identical']


def test_live_cohort_has_five_new_provider_inputs():
    assert len(COHORT['inputs']) == 5
    assert {s['evidence'] for s in COHORT['inputs']} == {'live'}
    assert len({s['reference'] for s in COHORT['inputs']}) == 3
