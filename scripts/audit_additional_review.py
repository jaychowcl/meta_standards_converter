# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Independent review of additional-cohort outputs; never reconverts evidence."""
from collections import Counter, defaultdict
import json
from pathlib import Path

from scripts.audit_authorship_report import read, build
from tests.support.authorship_capture import write_json
from tests.support.authorship_fidelity import EQUIVALENT, LOSS, UNRESOLVED, idf_cells, idf_comparison, pair_sdrf, scoped_sdrf_cells, table
from tests.support.deposited_fidelity import HEADER_ALIASES, without_derived_files, without_added_protocol_placeholders, assert_deposited_sdrf


def blank_projection(rows, accession):
    """Only observed whitespace literals; no trimming of scientific text."""
    normalized = [list(row) for row in rows]
    permitted = {'E-GEOD-13398': {'Characteristics[site]', 'Term Source REF', 'Term Accession Number'},
                 'E-GEOD-31356': {'Term Source REF'}}.get(accession, set())
    changes = []
    for i, row in enumerate(normalized[1:], 1):
        for j, value in enumerate(row):
            header = HEADER_ALIASES.get(rows[0][j], rows[0][j])
            if header in permitted and value == '  ':
                changes.append({'row': i, 'column': j + 1, 'field': header, 'source': value,
                    'classification': EQUIVALENT, 'reason': 'Exact observed whitespace-only source cell represents missing metadata.'})
                row[j] = ''
    omitted = [j for j,h in enumerate(normalized[0]) if accession == 'E-GEOD-31356' and h == 'Term Source REF'
               and all(row[j] == '' for row in normalized[1:])]
    # Omitting an empty qualifier cannot discard a populated qualifier or owner.
    keep = [j for j in range(len(normalized[0])) if j not in omitted]
    return [[row[j] for j in keep] for row in normalized], changes


def review_deposited(root):
    reports = read(root / 'deposited-results.json'); fields = []; totals = Counter()
    for report in reports:
        accession = report['accession']
        original = {kind: table(root / 'references' / f'{accession}.{kind}.txt') for kind in ('idf', 'sdrf')}
        expected, blank_changes = blank_projection(without_derived_files(original['sdrf']), accession)
        write_json(root / 'deposited' / accession / 'blank-transformations.json', blank_changes)
        totals['whitespace_only_source_cells'] += len(blank_changes)
        for route, detail in report['routes'].items():
            document = read(root / 'deposited' / accession / route / 'payload-0.json')[0]
            actual = without_added_protocol_placeholders(original['idf'], document)
            actual, _ = blank_projection(actual, accession)
            try: verdict = {'verdict': 'preserved', **assert_deposited_sdrf(expected, actual, accession=accession)}
            except AssertionError as exc: verdict = {'verdict': 'discrepancy', 'reason': str(exc)}
            for key, rows in [('source', original['sdrf']), ('projected', expected), ('generated', actual)]:
                verdict[key + '_populated_cells'] = len(scoped_sdrf_cells(rows, accession)[1])
            detail['sdrf'] = verdict
            pair = pair_sdrf(expected, actual, accession)
            for row in pair['comparisons']:
                row['source'] = row['reference']
                joined = 'reference_row' in row['identity'] and 'generated_row' in row['identity']
                if joined and row['source'] is not None:
                    row.update(classification=EQUIVALENT if row['source'] == row['generated'] else LOSS,
                        reason='Same deposited input; explicit row/channel match and attribute/qualifier ownership.')
                else:
                    row.update(classification=UNRESOLVED, reason='No unique row/owner correspondence; retain for review.')
                row.update(case=accession, route=route, mode='off', assessment='curated-input-sdrf', reference_assisted=False)
                fields.append(row)
            detail['sdrf']['cell_counts'] = dict(Counter(row['classification'] for row in pair['comparisons']))
            if route == 'via-json':
                totals['eligible_source_cells'] += len(scoped_sdrf_cells(expected, accession)[1])
                totals.update(Counter(row['classification'] for row in pair['comparisons'] if row['source'] is not None))
        write_json(root / 'deposited' / accession / 'comparison.json', report)
    write_json(root / 'deposited-results.json', reports)
    return fields, dict(totals)


def contact_review(root, cases):
    result = []
    for case in cases:
        name = case['provider'] + '-' + case['accession']
        for mode, detail in case['modes'].items():
            directory = root / 'cases' / name / mode
            if not (directory / 'source-json/payload-0.json').exists() or not (directory / 'replayed-magetab/payload-0.json').exists(): continue
            package = read(directory / 'source-json/payload-0.json')[0]
            document = read(directory / 'replayed-magetab/payload-0.json')[0]
            rows = {r[0]: r[1:] for r in document if r and isinstance(r[0], str)}
            contacts = []
            for index, last in enumerate(rows.get('Person Last Name', [])):
                contacts.append({key: rows.get(key, [])[index] if index < len(rows.get(key, [])) else None
                    for key in ('Person First Name', 'Person Mid Initials', 'Person Last Name', 'Person Email', 'Person Address', 'Person Affiliation', 'Person Roles')})
            references = [{'sample': s['iid'], 'contact_ref': s.get('contact_ref', []), 'contributor_ref': s.get('contributor_ref', [])} for s in package.get('sample', [])]
            counts = Counter((str(p['Person First Name'] or '').casefold(), str(p['Person Last Name'] or '').casefold()) for p in contacts)
            result.append({'case': name, 'mode': mode, 'json_contributors': package.get('contributor', []),
                'json_contributor_count': len(package.get('contributor', [])), 'sample_references': references,
                'idf_contacts': contacts, 'idf_contact_count': len(contacts),
                'repeated_first_last_names': [{'first': key[0], 'last': key[1], 'count': n} for key,n in counts.items() if n > 1]})
    write_json(root / 'contact-review.json', result)
    return result


def review(root, cohort):
    root = Path(root)
    fields, totals = review_deposited(root)
    summary = build(root, cohort=cohort)
    summary['curated_input_sdrf_cell_review'] = totals
    summary['planned_route_dispositions'] = []
    for spec in cohort['inputs']:
        case = spec['provider'] + '-' + spec['accession']
        replay = read(root / 'cases' / case / 'replay.json')
        for mode in ('off', 'curators', 'standard'):
            for route in ('json', 'magetab', 'saved-0'):
                value = replay.get(mode, {}).get(route)
                summary['planned_route_dispositions'].append({'case': case, 'mode': mode, 'route': route,
                    'status': value['status'] if value else 'unavailable',
                    'limitation': replay.get(mode, {}).get(route + '_limitation') if not value else None})
    contact_review(root, summary['cases'])
    # The generic report keeps source fidelity, reference agreement, and input
    # preservation separate. Add cell-level deposited evidence, not only a gate.
    with (root / 'field-discrepancies.jsonl').open('a') as stream:
        for row in fields: stream.write(json.dumps(row, ensure_ascii=False) + '\n')
    summary['field_report_records'] += len(fields)
    write_json(root / 'summary.json', summary)
    write_json(root / 'curated-input-sdrf-fields.json', fields)
    return summary
