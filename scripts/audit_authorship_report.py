# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Assemble the frozen audit without network access or changes to references."""
from collections import Counter, defaultdict
import csv
import hashlib
import io
import json
from pathlib import Path
import re

from tests.support.authorship_capture import write_json
from tests.support.authorship_fidelity import (
    EQUIVALENT, LOSS, CURATION, UNRESOLVED, ENRICHMENT, idf_cells, idf_comparison,
    table, pair_sdrf, row_identity, scoped_sdrf_cells,
)
from tests.support.authorship_sources import (
    documents, leaf_inventory, citation_facts, geo_checks, geo_sdrf_checks,
    publication_checks, native_coverage, geo_contact_idf_checks, native_contact_checks,
)
from tests.support.deposited_fidelity import assert_deposited_sdrf, without_derived_files


def read(path):
    return json.loads(Path(path).read_text())


def source_link_evidence(docs, accession, reference_idf):
    secondary = [v for row in reference_idf if row and row[0].replace(' ', '') == 'Comment[SecondaryAccession]' for v in row[1:] if v]
    source_values = []
    for doc in docs:
        if doc['reference_assisted']: continue
        if doc['kind'] == 'xml':
            for node in doc['xml'].iter():
                source_values.extend(node.attrib.values())
                if node.text: source_values.append(node.text)
        elif doc['kind'] == 'json': source_values.append(json.dumps(doc['data']))
    explicit = accession in secondary
    source_confirmed = any(re.search(r'(?<![A-Z0-9])' + re.escape(accession) + r'(?!\d)', value) for value in source_values)
    return {'requested_accession': accession, 'reference_secondary_accessions': secondary,
            'reference_explicitly_declares_input': explicit,
            'input_identity_present_in_raw_source': source_confirmed,
            'verified_study_link': explicit and source_confirmed}


def reference_assistance(case, mode, reference):
    replay_trace = case / mode / 'trace-json.json'
    if replay_trace.exists():
        trace = read(replay_trace)['events']
        consumed = [event for event in trace if '/biostudies/' in event['identity']['url']]
        hashes = {event.get('response') for event in consumed}
    else:
        trace = read(case / 'http/trace.json')['events']
        consumed = [event for event in trace if event['phase'] == mode and '/biostudies/' in event['identity']['url']]
        hashes = {event.get('attempt', {}).get('response') for event in consumed}
    return {'trace': str(replay_trace.relative_to(case)) if replay_trace.exists() else 'http/trace.json (acquisition)',
            'reference_requests': [event['identity']['url'] for event in consumed],
            'reference_assisted': bool(consumed),
            'same_comparison_reference_accession_consumed': bool(reference and any(reference in event['identity']['url'] for event in consumed)),
            'same_frozen_reference_bytes_consumed': [] if not reference else [kind for kind in ('idf', 'sdrf')
                if (case.parents[1] / 'references' / f'{reference}.{kind}.txt').exists()
                and hashlib.sha256((case.parents[1] / 'references' / f'{reference}.{kind}.txt').read_bytes()).hexdigest() in hashes]}


def annotate_literal_support(rows, facts):
    """Evidence pointers for candidate values; no cross-entity fidelity credit."""
    exact = defaultdict(list)
    for fact in facts:
        if isinstance(fact['value'], str): exact[fact['value']].append(fact['evidence'])
    for row in rows:
        source_candidates = []
        for key in ('generated', 'reference'):
            value = row.get(key)
            if isinstance(value, str) and value in exact:
                source_candidates.append({'matches': key, 'value': value, 'occurrences': len(exact[value]),
                    'evidence': exact[value][:3], 'scope_verified': False})
        if source_candidates:
            row['source_literal_candidates'] = source_candidates
            row['source'] = source_candidates[0]['value']
            row['source_assessment'] = 'Literal evidence only; ownership and transformation require the separate source checks.'
            if row['classification'] == CURATION:
                row['classification'] = UNRESOLVED
                row['reason'] = 'Reference annotation has a literal source candidate; curation-only origin is not established.'


def build(root, *, candidate=False, cohort=None):
    root = Path(root)
    cases, field_rows = [], []
    for case in sorted((root / 'cases').iterdir()):
        if not (case / 'acquisition.json').exists(): continue
        acquisition = read(case / 'acquisition.json')
        replay = read(case / 'replay.json') if (case / 'replay.json').exists() else {}
        reference = acquisition['reference']
        reference_idf = table(root / 'references' / f'{reference}.idf.txt') if reference else []
        reference_sdrf = table(root / 'references' / f'{reference}.sdrf.txt') if reference else []
        entry = {key: acquisition.get(key) for key in ('provider', 'accession', 'reference', 'physical_requests', 'elapsed_seconds')}
        entry['modes'] = {}
        for mode, phase in acquisition['phases'].items():
            output = case / mode
            detail = {'acquisition_status': phase.get('result', {}).get('status'),
                      'acquisition_limitation': phase.get('limitation'),
                      'replay_limitations': {key: value for key, value in replay.get(mode, {}).items() if key.endswith('_limitation')},
                      'assistance': reference_assistance(case, mode, reference)}
            detail['routes'] = {name: {key: result.get(key) for key in ('status', 'items')}
                                for name, result in replay.get(mode, {}).items()
                                if name in {'json', 'magetab', 'saved-0'} and isinstance(result, dict)}
            package_path = output / 'source-json/payload-0.json'
            direct_path = output / 'replayed-magetab/payload-0.json'
            if not package_path.exists() or not direct_path.exists() or not read(package_path) or not read(direct_path):
                detail['verification_limitation'] = 'No usable captured JSON or replayed direct MAGE-TAB payload.'
                entry['modes'][mode] = detail
                continue
            package, document = read(package_path)[0], read(direct_path)[0]
            docs = documents(case, mode)
            facts = leaf_inventory(docs)
            write_json(output / 'source-leaf-inventory.json', facts)
            write_json(output / 'source-publications.json', citation_facts(docs))
            sdrf = next(row[1] for row in document if row[0] == 'SDRF File')
            detail['population'] = {'samples': len(package.get('sample', [])), 'sdrf_rows': len(sdrf) - 1,
                'idf_populated_cells': len(idf_cells(document)), 'sdrf_populated_cells': len(scoped_sdrf_cells(sdrf, reference or acquisition['accession'])[1]),
                'series_publications': len(package.get('series', {}).get('pubmed_publication') or [])}
            if candidate:
                detail['candidate_json_equals_baseline_capture'] = read(package_path) == read(output / 'baseline-source-json/payload-0.json')
                detail['saved_json_origin'] = 'fresh candidate source-to-JSON conversion'
            elif cohort:
                detail['saved_json_origin'] = 'fresh source-to-JSON conversion'
                captured = output / 'captured-json/payload-0.json'
                detail['replayed_json_equals_capture'] = read(package_path) == read(captured) if captured.exists() else None
            else:
                detail['replayed_json_equals_capture'] = read(package_path) == read(output / 'replayed-json/payload-0.json')
            saved_path = output / 'saved-0-magetab/payload-0.json'
            if not saved_path.exists() or not read(saved_path):
                detail['verification_limitation'] = 'No usable freshly saved JSON export.'
                entry['modes'][mode] = detail
                continue
            saved = read(saved_path)[0]
            route_diff = idf_comparison(document, saved)
            detail['direct_vs_saved_json'] = {'idf_discrepancies': [r for r in route_diff if r['source'] != r['generated']]}
            saved_sdrf = next(row[1] for row in saved if row[0] == 'SDRF File')
            try:
                detail['direct_vs_saved_json']['sdrf'] = assert_deposited_sdrf(sdrf, saved_sdrf, accession=reference or acquisition['accession'])
                detail['direct_vs_saved_json']['sdrf_verdict'] = 'equivalent'
            except AssertionError as exc:
                detail['direct_vs_saved_json'].update(sdrf_verdict='discrepancy', reason=str(exc))
            checks = []
            if acquisition['provider'] == 'geo':
                checks.extend(geo_checks(docs, package))
                checks.extend(geo_sdrf_checks(docs, document))
                checks.extend(geo_contact_idf_checks(docs, document))
            pubs = publication_checks(docs, package, document, acquisition['provider'])
            # Existing authored AE fields outrank PubMed; linked-reference import
            # has the same policy in standard native conversion. This is a source
            # conflict assessment, not a false loss of the PubMed representation.
            if acquisition['provider'] == 'ae' or detail['assistance']['reference_assisted']:
                for row in pubs:
                    if row['classification'] == LOSS:
                        row.update(classification=UNRESOLVED,
                            reason='Authored/deposited citation has precedence over fetched PubMed fields; assess source conflict and reference quality, not automatic loss.')
            checks.extend(pubs)
            if acquisition['provider'] in {'sra', 'ena'}:
                detail['native_coverage'] = native_coverage(docs, package)
                contacts = native_contact_checks(docs, package, document)
                # Standard can prefer different authored contact details. Keep
                # those differences explicit without overriding source precedence.
                if detail['assistance']['reference_assisted']:
                    for row in contacts['checks']:
                        if row['classification'] == LOSS:
                            row.update(classification=UNRESOLVED, reason='Contact presentation after linked-reference enrichment requires precedence assessment.')
                checks.extend(contacts['checks'])
                detail['organization_only_evidence'] = contacts['organization_only_evidence']
            detail['source_checks'] = checks
            detail['source_check_counts'] = dict(Counter(row['classification'] for row in checks))
            detail['source_check_limitations'] = [
                'Independent semantic checks cover GEO study/sample/channel/contact scalars, characteristics, protocols and bound PubMed fields; native checks cover citation fields and study-scoped sample/run coverage.',
                'All other populated IDF/SDRF cells are inventoried below. Literal source support alone does not verify entity ownership or every native projection; unresolved cells are excluded from source-fidelity passes.']
            generated_rows = []
            for cell in idf_comparison([], document):
                cell.update(classification=UNRESOLVED, reason='Generated-field inventory; consult source checks and provenance for a semantic verdict.')
                generated_rows.append(cell)
            for cell in scoped_sdrf_cells(sdrf, reference or acquisition['accession'])[1]:
                generated_rows.append({'scope': 'SDRF', 'identity': {'generated_row': cell['row']}, 'owner': cell['owner'],
                    'field': cell['field'], 'source': None, 'generated': cell['value'], 'reference': None,
                    'classification': UNRESOLVED, 'reason': 'Generated-field inventory; no presumed source/reference correspondence.'})
            annotate_literal_support(generated_rows, facts)
            write_json(output / 'generated-field-inventory.json', generated_rows)
            if reference:
                link = source_link_evidence(docs, acquisition['accession'], reference_idf)
                if acquisition['provider'] == 'ae':
                    link.update(verified_study_link=True, reason='Direct requested AE accession and its deposited reference; assessed as curated input.')
                detail['study_pairing'] = link
                pair = pair_sdrf(without_derived_files(reference_sdrf), sdrf, reference)
                detail['reference_alignment'] = {key: value for key, value in pair.items()
                    if key not in {'reference_inventory', 'generated_inventory', 'comparisons'}}
                if acquisition['provider'] == 'ae':
                    authored = next((d for d in docs if d['kind'] == 'text' and '.idf.txt' in d['url']), None)
                    input_rows = list(csv.reader(io.StringIO(authored['text']), delimiter='\t')) if authored else []
                    ref_rows = idf_comparison(input_rows, document, reference_idf, same_input=bool(authored),
                                             documented=True, source_sdrf=reference_sdrf)
                else:
                    ref_rows = idf_comparison([], document, reference_idf, same_input=False)
                ref_rows.extend(pair['comparisons'])
                if not link['verified_study_link']:
                    for row in ref_rows:
                        row.update(classification=UNRESOLVED, reason='Study link unverified; no reference-agreement verdict.')
                annotate_literal_support(ref_rows, facts)
                detail['reference_comparison_counts'] = dict(Counter(row['classification'] for row in ref_rows))
                write_json(output / 'reference-comparison.json', {'pairing': link, 'sdrf': pair, 'fields': ref_rows})
            else:
                ref_rows = []
                detail['curated_reference_limitation'] = ('No verified curated counterpart for DRP000001; no curated agreement verdict.'
                    if acquisition['accession'] == 'DRP000001' else 'SRP002056 is the paired SRA/ENA source comparison; linked E-GEOD-18729 consumption is logged, not used as an independent declared reference.')
            for assessment, rows in [('source', checks), ('curated-reference', ref_rows), ('inventory', generated_rows)]:
                field_rows.extend(dict(row, case=case.name, mode=mode, assessment=assessment,
                                       reference_assisted=detail['assistance']['reference_assisted']) for row in rows)
            entry['modes'][mode] = detail
        cases.append(entry)
        write_json(case / 'assessment.json', entry)
    deposited = read(root / 'deposited-results.json')
    for result in deposited:
        original_idf = table(root / 'references' / (result['accession'] + '.idf.txt'))
        original_sdrf = table(root / 'references' / (result['accession'] + '.sdrf.txt'))
        for route, detail in result['routes'].items():
            document = read(root / 'deposited' / result['accession'] / route / 'payload-0.json')[0]
            detail['idf'] = idf_comparison(original_idf, document, documented=True, source_sdrf=original_sdrf)
            field_rows.extend(dict(row, case=result['accession'], mode='off', route=route,
                                   assessment='curated-input', reference_assisted=False) for row in detail['idf'])
        # Preserve every excluded source occurrence as an explicitly classified
        # export transformation instead of silently removing it from the audit.
        retained = without_derived_files([original_sdrf[0], list(range(len(original_sdrf[0])))])[1]
        excluded = []
        for index, row in enumerate(original_sdrf[1:], 1):
            for column, value in enumerate(row):
                if column not in retained and value not in ('', None):
                    excluded.append({'scope': 'SDRF', 'identity': {'reference_row': index, 'column': column + 1},
                        'field': original_sdrf[0][column], 'source': value, 'generated': None, 'reference': value,
                        'classification': EQUIVALENT, 'reason': 'Documented exclusion of derived-file nodes and their owned annotations/SAMPLE_FILE comments.',
                        'case': result['accession'], 'mode': 'off', 'assessment': 'curated-input-export-exclusion'})
        write_json(root / 'deposited' / result['accession'] / 'export-exclusions.json', excluded)
        field_rows.extend(excluded)
        write_json(root / 'deposited' / result['accession'] / 'comparison.json', result)
    write_json(root / 'deposited-results.json', deposited)
    source_only = [r for r in field_rows if r['assessment'] == 'source' and r['mode'] != 'standard']
    summary = {'cases': cases, 'source_checks_off_and_curators': dict(Counter(r['classification'] for r in source_only)),
        'source_check_denominator': len(source_only),
        'curated_input_studies': len(deposited), 'curated_input_idf_counts': dict(Counter(r['classification'] for d in deposited for r in d['routes']['via-json']['idf'])),
        'curated_input_idf_source_populated_cells': sum(len(idf_cells(table(root / 'references' / (d['accession'] + '.idf.txt')), imported_scalars=True)) for d in deposited),
        'curated_input_idf_source_verdicts': dict(Counter(r['classification'] for d in deposited for r in d['routes']['via-json']['idf'] if r['source'] is not None)),
        'curated_input_sdrf_preserved': sum(d['routes']['via-json']['sdrf']['verdict'] == 'preserved' for d in deposited),
        'curated_input_sdrf_populated_cells': sum(d['routes']['via-json']['sdrf']['projected_populated_cells'] for d in deposited),
        'curated_input_sdrf_excluded_cells': sum(d['routes']['via-json']['sdrf']['source_populated_cells'] - d['routes']['via-json']['sdrf']['projected_populated_cells'] for d in deposited),
        'field_report_records': len(field_rows)}
    routes = [route for case in cases for mode in case['modes'].values() for route in mode.get('routes', {}).values()]
    summary['live_routes'] = {'expected': len(cohort['inputs']) * 9 if cohort else 81, 'observed': len(routes),
        'execution': dict(Counter(item['execution'] for route in routes for item in route.get('items', []))),
        'validation': dict(Counter(item['validation'] for route in routes for item in route.get('items', []))),
        'status': dict(Counter(route['status'] for route in routes))}
    # Mirrored provider comparison: compare biological/run identity without
    # demanding provider-specific names, files or graph encodings be identical.
    mirror = {}
    mirror_accession = cohort.get('mirror_accession', 'SRP002056') if cohort else 'SRP002056'
    for mode in ('off', 'curators', 'standard'):
        pair = [next((c['modes'].get(mode, {}).get('native_coverage') for c in cases if c['provider'] == provider and c['accession'] == mirror_accession), None)
                for provider in ('sra', 'ena')]
        mirror[mode] = {'sra': pair[0], 'ena': pair[1], 'same_runs': bool(all(pair) and pair[0]['generated_runs'] == pair[1]['generated_runs']),
                        'same_samples': bool(all(pair) and pair[0]['generated_samples'] == pair[1]['generated_samples'])}
    summary[mirror_accession + '_provider_comparison'] = mirror
    write_json(root / 'summary.json', summary)
    with (root / 'field-discrepancies.jsonl').open('w') as stream:
        for row in field_rows: stream.write(json.dumps(row, ensure_ascii=False) + '\n')
    if not candidate and cohort is None:
        write_report(root, summary)
    return summary


def write_report(root, summary):
    cases = summary['cases']
    source = summary['source_checks_off_and_curators']
    total_requests = sum(case['physical_requests'] for case in cases)
    modes = [mode for case in cases for mode in case['modes'].values()]
    captured_equal = sum(mode.get('replayed_json_equals_capture') is True for mode in modes)
    routes_equal = sum(mode.get('direct_vs_saved_json', {}).get('sdrf_verdict') == 'equivalent'
                       and not mode['direct_vs_saved_json']['idf_discrepancies'] for mode in modes)
    routes = summary['live_routes']
    lines = [
        '# MSC conversion audit: source, publication authors and curated metadata', '',
        'Baseline: `942ef678e9625ed61b93438c02fd65b58079f466`. Acquisition: 2026-11-19 UTC. '
        'Production code and integration pins are unchanged. This is a bounded audit, not a blanket certification of every native metadata transformation.', '',
        'Three defects are independently reproduced: native PubMed names lose initials-only given names and suffixes; '
        'untagged GEO characteristics disappear; repeated recognized RelatedExperiment IDF rows retain only their first row. '
        'No production fixes were made.', '',
        '## Results with denominators', '',
        f'- **Execution and replay:** {routes["observed"]}/{routes["expected"]} planned live route outcomes '
        '(9 inputs × 3 presets × source→JSON/direct MAGE-TAB/saved JSON→MAGE-TAB). '
        f'{routes["execution"].get("succeeded", 0)} execute successfully; {routes["validation"].get("valid", 0)} validate. '
        'ERP000263 standard is partial because linked metadata conflict. '
        f'{captured_equal}/27 replayed JSON payloads equal their captures; {routes_equal}/27 direct/saved-JSON IDF/SDRF comparisons are equivalent. '
        'The frozen GSE328265 one-sample contract adds 9 offline route outcomes.',
        f'- **Independent source checks, off + curators:** {source.get(EQUIVALENT, 0)}/{summary["source_check_denominator"]} equivalent, '
        f'{source.get(LOSS, 0)} loss observations, {source.get(UNRESOLVED, 0)} unresolved source conflicts. '
        'These are field assertions, not independent datasets. The 32 losses repeat eight untagged GSE100 `none` annotations across JSON/SDRF and two presets. '
        'Their biological impact in that study is limited; a constructed two-channel exposure/control case proves meaningful text is lost too. '
        'The two unresolved assertions concern E-MTAB-6486’s supplied DOI versus PubMed, not converter loss.',
        f'- **Curated-input preservation:** {summary["curated_input_idf_source_verdicts"].get(EQUIVALENT, 0)}/'
        f'{summary["curated_input_idf_source_populated_cells"]} supplied populated IDF cells preserved or explicitly equivalent; '
        f'{summary["curated_input_idf_source_verdicts"].get(LOSS, 0)} lost RelatedExperiment cells in E-MTAB-6486. '
        'All 18 established panel studies preserve their populated IDF fields. E-MTAB-6486 is an additional study with the documented defect.',
        f'- **Curated-input SDRF:** {summary["curated_input_sdrf_preserved"]}/19 studies preserve all '
        f'{summary["curated_input_sdrf_populated_cells"]:,} eligible populated source cells and their row/owner occurrences. '
        f'{summary["curated_input_sdrf_excluded_cells"]:,} cells are excluded by the existing derived-file export policy; '
        'their original bytes and explicit exclusion inventories remain available. Named blank-literal rules, header aliases, identical ontology declaration consolidation and importer scalar trimming are listed in the comparison helpers.',
        f'- **Live acquisition:** {total_requests} physical requests across {len(cases)} inputs; maximum '
        f'{max(c["physical_requests"] for c in cases)} requests and {max(c["elapsed_seconds"] for c in cases):.2f} seconds for one input. '
        'The enforced ceiling is 200 requests/300 seconds per input across all presets. Prior-phase successful bodies are reused explicitly, preserving a consistent source version. '
        'Any acquisition/replay limitation is retained in the case assessment.', '',
        '## Cohort and pairing', '',
        '| Input | Reference | Samples | SDRF rows off / curators / standard | Reference rows joined off / standard |',
        '| --- | --- | ---: | --- | --- |',
    ]
    for case in cases:
        modes = case['modes']
        populations = [modes.get(m, {}).get('population', {}) for m in ('off', 'curators', 'standard')]
        joined = []
        for mode in ('off', 'standard'):
            alignment = modes.get(mode, {}).get('reference_alignment', {}).get('alignment', {})
            joined.append(f'{len(alignment.get("matches", []))}/{alignment["reference_rows"]}' if alignment else '—')
        lines.append(f'| {case["provider"]} {case["accession"]} | {case["reference"] or "No declared curated reference"} | '
                     f'{populations[0].get("samples", "?")} | ' + ' / '.join(str(p.get('sdrf_rows', '?')) for p in populations) +
                     ' | ' + ' / '.join(joined) + ' |')
    lines.extend([
        '', 'Study pairing uses explicit deposited secondary accessions and captured source identities. SDRF pairing uses sample/channel, assay, run and file identities; '
        'it retains duplicate occurrences and leaves ambiguous joins unresolved. GSE100’s two channels align through their source-name annotations, not row position. '
        'GSE52564 off/curators have 17 sample rows; standard supplies the raw-file detail needed for all 34 reference rows. '
        'SRA SRP045534 standard has 24 file rows versus 12 deposited rows, with 12 extra rows left explicitly visible.', '',
        'SRP002056 has the same four samples and four runs through SRA and ENA in all three presets. Provider-specific file/metadata representation is retained. '
        'ERP000263’s seven study-listed runs and seven samples are preserved; 14 additional run records reached through other links are recorded outside that study-list denominator. '
        'DRP000001 is handled by ENA, with one study-listed sample/run. No curated counterpart was verified.', '',
        '## Publication and contact findings', '',
        '- GEO GSE100, GSE60450 and GSE52564 preserve the ordered PubMed author lists when publication enrichment is enabled. '
        'The frozen GSE328265 source preserves all 11 names, including Börsch, Müller-Miny and Meyer Zu Hörste. Enrichment off intentionally does not fetch missing author details.',
        '- The live native papers checked have source-supported author identity/order and identifier alignment. The constructed Initials-only and Suffix boundaries fail for both SRA and ENA, '
        'so these successful live examples do not establish general author correctness.',
        '- All directly checked GEO contact names, initials, contact details and composed affiliations/addresses are preserved in IDF. '
        'The 12 sample-scoped BioSample contacts in SRP045534 preserve their names and affiliation; repeated display contacts may consolidate. '
        'Organization-only owners are recorded separately and are not treated as paper authors. Existing deterministic tests cover roles, conflicting/duplicate contacts and sample/platform scope.',
        '- The frozen E-GEOD-100 IDF contains `12524544` in Publication DOI; the source-linked PubMed record supplies `10.1038/ng1080`. '
        'The reference is defective at that field. E-MTAB-6486 supplies `UNKNOWN_1625850_E-MTAB-6483` while PubMed supplies `10.1126/scitranslmed.aao2565`. '
        'MSC preserves the supplied placeholder according to source precedence. Its author list differs from PubMed only by a final period.',
        '- SRP002056 standard consumes E-GEOD-18729 and changes full given names to the matching ordered surname/initials representation; this is an equivalent representation with reference assistance. '
        'DRP000001 standard obtains two papers through SRA/BioProject links; off/curators supply no study author list. This is peer-source enrichment, not evidence that DDBJ itself supplied those author names.', '',
        '## Confirmed defects and reproductions', '',
        '| ID | Supplied evidence | Observed output | Coverage |',
        '| --- | --- | --- | --- |',
        '| MSC-AUTH-001 | Initials=AB, LastName=Smith; ForeName=Alice, LastName=Jones, Suffix=Jr | Smith; Alice Jones | SRA/ENA × JSON/IDF × initials/suffix: 8 strict expected failures |',
        '| MSC-AUTH-002 | Untagged GEO Characteristics on separate treated/control channels | Both values absent | JSON and SDRF: 2 strict expected failures; eight live GSE100 literal occurrences corroborate the path |',
        '| MSC-AUTH-003 | Three repeated RelatedExperiment rows | Only the first relationship survives | JSON and IDF: 2 strict expected failures; real E-MTAB-6483 and E-MTAB-6432 links are lost from E-MTAB-6486 |', '',
        'Failures before marking are retained under `tests/*before-xfail.txt`. The focused deposited test budgets only the two exact E-MTAB-6486 losses; '
        'it continues to assert all other annotations and fails if the budget changes. The raw references are never rewritten from generated output.', '',
        '```sh',
        'cd ../source',
        'MSC_AUDIT_PY=/home/chow/projects/embl/runs/gsk_ThematicAtlases_adapter/.venv/bin/python',
        'PYTHONDONTWRITEBYTECODE=1 PYTHONPATH="$PWD/src:$PWD" "$MSC_AUDIT_PY" -m pytest -q -p no:cacheprovider --runxfail tests/archive/citations/test_authorship_audit_regressions.py',
        '```',
        '', 'Expected reproduction outcome at the baseline: 12 failed, 2 passed. A normal run retains those 12 as strict expected failures.', '',
        '## Curated-reference agreement and remaining curation', '',
        'Reference-only annotations remain visible and are excluded from converter-failure counts. The field ledger separates source loss, exact/equivalent agreement, '
        'reference-only annotation, evidence-backed enrichment and unresolved source/version/reference conflicts. “Reference-only” does not establish human authorship of the annotation.', '',
        'The exact comparison denominators differ by route because graph structure, file alternatives and enrichment change the populated cell inventory. '
        'Read `reference_comparison_counts` and `reference_alignment` together in each case’s `assessment.json`; no aggregate “percent correct” is claimed. '
        'Unmatched rows, repeated processing-step correspondence, ontology interpretations and unsupported literal joins remain explicit verification limitations. '
        'A raw value found somewhere in a source response is only a candidate, not proof of correct owner assignment.', '',
        'Source-only off and curators traces for GEO/native inputs contain no linked ArrayExpress retrieval. Standard explicitly consumes the same frozen E-GEOD-60450 '
        'and E-MTAB-308 IDF/SDRF bytes used for comparison; their agreement is reference-assisted. Direct E-MTAB-6486 necessarily starts from its deposited input and belongs to curated-input preservation. '
        'SRP002056 standard consumes E-GEOD-18729, which is logged but is outside the declared independent-reference cohort.', '',
        'Additional ontology labels, experimental designs/factors and interpreted protocol or biological annotations still require curation or better evidence. '
        'The audit inventories all populated IDF/SDRF cells but does not claim exhaustive independent semantic validation of every native projection. '
        'The full deposited-preservation oracle does cover every eligible cell and qualifier/occurrence, including the large E-GEOD-48034/48035 fixtures.', '',
        '## Evidence and rerun', '',
        '- `environment.json`, `environment-packages.json`, `environment-pip-check.txt`: actual runtime and capture revision. Explicit source import path isolates the baseline from the separately edited primary checkout. '
        'Python 3.12 and the existing environment were reused; pip check passes. Installed pytest 8.4.2 is below the source-declared test-extra minimum 9.1.1; no environment upgrade was made.',
        '- `cases/*/http`: original decoded HTTP application-body bytes, hashes, timestamps, phase reuse and redacted request identities. These are not compressed wire-level packets. '
        'Every offline manifest fails closed; no unrecorded network fallback is allowed.',
        '- `references`: original deposited bytes, their historical provenance and current hashes. E-MTAB-6486 uses the complete study fixture, not the reduced parser unit fixture.',
        '- `cases/*/{off,curators,standard}`: captured/replayed JSON, saved JSON, direct and round-trip IDF/SDRF, request traces, independent source checks and per-field source/generated/reference evidence.',
        '- `deposited-results.json`, `deposited/*/export-exclusions.json`: complete curated-input comparisons and intentional file exclusions.',
        '- `field-discrepancies.jsonl`: all material comparisons, classifications, source pointers and generated/reference inventories. `summary.json` supplies denominators and route limitations.',
        '- `tests/final-focused.txt` and XML: final focused validation (300 passed, 12 strict expected failures in this bundle); '
        '`tests/runxfail.txt`: unmasked confirmed defects (12 failed, 2 passed). Earlier harness setup failures are retained and excluded from metadata failure counts. '
        'They comprise AE handler selection before any network request, initial use of a reduced parser fixture, typed-payload serialization and initially unrecognized source aliases/blanks.', '',
        'From the audit source worktree, use the same Python executable:', '',
        '```sh',
        '"$MSC_AUDIT_PY" scripts/audit_authorship.py replay --bundle ../run-20261119',
        '"$MSC_AUDIT_PY" scripts/audit_authorship.py deposited --bundle ../run-20261119',
        '"$MSC_AUDIT_PY" scripts/audit_authorship.py frozen --bundle ../run-20261119',
        '"$MSC_AUDIT_PY" scripts/audit_authorship.py report --bundle ../run-20261119',
        '# New live evidence must use a fresh bundle; capture refuses an existing case.',
        '"$MSC_AUDIT_PY" scripts/audit_authorship.py capture --bundle ../new-capture',
        '```', '',
        'Semantic basis: [MAGE-TAB 1.1](https://www.ebi.ac.uk/biostudies/misc/MAGE-TABv1.1_2011_07_28.pdf) and '
        '[NLM PubMed Author model](https://dtd.nlm.nih.gov/ncbi/pubmed/doc/out/250101/el-Author.html).', '',
    ])
    validation = root / 'validation.json'
    if validation.exists():
        counts = read(validation)['latest_result_per_unique_test']
        lines[7:7] = [f'Final validation across deduplicated focused tests: **{counts["passed"]} passed, '
                      f'{counts["skipped_expected_failure"]} strict expected failures, {counts["failed"]} unexpected failures**. '
                      'The 300-pass full focused run and subsequent 49/39-pass targeted runs overlap; see `validation.json`.', '']
    (root / 'REPORT.md').write_text('\n'.join(lines))
