# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Scoped library fact selection and consistent run/workflow projection."""
from copy import deepcopy
import re

LIBRARY_FIELDS = frozenset({'library_strategy', 'library_source', 'library_selection', 'library_layout'})
GEOMETRY_FIELDS = frozenset({'nominal_length', 'nominal_sdev', 'orientation'})
FIELDS = LIBRARY_FIELDS | GEOMETRY_FIELDS


def field_name(value):
    name = re.sub(r'[\s-]+', '_', str(value).strip().lower())
    return name if name in FIELDS else None


def path_facts(steps):
    return [{key: c['value']} for step in steps for c in step.get('comments', [])
            if (key := field_name(c.get('name'))) and 'value' in c]


def resolve_library_facts(native, incoming, issues, identity):
    """Incoming records are one accepted source rank, never arbitrary neighbours."""
    from ..metadata.archive_enrichment import informative
    selected = {k: deepcopy(v) for k, v in native.items() if k in FIELDS}
    accepted, conflicts = {}, set()
    for key in FIELDS:
        values = [r[key] for r in incoming if informative(r.get(key))]
        signatures = {' '.join(str(v).split()).casefold() for v in values}
        if len(signatures) > 1:
            conflicts.add(key)
            message = f'{identity}: conflicting same-source {key} facts'
            if message not in issues:
                issues.append(message)
        elif values:
            accepted[key] = deepcopy(values[0])
    if 'library_layout' in conflicts:
        for key in GEOMETRY_FIELDS:
            accepted.pop(key, None)
    layout_changed = ('library_layout' in accepted and str(accepted['library_layout']).casefold()
                      != str(selected.get('library_layout', '')).casefold())
    if layout_changed or GEOMETRY_FIELDS & accepted.keys():
        for key in GEOMETRY_FIELDS:
            selected.pop(key, None)
    selected.update(accepted)
    if str(selected.get('library_layout', '')).upper() == 'SINGLE':
        for key in GEOMETRY_FIELDS:
            selected.pop(key, None)
    return selected


def apply_library_facts(run, steps, facts, *, preserve=()):
    from ..metadata.archive_enrichment import informative
    for key in FIELDS:
        run.pop(key, None)
    run.update(deepcopy(facts))
    for step in steps:
        comments = step.get('comments', [])
        rewritten = []
        for comment in comments:
            key = field_name(comment.get('name'))
            if key is None or key in preserve or (key in LIBRARY_FIELDS and step.get('kind') in {'extract', 'labeled_extract'}
                               and not informative(comment.get('value'))):
                # Missing authored preparation cells must not become native
                # fallback facts that saved-JSON recovery treats as AE evidence.
                rewritten.append(comment)
            elif key in facts:
                value = {**deepcopy(comment), 'value': str(facts[key])}
                if value not in rewritten:
                    rewritten.append(value)
        if 'comments' in step:
            step['comments'] = rewritten


def _shared_facts(choices, issues, identity):
    """Project a shared workflow only where every explicitly bound run agrees."""
    common, conflicts = {}, set()
    for key in FIELDS:
        values = [c.get(key) for c in choices]
        if values and len({str(v) for v in values}) == 1:
            if values[0] is not None:
                common[key] = deepcopy(values[0])
        elif values:
            conflicts.add(key)
    if 'library_layout' not in common:
        for key in GEOMETRY_FIELDS:
            common.pop(key, None)
    if conflicts:
        message = f'{identity}: shared workflow has conflicting run library facts: {sorted(conflicts)}'
        if message not in issues:
            issues.append(message)
    return common


def synchronize_library_facts(data, issues=None):
    """Complete scalar projections; recover saved authored paths by explicit IDs."""
    if data.get('source', {}).get('format') not in {'ENA', 'SRA'}:
        return
    from ..metadata.archive_workflows import path_ids, compatible
    from ..metadata.archive_enrichment import informative
    log_issues = issues is None
    issues = issues if issues is not None else []
    linked = {a['value'] for a in data.get('series', {}).get('accession', [])
              if re.fullmatch(r'E-[A-Z]+-\d+|GSE\d+', a.get('value', ''))}
    by_sample = {s['iid']: s for s in data.get('sample', [])}
    run_index = {(s['iid'], r['run']): r for s in by_sample.values() for r in s.get('sra_run', [])}
    experiment_index = {}
    for (sid, rid), run in run_index.items():
        experiment_index.setdefault((sid, run.get('experiment')), []).append(rid)
    paths_by_run = {}
    for path in data.get('series', {}).get('assay_paths', []):
        refs = {s['sample_ref'] for s in path['steps'] if s.get('sample_ref')}
        if len(refs) != 1:
            continue
        sid = next(iter(refs))
        scope = path_ids(path)
        candidates = {v for v in scope if re.fullmatch(r'[SED]RR\d+', v)}
        if not candidates:
            candidates = {r for exp in scope for r in experiment_index.get((sid, exp), [])}
        for rid in candidates:
            run = run_index.get((sid, rid))
            if run and compatible(scope, {rid, run.get('experiment')}):
                paths_by_run.setdefault((sid, rid), []).append(path)
    assignments = {}
    for sample in data.get('sample', []):
        for run in sample.get('sra_run', []):
            paths = paths_by_run.get((sample['iid'], run['run']), [])
            incoming = []
            for path in paths:
                # Older native JSON can contain a source-authored extract value
                # beside a generated contradictory assay comment. Registered,
                # linked workflow identity is required for this recovery.
                refs = {s.get('protocol_ref', '') for s in path['steps']}
                if any(ref.startswith(acc + ':') for acc in linked for ref in refs):
                    incoming.extend(path_facts([s for s in path['steps'] if s['kind'] in {'extract', 'labeled_extract'}]))
            facts = resolve_library_facts(run, incoming, issues, run['run'])
            apply_library_facts(run, [], facts)
            for path in paths:
                assignments.setdefault(id(path), (path, []))[1].append(facts)
        for key in LIBRARY_FIELDS:
            values = [r.get(key) for r in sample.get('sra_run', [])]
            if values:
                if all(informative(v) for v in values) and len({str(v) for v in values}) == 1:
                    sample[key] = deepcopy(values[0])
                else:
                    sample.pop(key, None)
    for path, choices in assignments.values():
        common = _shared_facts(choices, issues, str(sorted(path_ids(path))))
        apply_library_facts({}, path['steps'], common)
    if log_issues:
        import logging
        for message in issues:
            logging.getLogger(__name__).warning(message)
