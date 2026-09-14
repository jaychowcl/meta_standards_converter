"""Submission presentation policies for explicitly bound native file inventories."""
from copy import deepcopy
import json
import logging

logger = logging.getLogger(__name__)


def source_annotations(paths, sample_ids, record, comments):
    """Source-only results are annotations, never attached to an arbitrary run."""
    from ..miniml.archive_residuals import contains
    kept, annotations = [], {}
    for path in paths:
        steps = path.get('steps', [])
        bound = {s['sample_ref'] for s in steps if s.get('sample_ref')}
        files = [s for s in steps if s['kind'].startswith('derived_')]
        if (len(bound) != 1 or not bound <= sample_ids or not files
                or any(s['kind'] not in {'source', 'sample', 'protocol_application',
                                       'derived_array_data_file', 'derived_array_data_matrix_file'} for s in steps)):
            kept.append(path)
            continue
        sid = next(iter(bound))
        source = next((s for s in steps if s['kind'] in ('source', 'sample')), None)
        if source is None:
            kept.append(path)
            continue
        if sum(s['kind'] in ('source', 'sample') for s in steps) != 1:
            kept.append(path)
            continue
        binding = json.dumps(source, sort_keys=True, ensure_ascii=False)
        entry = annotations.setdefault(binding, (sid, deepcopy(source), []))
        for index, node in enumerate(steps):
            if node not in files:
                continue
            value = record(node)
            methods = [s for s in steps[:index] if s['kind'] == 'protocol_application']
            if methods:
                value['PROCESSING'] = json.dumps(methods, sort_keys=True, ensure_ascii=False)
            inputs = [record(s) for s in steps[:index] if s['kind'].startswith('derived_')]
            if inputs:
                value['INPUTS'] = json.dumps(inputs, sort_keys=True, ensure_ascii=False)
            context = {k: v for k, v in path.items() if k != 'steps'}
            if context:
                value['CONTEXT'] = json.dumps(context, sort_keys=True, ensure_ascii=False)
            values = entry[2]
            if value not in values:
                values.append(value)
    attached = set()
    for path in kept:
        for step in path.get('steps', []):
            sid = step.get('sample_ref')
            if step['kind'] not in ('source', 'sample'):
                continue
            for binding, (bound, source, values) in annotations.items():
                if bound == sid and contains(source, step):
                    for value in values:
                        step.setdefault('comments', []).extend(comments(value, 'SAMPLE_FILE_', fixed=True))
                    attached.add(binding)
    # A sample with no acquisition retains a biological row, not an invented run.
    for binding in annotations.keys() - attached:
        _, source, values = annotations[binding]
        for value in values:
            source.setdefault('comments', []).extend(comments(value, 'SAMPLE_FILE_', fixed=True))
        kept.append({'steps': [source]})
    return kept


def _set(record):
    repository = record.get('_repository', '')
    role = str(record.get('ROLE', '')).upper()
    if repository == 'ArrayExpress': rank = 0
    elif repository == 'GEO': rank = 1
    elif role == 'GENERATED_FILE': rank = 2
    elif not role: rank = 3  # Legacy unlabelled linked inventory; no asserted repository.
    elif role in ('ORIGINAL', 'SUBMISSION_FILE'): rank = 4
    else: rank = 5
    return rank, repository, record.get('_source', '')


def primary_sets(groups, key, is_fastq, combine, archives):
    families = {}
    for identity, (path, scope, files) in groups.items():
        scan = next(i for i, s in enumerate(path['steps']) if s['kind'] == 'scan')
        acquisition = key({**path, 'steps': path['steps'][:scan+1]})
        families.setdefault(acquisition, []).append((identity, path, scope, files))
    selected = {}
    for family in families.values():
        sets = {}
        for _, _, _, files in family:
            for value in files:
                if str(value.get('FORMAT', '')).lower() in ('fastq', 'fastq.gz'):
                    sets.setdefault(_set(value), []).append(value)
        def roles(values):
            return {(str(f.get('LANE', '')), str(f['READ_INDEX'])) for f in values if f.get('READ_INDEX')}
        required = set().union(*(roles(v) for v in sets.values())) if sets else set()
        complete = {k: v for k, v in sets.items() if v and all(is_fastq(f) for f in v)
                    and (not required or roles(v) >= required)}
        if not complete:
            continue
        priority = min(k[0] for k in complete)
        choices = [k for k in complete if k[0] == priority]
        if len(choices) != 1:
            logger.warning('Ambiguous primary read sets; retaining explicitly distinct inventories')
            continue
        chosen = choices[0]
        if any(k[0] < priority for k in sets):
            logger.warning('Incomplete preferred read set; using complete %s inventory', chosen)
        logger.debug('Primary read set %s for %s', chosen, family[0][2])
        for identity, path, scope, files in family:
            reads = [f for f in files if is_fastq(f) and _set(f) == chosen]
            for value in files:
                if is_fastq(value) and _set(value) != chosen:
                    combine(archives[scope], value)
            if reads:
                selected[identity] = reads
            elif any(s['kind'].startswith('derived_') for s in path['steps']):
                # A genuinely distinct processing branch must not be discarded
                # or rebound to different raw inputs by presentation preference.
                selected[identity] = [f for f in files if is_fastq(f)]
            else:
                selected[identity] = None
    return selected


def prune_empty_file_columns(rows):
    """Keep all occurrences of a populated field so record columns stay aligned."""
    if not rows:
        return rows
    prefixes = ('Comment[ARCHIVE_FILE_', 'Comment[SUBMITTED_FILE_', 'Comment[SAMPLE_FILE_',
                'Comment[FASTQ_ALTERNATIVE_')
    populated = {label for i, label in enumerate(rows[0]) if any(row[i] not in ('', None) for row in rows[1:])}
    keep = [i for i, label in enumerate(rows[0]) if not label.startswith(prefixes) or label in populated]
    return [[row[i] for i in keep] for row in rows]
