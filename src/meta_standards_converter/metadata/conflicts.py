# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Portable serialization of existing scientific conflict evidence."""
from copy import deepcopy
import json


def _json(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(',', ':'))


def _valid_record(record):
    if not isinstance(record, dict) or record.get('kind') != 'source_conflict': return False
    fact = record.get('metadata')
    if not isinstance(fact, dict): return False
    if not all(isinstance(fact.get(k), str) and fact[k] for k in ('field', 'path', 'reason')): return False
    statements = fact.get('statements')
    if not isinstance(statements, list) or len(statements) < 2: return False
    if not all(isinstance(s, dict) and all(k in s for k in ('provider', 'accession', 'value')) for s in statements): return False
    selected = fact.get('selected')
    return isinstance(selected, dict) and any(all(selected.get(k) == s[k] for k in ('provider', 'accession', 'value')) for s in statements)


def _conflicts(data, sample=None):
    from .preparation import package_mapping
    data = package_mapping(data)
    records = data.get('extensions', {}).get('insdc', {}).get('records', [])
    result, seen = [], set()
    ids = ({str(sample.get('iid'))} | {str(a.get('value')) for a in sample.get('accession', [])}) if sample else set()
    runs = {str(r.get(k)) for r in sample.get('sra_run', []) for k in ('run', 'experiment')} if sample else set()
    for record in records:
        if not _valid_record(record): continue
        path = record['metadata']['path']
        if sample and path.startswith('/sample/') and not any(path.startswith('/sample/' + sid + '/') for sid in ids) and not any('/' + rid + '/' in path for rid in runs):
            continue
        key = _json(record)
        if key not in seen:
            seen.add(key); result.append(deepcopy(record))
    return result


def _decode_comment(value):
    try:
        envelope = json.loads(value)
        if isinstance(envelope, dict) and envelope.get('schema_version') == '1.0' and _valid_record(envelope.get('record')):
            return envelope['record']
    except (ValueError, TypeError):
        pass
    return None


def _insert_conflict_comments(data, rows):
    records = _conflicts(data)
    retained = []
    for row in rows:
        if row and row[0] == 'Comment[MSCConflict]':
            values = [v for v in row[1:] if _decode_comment(v) is None]
            if values: retained.append([row[0], *values])
        elif not row or row[0] != 'Comment[MSCConflictSummary]': retained.append(row)
    for record in records:
        retained.append(['Comment[MSCConflict]', _json({'schema_version': '1.0', 'record': record})])
    if records:
        retained.append(['Comment[MSCConflictSummary]', *[
            f"{r['metadata']['path']}: selected {r['metadata']['selected']['provider']} {r['metadata']['selected']['value']!r}; {r['metadata']['reason']}"
            for r in records]])
    rows[:] = retained


def _restore_conflict_comments(data, warn):
    kept, found = [], []
    for comment in data.get('series', {}).get('comments', []):
        if comment.get('name') == 'MSCConflict':
            record = _decode_comment(comment.get('value'))
            if record is not None:
                found.append(record); continue
            warn('Unrecognized or malformed MSCConflict comment retained without interpretation')
        if comment.get('name') != 'MSCConflictSummary': kept.append(comment)
    data['series']['comments'] = kept
    if found:
        extension = data.setdefault('extensions', {}).setdefault('insdc', {'version': '1.0', 'records': []})
        records = extension.setdefault('records', [])
        for record in found:
            if record not in records: records.append(record)
