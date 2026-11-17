# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Conservative source-text diagnostics and corroborated reversible decoding."""
import json
from copy import deepcopy

_PROSE = {'title', 'description', 'summary', 'overall_design', 'source',
          'growth_protocol', 'treatment_protocol', 'extract_protocol',
          'label_protocol', 'hyb_protocol', 'scan_protocol', 'data_processing'}


def _alternatives(value):
    if not isinstance(value, str): return
    for encoding in ('cp1252', 'latin1'):
        try:
            clean = value.encode(encoding).decode('utf-8')
            if clean != value and clean.encode('utf-8').decode(encoding) == value:
                yield clean, encoding
        except UnicodeError:
            continue


def _text_choice(left, right, path):
    if path.rsplit('/', 1)[-1] not in _PROSE: return None
    for original, selected, chosen in ((left, right, True), (right, left, False)):
        for clean, encoding in _alternatives(original.get('value')):
            if clean == selected.get('value'):
                return chosen, {'provider': selected['provider'], 'kind': 'text_normalization',
                    'accession': selected['accession'], 'metadata': {
                        'path': path, 'statements': [deepcopy(left), deepcopy(right)],
                        'selected': deepcopy(selected), 'original': deepcopy(original),
                        'encoding': encoding, 'reason': 'same scoped statement corroborates reversible UTF-8 decoding'}}
    return None


def _text_issues(package):
    from .preparation import package_mapping
    data = package_mapping(package)
    def walk(value, path):
        if isinstance(value, dict):
            for key, literal in value.items():
                if key in _PROSE and isinstance(literal, str):
                    if next(_alternatives(literal), None) is not None or '\ufffd' in literal:
                        yield f'possible source-text encoding corruption at {path}/{key}; original retained without same-field corroboration'
                elif key in {'channel', 'protocols'}:
                    yield from walk(literal, path + '/' + key)
        elif isinstance(value, (list, tuple)):
            for i, child in enumerate(value): yield from walk(child, path + '/' + str(i))
    yield from walk(data.get('series', {}), '/series')
    for sample in data.get('sample', []):
        yield from walk(sample, '/sample/' + sample['iid'])


def _decode_normalization(value):
    try:
        envelope = json.loads(value)
        if not isinstance(envelope, dict) or envelope.get('schema_version') != '1.0': return None
        record = envelope['record']; fact = record['metadata']
        if record.get('kind') not in {'text_normalization', 'text_normalization_suggestion'} or not isinstance(fact['path'], str): return None
        statements = fact['statements']
        if not isinstance(statements, list) or len(statements) != 2: return None
        if not all(isinstance(s, dict) and all(k in s for k in ('provider', 'accession', 'value')) for s in statements): return None
        choice = _text_choice(*statements, fact['path'])
        if choice and {**choice[1], 'kind': record['kind']} == record: return record
    except (ValueError, TypeError, KeyError): pass
    return None


def _normalization_comments(data, rows=None, warn=None):
    from .preparation import package_mapping
    from .conflicts import _json
    data = package_mapping(data)
    if rows is not None:
        rows[:] = [r for r in rows if not (r and r[0] == 'Comment[MSCTextNormalization]' and all(_decode_normalization(v) for v in r[1:]))]
        values = dict.fromkeys(_json({'schema_version': '1.0', 'record': r}) for r in data.get('extensions', {}).get('insdc', {}).get('records', []) if r.get('kind') in {'text_normalization', 'text_normalization_suggestion'})
        rows.extend(['Comment[MSCTextNormalization]', value] for value in values)
    else:
        kept = []
        for comment in data.get('series', {}).get('comments', []):
            if comment.get('name') == 'MSCTextNormalization':
                record = _decode_normalization(comment.get('value'))
                if record:
                    records = data.setdefault('extensions', {}).setdefault('insdc', {'version': '1.0', 'records': []})['records']
                    if record not in records: records.append(record)
                    continue
                warn('Unrecognized or malformed MSCTextNormalization comment retained without interpretation')
            kept.append(comment)
        data['series']['comments'] = kept
