# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Select scoped IDF contacts and consolidate presentation facts, not identities."""
from copy import deepcopy
import re
import json


def _refs(entity):
    return {r['ref'] for field in ('contact_ref', 'contributor_ref')
            for r in entity.get(field, []) if r.get('ref')}


def contributors(data):
    series = data.get('series') or {}
    eligible = _refs(series)
    for sample in data.get('sample', []):
        eligible.update(_refs(sample))
    for protocol in series.get('protocols', []):
        eligible.update(protocol.get('contacts', []))
        eligible.update(protocol.get('performers', []))
    for path in series.get('assay_paths', []):
        eligible.update(s['performer'] for s in path.get('steps', []) if s.get('performer'))
    platform_only = set().union(*(_refs(p) for p in data.get('platform', []))) - eligible
    native = data.get('source', {}).get('format') in ('SRA', 'ENA')
    result = []
    nested = [*series.get('contributor', []), *series.get('contact', []),
              *(c for sample in data.get('sample', []) for c in sample.get('contact', []))]
    for c in [*data.get('contributor', []), *nested]:
        iid = c.get('iid')
        if iid and iid in platform_only and c not in nested:
            continue
        if native and not any(c.get(k) for k in ('person', 'email', 'phone', 'fax')):
            continue
        result.append(deepcopy(c))
    return result


def _normal(value):
    return ' '.join(str(value or '').split()).casefold()


def identity_evidence(contributor):
    tokens = set()
    if contributor.get('iid'):
        tokens.add(('iid', contributor['iid']))
    for accession in contributor.get('accession', []):
        if accession.get('database') and accession.get('value'):
            tokens.add(('accession', accession['database'], accession['value']))
    for relation in contributor.get('relation', []):
        if str(relation.get('type', '')).casefold() in {'same as', 'same_as', 'equivalent'} and relation.get('target'):
            tokens.add(('iid', relation['target']))
    return tokens


def _facts_compatible(a, b):
    if not a or not b:
        return True
    if isinstance(a, dict) and isinstance(b, dict):
        return all(_facts_compatible(a.get(k), b.get(k)) for k in a.keys() | b.keys())
    return a == b


def _compatible(a, b):
    for key in a.keys() | b.keys():
        if key == '_identity':
            continue
        if key == '_facts':
            if not _facts_compatible(json.loads(a.get(key) or '{}'), json.loads(b.get(key) or '{}')):
                return False
        elif a.get(key) and b.get(key) and _normal(a[key]) != _normal(b[key]):
            return False
    return True


def unique_people(columns):
    """Exact exported facts first; reliable identity joins require compatible facts."""
    unique = []
    for values in columns:
        if tuple(_normal(values.get(k)) for k in ('First Name', 'Last Name', 'Email', 'Affiliation')) == (
                'geo', 'curators', 'geo-group@ncbi.nlm.nih.gov', 'ncbi'):
            continue
        unique.append(deepcopy(values))
    groups, owner = {}, {}
    for i, values in enumerate(unique):
        tokens = set(values.get('_identity', ()))
        email, first, last = (_normal(values.get(k)) for k in ('Email', 'First Name', 'Last Name'))
        if (first and last and re.fullmatch(r'[^\s@]+@[^\s@]+\.[^\s@]+', email)
                and email.split('@')[0] not in {'info', 'contact', 'support', 'admin', 'office', 'help', 'team'}):
            tokens.add(('email', email, first, last))
        connected = {owner[token] for token in tokens if token in owner}
        leader = min(connected | {i})
        members = {i}
        for previous in connected:
            members.update(groups.pop(previous))
        groups.setdefault(leader, set()).update(members)
        for token, previous in list(owner.items()):
            if previous in connected:
                owner[token] = leader
        for token in tokens:
            owner[token] = leader
    removed = set()
    for group in groups.values():
        indices = sorted(group)
        # Check the entire component; sparse facts cannot bridge conflicts.
        if not all(_compatible(unique[i], unique[j]) for i in indices for j in indices):
            continue
        target = unique[indices[0]]
        for i in indices[1:]:
            target.update({k:v for k,v in unique[i].items() if v and not target.get(k)})
            removed.add(i)
    result, seen = [], set()
    for i, values in enumerate(unique):
        if i in removed:
            continue
        key = tuple(sorted((k, v if k == '_facts' else _normal(v))
                           for k, v in values.items() if k != '_identity'))
        if key not in seen:
            seen.add(key)
            result.append(values)
    return result
