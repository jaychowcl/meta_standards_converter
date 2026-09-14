"""Select scoped IDF contacts and consolidate presentation facts, not identities."""
from copy import deepcopy
import re


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


def _compatible(a, b):
    return all(not a.get(k) or not b.get(k) or
               (a[k] == b[k] if k == '_facts' else _normal(a[k]) == _normal(b[k]))
               for k in a.keys() | b.keys())


def unique_people(columns):
    """Exact facts first; incomplete identities cannot bridge conflicting people."""
    unique = []
    seen = set()
    for values in columns:
        # This explicit repository contact is supplied as a named group by SRA.
        # Do not classify arbitrary unknown names or a generic 'curator' role.
        if tuple(_normal(values.get(k)) for k in ('First Name', 'Last Name', 'Email', 'Affiliation')) == (
                'geo', 'curators', 'geo-group@ncbi.nlm.nih.gov', 'ncbi'):
            continue
        key = tuple((k, v if k == '_facts' else _normal(v)) for k, v in values.items())
        if key not in seen:
            seen.add(key)
            unique.append(deepcopy(values))
    groups = {}
    for i, values in enumerate(unique):
        email, first, last = (_normal(values.get(k)) for k in ('Email', 'First Name', 'Last Name'))
        if not first or not last or not re.fullmatch(r'[^\s@]+@[^\s@]+\.[^\s@]+', email):
            continue
        if email.split('@')[0] in {'info', 'contact', 'support', 'admin', 'office', 'help', 'team'}:
            continue
        groups.setdefault((email, first, last), []).append(i)
    removed = set()
    for indices in groups.values():
        if not all(_compatible(unique[i], unique[j]) for i in indices for j in indices):
            continue
        target = unique[indices[0]]
        for i in indices[1:]:
            target.update({k: v for k, v in unique[i].items() if v and not target.get(k)})
            removed.add(i)
    return [v for i, v in enumerate(unique) if i not in removed]
