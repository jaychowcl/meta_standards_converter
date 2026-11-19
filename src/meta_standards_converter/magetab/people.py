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
from itertools import product
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


def contact_facts(contributor, organization):
    """Keep source components on the private export copy until consolidation."""
    return deepcopy({
        'organization': contributor.get('organization') or organization.get('name'),
        **{key: contributor.get(key) for key in ('laboratory', 'department', 'company')},
        'address': contributor.get('address') or organization.get('address'),
    })


def render_contact(facts):
    from ..helpers.json_helper import JSONHandler
    affiliation = ', '.join(dict.fromkeys(str(facts[key]) for key in (
        'laboratory', 'department', 'organization', 'company') if facts.get(key))) or None
    address = facts.get('address')
    parts = JSONHandler()._flatten_values(address) if address else []
    if address and not isinstance(address, str):
        parts = [affiliation, *parts]
    return {'Affiliation': affiliation, 'Address': ', '.join(str(x) for x in parts if x) or None}


def _structured_compatible(a, b):
    if not a or not b:
        return True
    if isinstance(a, dict) and isinstance(b, dict):
        return all(_structured_compatible(a.get(k), b.get(k)) for k in a.keys() | b.keys())
    if isinstance(a, (list, tuple)) and isinstance(b, (list, tuple)):
        return len(a) == len(b) and all(_structured_compatible(x, y) for x, y in zip(a, b))
    return type(a) is type(b) and _normal(a) == _normal(b)


def _opaque_address_matches(opaque, structured):
    """Exact render variants only: retain every street/locality/country component."""
    address = structured['address']
    if not isinstance(opaque, str) or not address.get('line') or set(address) - {
            'line', 'city', 'state', 'province', 'zip_code', 'postal_code', 'country'}:
        return False
    prefixes = [key for key in ('laboratory', 'department', 'organization', 'company') if structured.get(key)]
    postcodes = [key for key in ('zip_code', 'postal_code') if address.get(key)]
    optional = [(False, key) for key in prefixes] + [(True, key) for key in postcodes]
    for omissions in product((False, True), repeat=len(optional)):
        candidate = deepcopy(structured)
        for omit, (postal, key) in zip(omissions, optional):
            if omit:
                (candidate['address'] if postal else candidate).pop(key)
        if _normal(opaque) == _normal(render_contact(candidate)['Address']):
            return True
    return False


def _contacts_compatible(a, b):
    if not all(_structured_compatible(a.get(k), b.get(k)) for k in (
            'organization', 'laboratory', 'department', 'company')):
        return False
    left, right = a.get('address'), b.get('address')
    if not left or not right or _structured_compatible(left, right):
        return True
    if isinstance(left, dict) and isinstance(right, str):
        return _opaque_address_matches(right, a)
    if isinstance(right, dict) and isinstance(left, str):
        return _opaque_address_matches(left, b)
    return False


def _complete(a, b):
    result = deepcopy(a)
    for key, value in b.items():
        if not result.get(key):
            result[key] = deepcopy(value)
        elif isinstance(result[key], dict) and isinstance(value, dict):
            result[key] = _complete(result[key], value)
    return result


def _complete_contact(a, b):
    result = _complete(a, b)
    if isinstance(b.get('address'), dict) and not isinstance(a.get('address'), dict):
        result['address'] = deepcopy(b['address'])
    return result


def _persistent_ids(values):
    result = {}
    for token in values.get('_identity', ()):
        if len(token) == 3 and token[0] == 'accession':
            result.setdefault(_normal(token[1]), set()).add(_normal(token[2]))
    return result


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
    left, right = _persistent_ids(a), _persistent_ids(b)
    if any(len(ids) > 1 for ids in [*left.values(), *right.values()]):
        return False
    if any(left[key] != right[key] for key in left.keys() & right.keys()):
        return False
    structured = '_contact' in a and '_contact' in b
    if structured and not _contacts_compatible(a['_contact'], b['_contact']):
        return False
    for key in a.keys() | b.keys():
        if key in {'_identity', '_contact'} or (structured and key in {'Affiliation', 'Address'}):
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
            incoming = unique[i]
            if '_contact' in target and '_contact' in incoming:
                target['_contact'] = _complete_contact(target['_contact'], incoming['_contact'])
                target.update(render_contact(target['_contact']))
            target.update({k:v for k,v in incoming.items() if v and not target.get(k)})
            target['_identity'] = set(target.get('_identity', ())) | set(incoming.get('_identity', ()))
            removed.add(i)
    result, seen = [], {}
    for i, values in enumerate(unique):
        if i in removed:
            continue
        persistent = tuple(sorted((kind, tuple(sorted(ids))) for kind, ids in _persistent_ids(values).items()))
        key = (persistent, tuple(sorted((k, v if k == '_facts' else _normal(v))
                           for k, v in values.items() if k not in {'_identity', '_contact'})))
        previous = seen.setdefault(key, [])
        if not previous or not all(_compatible(other, values) for other in previous):
            result.append(values)
        previous.append(values)
    return result
