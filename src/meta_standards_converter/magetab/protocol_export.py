# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Protocol identity allocation on the operation-local MAGE-TAB export copy."""
from copy import deepcopy
from ._performance import measured_phase, work
import logging
import re

logger = logging.getLogger(__name__)


from ..miniml.protocol_text import comparable_protocol_text as _text


def _value(value, description=False):
    if isinstance(value, str):
        return _text(value) if description else " ".join(value.split())
    if isinstance(value, list):
        return [_value(v) for v in value]
    if isinstance(value, dict):
        return {k: _value(v) for k, v in value.items() if v not in (None, '', [], {})}
    return value


def _compatible(a, b):
    if not a.get('description') or not b.get('description'):
        return False
    return all((_type_compatible(a[k], b[k]) if k == 'type' else
                _value(a[k], k == 'description') == _value(b[k], k == 'description'))
               for k in a.keys() & b.keys() if k != 'name'
               and a[k] not in (None, '', [], {}) and b[k] not in (None, '', [], {}))


def _type_compatible(a, b):
    aliases = {'nucleic acid library construction protocol': 'library construction protocol',
               'sample treatment protocol': 'treatment protocol'}
    a, b = (deepcopy(v) if isinstance(v, dict) else {'value': v} for v in (a, b))
    for value in (a, b):
        label = ' '.join(str(value.get('value', '')).split()).casefold()
        value['value'] = aliases.get(label, label)
    return all(_value(a[k]) == _value(b[k]) for k in a.keys() & b.keys()
               if a[k] not in (None, '') and b[k] not in (None, ''))


def _registered(name):
    match = re.search(r'(?:^|:)(P-MTAB-\d+)$', name)
    return match[1] if match else None


def _prepared(source, *, verbatim=False):
    """Normalize one definition without changing its authored representation."""
    fields = {}
    for key, value in source.items():
        if key == 'name' or value in (None, '', [], {}):
            continue
        if key == 'type':
            aliases = {'nucleic acid library construction protocol': 'library construction protocol',
                       'sample treatment protocol': 'treatment protocol'}
            values = dict(value) if isinstance(value, dict) else {'value': value}
            label = ' '.join(str(values.get('value', '')).split()).casefold()
            values['value'] = aliases.get(label, label)
            fields[key] = {k: _value(v) for k, v in values.items() if v not in (None, '')}
        elif key == 'description' and verbatim:
            fields[key] = value
        else:
            fields[key] = _value(value, key == 'description')
    return bool(source.get('description')), fields


def _prepared_compatible(a, b):
    if not a[0] or not b[0]:
        return False
    left, right = a[1], b[1]
    for key in left.keys() & right.keys():
        if key == 'type':
            if any(left[key][k] != right[key][k] for k in left[key].keys() & right[key].keys()):
                return False
        elif left[key] != right[key]:
            return False
    return True


def _comparison_key(value):
    """Hash equality-preserving JSON comparison values, including sparse fields."""
    if isinstance(value, dict):
        return ('mapping', frozenset((k, _comparison_key(v)) for k, v in value.items()))
    if isinstance(value, list):
        return ('list', tuple(_comparison_key(v) for v in value))
    if isinstance(value, tuple):
        return ('tuple', tuple(_comparison_key(v) for v in value))
    return ('value', value)


def _neighborhoods(protocols, prepared):
    # Exact comparison-equivalent definitions share one neighborhood. Names and
    # registered identities still belong to the original ordered occurrences.
    classes = {}
    membership = []
    buckets = {}
    for position, item in enumerate(prepared):
        key = _comparison_key(item)
        if key not in classes:
            classes[key] = {'value': item, 'positions': [], 'neighbors': frozenset()}
            if item[0]:
                description = _comparison_key(item[1].get('description'))
                buckets.setdefault(description, []).append(key)
        classes[key]['positions'].append(position)
        membership.append(key)
    comparisons = 0
    for bucket in buckets.values():
        for key in bucket:
            members = []
            for other in bucket:
                comparisons += 1
                if _prepared_compatible(classes[key]['value'], classes[other]['value']):
                    members.extend(classes[other]['positions'])
            classes[key]['neighbors'] = frozenset(members)
    registered = [_registered(p['name']) for p in protocols]
    for definition in classes.values():
        definition['registered_choices'] = frozenset(
            registered[j] for j in definition['neighbors'] if registered[j])
    work(protocol_candidate_pairs=comparisons, protocol_definition_classes=len(classes))
    return (
        [classes[key]['neighbors'] for key in membership],
        [classes[key]['registered_choices'] for key in membership],
    )


@measured_phase("magetab_protocol_matching")
def prepare_protocols(data):
    """Retain all definitions, coalesce compatible equivalents and rewrite references."""
    from meta_standards_converter.miniml.archive_protocols import prepare_native_protocols
    prepare_native_protocols(data)
    series = data.get('series', {})
    protocols = series.get('protocols', [])
    work(protocol_definitions=len(protocols))
    if not protocols:
        return
    reserved = {_registered(p['name']) or p['name'] for p in protocols
                if _registered(p['name']) or p['name'].startswith('P-')}
    groups = []
    groups_by_neighbors = {}
    bindings = {}
    # Compatibility is not transitive: an incomplete record must not choose
    # arbitrarily between conflicting richer definitions.
    from .text import geo_origin
    verbatim = geo_origin(data)
    prepared = [_prepared(source, verbatim=verbatim) for source in protocols]
    neighbors, registered_choices = _neighborhoods(protocols, prepared)
    # Registered definitions have naming priority, regardless of source order.
    for position, source in sorted(enumerate(protocols), key=lambda item: not bool(_registered(item[1]['name']))):
        registered = _registered(source['name'])
        group = next((g for g in groups_by_neighbors.get(neighbors[position], [])
                      if not (registered and g['registered'] and registered != g['registered'])
                      and not (g['registered'] and not registered and len(registered_choices[position]) > 1)
                      and _prepared_compatible(g['prepared'], prepared[position])), None)
        if group is None:
            group = {'definition': deepcopy(source), 'registered': registered, 'position': position,
                     'prepared': prepared[position]}
            groups.append(group)
            groups_by_neighbors.setdefault(neighbors[position], []).append(group)
        else:
            group['position'] = min(group['position'], position)
            changed = False
            for key, value in source.items():
                if key != 'name' and group['definition'].get(key) in (None, '', [], {}):
                    group['definition'][key] = deepcopy(value)
                    changed = True
                elif key == 'type' and isinstance(value, dict) and isinstance(group['definition'].get(key), dict):
                    additions = {k: deepcopy(v) for k, v in value.items()
                                 if v and not group['definition'][key].get(k)}
                    group['definition'][key].update(additions)
                    changed = changed or bool(additions)
            if changed:
                group['prepared'] = _prepared(group['definition'], verbatim=verbatim)
        bindings[source['name']] = group
    used = set()
    counter = 1
    for group in sorted(groups, key=lambda g: g['position']):
        old = group['definition']['name']
        name = group['registered'] or (old if old.startswith('P-') else None)
        if name in used:
            logger.warning('Conflicting protocol definitions for %s; allocating a separate local identifier', name)
            name = None
        if name is None:
            while f"P-{series.get('iid') or 'study'}-{counter}" in reserved | used:
                counter += 1
            name = f"P-{series.get('iid') or 'study'}-{counter}"
            counter += 1
        group['definition']['name'] = name
        used.add(name)
    series['protocols'] = [g['definition'] for g in sorted(groups, key=lambda g: g['position'])]
    for path in series.get('assay_paths', []):
        for step in path.get('steps', []):
            if step.get('protocol_ref') in bindings:
                step['protocol_ref'] = bindings[step['protocol_ref']]['definition']['name']
