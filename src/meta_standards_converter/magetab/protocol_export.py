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


def prepare_protocols(data):
    """Retain all definitions, coalesce compatible equivalents and rewrite references."""
    from meta_standards_converter.miniml.archive_protocols import prepare_native_protocols
    prepare_native_protocols(data)
    series = data.get('series', {})
    protocols = series.get('protocols', [])
    if not protocols:
        return
    reserved = {_registered(p['name']) or p['name'] for p in protocols
                if _registered(p['name']) or p['name'].startswith('P-')}
    groups = []
    bindings = {}
    # Compatibility is not transitive: an incomplete record must not choose
    # arbitrarily between conflicting richer definitions.
    neighbors = [{j for j, other in enumerate(protocols) if _compatible(source, other)}
                 for source in protocols]
    registered_choices = [{_registered(protocols[j]['name']) for j in neighborhood
                           if _registered(protocols[j]['name'])} for neighborhood in neighbors]
    # Registered definitions have naming priority, regardless of source order.
    for position, source in sorted(enumerate(protocols), key=lambda item: not bool(_registered(item[1]['name']))):
        registered = _registered(source['name'])
        group = next((g for g in groups if not (registered and g['registered'] and registered != g['registered'])
                      and not (g['registered'] and not registered and len(registered_choices[position]) > 1)
                      and neighbors[position] == g['neighbors'] and _compatible(g['definition'], source)), None)
        if group is None:
            group = {'definition': deepcopy(source), 'registered': registered, 'position': position,
                     'neighbors': neighbors[position]}
            groups.append(group)
        else:
            group['position'] = min(group['position'], position)
            for key, value in source.items():
                if key != 'name' and group['definition'].get(key) in (None, '', [], {}):
                    group['definition'][key] = deepcopy(value)
                elif key == 'type' and isinstance(value, dict) and isinstance(group['definition'].get(key), dict):
                    group['definition'][key].update({k: deepcopy(v) for k, v in value.items()
                                                     if v and not group['definition'][key].get(k)})
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
