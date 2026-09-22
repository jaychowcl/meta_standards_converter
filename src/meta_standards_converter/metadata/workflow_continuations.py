# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Absorb only proven generated acquisition prefixes into authored continuations."""
from copy import deepcopy
import re


def _node_identity(step):
    if step['kind'].endswith('_file'):
        return step['kind'], step.get('link', {}).get('value'), step.get('name')
    return tuple(step.get(k) for k in ('kind', 'name', 'sample_ref', 'protocol_ref'))


def _complete(old, new):
    """Complete new evidence without overwriting any incompatible assertion."""
    if isinstance(old, dict) and isinstance(new, dict):
        for key, value in old.items():
            if key not in new:
                new[key] = deepcopy(value)
            elif key == 'type' and isinstance(value, str) and isinstance(new[key], str) and value.casefold() == new[key].casefold():
                continue
            elif not _complete(value, new[key]):
                return False
        return True
    if isinstance(old, list) and isinstance(new, list):
        for value in old:
            if value in new:
                continue
            if isinstance(value, dict) and value.get('name'):
                matches = [v for v in new if isinstance(v, dict) and v.get('name') == value['name']]
                if matches:
                    if len(matches) != 1 or not _complete(value, matches[0]):
                        return False
                    continue
            new.append(deepcopy(value))
        return True
    return old == new


def _bound_acquisition(steps):
    sources = {s.get('sample_ref') for s in steps if s.get('sample_ref')}
    assays = [s.get('name', '') for s in steps if s['kind'] == 'assay']
    runs = [s.get('name', '') for s in steps if s['kind'] == 'scan']
    return (len(sources) == len(assays) == len(runs) == 1
            and re.fullmatch(r'[SED]RX\d+', assays[0])
            and re.fullmatch(r'[SED]RR\d+', runs[0]))


def _reconcile_continuations(paths, generated, issues):
    # Complementary facts on an identical acquisition belong to every file
    # alternative of that acquisition. Result/factor nodes are excluded.
    acquisitions = {}
    for path in paths:
        boundary = next((i for i, node in enumerate(path['steps']) if node['kind'].endswith('_file')), None)
        if boundary is not None and _bound_acquisition(path['steps'][:boundary]):
            prefix = path['steps'][:boundary]
            acquisitions.setdefault(tuple(_node_identity(s) for s in prefix), []).append((path, boundary))
    for group in acquisitions.values():
        if not any(id(p) in generated for p, _ in group):
            continue
        combined = deepcopy(group[0][0]['steps'][:group[0][1]])
        if all(all(_complete(old, new) for old, new in zip(path['steps'][:end], combined)) for path, end in group):
            for path, end in group:
                path['steps'][:end] = deepcopy(combined)
        else:
            message = 'enrichment workflow: unresolved annotation scope on shared acquisition'
            if message not in issues: issues.append(message)
    # Index every raw-file prefix once. Terminal result branches are never
    # candidates for removal, and authored occurrences remain distinct.
    index = {}
    for path in paths:
        identity = []
        for position, node in enumerate(path['steps']):
            identity.append(_node_identity(node))
            if node['kind'] == 'array_data_file' and position < len(path['steps']) - 1:
                index.setdefault(tuple(identity), []).append(path)
    removed = set()
    for path in paths:
        if (id(path) not in generated or path['steps'][-1]['kind'] != 'array_data_file'
                or not _bound_acquisition(path['steps'])):
            continue
        key = tuple(_node_identity(s) for s in path['steps'])
        choices = index.get(key, [])
        if not choices:
            continue
        prepared = []
        for choice in choices:
            prefix = deepcopy(choice['steps'][:len(key)])
            metadata = {k: deepcopy(v) for k, v in choice.items() if k != 'steps'}
            if (all(_complete(old, new) for old, new in zip(path['steps'], prefix))
                    and _complete({k: v for k, v in path.items() if k != 'steps'}, metadata)):
                prepared.append((choice, prefix, metadata))
        if len(prepared) != len(choices):
            message = 'enrichment workflow: unresolved annotation scope or conflicting evidence on acquisition continuation'
            if message not in issues: issues.append(message)
            continue
        for choice, prefix, metadata in prepared:
            choice.update(metadata)
            choice['steps'][:len(key)] = prefix
        removed.add(id(path))
    paths[:] = [path for path in paths if id(path) not in removed]
