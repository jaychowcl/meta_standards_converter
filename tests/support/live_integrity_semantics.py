# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Independent per-row SDRF oracle for sparse protocol-column presentation.

No MSC imports. Ordered nonempty applications, duplicate rows/attributes and
qualifier owners remain significant; physical empty column slots do not.
"""
from collections import Counter, defaultdict, deque
import json
import re

from tests.support.deposited_fidelity import HEADER_ALIASES, NODE_HEADERS, QUALIFIER_HEADERS


def header(value):
    value = HEADER_ALIASES.get(value, value)
    match = re.fullmatch(r'(Characteristics|Comment|Parameter\s+Value|Factor\s*Value|Unit)\s*\[(.*)\]', value)
    if match:
        prefix = match[1]
        if prefix.replace(' ', '') == 'FactorValue':
            prefix = 'Factor Value'
        return prefix + '[' + match[2] + ']'
    return value


def nonempty(value):
    return value is not None and (not isinstance(value, str) or bool(value.strip()))


def stable(value):
    return json.dumps(value, sort_keys=True, ensure_ascii=False)


def normalize(table):
    return [[header(v) for v in table[0]], *[list(row) for row in table[1:]]]


def organism_taxids(headers, row):
    """Read only identifiers attached to an organism characteristic occurrence."""
    result = set()
    for i, label in enumerate(headers):
        if header(label).casefold() != 'characteristics[organism]': continue
        qualifiers = []
        for h, value in zip(headers[i+1:], row[i+1:]):
            h = header(h)
            if h not in QUALIFIER_HEADERS and not h.startswith('Comment['): break
            qualifiers.append((h.casefold(), value))
        for h, value in qualifiers:
            if h == 'term accession number':
                match = re.search(r'NCBITaxon[:_/#](\d+)$', str(value), re.I)
                if match: result.add(match[1])
                elif str(value).isdigit() and ('term source ref','NCBITaxon') in qualifiers:
                    result.add(str(value))
            elif h in ('comment[taxid]', 'comment[ncbi taxon id]') and str(value).isdigit():
                result.add(str(value))
    return result


def row_graph(headers, row):
    assert len(headers) == len(row)
    graph, factors = [], []
    node, owner = None, None
    for label, value in zip(headers, row):
        label = header(label)
        value = value if nonempty(value) else ''
        if label in NODE_HEADERS:
            node = {'kind': label, 'value': value, 'attributes': [], 'qualifiers': []}
            graph.append(node); owner = node
        elif label.startswith('Factor Value['):
            owner = {'kind': label, 'value': value, 'qualifiers': []}
            factors.append(owner)
        elif label.startswith('Comment[') or label in QUALIFIER_HEADERS:
            assert owner is not None, f'orphan qualifier: {label}'
            owner['qualifiers'].append({'kind': label, 'value': value, 'qualifiers': []})
        elif label == 'Unit' or label.startswith('Unit['):
            assert owner is not None, 'orphan unit'
            unit = {'kind': label, 'value': value, 'qualifiers': []}
            owner['qualifiers'].append(unit); owner = unit
        else:
            assert node is not None, f'attribute without node: {label}'
            owner = {'kind': label, 'value': value, 'qualifiers': []}
            node['attributes'].append(owner)
    def clean(item):
        qualifiers = [x for q in item['qualifiers'] if (x := clean(q)) is not None]
        attributes = [x for a in item.get('attributes', []) if (x := clean(a)) is not None]
        if not nonempty(item['value']) and not qualifiers and not attributes:
            return None
        result = {'kind': item['kind'], 'value': item['value'],
            'qualifiers': sorted(qualifiers, key=lambda q: q['kind'])}
        if 'attributes' in item:
            result['attributes'] = sorted(attributes, key=lambda a: a['kind'])
        return result
    return {'graph': [x for n in graph if (x := clean(n)) is not None],
        'factors': sorted([x for f in factors if (x := clean(f)) is not None], key=lambda f: f['kind'])}


def compare_rows(source, generated):
    """Match entire meaningful row graphs, including multiplicity and ownership."""
    left = [row_graph(source[0], row) for row in source[1:]]
    right = [row_graph(generated[0], row) for row in generated[1:]]
    available = defaultdict(deque)
    for index, value in enumerate(right, 1):
        available[stable(value)].append(index)
    matches, missing = [], []
    for index, value in enumerate(left, 1):
        if available[stable(value)]:
            matches.append({'source_row': index, 'generated_row': available[stable(value)].popleft()})
        else:
            missing.append({'source_row': index, 'graph': value})
    added = [{'generated_row': index, 'graph': right[index-1]} for rows in available.values() for index in rows]
    return {'equivalent': not missing and not added, 'matches': matches,
        'missing': missing, 'added': added, 'source_rows': len(left), 'generated_rows': len(right),
        'meaningful_source_cells': sum(nonempty(v) for row in source[1:] for v in row),
        'whitespace_only_source_cells': sum(isinstance(v, str) and bool(v) and not v.strip() for row in source[1:] for v in row)}
