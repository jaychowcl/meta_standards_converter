# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Bounded human presentation; structured evidence is never reduced here."""
from collections import OrderedDict
import hashlib
import json


def _priority(diagnostic):
    severity = {'error': 0, 'warning': 1, 'info': 2, 'debug': 3}.get(diagnostic.get('severity'), 1)
    code = diagnostic.get('code')
    detail = 2 if code == 'source_partial' else 1 if code == 'source_fallback' else 0
    return severity, detail


class _Groups:
    def __init__(self, verbose=False, detailed=False):
        self.verbose = verbose
        self.detailed = detailed
        self.seen = set()
        self.groups = OrderedDict()
        self.counts = {}
        self.curation = OrderedDict()

    def add(self, diagnostic):
        d = dict(diagnostic)
        if d.get('curation'):
            if self.verbose:
                c = d['curation']
                key = (c['dataset_id'], c['field'], c['reason'])
                locations = self.curation.setdefault(key, OrderedDict())
                occurrence = (c.get('package_index', 0), c['entity_id'], d.get('location'))
                locations.setdefault(occurrence, d)
            return
        # A shared cause may have several consequences. Uncorrelated diagnostics
        # are grouped only when all their affected context and text agree.
        key = d.get('correlation_id') or tuple(d.get(k) for k in (
            'code', 'provider', 'operation', 'stage', 'location', 'message'))
        if not d.get('correlation_id'):
            key = hashlib.sha256(json.dumps(key, ensure_ascii=True).encode()).digest()
        self.seen.add(key)
        if key not in self.groups:
            if not self.verbose and len(self.groups) >= 5:
                worst = max(self.groups, key=lambda k: min(_priority(v) for v in self.groups[k]))
                if _priority(d) >= min(_priority(v) for v in self.groups[worst]):
                    return
                del self.groups[worst]
                self.counts.pop(worst, None)
            self.groups[key] = []
        self.counts[key] = self.counts.get(key, 0) + 1
        group = self.groups[key]
        if not any((v.get('code'), v.get('message')) == (d.get('code'), d.get('message')) for v in group):
            group.append(d)
            group.sort(key=_priority)
            if not self.verbose:
                # Keep the actionable cause and the fallback consequence.
                fallback = next((v for v in group if v.get('code') == 'source_fallback'), None)
                group[:] = [group[0]] + ([fallback] if fallback is not None and fallback is not group[0] else [])

    def lines(self):
        for key, group in sorted(self.groups.items(), key=lambda item: _priority(item[1][0])):
            for d in group:
                context = ' '.join(str(d[k]) for k in ('provider', 'location') if d.get(k))
                yield f"  {d.get('stage') or 'conversion'} [{d.get('code')}]{' '+context if context else ''}: {d.get('message')}"
            action = next((d.get('action') for d in reversed(group) if d.get('action')), None)
            if action and not any(action in d.get('message', '') for d in group):
                yield f'    Action: {action}'
            if group[0].get('severity') == 'error' and group[0].get('correlation_id'):
                yield f"    Correlation: {group[0]['correlation_id']}"
            if self.counts[key] > 1 and self.verbose:
                yield f'    {self.counts[key]} occurrences (full evidence in log/report)'
        extra = len(self.seen) - len(self.groups)
        if extra:
            yield f'  {extra} additional diagnostic groups; use -v or a DEBUG log for details.'
        for (dataset, field, reason), locations in self.curation.items():
            values = list(locations.values())
            first = values[0]
            types = {d['curation']['entity_type'] for d in values}
            kind = next(iter(types)) if len(types) == 1 else 'record'
            entities = {(d['curation'].get('package_index', 0), d['curation']['entity_id']) for d in values}
            count = len(entities)
            noun = {'study': 'studies'}.get(kind, kind + 's') if count != 1 else kind
            label = field.split('.', 1)[-1].replace('_', ' ')
            condition = {'missing': 'missing', 'reported_unknown': 'reported as unknown',
                         'lookup_no_match': 'not resolved by lookup',
                         'lookup_incomplete': 'incomplete in the lookup response'}.get(reason, reason.replace('_', ' '))
            yield f"  WARNING [{first['code']}] {dataset}: {label} {condition} in {count} {noun}."
            lookup_states = sorted({(d['curation'].get('lookup_provider') or '', d['curation']['lookup_status'])
                                    for d in values if d['curation'].get('lookup_status')})
            if lookup_states:
                yield '    Lookup: ' + '; '.join(f'{provider} {status}'.strip().replace('_', ' ')
                                                for provider, status in lookup_states) + '.'
            chosen = values if self.detailed else values[:3]
            heading = 'Locations' if self.detailed else 'Examples'
            examples = [f"{d['curation']['entity_id']} (package {d['curation'].get('package_index', 0)} {d['location']})" for d in chosen]
            yield f'    {heading}: ' + '; '.join(examples)
            yield f"    Action: {first['action']}"
            if not self.detailed and len(values) > 3:
                yield '    Use -vv or --report-json for all affected locations.'


def _item_lines(source, identity, status, artifacts):
    visible = {k: v for k, v in artifacts.items() if 'bundle_pointer' not in k}
    effect = ('partial output written' if status == 'partial' else
              'failure with retained artifacts' if status in {'failed', 'interrupted'} else 'output written') if visible else 'no output produced'
    label = f'{identity}: {source}' if identity and source and identity != source else source or identity or 'input'
    yield f'{label}: {status}; {effect}'
    for name, path in visible.items():
        yield f'  {name}: {path}'
