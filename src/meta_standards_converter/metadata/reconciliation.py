# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Source value selection is independent of archive identity and membership joins."""
from copy import deepcopy

_RANK = {'ArrayExpress': 3, 'GEO': 2, 'SRA': 1, 'ENA': 1, 'DRA': 1}


def provider(document):
    value = str(document.get('source', {}).get('format', '')).casefold()
    return {'mage-tab': 'ArrayExpress', 'magetab': 'ArrayExpress', 'arrayexpress': 'ArrayExpress',
            'biostudies': 'ArrayExpress', 'ae': 'ArrayExpress', 'geo miniml': 'GEO', 'geo': 'GEO',
            'sra': 'SRA', 'ena': 'ENA', 'dra': 'DRA'}.get(value, value)



def _meaning(value):
    if isinstance(value, str): return ' '.join(value.split()).casefold()
    if isinstance(value, dict):
        if 'value' in value: return _meaning(value['value'])
        return tuple(sorted((k, _meaning(v)) for k, v in value.items()))
    if isinstance(value, (tuple, list)): return tuple(_meaning(v) for v in value)
    return value


class SourceSelection:
    def __init__(self, data, extra, issues, *, mode=None):
        self.mode = mode
        self.current = {'provider': provider(data), 'accession': data['series']['iid']}
        self.geo_primary = self.current['provider'] == 'GEO'
        self.incoming = {'provider': provider(extra), 'accession': extra['series']['iid']}
        self.issues, self.records, self.selections = issues, [], {}
        for record in data.get('extensions', {}).get('insdc', {}).get('records', []):
            if record.get('kind') == 'source_selection':
                self.selections.update(deepcopy(record['metadata'].get('selections', {})))

    def prefer(self, old, new, path, default, *, scientific=False):
        if self.mode not in {'standard', 'curators'} and not self.geo_primary:
            return default
        from .archive_enrichment import informative
        previous = self.selections.get(path, {})
        origin = previous if previous.get('value') == old else self.current
        left = {**origin, 'value': deepcopy(old)}
        right = {**self.incoming, 'value': deepcopy(new)}
        if self.mode == 'standard':
            from .source_text import _text_choice
            text_choice = _text_choice(left, right, path)
            if text_choice is not None:
                chosen, record = text_choice
                if self.geo_primary:
                    record = {**record, 'kind': 'text_normalization_suggestion'}
                if record not in self.records: self.records.append(record)
                if not self.geo_primary:
                    self.selections[path] = right if chosen else left
                    return chosen
        supplied = old not in (None, '', [], {}) if self.geo_primary else informative(old)
        chosen = informative(new) and (not supplied or (
            self.mode == 'standard' and not self.geo_primary
            and _RANK.get(right['provider'], 0) > _RANK.get(left['provider'], 0)))
        selected = right if chosen else left
        reason = 'primary input value retained' if self.mode == 'curators' else 'AE > GEO > native; retain current value for equal source rank'
        if self.geo_primary: reason = 'GEO primary input value retained; enrichment is additive'
        if not supplied and informative(new): reason = 'missing value completed'
        if informative(selected['value']): self.selections[path] = selected
        different = old != new if self.geo_primary else _meaning(old) != _meaning(new)
        if (scientific or self.geo_primary) and supplied and informative(new) and different:
            metadata = {'field': path.rsplit('/', 1)[-1], 'path': path,
                        'statements': [left, right], 'selected': selected, 'reason': reason, 'mode': self.mode}
            record = {'provider': right['provider'], 'kind': 'source_conflict', 'accession': right['accession'], 'metadata': metadata}
            if record not in self.records: self.records.append(record)
            message = (f"source conflict at {path}: {left['provider']} {left['accession']} {old!r} versus "
                       f"{right['provider']} {right['accession']} {new!r}; selected {selected['provider']} "
                       f"{selected['value']!r} ({reason})")
            if message not in self.issues: self.issues.append(message)
        return chosen

    def prefer_workflow(self, sample, identifiers):
        path = f"/sample/{sample}/workflow/{','.join(identifiers)}/source"
        previous = self.selections.get(path, {}).get('value', self.current['provider'])
        return self.prefer(previous, self.incoming['provider'], path, True)

    def evidence(self):
        if self.mode not in {'standard', 'curators'}: return []
        return self.records + [{'provider': self.current['provider'], 'kind': 'source_selection',
            'accession': self.current['accession'], 'metadata': {'selections': self.selections}}]


def inspect_recorded_conflicts(package, mode):
    """Inspect supplied alternatives only; this operation never retrieves metadata."""
    from ..miniml import MINiMLCodec
    from ..miniml.archive_libraries import FIELDS, path_facts
    from .archive_workflows import path_ids, compatible
    import re
    from .archive_enrichment import entity_ids
    data = package.to_mapping()
    records = data.get('extensions', {}).get('insdc', {}).get('records', [])
    issues, additions = [], []
    for record in list(records):
        if record.get('kind') == 'source_conflict':
            fact = record['metadata']
            issues.append(f"recorded source conflict at {fact['path']}: {fact['statements']!r}; selected {fact['selected']!r}")
        if mode != 'curators' or record.get('kind') != 'MINiML' or not isinstance(record.get('metadata'), dict): continue
        extra = record['metadata']
        if not extra.get('sample'): continue
        extra = {**extra, 'source': extra.get('source') or {'format': record['provider']},
                 'series': {**extra.get('series', {}), 'iid': record.get('accession') or extra.get('series', {}).get('iid', '')}}
        selection = SourceSelection(data, extra, issues, mode='curators')
        bound_facts = {}
        for assay_path in extra.get('series', {}).get('assay_paths', []):
            refs = {step.get('sample_ref') for step in assay_path.get('steps', [])} - {None}
            scope = path_ids(assay_path)
            # Only explicit same-sample, run/library bindings support comparison.
            if len(refs) != 1 or not scope:
                continue
            identifiers = {v for v in scope if re.fullmatch(r'[SED]RR\d+', v)} or scope
            for identifier in identifiers:
                bound_facts.setdefault((next(iter(refs)), identifier), []).append(
                    (scope, path_facts(assay_path['steps'])))
        for sample in data.get('sample', []):
            matches = [s for s in extra['sample'] if entity_ids(sample, sample=True) & entity_ids(s, sample=True)]
            if len(matches) != 1: continue
            reverse = [s for s in data['sample'] if entity_ids(s, sample=True) & entity_ids(matches[0], sample=True)]
            if len(reverse) != 1: continue
            channels, other_channels = sample.get('channel', []), matches[0].get('channel', [])
            if len(channels) == len(other_channels) == 1:
                for field in ('organism', 'molecule'):
                    if field in other_channels[0]:
                        selection.prefer(channels[0].get(field), other_channels[0][field],
                                         f"/sample/{sample['iid']}/channel/0/{field}", False, scientific=True)
            incoming = {r.get('run'): r for r in matches[0].get('sra_run', []) if r.get('run')}
            for run in sample.get('sra_run', []):
                other = incoming.get(run.get('run'), {})
                if other.get('experiment') and run.get('experiment') and other['experiment'] != run['experiment']: continue
                supplied = [other]
                native_scope = {run.get('run'), run.get('experiment')} - {None}
                for identifier in dict.fromkeys((run.get('run'), run.get('experiment'))):
                    if identifier is None:
                        continue
                    for scope, facts in bound_facts.get((matches[0]['iid'], identifier), []):
                        if compatible(scope, native_scope):
                            supplied.extend(facts)
                for field in sorted(FIELDS):
                    seen = []
                    for statement in supplied:
                        if field in statement and statement[field] not in seen:
                            seen.append(statement[field])
                            selection.prefer(run.get(field), statement[field], f"/sample/{sample['iid']}/sra_run/{run['run']}/{field}", False, scientific=True)
        additions.extend(selection.records)
    if additions:
        records.extend(r for r in additions if r not in records)
    return MINiMLCodec().decode(data, strict=True).package, list(dict.fromkeys(issues))
