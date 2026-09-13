"""Accession-bound enrichment supporting native archive converters."""
import re
from urllib.parse import parse_qs, urlsplit


def linked_accessions(package):
    found = set()
    def accept(value):
        if not isinstance(value, str):
            return
        if re.fullmatch(r'GSE\d+|E-[A-Z]+-\d+', value):
            found.add(value)
        if value.startswith(('https://', 'http://')):
            for acc in parse_qs(urlsplit(value).query).get('acc', []):
                if re.fullmatch(r'GSE\d+', acc):
                    found.add(acc)
    data = package.to_mapping() if hasattr(package, 'to_mapping') else package
    for entity in [data['series'], *data.get('sample', [])]:
        for item in entity.get('accession', []):
            accept(item.get('value'))
        for item in entity.get('relation', []):
            accept(item.get('target'))
    def visit(node):
        if not isinstance(node, dict):
            return
        tag = node.get('tag')
        if tag in ('EXTERNAL_ID', 'SECONDARY_ID'):
            accept((node.get('text') or '').strip())
        if tag == 'XREF_LINK':
            fields = {c.get('tag'): (c.get('text') or '').strip() for c in node.get('children', [])}
            if fields.get('DB', '').upper() in ('GEO', 'ARRAYEXPRESS', 'GDS'):
                accept(fields.get('ID'))
        for c in node.get('children', []):
            visit(c)
    for record in data.get('extensions', {}).get('insdc', {}).get('records', []):
        visit(record.get('metadata'))
    return sorted(found)


# These are absence markers for precedence only; the source literals remain retained.
_MISSING = {'', 'na', 'n/a', 'null', 'none', 'unknown', 'not applicable',
            'not available', 'not provided', 'not collected', 'missing', 'unspecified'}
_ID = re.compile(r'(?<![A-Z0-9])(?:[SED]R[PSXR]\d+|SAM(?:N|EA|D)\d+|PRJ(?:NA|EB|DB)\d+|GS[EM]\d+|E-[A-Z]+-\d+)(?![A-Z0-9])')


def informative(value):
    if isinstance(value, str):
        return value.strip().lower() not in _MISSING
    if isinstance(value, dict):
        return informative(value.get('value')) if 'value' in value else any(informative(v) for v in value.values())
    if isinstance(value, list):
        return any(informative(v) for v in value)
    return value is not None


def entity_ids(entity, *, sample=False):
    """Only identifier-bearing fields participate in joins, never descriptions."""
    values = [entity.get('iid', '')] + [a.get('value', '') for a in entity.get('accession', [])]
    values += [a.get('target', '') for a in entity.get('relation', [])]
    if sample:
        values += entity.get('sra_accession', [])
        for run in entity.get('sra_run', []):
            values += [run.get(k, '') for k in ('run', 'experiment', 'sample', 'biosample')]
    result = {acc for value in values for acc in _ID.findall(str(value))}
    if sample:
        # A common study/project does not establish a sample identity.
        result = {a for a in result if not re.match(r'(?:[SED]RP|PRJ|GSE|E-)', a)}
    return result


def _merge_entity(target, extra, prefer, protected=()):
    from copy import deepcopy
    additive = {'raw_data', 'supplementary_data', 'relation', 'accession'}
    for key, value in extra.items():
        if key in protected:
            continue
        if key in additive:
            # Files retain all alternatives. Identity declarations need no duplicates.
            target.setdefault(key, [])
            target[key].extend(deepcopy(v) for v in value if key in {'raw_data', 'supplementary_data'} or v not in target[key])
        elif key == 'characteristics':
            groups = {}
            for item in value:
                groups.setdefault(item['name'], []).append(item)
            current = target.setdefault(key, [])
            for name, group in groups.items():
                old = [v for v in current if v['name'] == name]
                if informative(group) and (prefer or not informative(old)):
                    current[:] = [v for v in current if v['name'] != name] + deepcopy(group)
                elif not old:
                    current.extend(deepcopy(group))
        elif informative(value) and (prefer or not informative(target.get(key))):
            # Replace complete typed values (including units / term references) together.
            target[key] = deepcopy(value)


def merge_archive_metadata(package, other, *, prefer=False, linked_accession=None):
    """Merge an explicitly related package; peer additions require a shared read study."""
    from copy import deepcopy
    from ..miniml import MINiMLCodec
    data, extra = package.to_mapping(), other.to_mapping()
    issues = []
    native_ids, extra_ids = entity_ids(data['series']), entity_ids(extra['series'])
    shared_read = {a for a in native_ids & extra_ids if re.fullmatch(r'[SED]RP\d+', a)}
    linked_match = linked_accession and linked_accession in linked_accessions(package) and linked_accession in extra_ids
    if (not prefer and not shared_read) or (prefer and not (native_ids & extra_ids or linked_match)):
        return package, ['enrichment: unresolved study identity']
    extension = data.setdefault('extensions', {}).setdefault('insdc', {'version': '1.0', 'records': []})
    if 'insdc' in extra.get('extensions', {}) and not prefer:
        extension['records'].extend(deepcopy(extra['extensions']['insdc']['records']))
    else:
        extension['records'].append({'provider': extra.get('source', {}).get('format', 'linked'),
            'kind': 'MINiML', 'accession': extra['series'].get('iid'), 'metadata': deepcopy(extra)})
    _merge_entity(data['series'], extra['series'], prefer,
                  {'iid', 'sample_ref', 'assay_paths', 'protocols', 'contact_ref'})
    native_samples, extra_samples = data.get('sample', []), extra.get('sample', [])
    joins = {i: [j for j, s in enumerate(extra_samples) if entity_ids(n, sample=True) & entity_ids(s, sample=True)]
             for i, n in enumerate(native_samples)}
    reverse = {j: [i for i, js in joins.items() if j in js] for j in range(len(extra_samples))}
    matched = {}
    for i, candidates in joins.items():
        if len(candidates) != 1 or len(reverse[candidates[0]]) != 1:
            if candidates:
                issues.append(f"{native_samples[i]['iid']}: ambiguous enrichment sample identity")
            continue
        j = candidates[0]
        target, source = native_samples[i], extra_samples[j]
        matched[source['iid']] = target['iid']
        _merge_entity(target, source, prefer, {'iid', 'channel', 'channel_count', 'sra_run', 'ena_accession', 'sra_accession', 'contact_ref', 'platform_ref'})
        channels, incoming = target.get('channel', []), source.get('channel', [])
        # A single channel is unambiguous; multi-channel arrays have no generic join.
        if len(channels) == len(incoming) == 1:
            _merge_entity(channels[0], incoming[0], prefer)
        elif incoming:
            issues.append(f"{target['iid']}: ambiguous enrichment channels")
    for j, sample in enumerate(extra_samples):
        if not reverse[j]:
            if prefer:
                issues.append(f"{sample['iid']}: unmatched linked sample retained outside native membership")
            else:
                native_samples.append(deepcopy(sample))
                matched[sample['iid']] = sample['iid']
    # Namespace linked protocol definitions and references as a coupled group.
    namespace = str(extra['series'].get('iid') or 'linked') + ':'
    proto_names = {}
    for protocol in extra['series'].get('protocols', []):
        protocol = deepcopy(protocol)
        old = protocol['name']
        protocol['name'] = namespace + old
        proto_names[old] = protocol['name']
        if protocol not in data['series'].setdefault('protocols', []):
            data['series']['protocols'].append(protocol)
    sample_map = {s['iid']: s for s in native_samples}
    paths = data['series'].setdefault('assay_paths', [])
    # Attach linked sample protocols/factors to native paths without replacing run IDs.
    additions = {}
    for path in extra['series'].get('assay_paths', []):
        refs = {s.get('sample_ref') for s in path.get('steps', []) if s.get('sample_ref')}
        targets = {matched[r] for r in refs if r in matched}
        if len(targets) != 1 or any(r not in matched for r in refs):
            continue
        target = next(iter(targets))
        for step in path['steps']:
            step = deepcopy(step)
            if step.get('protocol_ref') in proto_names:
                step['protocol_ref'] = proto_names[step['protocol_ref']]
            if step.get('sample_ref'):
                step['sample_ref'] = target
            if step.get('kind') == 'protocol_application' or step.get('factor_values'):
                additions.setdefault(target, []).append(step)
        if not prefer and not any(any(s.get('sample_ref') == target for s in p['steps']) for p in paths):
            cp = deepcopy(path)
            for step in cp['steps']:
                if step.get('sample_ref'):
                    step['sample_ref'] = target
                if step.get('protocol_ref') in proto_names:
                    step['protocol_ref'] = proto_names[step['protocol_ref']]
            paths.append(cp)
    for path in paths:
        target = next((s.get('sample_ref') for s in path['steps'] if s.get('sample_ref')), None)
        if target not in sample_map:
            continue
        channels = sample_map[target].get('channel', [])
        if len(channels) == 1:
            for step in path['steps']:
                if step['kind'] == 'source':
                    organisms = [{'name': 'organism', 'value': o['value'], 'term_source_ref': 'NCBITaxon',
                                  'term_accession_number': o.get('taxid')} for o in channels[0].get('organism', [])]
                    step['characteristics'] = deepcopy(channels[0].get('characteristics', [])) + organisms
                    for addition in additions.get(target, []):
                        if addition.get('factor_values'):
                            step['factor_values'] = deepcopy(addition['factor_values'])
            protocols = [s for s in additions.get(target, []) if s.get('kind') == 'protocol_application']
            path['steps'][1:1] = [deepcopy(s) for s in protocols if s not in path['steps']]
    data['series']['sample_ref'] = [{'ref': s['iid']} for s in native_samples]
    for db in extra.get('database', []):
        if db not in data.setdefault('database', []):
            data['database'].append(deepcopy(db))
    return MINiMLCodec().decode(data).package, issues


class LinkedArchiveEnricher:
    def __init__(self, resource_profile='standard', geo_converter=None, ae_converter=None):
        self.profile, self.geo, self.ae = resource_profile, geo_converter, ae_converter

    def enrich(self, package):
        issues = []
        # Stable ascending priority: informative AE values are applied last.
        links = sorted(linked_accessions(package), key=lambda a: (a.startswith('E-'), a))
        for accession in links:
            try:
                if accession.startswith('GSE'):
                    if self.geo is None:
                        from ..converters.geo2json import GEO2JSONConverter
                        self.geo = GEO2JSONConverter(resource_profile=self.profile)
                    packages = self.geo.convert(accession, enrich=False, related_series=False)
                else:
                    if self.ae is None:
                        from ..converters.ae2json import AE2JSONConverter
                        self.ae = AE2JSONConverter(resource_profile=self.profile)
                    packages = self.ae.convert(accession)
                candidates = [p for p in packages if accession in entity_ids(p.to_mapping()['series'])]
                if len(candidates) != 1:
                    issues.append(f'{accession}: ambiguous or unavailable enrichment study')
                    continue
                package, merge_issues = merge_archive_metadata(package, candidates[0], prefer=True, linked_accession=accession)
                issues.extend(merge_issues)
            except Exception as error:
                issues.append(f'{accession}: enrichment unavailable ({type(error).__name__})')
        return package, issues
