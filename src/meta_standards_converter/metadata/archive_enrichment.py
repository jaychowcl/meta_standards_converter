# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
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
_ID = re.compile(r'(?<![A-Z0-9])(?:[SED]R[PSXR]\d+|SAM(?:N|EA|D)\d+|PRJ(?:NA|EB|DB|DA)\d+|GS[EM]\d+|E-[A-Z]+-\d+)(?![A-Z0-9])')


def informative(value):
    if isinstance(value, str):
        return value.strip().lower() not in _MISSING and not value.strip().lower().startswith('missing:')
    if isinstance(value, dict):
        return informative(value.get('value')) if 'value' in value else any(informative(v) for v in value.values())
    if isinstance(value, list):
        return any(informative(v) for v in value)
    return value is not None


def entity_ids(entity, *, sample=False):
    """Only identifier-bearing fields participate in joins, never descriptions."""
    values = [entity.get('iid', '')] + [a.get('value', '') for a in entity.get('accession', [])]
    identity_relations = {'geo', 'arrayexpress', 'biosample', 'bioproject', 'sra', 'ena', 'dra', 'ddbj'}
    identity_relations.update(provider + '-' + kind for provider in ('ena','sra','dra','ddbj')
                              for kind in ('study','sample','experiment','run'))
    values += [a.get('target', '') for a in entity.get('relation', [])
               if str(a.get('type', '')).casefold() in identity_relations]
    if sample:
        values += entity.get('sra_accession', [])
        for run in entity.get('sra_run', []):
            values += [run.get(k, '') for k in ('run', 'experiment', 'sample', 'biosample')]
    result = {acc for value in values for acc in _ID.findall(str(value))}
    if sample:
        # A common study/project does not establish a sample identity.
        result = {a for a in result if not re.match(r'(?:[SED]RP|PRJ|GSE|E-)', a)}
    return result


def complete_characteristic_groups(preferred, fallback):
    """Complete annotations only on unique, exact, compatible value occurrences."""
    from copy import deepcopy
    result = deepcopy(preferred)
    term_keys = {'term_source_ref', 'term_accession_number'}

    def compatible(a, b):
        if isinstance(a, dict) and isinstance(b, dict):
            left = {k for k in term_keys if informative(a.get(k))}
            right = {k for k in term_keys if informative(b.get(k))}
            if left and right and not left & right:
                return False
            return all(compatible(a[k], b[k]) for k in a.keys() & b.keys())
        return not informative(a) or not informative(b) or a == b

    def complete(a, b):
        for key, value in b.items():
            if isinstance(a.get(key), dict) and isinstance(value, dict):
                complete(a[key], value)
            elif not informative(a.get(key)) and informative(value):
                a[key] = deepcopy(value)

    for value in result:
        same = lambda v: v.get('name') == value.get('name') and v.get('value') == value.get('value')
        candidates = [v for v in fallback if same(v)]
        if len(candidates) != 1 or sum(same(v) for v in result) != 1:
            continue
        if not informative(value.get('value')):
            continue
        group = {k:deepcopy(v) for k,v in value.items() if k in term_keys | {'unit'}}
        old = {k:v for k,v in candidates[0].items() if k in term_keys | {'unit'}}
        if compatible(group, old):
            complete(group, old)
            value.update(group)
    return result


def _merge_publications(current, incoming, prefer, *, issues=None, scope='entity'):
    """Retain scoped papers and complete only compatible identifier-bound groups."""
    from collections import Counter
    from copy import deepcopy
    from itertools import combinations
    import json
    import logging
    from ..miniml.publication_identifiers import publication_identity

    result = deepcopy(current)
    status_fields = {'status', 'status_term_source_ref', 'status_term_accession_number'}
    signature = lambda value: json.dumps(value, sort_keys=True)
    counts = Counter(signature(v) for v in result)
    occurrences = Counter()
    inventory = [publication_identity(v) for v in [*current, *incoming]]
    identities = inventory[len(current):]

    def compatible(left, right):
        return all(left[k] == right[k] for k in left.keys() & right.keys())

    def retain(value, component):
        key = signature(value)
        occurrences[key] += 1
        if counts[key] < occurrences[key]:
            result.append(deepcopy(value))
            result_components.append(component)
            counts[key] += 1

    def conflict(identifier):
        message = f'{scope}: conflicting publication evidence for {identifier}; retaining separate source citations'
        if issues is not None and message not in issues:
            issues.append(message)
        logging.getLogger(__name__).warning(message)

    # Validate the original inventory before accepting any new identity. Two
    # individually compatible rows can otherwise bridge different native papers.
    parents = list(range(len(inventory)))

    def root(index):
        while parents[index] != index:
            parents[index] = parents[parents[index]]
            index = parents[index]
        return index

    owners = {}
    for index, identity in enumerate(inventory):
        for identifier in identity.items():
            if identifier in owners:
                parents[root(index)] = root(owners[identifier])
            else:
                owners[identifier] = index
    components = {}
    for index in range(len(inventory)):
        components.setdefault(root(index), []).append(index)

    def details_conflict(i, j):
        left, right = incoming[i], incoming[j]
        differing = any(informative(left[k]) and informative(right[k]) and left[k] != right[k]
                        for k in left.keys() & right.keys() - {'pubmed_id', 'doi', 'pmcid'})
        terms = status_fields - {'status'}
        a = {k for k in terms if informative(left.get(k))}
        b = {k for k in terms if informative(right.get(k))}
        return differing or bool(a and b and not a & b)

    conflicted = set()
    for members in components.values():
        added = [i - len(current) for i in members if i >= len(current)]
        identifiers = {}
        for i in members:
            for kind, literal in inventory[i].items():
                identifiers.setdefault(kind, set()).add(literal)
        if (any(len(values) > 1 for values in identifiers.values())
                or any(details_conflict(i, j) for i, j in combinations(added, 2))):
            conflicted.update(added)
    incoming_counts = Counter()
    result_components = [root(index) for index in range(len(current))]
    for index, value in enumerate(incoming):
        identity = identities[index]
        component = root(len(current) + index)
        incoming_counts[component] += 1
        matches = [(i, old, publication_identity(old)) for i, old in enumerate(result)
                   if identity and result_components[i] == component]
        groups = [identity, *[ids for _, _, ids in matches]]
        if index in conflicted or any(not compatible(a, b) for a, b in combinations(groups, 2)):
            conflict(identity)
            retain(value, component)
            continue
        if not matches:
            retain(value, component)
            continue
        for position, old, old_ids in matches:
            primary, secondary = (value, old) if prefer else (old, value)
            primary_ids, secondary_ids = (identity, old_ids) if prefer else (old_ids, identity)
            merged = deepcopy(primary)
            for field, fallback in secondary.items():
                if field in status_fields:
                    continue
                present = field in primary_ids if field in {'pubmed_id', 'doi', 'pmcid'} else informative(primary.get(field))
                supplied = field in secondary_ids if field in {'pubmed_id', 'doi', 'pmcid'} else informative(fallback)
                if not present and supplied:
                    merged[field] = deepcopy(fallback)
            selected = {k: deepcopy(primary[k]) for k in status_fields if informative(primary.get(k))}
            fallback = {k: deepcopy(secondary[k]) for k in status_fields if informative(secondary.get(k))}
            if not informative(selected.get('status')) and informative(fallback.get('status')):
                selected, fallback = fallback, selected
            terms = status_fields - {'status'}
            a, b = terms & selected.keys(), terms & fallback.keys()
            if compatible(selected, fallback) and not (a and b and not a & b):
                selected.update({k: v for k, v in fallback.items() if k not in selected})
            for field in status_fields:
                merged.pop(field, None)
            merged.update(selected)
            result[position] = merged
        # Matching an earlier incoming row does not consume another supplied
        # occurrence. Retain the larger multiplicity across source lists. Only
        # unanimous completed facts may fill a newly retained occurrence.
        if len(matches) < incoming_counts[component]:
            completed = [result[position] for position, _, _ in matches]
            added = deepcopy(value)
            common = {k: deepcopy(v) for k, v in completed[0].items()
                      if k not in status_fields and all(k in other and other[k] == v for other in completed[1:])}
            statuses = [{k: v for k, v in item.items() if k in status_fields} for item in completed]
            if all(group == statuses[0] for group in statuses[1:]):
                for field in status_fields:
                    added.pop(field, None)
                common.update(deepcopy(statuses[0]))
            added.update(common)
            result.append(added)
            result_components.append(component)
    return result


def _merge_entity(target, extra, prefer, protected=(), *, issues=None):
    from copy import deepcopy
    additive = {'raw_data', 'supplementary_data', 'relation', 'accession', 'contact_ref', 'contributor_ref', 'pubmed_id'}
    for key, value in extra.items():
        if key in protected:
            continue
        if key == 'pubmed_publication':
            target[key] = _merge_publications(target.get(key, []), value, prefer, issues=issues,
                                             scope=target.get('iid') or target.get('run') or 'entity')
        elif key in additive:
            # Files retain all alternatives. Identity declarations need no duplicates.
            target.setdefault(key, [])
            target[key].extend(deepcopy(v) for v in value if key in {'raw_data', 'supplementary_data'} or v not in target[key])
        elif key == 'status':
            current = target.setdefault(key, [])
            for status in value:
                same = next((v for v in current if v.get('database') == status.get('database') and v.get('accession') == status.get('accession')
                             and all(not informative(v.get(k)) or not informative(literal) or v[k] == literal
                                     for k, literal in status.items())), None)
                if same is None: current.append(deepcopy(status))
                else:
                    for field, literal in status.items():
                        if informative(literal) and (prefer or not informative(same.get(field))): same[field] = deepcopy(literal)
        elif key == 'characteristics':
            groups = {}
            for item in value:
                groups.setdefault(item['name'], []).append(item)
            current = target.setdefault(key, [])
            for name, group in groups.items():
                old = [v for v in current if v['name'] == name]
                if informative(group) and (prefer or not informative(old)):
                    current[:] = [v for v in current if v['name'] != name] + complete_characteristic_groups(group, old)
                elif not old:
                    current.extend(deepcopy(group))
        elif key == 'organism' and informative(value) and (prefer or not informative(target.get(key))):
            organisms = deepcopy(value)
            for organism in organisms:
                same = [o for o in target.get(key, []) if o.get('value', '').strip().casefold() == organism.get('value', '').strip().casefold()]
                if len(same) == 1 and not organism.get('taxid') and same[0].get('taxid'):
                    organism['taxid'] = same[0]['taxid']
            target[key] = organisms
        elif informative(value) and (prefer or not informative(target.get(key))):
            # Replace complete typed values (including units / term references) together.
            target[key] = deepcopy(value)


def merge_archive_metadata(package, other, *, prefer=False, linked_accession=None):
    """Merge an explicitly related package; peer additions require a shared read study."""
    from copy import deepcopy
    from ..miniml import MINiMLCodec
    data, extra = package.to_mapping(), other.to_mapping()
    from ..miniml.archive_residuals import source_records, finalize
    evidence = source_records(package)
    if not hasattr(package, '_archive_source_records'):
        evidence.append({'provider':data.get('source',{}).get('format'), 'kind':'MINiML', 'accession':data['series']['iid'],
                         'metadata':{k:deepcopy(v) for k,v in data.items() if k not in ('extensions','source','miniml_schema_version')}})
    issues = []
    for document in (data, extra):
        source = document.get('source', {}).get('format', '')
        database = {'MAGE-TAB': 'ArrayExpress', 'GEO MINiML': 'GEO'}.get(source, source)
        for entity in [document['series'], *document.get('sample', [])]:
            for status in entity.get('status', []):
                status.setdefault('database', database)
    native_ids, extra_ids = entity_ids(data['series']), entity_ids(extra['series'])
    shared_read = {a for a in native_ids & extra_ids if re.fullmatch(r'[SED]RP\d+', a)}
    linked_match = linked_accession and linked_accession in linked_accessions(package) and linked_accession in extra_ids
    if (not prefer and not shared_read) or (prefer and not (native_ids & extra_ids or linked_match)):
        return package, ['enrichment: unresolved study identity']
    extension = data.setdefault('extensions', {}).setdefault('insdc', {'version': '1.0', 'records': []})
    if 'insdc' not in extra.get('extensions', {}) or prefer:
        extension['records'].append({'provider': extra.get('source', {}).get('format', 'linked'),
            'kind': 'MINiML', 'accession': extra['series'].get('iid'), 'metadata': deepcopy(extra)})
    from .archive_workflows import mark_file_origins
    if prefer:
        mark_file_origins(extra)
    from .archive_workflows import merge_declarations, merge_workflows
    namespace = str(extra['series'].get('iid') or 'linked') + ':'
    declaration_mappings = merge_declarations(data, extra, namespace)
    _merge_entity(data['series'], extra['series'], prefer,
                  {'iid', 'sample_ref', 'assay_paths', 'protocols'}, issues=issues)
    native_samples, extra_samples = data.get('sample', []), extra.get('sample', [])
    joins = {i: [j for j, s in enumerate(extra_samples) if entity_ids(n, sample=True) & entity_ids(s, sample=True)]
             for i, n in enumerate(native_samples)}
    reverse = {j: [i for i, js in joins.items() if j in js] for j in range(len(extra_samples))}
    matched = {}
    from ..miniml.archive_libraries import FIELDS as library_fields
    original_channels = {s['iid']: deepcopy(s.get('channel', [])) for s in native_samples}
    workflow_fields = ('scan_protocol', 'hybridization_protocol', 'data_processing')
    original_protocols = {s['iid']: {k:deepcopy(s[k]) for k in workflow_fields if k in s} for s in native_samples}
    for i, candidates in joins.items():
        if len(candidates) != 1 or len(reverse[candidates[0]]) != 1:
            if candidates:
                issues.append(f"{native_samples[i]['iid']}: ambiguous enrichment sample identity")
            continue
        j = candidates[0]
        target, source = native_samples[i], extra_samples[j]
        matched[source['iid']] = target['iid']
        _merge_entity(target, source, prefer, {'iid', 'channel', 'channel_count', 'sra_run', 'ena_accession', 'sra_accession', 'platform_ref'} | library_fields, issues=issues)
        runs = {r['run']: r for r in target.get('sra_run', [])}
        for run in source.get('sra_run', []):
            if run['run'] in runs:
                current = runs[run['run']]
                if run.get('experiment') and current.get('experiment') and run['experiment'] != current['experiment']:
                    issues.append(f"{run['run']}: conflicting enrichment experiment identity")
                    continue
                publication_fields = {'pubmed_id', 'pubmed_publication'}
                _merge_entity(current, {k: v for k, v in run.items() if k in publication_fields}, prefer, issues=issues)
                _merge_entity(current, run, False, {'run', 'sample', 'study', 'files', 'fastq_files'} | publication_fields | (library_fields if prefer else set()), issues=issues)
                for key in ('files', 'fastq_files'):
                    current.setdefault(key, []).extend(deepcopy(run.get(key, [])))
            elif not prefer:
                target.setdefault('sra_run', []).append(deepcopy(run))
        if source.get('platform_ref') and (prefer or not target.get('platform_ref')):
            target['platform_ref'] = deepcopy(source['platform_ref'])
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
    explicit = merge_workflows(data, extra, matched, proto_names, prefer, issues, original_channels=original_channels)
    from ..miniml.archive_libraries import synchronize_library_facts, resolve_library_facts, apply_library_facts
    if prefer and not extra['series'].get('assay_paths') and extra.get('source', {}).get('format') in {'GEO', 'GEO MINiML'}:
        for source in extra_samples:
            target = sample_map.get(matched.get(source['iid']))
            if target is None:
                continue
            for run in target.get('sra_run', []):
                incoming = [r for r in source.get('sra_run', []) if r.get('run') == run['run']
                            and (not r.get('experiment') or r.get('experiment') == run.get('experiment'))]
                if incoming:
                    apply_library_facts(run, [], resolve_library_facts(run, incoming, issues, run['run']))
    synchronize_library_facts(data, issues)
    if prefer and extra.get('source', {}).get('format') == 'MAGE-TAB':
        uncovered = [s['iid'] for s in native_samples if s['iid'] in set(matched.values()) - explicit]
        # AE sample scalars are projections of its workflows. Rejected paths
        # cannot supply a material/protocol for an unrelated native acquisition.
        from ..miniml.archive_paths import complete_native_paths
        # Supplied methods remain unused definitions even when their application
        # cannot be projected onto all acquisitions of the matched sample.
        complete_native_paths({'source': data.get('source', {}),
            'series': {'iid': data['series']['iid'], 'protocols': data['series'].setdefault('protocols', []), 'assay_paths': []},
            'sample': [sample_map[target] for target in uncovered]})
        scoped_fields = ('molecule', 'growth_protocol', 'treatment_protocol', 'extract_protocol', 'label_protocol')
        for target in uncovered:
            for field in workflow_fields:
                if field in original_protocols[target]: sample_map[target][field] = deepcopy(original_protocols[target][field])
                else: sample_map[target].pop(field, None)
            current = sample_map[target].get('channel', [])
            previous = original_channels.get(target, [])
            if len(current) == len(previous) == 1:
                for field in scoped_fields:
                    if field in previous[0]: current[0][field] = deepcopy(previous[0][field])
                    else: current[0].pop(field, None)
    from ..miniml.archive_paths import complete_native_paths
    complete_native_paths(data)
    for sample in native_samples:
        for key in ('library_strategy', 'library_selection', 'library_source'):
            values = {r.get(key) for r in sample.get('sra_run', [])}
            if len(values) > 1:
                sample.pop(key, None)
    data['series']['sample_ref'] = [{'ref': s['iid']} for s in native_samples]
    from ..miniml.archive_entities import declare_ontologies
    declare_ontologies(data, issues)
    # Keep only unrepresented incoming fields. Rewrite the snapshot's references
    # using the same verified joins and protocol namespace as the core merge.
    snapshot = deepcopy(extra)
    def rewrite(value):
        if isinstance(value, list):
            for v in value: rewrite(v)
        elif isinstance(value, dict):
            for key, v in list(value.items()):
                if key in ('iid','sample_ref','ref') and isinstance(v,str): value[key] = matched.get(v,v)
                elif key == 'protocol_ref' and isinstance(v,str) and v in proto_names: value[key] = proto_names[v]
                elif key != 'extensions': rewrite(v)
    rewrite(snapshot)
    for protocol in snapshot['series'].get('protocols', []):
        protocol['name'] = proto_names.get(protocol['name'], protocol['name'])
    snapshot['series']['iid'] = data['series']['iid']
    evidence.append({'provider':extra.get('source',{}).get('format'), 'kind':'MINiML', 'accession':extra['series']['iid'],
                     'metadata':{k:v for k,v in snapshot.items() if k not in ('extensions','source','miniml_schema_version')}})
    for occurrence in extra.get('extensions', {}).get('magetab', {}).get('unbound_annotations', []):
        retained = deepcopy(occurrence)
        rewrite(retained)
        evidence.append({'provider':'ArrayExpress', 'kind':'magetab_unbound_annotation',
                         'accession':extra['series']['iid'], 'metadata':retained})
        issues.append(f"Unbound {occurrence.get('header')} at {occurrence.get('sdrf')} row "
                      f"{occurrence.get('row_index')}, column {occurrence.get('column_index')}; occurrence retained in residual metadata.")
    if not prefer:
        for record in source_records(other):
            if record['kind'] == 'term_source_declaration':
                old = record['metadata']['iid']
                record['metadata']['iid'] = declaration_mappings.get(old, old)
            evidence.append(record)
    evidence.extend(r for r in data['extensions']['insdc']['records'] if r['kind']=='annotation' and r not in evidence)
    return finalize(data, evidence), issues


class LinkedArchiveEnricher:
    def __init__(self, resource_profile='standard', geo_converter=None, ae_converter=None):
        self.profile, self.geo, self.ae = resource_profile, geo_converter, ae_converter

    def enrich(self, package):
        return self.enrich_selected(package)

    def enrich_selected(self, package, *, accessions=None):
        from .preparation_scope import convert_source
        issues = []
        # Stable ascending priority: informative AE values are applied last.
        discovered = set(linked_accessions(package) if accessions is None else accessions)
        # The legacy E-GEOD namespace is an exact GEO accession lookup, not an
        # asserted identity. Accept it only after retrieving matching metadata.
        candidates = {'E-GEOD-' + acc[3:] for acc in discovered if acc.startswith('GSE')}
        probes = candidates - discovered
        links = sorted(discovered | candidates, key=lambda a: (a.startswith('E-'), a))
        for accession in links:
            try:
                if accession.startswith('GSE'):
                    if self.geo is None:
                        from ..converters.geo2json import GEO2JSONConverter
                        self.geo = GEO2JSONConverter(resource_profile=self.profile)
                    packages = convert_source(self.geo, accession, enrich=False, related_series=False)
                else:
                    if self.ae is None:
                        from ..converters.ae2json import AE2JSONConverter
                        self.ae = AE2JSONConverter(resource_profile=self.profile)
                    packages = convert_source(self.ae, accession)
                candidates = [p for p in packages if accession in entity_ids(p.to_mapping()['series'])]
                if len(candidates) != 1:
                    issues.append(f'{accession}: ambiguous or unavailable enrichment study')
                    continue
                if accession in probes:
                    common = entity_ids(package.to_mapping()['series']) & entity_ids(candidates[0].to_mapping()['series'])
                    if 'GSE' + accession.removeprefix('E-GEOD-') not in common and not any(re.fullmatch(r'[SED]RP\d+', a) for a in common):
                        issues.append(f'{accession}: unresolved enrichment study identity')
                        continue
                from .archive_diagnostics import consistency_issues
                issues.extend(f'{accession}: {issue}' for issue in consistency_issues(candidates[0].to_mapping()))
                package, merge_issues = merge_archive_metadata(package, candidates[0], prefer=True, linked_accession=accession)
                issues.extend(merge_issues)
            except Exception as error:
                issues.append(f'{accession}: enrichment unavailable ({type(error).__name__})')
        return package, issues
