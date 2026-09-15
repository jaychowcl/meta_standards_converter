# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Accession-bound workflow and declaration merging; no retrieval or orchestration."""
from copy import deepcopy
import re
from pathlib import PurePosixPath
from urllib.parse import urlsplit, unquote, quote


def remap_references(value, mappings):
    if isinstance(value, list):
        for item in value:
            remap_references(item, mappings)
    elif isinstance(value, dict):
        for key, item in value.items():
            if key == 'extensions':
                continue
            if key in {'database', 'term_source_ref', 'ref', 'organization_ref', 'contributor_ref'} and isinstance(item, str):
                value[key] = mappings.get(item, item)
            else:
                remap_references(item, mappings)


def merge_declarations(data, extra, namespace):
    """Global database IDs agree by identity; conflicting local IDs are scoped."""
    mappings = {}
    global_ids = {'GEO', 'ArrayExpress', 'BioProject', 'BioSample', 'ENA', 'SRA', 'DRA', 'INSDC', 'NCBITaxon'}
    for document in (data, extra):
        known = {d['iid'] for d in document.get('database', [])}
        for record in document.get('extensions', {}).get('insdc', {}).get('records', []):
            if record['kind'] == 'term_source_declaration':
                declaration = record['metadata']
                if declaration['iid'] not in known:
                    document.setdefault('database', []).append(deepcopy(declaration))
                    known.add(declaration['iid'])
    existing = {d['iid']: d for d in data.get('database', [])}
    def conflicts(old, incoming):
        return old and any(old.get(k) and incoming.get(k) and old[k] != incoming[k]
                           for k in ('name', 'url', 'uri', 'version'))
    for incoming in extra.get('database', []):
        iid = incoming['iid']
        old = existing.get(iid)
        if conflicts(old, incoming) and iid not in global_ids:
            candidate, suffix = namespace + iid, 2
            while conflicts(existing.get(candidate), incoming):
                candidate = namespace + iid + ':' + str(suffix)
                suffix += 1
            incoming['iid'] = mappings[iid] = candidate
            old = existing.get(candidate)
        if old is None:
            new = deepcopy(incoming)
            data.setdefault('database', []).append(new)
            existing[new['iid']] = new
        else:
            for key, value in incoming.items():
                if value and not old.get(key):
                    old[key] = deepcopy(value)
    for kind in ('organization', 'contributor', 'platform'):
        for incoming in extra.get(kind, []):
            iid = incoming['iid']
            mappings[iid] = namespace + iid
            incoming['iid'] = mappings[iid]
    remap_references(extra, mappings)
    for kind in ('organization', 'contributor', 'platform'):
        current = {v['iid']: v for v in data.setdefault(kind, [])}
        for incoming in extra.get(kind, []):
            if incoming['iid'] in current:
                current[incoming['iid']].update(deepcopy(incoming))
            else:
                data[kind].append(deepcopy(incoming))
    # IDF contributors are study-scoped even when the input omitted explicit refs.
    if (extra.get('source', {}).get('format') == 'MAGE-TAB'
            and extra.get('contributor') and not extra['series'].get('contributor_ref')):
        extra['series']['contributor_ref'] = [{'ref': p['iid']} for p in extra['contributor']]
    return mappings


def path_ids(path):
    values = []
    for step in path.get('steps', []):
        if step.get('kind') in ('assay', 'scan'):
            values.append(step.get('name', ''))
            values.extend(c.get('value', '') for c in step.get('comments', [])
                          if c.get('name', '').upper() in {'ENA_RUN', 'ENA_EXPERIMENT', 'SRA_RUN', 'SRA_EXPERIMENT'})
    return {acc for value in values for acc in re.findall(r'\b[SED]R[RX]\d+\b', str(value))}


def compatible(scope, native):
    for pattern in (r'[SED]RR\d+', r'[SED]RX\d+'):
        wanted = {v for v in scope if re.fullmatch(pattern, v)}
        actual = {v for v in native if re.fullmatch(pattern, v)}
        if wanted and not wanted & actual:
            return False
    return True


def is_file(step):
    return step.get('kind', '').endswith('_file') or step.get('kind') == 'image_file'


def mark_file_origins(data):
    """Record-level source identity supports read-set selection without URL guesses."""
    repository = {'MAGE-TAB': 'ArrayExpress', 'GEO MINiML': 'GEO', 'GEO': 'GEO'}.get(data.get('source', {}).get('format'))
    if not repository:
        return
    origin = {'repository': repository, 'source_accession': data['series']['iid']}
    for path in data['series'].get('assay_paths', []):
        for node in path.get('steps', []):
            if is_file(node) and node.get('link'):
                node['link'].update(origin)
    for entity in [data['series'], *data.get('sample', [])]:
        for field in ('raw_data', 'supplementary_data'):
            for file in entity.get(field, []):
                file.update(origin)
        for run in entity.get('sra_run', []):
            for field in ('files', 'fastq_files'):
                for file in run.get(field, []):
                    file.update(origin)


def file_node(file, kind='array_data_file'):
    uri = _file_uri(file.get('uri') or file.get('value'))
    node = {'kind': kind, 'name': file.get('filename') or unquote(PurePosixPath(urlsplit(uri or '').path).name) or uri or ''}
    if uri:
        node['link'] = {'value': uri}
        if file.get('format') or file.get('type'):
            node['link']['type'] = file.get('format') or file['type']
        node['link'].update({k: deepcopy(file[k]) for k in ('repository', 'source_accession', 'companion_files') if k in file})
    node['comments'] = [{'name': label, 'value': str(file[key])} for key, label in
                        [('md5', 'MD5'), ('bytes', 'File size'), ('format', 'File format'),
                         ('checksum', 'Checksum'), ('checksum_method', 'Checksum method')]
                        if file.get(key)]
    return node


def _file_uri(value):
    # ArrayExpress FASTQ_URI historically supplies this explicit FTP location
    # without a scheme. Treat the remainder as a literal path, including '#'.
    if isinstance(value, str) and value.startswith('ftp.sra.ebi.ac.uk/'):
        return 'ftp://ftp.sra.ebi.ac.uk/' + quote(value.split('/', 1)[1], safe="/%:@!$&'()*+,;=-._~")
    return value


def _file_projection_facts(node):
    """Comparable file annotations, independent of link/comment spelling."""
    aliases = {'file format': 'type', 'format': 'type', 'file uri': 'value',
               'uri': 'value', 'fastq uri': 'value', 'file size': 'bytes',
               'size': 'bytes', 'fastq bytes': 'bytes', 'fastq md5': 'md5',
               'checksum algorithm': 'checksum method'}
    facts = {}
    entries = list(node.get('link', {}).items()) + [
        (c.get('name', ''), c.get('value')) for c in node.get('comments', [])]
    for name, value in entries:
        if value in (None, ''):
            continue
        key = name.lower().replace('_', ' ')
        key = aliases.get(key, key)
        if key == 'md5':
            facts.setdefault('checksum method', set()).add('md5')
            key = 'checksum'
        text = str(value)
        if key == 'checksum method':
            text = text.lower()
        facts.setdefault(key, set()).add(text)
    return facts


def _compatible_file_projection(existing, incoming):
    """Decline consolidation when any shared annotation conflicts."""
    old, new = _file_projection_facts(existing), _file_projection_facts(incoming)
    return all(old[k] == new[k] and len(old[k]) == 1 for k in old.keys() & new.keys())


def _prepare_path(path, target, proto_names):
    result = deepcopy(path)
    # Source document names are evidence, not separate native output datasets.
    result.pop('document', None)
    steps = result['steps']
    scope = path_ids(path)
    for pattern in (r'[SED]RR\d+', r'[SED]RX\d+'):
        if len({v for v in scope if re.fullmatch(pattern, v)}) > 1:
            raise ValueError('conflicting enrichment run or experiment identifiers')
    file_names = {s.get('name') for s in steps if is_file(s)}
    for step in steps:
        file_names.update(unquote(PurePosixPath(urlsplit(_file_uri(c['value'])).path).name)
                          for c in step.get('comments', [])
                          if c.get('name', '').upper() == 'FASTQ_URI' and c.get('value'))
    expanded = []
    for step in steps:
        if step.get('sample_ref') or step.get('kind') in ('source', 'sample'):
            step['sample_ref'] = target
        if step.get('protocol_ref') in proto_names:
            step['protocol_ref'] = proto_names[step['protocol_ref']]
        comments = step.get('comments', [])
        submitted_alias = (any(c.get('name', '').upper() == 'SUBMITTED_FILE_NAME' and c.get('value') == step.get('name') for c in comments)
                           and any(c.get('name', '').upper() == 'FASTQ_URI' and c.get('value') for c in comments))
        if step.get('kind') == 'scan' and (step.get('name') in file_names or submitted_alias):
            runs = {v for v in scope if re.fullmatch(r'[SED]RR\d+', v)}
            if len(runs) == 1:
                step['name'] = next(iter(runs))
        if step.get('link', {}).get('value'):
            step['link']['value'] = _file_uri(step['link']['value'])
        comments = step.get('comments', [])
        reads = [file_node({'uri': c['value'], 'format': 'fastq'}) for c in comments
                 if c.get('name', '').upper() == 'FASTQ_URI' and c.get('value')]
        fields = {'FASTQ_FILE_NAME': 'filename', 'FASTQ_MD5': 'MD5', 'FASTQ_BYTES': 'File size',
                  'FASTQ_FORMAT': 'File format', 'FASTQ_CHECKSUM': 'Checksum',
                  'FASTQ_CHECKSUM_METHOD': 'Checksum method', 'READ_TYPE': 'READ_TYPE', 'READ_INDEX': 'READ_INDEX',
                  'SUBMITTED_FILE_NAME': 'SUBMITTED_FILE_NAME'}
        moved = set()
        for label, destination in fields.items():
            values = [(i, c) for i, c in enumerate(comments) if c.get('name', '').upper() == label]
            if reads and len(values) == len(reads):
                for read, (i, comment) in zip(reads, values):
                    if destination == 'filename':
                        if read['name'] and read['name'] != comment['value']:
                            read['comments'].append(deepcopy(comment))
                        else:
                            read['name'] = comment['value']
                    else:
                        read['comments'].append({**deepcopy(comment), 'name': destination})
                    moved.add(i)
        kept = [c for i, c in enumerate(comments) if i not in moved
                and not (c.get('name', '').upper() == 'FASTQ_URI' and c.get('value'))]
        if 'comments' in step:
            step['comments'] = kept
        expanded.append(step)
        # A supplied run FASTQ is an acquisition output, preceding processing.
        expanded.extend(reads)
    result['steps'] = expanded
    return result


def _acquisition(steps):
    end = next((i for i, s in enumerate(steps) if is_file(s)), len(steps))
    prefix = steps[:end]
    # Without an explicit raw node, operations after the scan belong to results.
    if end == len(steps) or steps[end]['kind'] != 'array_data_file':
        scan = next((i for i, s in enumerate(prefix) if s['kind'] == 'scan'), None)
        if scan is not None:
            prefix = prefix[:scan + 1]
    return deepcopy(prefix)


def _biological_characteristics(preferred, fallback):
    from .archive_enrichment import _merge_entity
    result = {'characteristics': deepcopy(fallback)}
    _merge_entity(result, {'characteristics': preferred}, True)
    return result['characteristics']


def _bind_native(steps, target, native):
    for step in steps:
        if step.get('kind') in ('source', 'sample', 'assay', 'scan'):
            candidates = [s for s in native if s.get('kind') == step['kind']]
            old = next(iter(candidates), None)
            if len(candidates) == 1 and step['kind'] in ('source', 'sample'):
                step['characteristics'] = _biological_characteristics(
                    step.get('characteristics', []), old.get('characteristics', []))
            if old:
                step['name'] = old['name']
                step['sample_ref'] = target
                step.setdefault('comments', []).extend(c for c in deepcopy(old.get('comments', []))
                                                        if c not in step.get('comments', []))
    return steps


def _biological_node(sample):
    node = {'kind': 'source', 'name': sample['iid'], 'sample_ref': sample['iid']}
    if len(sample.get('channel', [])) == 1:
        channel = sample['channel'][0]
        node['characteristics'] = deepcopy(channel.get('characteristics', []))
        node['characteristics'] += [{'name': 'organism', 'value': o['value'],
            **({'term_source_ref': 'NCBITaxon', 'term_accession_number': o['taxid']} if o.get('taxid') else {})}
            for o in channel.get('organism', [])]
    return node


def merge_workflows(data, extra, matched, proto_names, prefer, issues, *, original_channels=None):
    """Bind complete ordered branches and acquisition prefixes by explicit identity."""
    paths = data['series'].setdefault('assay_paths', [])
    samples = {s['iid']: s for s in data['sample']}
    source_samples = {matched[s['iid']]: s for s in extra.get('sample', []) if s['iid'] in matched}
    # Refresh native biological projections before accepting authored workflows.
    # The actual incoming channel has priority; a scoped path can retain its own
    # compatible annotations when that channel supplies only a plain literal.
    from .archive_enrichment import _merge_entity
    for path in paths:
        for step in path['steps']:
            target = step.get('sample_ref')
            sample = samples.get(target)
            if sample and step['kind'] in ('source', 'sample'):
                old_channels = (original_channels or {}).get(target, [])
                if len(old_channels) == len(sample.get('channel', [])) == 1:
                    old = _biological_node({'iid': target, 'channel': old_channels})
                    if step.get('characteristics', []) == old.get('characteristics', []):
                        # Only a complete generated projection establishes which
                        # occurrence is taxonomy versus a supplied characteristic.
                        step['characteristics'] = _biological_node(sample).get('characteristics', [])
                        continue
                step['characteristics'] = _biological_characteristics(
                    step.get('characteristics', []), _biological_node(sample).get('characteristics', []))
                channels = source_samples.get(target, {}).get('channel', [])
                if len(channels) == 1:
                    _merge_entity(step, {'characteristics': channels[0].get('characteristics', [])}, prefer)
    originals = deepcopy(paths)
    templates, incoming, explicit = {}, [], set()
    for path in extra['series'].get('assay_paths', []):
        refs = {s['sample_ref'] for s in path['steps'] if s.get('sample_ref')}
        scope = path_ids(path)
        targets = {matched[r] for r in refs if r in matched}
        if not refs and scope:
            targets = {target for target in matched.values() if any(
                compatible(scope, {r.get('run'), r.get('experiment')})
                for r in samples[target].get('sra_run', []))}
        if len(targets) != 1 or (refs - matched.keys()):
            issues.append('enrichment workflow: ambiguous or unresolved sample binding')
            continue
        target = next(iter(targets))
        try:
            prepared = _prepare_path(path, target, proto_names)
            channels = source_samples.get(target, {}).get('channel', [])
            for step in prepared['steps']:
                if step['kind'] in ('source', 'sample') and len(channels) == 1:
                    step['characteristics'] = _biological_characteristics(
                        step.get('characteristics', []), channels[0].get('characteristics', []))
            if prefer:
                mark_file_origins({'source': extra.get('source', {}),
                                   'series': {'iid': extra['series']['iid'], 'assay_paths': [prepared]}})
        except ValueError as error:
            issues.append(f'{target}: {error}')
            continue
        candidates = [p for p in originals if any(s.get('sample_ref') == target for s in p['steps'])
                      and compatible(scope, path_ids(p))]
        run_scopes = {tuple(sorted(path_ids(p))) for p in candidates}
        if not prefer:
            # Peer runs may be new to the native inventory but must be present
            # in the already validated merged sample/run registry.
            if scope and not any(compatible(scope, {r.get('run'), r.get('experiment')})
                                 for r in samples[target].get('sra_run', [])):
                issues.append(f'{target}: incompatible peer workflow scope')
                continue
            paths.append(prepared)
            continue
        if not scope and not any(s['kind'] in ('assay', 'scan', 'extract', 'labeled_extract') for s in prepared['steps']):
            if any(is_file(s) for s in prepared['steps']):
                paths.append(prepared)
            continue
        if not candidates or (not scope and len(run_scopes) != 1):
            issues.append(f'{target}: ambiguous or incompatible enrichment workflow scope')
            continue
        prefix = _acquisition(prepared['steps'])
        for key in run_scopes:
            identity = (target, key)
            if prefix not in templates.setdefault(identity, []):
                templates[identity].append(prefix)
            incoming.append((identity, prepared))
    if prefer:
        from ..miniml.archive_libraries import path_facts, resolve_library_facts, apply_library_facts, _shared_facts
        selections = {}
        for (target, key), choices in templates.items():
            if len(choices) != 1:
                continue
            selected_runs = []
            for run in samples[target].get('sra_run', []):
                if not compatible(set(key), {run.get('run'), run.get('experiment')}):
                    continue
                supplied = path_facts(choices[0])
                for source in extra.get('sample', []):
                    if matched.get(source['iid']) == target:
                        supplied.extend(r for r in source.get('sra_run', []) if r.get('run') == run['run']
                                        and (not r.get('experiment') or r['experiment'] == run.get('experiment')))
                facts = resolve_library_facts(run, supplied, issues, run['run'])
                selected_runs.append(facts)
                apply_library_facts(run, [], facts)
            selections[(target, key)] = _shared_facts(selected_runs, issues, f'{target} {list(key)}')
        for path in paths:
            target = next((s.get('sample_ref') for s in path['steps'] if s.get('sample_ref')), None)
            choices = templates.get((target, tuple(sorted(path_ids(path)))), [])
            if len(choices) > 1:
                message = f'{target}: ambiguous enrichment workflows retained outside core'
                if message not in issues:
                    issues.append(message)
                continue
            if not choices:
                continue
            prefix = _bind_native(deepcopy(choices[0]), target, path['steps'])
            facts = selections.get((target, tuple(sorted(path_ids(path)))))
            if facts is not None:
                apply_library_facts({}, prefix, facts)
            # An experiment-level workflow does not delete the native run.
            for kind in ('assay', 'scan'):
                if not any(s['kind'] == kind for s in prefix):
                    native = next((s for s in path['steps'] if s['kind'] == kind), None)
                    if native:
                        position = next((i for i, s in enumerate(prefix) if s['kind'] == 'scan'), len(prefix)) if kind == 'assay' else len(prefix)
                        prefix.insert(position, deepcopy(native))
            if any(s['kind'] in ('extract', 'labeled_extract') for s in prefix):
                explicit.add((target, tuple(sorted(path_ids(path)))))
            # Preserve the complete original branch after raw acquisition.
            end = next((i for i, s in enumerate(path['steps']) if is_file(s)), len(path['steps']))
            path['steps'] = prefix + deepcopy(path['steps'][end:])
        for (target, key), prepared in incoming:
            if len(templates[(target, key)]) != 1 or not any(is_file(s) for s in prepared['steps']):
                continue
            candidates = [p for p in originals if any(s.get('sample_ref') == target for s in p['steps'])
                          and tuple(sorted(path_ids(p))) == key]
            # All candidates have the same verified sample/run identity.
            steps = _bind_native(deepcopy(prepared['steps']), target, candidates[0]['steps'])
            facts = selections.get((target, key))
            if facts is not None:
                apply_library_facts({}, steps, facts)
            paths.append({**deepcopy(prepared), 'steps': steps})
        for source in extra.get('sample', []):
            target = matched.get(source['iid'])
            if target is None:
                continue
            for run in source.get('sra_run', []):
                scope = {run[k] for k in ('run', 'experiment') if run.get(k)}
                candidates = [p for p in paths if any(s.get('sample_ref') == target for s in p['steps'])
                              and compatible(scope, path_ids(p))]
                prefixes = []
                for path in candidates:
                    prefix = _acquisition(path['steps'])
                    if prefix not in prefixes:
                        prefixes.append(prefix)
                if scope and len(prefixes) == 1:
                    for file in run.get('files') or run.get('fastq_files', []):
                        if file.get('uri') or file.get('filename'):
                            paths.append({'steps': deepcopy(prefixes[0]) + [file_node(file)]})
            for field, kind in [('raw_data', 'array_data_file'), ('supplementary_data', 'derived_array_data_file')]:
                for link in source.get(field, []):
                    if link.get('value'):
                        node = file_node(link, kind)
                        matches = [s for p in paths if {n.get('sample_ref') for n in p['steps'] if n.get('sample_ref')} == {target}
                                   for s in p['steps'] if s['kind'] == kind and (
                                       s.get('link', {}).get('value') == _file_uri(link['value'])
                                       or (s.get('name') == node['name'] == link['value']
                                           and not s.get('link', {}).get('value')))]
                        compatible_files = matches and all(_compatible_file_projection(s, node) for s in matches)
                        if compatible_files:
                            for s in matches:
                                s.setdefault('link', {}).update(node.get('link', {}))
                                for comment in node.get('comments', []):
                                    if comment not in s.setdefault('comments', []): s['comments'].append(deepcopy(comment))
                        else:
                            paths.append({'steps': [_biological_node(samples[target]), node]})
    # Complete channel fallbacks only after binding scoped workflow evidence.
    # Never replace an authored biological-node annotation with a sample copy.
    for path in paths:
        for step in path['steps']:
            sample = samples.get(step.get('sample_ref'))
            if sample and step['kind'] in ('source', 'sample'):
                step['characteristics'] = _biological_characteristics(
                    step.get('characteristics', []), _biological_node(sample).get('characteristics', []))
    # Scalar material fields are safe only when every native acquisition for
    # that sample received the explicit incoming workflow.
    return {target for target in samples if any(t == target for t, _ in explicit)
            and all((target, tuple(sorted(path_ids(p)))) in explicit for p in originals
                    if path_ids(p) and any(s.get('sample_ref') == target for s in p['steps']))}
