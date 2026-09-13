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
from urllib.parse import urlsplit, unquote


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
    existing = {d['iid']: d for d in data.get('database', [])}
    for incoming in extra.get('database', []):
        iid = incoming['iid']
        old = existing.get(iid)
        conflicts = old and any(old.get(k) and incoming.get(k) and old[k] != incoming[k]
                                for k in ('name', 'url', 'uri', 'version'))
        if conflicts and iid not in global_ids:
            incoming['iid'] = mappings[iid] = namespace + iid
            old = None
        if old is None:
            data.setdefault('database', []).append(deepcopy(incoming))
        else:
            for key, value in incoming.items():
                if value and not old.get(key):
                    old[key] = deepcopy(value)
    for kind in ('organization', 'contributor'):
        for incoming in extra.get(kind, []):
            iid = incoming['iid']
            mappings[iid] = namespace + iid
            incoming['iid'] = mappings[iid]
    remap_references(extra, mappings)
    for kind in ('organization', 'contributor'):
        data.setdefault(kind, []).extend(deepcopy(extra.get(kind, [])))
    # IDF contributors are study-scoped even when the input omitted explicit refs.
    if extra.get('contributor') and not extra['series'].get('contributor_ref'):
        extra['series']['contributor_ref'] = [{'ref': p['iid']} for p in extra['contributor']]


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


def file_node(file, kind='array_data_file'):
    uri = file.get('uri') or file.get('value')
    node = {'kind': kind, 'name': file.get('filename') or unquote(PurePosixPath(urlsplit(uri or '').path).name) or uri or ''}
    if uri:
        node['link'] = {'value': uri}
        if file.get('format') or file.get('type'):
            node['link']['type'] = file.get('format') or file['type']
    node['comments'] = [{'name': label, 'value': str(file[key])} for key, label in
                        [('md5', 'MD5'), ('bytes', 'File size'), ('format', 'File format'),
                         ('checksum', 'Checksum'), ('checksum_method', 'Checksum method')]
                        if file.get(key)]
    return node


def merge_workflows(data, extra, matched, proto_names, prefer, issues):
    """Replace uniquely matched complete workflows while retaining every file branch."""
    paths = data['series'].setdefault('assay_paths', [])
    originals = deepcopy(paths)
    native_samples = {s['iid']: s for s in data['sample']}
    templates, incoming_files, explicit = {}, [], set()
    for path in extra['series'].get('assay_paths', []):
        refs = {s['sample_ref'] for s in path['steps'] if s.get('sample_ref')}
        scope = path_ids(path)
        targets = {matched[r] for r in refs if r in matched}
        # Older packages omitted sample_ref; recover only via verified accessions.
        if not targets and scope:
            targets = {target for target in matched.values() if any(
                compatible(scope, {r.get('run'), r.get('experiment')})
                for r in native_samples[target].get('sra_run', []))}
        if len(targets) != 1 or (refs & matched.keys() and refs - matched.keys()):
            issues.append('enrichment workflow: ambiguous or unresolved sample binding')
            continue
        target = next(iter(targets))
        candidates = [p for p in originals if any(s.get('sample_ref') == target for s in p['steps'])
                      and compatible(scope, path_ids(p))]
        run_scopes = {tuple(sorted(path_ids(p))) for p in candidates}
        if not candidates or (not scope and len(run_scopes) != 1):
            issues.append(f'{target}: ambiguous or incompatible enrichment workflow scope')
            continue
        steps = deepcopy(path['steps'])
        for step in steps:
            if step.get('sample_ref') or step.get('kind') in ('source', 'sample'):
                step['sample_ref'] = target
            if step.get('protocol_ref') in proto_names:
                step['protocol_ref'] = proto_names[step['protocol_ref']]
        for step in steps:
            kept = []
            for comment in step.get('comments', []):
                if comment.get('name', '').upper() == 'FASTQ_URI' and comment.get('value'):
                    node = file_node({'uri': comment['value'], 'format': 'fastq'})
                    node['comments'].append(deepcopy(comment))
                    incoming_files.append((target, scope, node))
                else:
                    kept.append(comment)
            if 'comments' in step:
                step['comments'] = kept
        workflow = [s for s in steps if not is_file(s)]
        explicit_workflow = any(s.get('kind') in ('extract', 'labeled_extract') for s in workflow)
        if explicit_workflow:
            explicit.add(target)
        for key in run_scopes:
            identity = (target, key)
            if workflow not in templates.setdefault(identity, []):
                templates[identity].append(workflow)
        for step in steps:
            if is_file(step):
                incoming_files.append((target, scope, step))
        if not prefer:
            paths.append({'steps': steps})
    if prefer and templates:
        data['extensions']['insdc']['records'].append({'provider': data.get('source', {}).get('format', 'native'),
            'kind': 'MINiML_workflows', 'accession': data['series']['iid'], 'metadata': {'assay_paths': originals}})
        for path in paths:
            target = next((s.get('sample_ref') for s in path['steps'] if s.get('sample_ref')), None)
            choices = templates.get((target, tuple(sorted(path_ids(path)))), [])
            if len(choices) > 1:
                issues.append(f'{target}: ambiguous enrichment workflows retained outside core')
                continue
            if not choices:
                continue
            steps = deepcopy(choices[0])
            for step in steps:
                if step.get('kind') in ('source', 'sample', 'assay', 'scan'):
                    native = next((s for s in path['steps'] if s.get('kind') == step['kind']), None)
                    if native:
                        step['name'] = native['name']
                        step['sample_ref'] = target
                        step.setdefault('comments', []).extend(c for c in deepcopy(native.get('comments', [])) if c not in step.get('comments', []))
            path['steps'] = steps + deepcopy([s for s in path['steps'] if is_file(s)])
    # Update biological projections only on biological nodes, never on extracts.
    for path in paths:
        for step in path['steps']:
            sample = native_samples.get(step.get('sample_ref'))
            if sample and step.get('kind') in ('source', 'sample') and len(sample.get('channel', [])) == 1:
                channel = sample['channel'][0]
                values = deepcopy(channel.get('characteristics', []))
                values += [{'name': 'organism', 'value': o['value'], **({'term_source_ref': 'NCBITaxon', 'term_accession_number': o['taxid']} if o.get('taxid') else {})} for o in channel.get('organism', [])]
                step['characteristics'] = values
    if prefer:
        for target, scope, node in incoming_files:
            candidate = next((p for p in paths if any(s.get('sample_ref') == target for s in p['steps']) and compatible(scope, path_ids(p))), None)
            if candidate:
                paths.append({'steps': deepcopy([s for s in candidate['steps'] if not is_file(s)]) + [node]})
        for source in extra.get('sample', []):
            target = matched.get(source['iid'])
            if target is None:
                continue
            for run in source.get('sra_run', []):
                scope = {run[k] for k in ('run', 'experiment') if run.get(k)}
                candidate = next((p for p in paths if any(s.get('sample_ref') == target for s in p['steps']) and compatible(scope, path_ids(p))), None)
                if candidate and scope:
                    files = run.get('files') or run.get('fastq_files', [])
                    for file in files:
                        if file.get('uri') or file.get('filename'):
                            paths.append({'steps': deepcopy([s for s in candidate['steps'] if not is_file(s)]) + [file_node(file)]})
            represented = {s.get('link', {}).get('value') for p in extra['series'].get('assay_paths', []) for s in p['steps']}
            represented.update(f.get('uri') for run in source.get('sra_run', []) for f in (run.get('files') or run.get('fastq_files', [])))
            for field, kind in [('raw_data', 'array_data_file'), ('supplementary_data', 'derived_array_data_file')]:
                for link in source.get(field, []):
                    if link.get('value') and link['value'] not in represented:
                        # This is a sample result; no run is implied by list order.
                        paths.append({'steps': [{'kind': 'source', 'name': target, 'sample_ref': target,
                                                 'characteristics': deepcopy(native_samples[target].get('channel', [{}])[0].get('characteristics', []))}, file_node(link, kind)]})
    return explicit
