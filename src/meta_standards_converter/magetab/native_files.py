# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Project native archive download alternatives onto experimental run rows.

This module operates only on the constructor's private export copy. Stored
MINiML keeps every file occurrence. Derived files remain explicit workflow nodes.
"""
from copy import deepcopy
import json
import re


_FIELDS = ('NAME', 'URI', 'FORMAT', 'ROLE', 'BYTES', 'MD5', 'CHECKSUM', 'CHECKSUM_METHOD')
_ALIASES = {'FILE FORMAT': 'FORMAT', 'FILE SIZE': 'BYTES', 'CHECKSUM METHOD': 'CHECKSUM_METHOD',
            'FILE URI': 'URI', 'FASTQ_URI': 'URI', 'FASTQ_FILE_NAME': 'NAME',
            'FASTQ_MD5': 'MD5', 'FASTQ_BYTES': 'BYTES'}


def _key(value):
    return json.dumps(value, sort_keys=True, ensure_ascii=False)


def _record(step):
    link = step.get('link') or {}
    result = {'NAME': step.get('name', ''), 'URI': link.get('value', '')}
    interpreted = {'value', 'type'}
    for field, target in (('bytes', 'BYTES'), ('role', 'ROLE'), ('checksum_method', 'CHECKSUM_METHOD'),
                          ('file_checksum', 'CHECKSUM'), ('aspera', 'ASPERA'), ('galaxy', 'GALAXY')):
        if link.get(field) not in (None, ''):
            result[target] = str(link[field])
            interpreted.add(field)
    if link.get('checksum'):
        result['MD5'] = str(link['checksum'])
        interpreted.add('checksum')
        if 'CHECKSUM' not in result:
            result['CHECKSUM'] = result['MD5']
            if result.get('CHECKSUM_METHOD', 'MD5').upper() != 'MD5':
                result.setdefault('_extra', []).append({'name':'CHECKSUM_METHOD', 'value':result['CHECKSUM_METHOD']})
            result['CHECKSUM_METHOD'] = 'MD5'
    for field, target in (('repository', '_repository'), ('source_accession', '_source')):
        if isinstance(link.get(field), str):
            result[target] = link[field]
            interpreted.add(field)
    for comment in step.get('comments', []):
        name = str(comment.get('name', '')).upper()
        name = _ALIASES.get(name, name)
        value = str(comment.get('value', ''))
        if name in result and result[name] != value:
            # Preserve repeated conflicting annotations rather than choosing one.
            result.setdefault('_extra', []).append(deepcopy(comment))
        else:
            result[name] = value
    if not result.get('FORMAT'):
        if any(c.get('name', '').upper() == 'FASTQ_URI' for c in step.get('comments', [])):
            result['FORMAT'] = 'fastq'
        elif link.get('type'):
            result['FORMAT'] = link['type']
    extra = {k: v for k, v in step.items() if k not in {'kind', 'name', 'link', 'comments'}}
    if extra:
        result['_node'] = extra
    if set(link) - interpreted:
        result['_link'] = {k: v for k, v in link.items() if k not in interpreted}
    return result


def file_metadata_columns(step):
    """Render structured link facts while retaining supplied comment spelling."""
    record = _record(step)
    comments = step.get('comments', [])
    result = [('Comment[File URI]', record['URI'])] if record.get('URI') else []
    # Preserve existing GEO/AE columns. Add represented link fields only where
    # a supplied comment has not already carried the same semantic value.
    for key in ('FORMAT', 'ROLE', 'BYTES', 'CHECKSUM', 'CHECKSUM_METHOD', 'ASPERA', 'GALAXY'):
        value = record.get(key)
        aliases = {key, 'MD5'} if key == 'CHECKSUM' and record.get('CHECKSUM_METHOD') == 'MD5' else {key}
        if value and not any(_ALIASES.get(str(c.get('name', '')).upper(), str(c.get('name', '')).upper()) in aliases
                             and str(c.get('value', '')) == str(value) for c in comments):
            result.append((f'Comment[{key}]', value))
    if record.get('MD5') and record.get('CHECKSUM_METHOD') != 'MD5':
        if not any(c.get('name', '').upper() == 'MD5' and c.get('value') == record['MD5'] for c in comments):
            result.append(('Comment[MD5]', record['MD5']))
    return result


def _checksums(record):
    if any(str(c.get('name', '')).upper() in {'MD5','CHECKSUM','CHECKSUM_METHOD'} for c in record.get('_extra', [])):
        return {}
    values = {}
    pairs = [('MD5', record.get('MD5'))]
    if record.get('CHECKSUM'):
        pairs.append((record.get('CHECKSUM_METHOD', ''), record['CHECKSUM']))
    for algorithm, value in pairs:
        if not value:
            continue
        algorithm = str(algorithm).upper().replace('-', '')
        length = {'MD5': 32, 'SHA1': 40, 'SHA256': 64, 'SHA512': 128}.get(algorithm)
        if not length or not re.fullmatch(r'[0-9a-fA-F]{%d}' % length, str(value)):
            return {}
        value = str(value).lower()
        if algorithm in values and values[algorithm] != value:
            return {}
        values[algorithm] = value
    return values


def _combine(records, candidate, *, preserve_primary=False):
    """Complete explicit identities or checksum-verified, compatible mirrors."""
    for record in records:
        locations = {r.get('URI') for r in [record, *record.get('_alternatives', [])]} - {None, ''}
        identity = bool(candidate.get('URI') and candidate['URI'] in locations)
        same_repository = not (candidate.get('_repository') and record.get('_repository')
                               and candidate['_repository'] != record['_repository'])
        identity = identity or (same_repository and candidate.get('NAME') and candidate.get('NAME') == record.get('NAME')
                                and (not candidate.get('URI') or not record.get('URI')))
        left, right = _checksums(record), _checksums(candidate)
        common = left.keys() & right.keys()
        mirror = (_fastq(record) and _fastq(candidate) and bool(common) and all(left[k] == right[k] for k in common)
                  and candidate.get('NAME') and candidate['NAME'] == record.get('NAME'))
        ignored = {'URI', '_alternatives', '_repository', '_source'}
        if mirror:
            ignored |= {'MD5', 'CHECKSUM', 'CHECKSUM_METHOD'}
        compatible = all(k in ignored or not record.get(k) or not v or record[k] == v
                         for k, v in candidate.items())
        if (identity or mirror) and compatible:
            from .native_file_selection import _set
            if not preserve_primary and _set(candidate) < _set(record) and candidate.get('URI'):
                previous = deepcopy(record)
                record.clear()
                record.update(deepcopy(candidate))
                for k, v in previous.items():
                    if k != '_alternatives' and v and not record.get(k): record[k] = v
                alternatives = previous.pop('_alternatives', [])
                if previous.get('URI') != record.get('URI'):
                    alternatives.insert(0, previous)
                if alternatives:
                    record.setdefault('_alternatives', []).extend(alternatives)
                return
            if candidate.get('URI') and locations and candidate['URI'] not in locations:
                record.setdefault('_alternatives', []).append(deepcopy(candidate))
            record.update({k: v for k, v in candidate.items()
                           if k != '_alternatives' and v and not record.get(k)})
            return
    records.append(deepcopy(candidate))


def _fastq(record):
    return str(record.get('FORMAT', '')).lower() in {'fastq', 'fastq.gz'} and bool(record.get('URI'))


def _scan_records(scan):
    records, kept = [], []
    current = {}
    for comment in scan.get('comments', []):
        name = comment.get('name', '')
        if name.startswith('SUBMITTED_FILE_') and name[15:] in _FIELDS:
            field = name[15:]
            if field in current:
                records.append(current)
                current = {}
            current[field] = comment.get('value', '')
        else:
            kept.append(comment)
    if current:
        records.append(current)
    scan['comments'] = kept
    return records


def _comments(record, prefix, *, fixed=False):
    fields = list(_FIELDS) + sorted(k for k in record if k not in _FIELDS and not k.startswith('_'))
    result = [{'name': prefix + k, 'value': record.get(k, '')}
              for k in fields if fixed or record.get(k)]
    for field in ('_extra', '_node', '_link'):
        if record.get(field):
            result.append({'name': prefix + field[1:].upper(), 'value': _key(record[field])})
    if prefix in ('ARCHIVE_FILE_', 'SUBMITTED_FILE_'):
        for alternative in record.get('_alternatives', []):
            result.extend(_comments(alternative, prefix, fixed=fixed))
    return result


def _complete_file(long, short):
    """Combine an exact file occurrence; a bare member name is not a URI."""
    left, right = _record(long), _record(short)
    if left.get('_extra') or right.get('_extra'):
        return None
    for record in (left, right):
        if record.get('URI') == record.get('NAME'):
            record['URI'] = ''
    # Mirrors are handled by run grouping, not by sparse workflow inference.
    if left.get('URI') and right.get('URI') and left['URI'] != right['URI']:
        return None
    combined = [deepcopy(left)]
    _combine(combined, right)
    if len(combined) != 1:
        return None
    node = deepcopy(long)
    for key, value in short.get('link', {}).items():
        if not node.setdefault('link', {}).get(key) or (key == 'value' and node['link'][key] == node.get('name')):
            node['link'][key] = value
    for comment in short.get('comments', []):
        key = str(comment.get('name', '')).upper()
        key = _ALIASES.get(key, key)
        if not left.get(key):
            empty = next((c for c in node.get('comments', []) if not c.get('value')
                          and _ALIASES.get(str(c.get('name', '')).upper(), str(c.get('name', '')).upper()) == key), None)
            if empty is not None:
                empty['value'] = comment['value']
                continue
            node.setdefault('comments', []).append(deepcopy(comment))
    for key, value in short.items():
        if key not in ('link', 'comments') and key not in node: node[key] = deepcopy(value)
    return node


def _complete_branches(paths, sample_ids, raw_inventories=None):
    """Suppress only fully represented sparse prefixes or result suffixes."""
    paths = deepcopy(paths)
    raw_records = [[_record(s) for s in p['steps'] if s['kind'] == 'array_data_file'] for p in paths]
    buckets = {}
    bound = {}
    for i, path in enumerate(paths):
        refs = {s['sample_ref'] for s in path.get('steps', []) if s.get('sample_ref')}
        if len(refs) != 1 or not refs <= sample_ids: continue
        bound[i] = next(iter(refs))
        for j, step in enumerate(path['steps']):
            if step['kind'] in ('array_data_file', 'derived_array_data_file'):
                buckets.setdefault((bound[i], step['kind'], step.get('name')), []).append((i, j))
    removed = set()
    for i in sorted(bound, key=lambda i: len(paths[i]['steps']), reverse=True):
        short = paths[i]; steps = short['steps']
        if not steps or steps[-1]['kind'] not in ('array_data_file', 'derived_array_data_file'): continue
        if sum(s['kind'] in ('array_data_file', 'derived_array_data_file') for s in steps) != 1: continue
        candidates = []
        for k, j in buckets.get((bound[i], steps[-1]['kind'], steps[-1].get('name')), []):
            long = paths[k]; full = long['steps']
            if k in removed or len(full) <= len(steps): continue
            if {key:v for key,v in short.items() if key != 'steps'} != {key:v for key,v in long.items() if key != 'steps'}: continue
            prefix = j == len(steps)-1 and steps[:-1] == full[:j]
            suffix = (steps[-1]['kind'] == 'derived_array_data_file' and len(steps) > 2 and j == len(full)-1
                      and steps[0] == full[0] and steps[1:-1] == full[j-len(steps)+2:j])
            if not (prefix or suffix): continue
            node = _complete_file(full[j], steps[-1])
            if node is not None: candidates.append((k,j,node))
        identities = {paths[k]['steps'][j].get('link', {}).get('value') for k,j,_ in candidates}
        identities -= {None, '', steps[-1].get('name')}
        # An unlocated filename cannot choose between distinct file versions.
        if len(identities) > 1: continue
        if any(_complete_file(a, b) is None for _,_,a in candidates for _,_,b in candidates): continue
        if candidates:
            for k,j,node in candidates:
                previous = _record(paths[k]['steps'][j])
                originals = [r for r in raw_records[k] if r == previous]
                if len(originals) == 1 and previous.get('URI') == previous.get('NAME') and node.get('link', {}).get('value'):
                    originals[0]['URI'] = node['link']['value']
                paths[k]['steps'][j] = node
                raw_records[k].extend(deepcopy(raw_records[i]))
            removed.add(i)
    retained = [p for i,p in enumerate(paths) if i not in removed]
    if raw_inventories is not None:
        raw_inventories.update({id(p): raw_records[i] for i,p in enumerate(paths) if i not in removed})
    return retained


def project_native_files(data):
    """Consolidate native raw-file paths, keeping run-scoped archive annotations."""
    if str(data.get('source', {}).get('format', '')).upper() not in {'ENA', 'SRA'}:
        return
    series = data.get('series', {})
    original_inventory = {}
    paths = _complete_branches(series.get('assay_paths', []), {s['iid'] for s in data.get('sample', [])}, original_inventory)
    from urllib.parse import unquote, urlsplit
    for sample in data.get('sample', []):
        sources = [s for p in paths for s in p['steps'] if s['kind'] in ('source', 'sample') and s.get('sample_ref') == sample['iid']]
        for link in sample.get('supplementary_data', []):
            if link.get('type') != 'assembly report' or not link.get('value') or not sources:
                continue
            if any(s.get('link', {}).get('value') == link['value'] for p in paths for s in p['steps']
                   if any(n.get('sample_ref') == sample['iid'] for n in p['steps'])):
                continue
            node = {'kind': 'derived_array_data_file', 'name': unquote(urlsplit(link['value']).path.rsplit('/', 1)[-1]), 'link': deepcopy(link)}
            annotation = _comments(_record(node), 'SAMPLE_FILE_', fixed=True)
            if any(any(s.get('comments', [])[i:i+len(annotation)] == annotation for i in range(len(s.get('comments', [])))) for s in sources):
                continue
            # A sample-level report cannot inherit an arbitrary aliquot's facts.
            source = deepcopy(sources[0]) if all(s == sources[0] for s in sources) else {
                'kind': 'source', 'name': sample['iid'], 'sample_ref': sample['iid']}
            paths.append({'steps': [source, node]})
    from .native_file_selection import source_annotations, primary_sets, _coalesced
    paths = source_annotations(paths, {s['iid'] for s in data.get('sample', [])}, _record, _comments)
    groups, order, archives, inventories = {}, [], {}, {}
    for original in paths:
        path = deepcopy(original)
        steps = path.get('steps', [])
        scans = [s for s in steps if s.get('kind') == 'scan']
        raw = [s for s in steps if s.get('kind') == 'array_data_file']
        if len(scans) != 1 or not raw or not scans[0].get('name'):
            # Already projected, sample-scoped, or ambiguous run association.
            order.append((None, original))
            continue
        scan = scans[0]
        scope = _key([(s.get('kind'), s.get('sample_ref'), s.get('name')) for s in steps
                      if s.get('kind') in {'source', 'sample', 'assay', 'scan'}])
        pool = archives.setdefault(scope, [])
        for record in _scan_records(scan):
            record['ROLE'] = 'SUBMISSION_FILE'
            pool.append(deepcopy(record))
        path['steps'] = [s for s in steps if s.get('kind') != 'array_data_file']
        identity = _key(path)
        if identity not in groups:
            groups[identity] = (path, scope, [])
            order.append((identity, None))
        files = groups[identity][2]
        for step in raw:
            record = _record(step)
            files.append(deepcopy(record))
        inventories.setdefault(identity, []).extend(deepcopy(original_inventory.get(id(original), [_record(s) for s in raw])))
    # Complete aliases before deciding their representation: an untyped copy of
    # an explicitly identified FASTQ is not a separate archival alternative.
    for path, scope, files in groups.values():
        files[:] = _coalesced(files, _combine)
        for record in files:
            if not _fastq(record):
                archives[scope].append(deepcopy(record))
    selected = primary_sets(groups, _key, _fastq, _combine, archives, inventories)
    for scope, values in archives.items():
        archives[scope] = _coalesced(values, _combine)
    result = []
    for identity, original in order:
        if identity is None:
            result.append(original)
            continue
        path, scope, files = groups[identity]
        if identity in selected and selected[identity] is None:
            continue
        reads = selected.get(identity, [record for record in files if _fastq(record)])
        for read in reads or [{}]:
            row = deepcopy(path)
            scan = next(s for s in row['steps'] if s.get('kind') == 'scan')
            comments = scan.setdefault('comments', [])
            comments.extend(_comments(read, 'FASTQ_'))
            for alternative in read.get('_alternatives', []):
                comments.extend(_comments(alternative, 'FASTQ_ALTERNATIVE_', fixed=True))
            # Use the established URI header and an explicit filename label.
            for comment in comments:
                if comment['name'] == 'FASTQ_NAME':
                    comment['name'] = 'FASTQ_FILE_NAME'
            if not read:
                comments.append({'name': 'FASTQ_URI', 'value': ''})
            for record in archives[scope]:
                prefix = 'SUBMITTED_FILE_' if record.get('ROLE', '').upper() in ('SUBMISSION_FILE', 'ORIGINAL') else 'ARCHIVE_FILE_'
                comments.extend(_comments(record, prefix, fixed=True))
            result.append(row)
    if paths:
        # Repeated, explicitly bound sample-result relationships do not imply
        # additional experiments. Do not match unbound sources by their names.
        seen = set()
        sample_ids = {s['iid'] for s in data.get('sample', [])}
        unique = []
        for path in result:
            bound = any(s.get('sample_ref') in sample_ids for s in path.get('steps', []))
            identity = _key(path)
            if not bound or identity not in seen:
                unique.append(path)
            if bound:
                seen.add(identity)
        series['assay_paths'] = unique
