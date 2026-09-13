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
    if set(link) - {'value', 'type'}:
        result['_link'] = {k: v for k, v in link.items() if k not in {'value', 'type'}}
    return result


def _checksums(record):
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


def _combine(records, candidate):
    """Complete explicit identities or checksum-verified, compatible mirrors."""
    for record in records:
        locations = {r.get('URI') for r in [record, *record.get('_alternatives', [])]} - {None, ''}
        identity = bool(candidate.get('URI') and candidate['URI'] in locations)
        identity = identity or (candidate.get('NAME') and candidate.get('NAME') == record.get('NAME')
                                and (not candidate.get('URI') or not record.get('URI')))
        left, right = _checksums(record), _checksums(candidate)
        common = left.keys() & right.keys()
        mirror = (_fastq(record) and _fastq(candidate) and bool(common) and all(left[k] == right[k] for k in common)
                  and candidate.get('NAME') and candidate['NAME'] == record.get('NAME'))
        ignored = {'URI', '_alternatives'}
        if mirror:
            ignored |= {'MD5', 'CHECKSUM', 'CHECKSUM_METHOD'}
        compatible = all(k in ignored or not record.get(k) or not v or record[k] == v
                         for k, v in candidate.items())
        if (identity or mirror) and compatible:
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
    return result


def project_native_files(data):
    """Consolidate native raw-file paths, keeping run-scoped archive annotations."""
    if str(data.get('source', {}).get('format', '')).upper() not in {'ENA', 'SRA'}:
        return
    series = data.get('series', {})
    paths = series.get('assay_paths', [])
    groups, order, archives = {}, [], {}
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
            _combine(pool, record)
        path['steps'] = [s for s in steps if s.get('kind') != 'array_data_file']
        identity = _key(path)
        if identity not in groups:
            groups[identity] = (path, scope, [])
            order.append((identity, None))
        files = groups[identity][2]
        for step in raw:
            _combine(files, _record(step))
    # Complete aliases before deciding their representation: an untyped copy of
    # an explicitly identified FASTQ is not a separate archival alternative.
    for path, scope, files in groups.values():
        for record in files:
            if not _fastq(record):
                _combine(archives[scope], record)
    result = []
    for identity, original in order:
        if identity is None:
            result.append(original)
            continue
        path, scope, files = groups[identity]
        reads = [record for record in files if _fastq(record)]
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
                prefix = 'SUBMITTED_FILE_' if record.get('ROLE') == 'SUBMISSION_FILE' else 'ARCHIVE_FILE_'
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
