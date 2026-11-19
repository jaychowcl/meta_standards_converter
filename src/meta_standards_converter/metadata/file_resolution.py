# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Pure file interpretation shared by source adapters, reconciliation and export.

These functions interpret detached records. Original source occurrences remain
in MINiML; a resolved download view never replaces the archival inventory.
"""
from collections.abc import Mapping
from copy import deepcopy
from dataclasses import dataclass
import os
import re
import json
from urllib.parse import urlsplit


def structural_uri(value):
    """Remove XML/TSV framing whitespace without changing literal URI contents."""
    return str(value).strip(' \t\r\n') if value is not None else ''


def extension(value):
    value = structural_uri(value)
    name = os.path.basename(urlsplit(value).path or value).lower()
    for suffix in ('.gz', '.zip', '.bz2', '.xz'):
        if name.endswith(suffix):
            name = name[:-len(suffix)]
            break
    return os.path.splitext(name)[1]


def format_name(value):
    name = str(value or '').strip().casefold().lstrip('.')
    aliases = {'fq': 'fastq', 'fq.gz': 'fastq', 'fastq.gz': 'fastq',
               'cel.gz': 'cel', 'gpr.gz': 'gpr', 'tif': 'tiff'}
    return aliases.get(name, name)


def field_value(key, value):
    if key in {'FORMAT', 'type', 'format'}:
        return format_name(value)
    if key in {'URI', 'uri', 'value'}:
        return structural_uri(value)
    if key in {'CHECKSUM_METHOD', 'checksum_method'}:
        return str(value).upper().replace('-', '')
    if key in {'MD5', 'CHECKSUM', 'md5', 'checksum'}:
        return str(value).casefold()
    return value


def classify(value):
    """Explicit role/format precedes a conservative filename fallback."""
    record = value if isinstance(value, Mapping) else {'value': value}
    role = str(record.get('file_role') or record.get('role') or '').casefold()
    kind = format_name(record.get('type') or record.get('format'))
    if role in {'derived', 'processed', 'matrix'} or kind in {'derived', 'processed', 'matrix', 'assembly report'}:
        return 'matrix_or_derived'
    suffix = extension(record.get('value') or record.get('uri') or record.get('filename'))
    known = kind or format_name(suffix)
    if known in {'fastq', 'bam', 'sam', 'cram', 'srf', 'sra'}:
        return 'sequencing_raw'
    if known in {'cel', 'gpr', 'idat', 'exp', 'rpt', 'cab', 'tiff'}:
        return 'array_raw'
    if known in {'mtx', 'h5ad', 'chp'}:
        return 'matrix_or_derived'
    if kind and kind not in {'txt', 'tsv', 'csv', 'h5', 'hdf5', 'supplementary', 'raw'}:
        return 'supplementary'
    # A generic text/container suffix does not establish that a file is derived.
    if suffix.lstrip('.') != known:
        return classify({'format': suffix.lstrip('.')})
    return 'supplementary'


def origin(record):
    role = str(record.get('role') or record.get('ROLE') or '').upper()
    if role in {'ORIGINAL', 'SUBMITTED', 'SUBMISSION_FILE'}:
        return 'submitted'
    if role in {'GENERATED_FILE', 'ARCHIVE_GENERATED'}:
        return 'archive_generated'
    return 'unknown'


def ena_read_role(record, run, layout):
    """Interpret only an ENA-generated filename with its known run and layout."""
    if origin(record) != 'archive_generated' or str(layout or '').upper() != 'PAIRED':
        return {}
    name = record.get('filename') or os.path.basename(urlsplit(structural_uri(record.get('uri'))).path)
    if name == f'{run}.fastq.gz':
        return {'read_role': 'unpaired'}
    for index in (1, 2):
        if name == f'{run}_{index}.fastq.gz':
            return {'read_role': 'paired', 'read_index': str(index)}
    return {}


def checksums(record):
    if any(str(c.get('name', '')).upper() in {'MD5', 'CHECKSUM', 'CHECKSUM_METHOD'} for c in record.get('_extra', [])):
        return {}
    values = {}
    pairs = [('MD5', record.get('MD5'))]
    if record.get('CHECKSUM'):
        pairs.append((record.get('CHECKSUM_METHOD', ''), record['CHECKSUM']))
    for algorithm, value in pairs:
        if not value:
            continue
        algorithm = str(algorithm).upper().replace('-', '')
        size = {'MD5': 32, 'SHA1': 40, 'SHA256': 64, 'SHA512': 128}.get(algorithm)
        if not size or not re.fullmatch(r'[0-9a-fA-F]{%d}' % size, str(value)):
            return {}
        value = str(value).lower()
        if algorithm in values and values[algorithm] != value:
            return {}
        values[algorithm] = value
    return values


def combine(records, candidate, *, preserve_primary=False, rank=None):
    """Complete compatible files within a verified sample/run scope.

An exact bare filename can complete a sparse descriptor in this scope; distinct
URIs require checksum evidence. Callers check pairwise compatibility before
joining components, so a sparse descriptor cannot bridge conflicting files.
"""
    for record in records:
        locations = {structural_uri(r.get('URI')) for r in [record, *record.get('_alternatives', [])]} - {''}
        identity = bool(candidate.get('URI') and structural_uri(candidate['URI']) in locations)
        repository_matches = not (candidate.get('_repository') and record.get('_repository')
                                  and candidate['_repository'] != record['_repository'])
        identity = identity or (repository_matches and candidate.get('NAME') and candidate['NAME'] == record.get('NAME')
                                and (not candidate.get('URI') or not record.get('URI')))
        left, right = checksums(record), checksums(candidate)
        common = left.keys() & right.keys()
        mirror = (format_name(record.get('FORMAT')) == format_name(candidate.get('FORMAT')) == 'fastq'
                  and bool(common) and all(left[k] == right[k] for k in common)
                  and candidate.get('NAME') and candidate['NAME'] == record.get('NAME'))
        ignored = {'URI', '_alternatives', '_repository', '_source'}
        if mirror:
            ignored |= {'MD5', 'CHECKSUM', 'CHECKSUM_METHOD'}
        compatible = not any(left[k] != right[k] for k in common) and all(
            k in ignored or not record.get(k) or not v or field_value(k, record[k]) == field_value(k, v)
            for k, v in candidate.items())
        if not (identity or mirror) or not compatible:
            continue
        if rank and not preserve_primary and rank(candidate) < rank(record) and candidate.get('URI'):
            previous = deepcopy(record)
            record.clear(); record.update(deepcopy(candidate))
            record.update({k: v for k, v in previous.items() if k != '_alternatives' and v and not record.get(k)})
            alternatives = previous.pop('_alternatives', [])
            if previous.get('URI') != record.get('URI'):
                alternatives.insert(0, previous)
            if alternatives:
                record.setdefault('_alternatives', []).extend(alternatives)
            return
        if candidate.get('URI') and locations and structural_uri(candidate['URI']) not in locations:
            record.setdefault('_alternatives', []).append(deepcopy(candidate))
        record.update({k: deepcopy(v) for k, v in candidate.items() if k != '_alternatives' and v and not record.get(k)})
        return
    records.append(deepcopy(candidate))


@dataclass(frozen=True)
class FileAsset:
    """Exact download identity and indices into the unchanged source inventory."""
    uri: str
    checksum: str
    occurrences: tuple[int, ...]


def download_assets(records):
    """Group exact URI/checksum occurrences without merging conflicting evidence."""
    grouped = {}
    for index, record in enumerate(records):
        uri = structural_uri(record.get('uri') or record.get('value'))
        if uri:
            checksum = str(record.get('md5') or record.get('checksum') or '').casefold()
            # Equal URLs do not authorize dropping contradictory size/role or
            # other supplied metadata. Only equivalent complete records share
            # one presentation occurrence.
            facts = json.dumps({k: field_value(k, v) for k, v in record.items()}, sort_keys=True, default=str)
            identity = (uri, checksum, facts)
            grouped.setdefault(identity, []).append(index)
    return tuple(FileAsset(uri, checksum, tuple(indices)) for (uri, checksum, _), indices in grouped.items())


def complete_path_files(data):
    """Complete exact sample/run/URI matches without replacing source occurrences.

    Supplied file annotations keep their order and ownership. All matching
    records must agree on populated facts; neither filename resemblance nor a
    sparse record may bridge conflicting identities. This works for older JSON
    when its original run inventory still contains the missing information.
    """
    from collections import defaultdict
    aliases = {'FILE FORMAT': 'FORMAT', 'FILE SIZE': 'BYTES', 'FASTQ_MD5': 'MD5',
               'FASTQ_BYTES': 'BYTES', 'FILE URI': 'URI', 'FASTQ_URI': 'URI',
               'FASTQ_FILE_NAME': 'NAME', 'CHECKSUM METHOD': 'CHECKSUM_METHOD'}
    fields = {'md5': 'MD5', 'bytes': 'BYTES', 'format': 'FORMAT',
              'read_role': 'READ_ROLE', 'read_index': 'READ_INDEX',
              'checksum': 'CHECKSUM', 'checksum_method': 'CHECKSUM_METHOD'}
    samples = {s['iid']: s for s in data.get('sample', [])}
    inventory = defaultdict(list)
    for sid, sample in samples.items():
        for run in sample.get('sra_run', []):
            for record in [*run.get('files', []), *run.get('fastq_files', [])]:
                uri = structural_uri(record.get('uri') or record.get('url'))
                if not uri:
                    continue
                facts = {target: record[key] for key, target in fields.items() if record.get(key)}
                facts['NAME'] = record.get('filename')
                inventory[sid, run.get('run'), uri].append(facts)
    for path in data.get('series', {}).get('assay_paths', []):
        refs = {s['sample_ref'] for s in path['steps'] if s.get('sample_ref')}
        if len(refs) != 1:
            continue
        sid = next(iter(refs))
        runs = {s.get('name') for s in path['steps'] if s['kind'] == 'scan'}
        runs.update(c['value'] for s in path['steps'] for c in s.get('comments', [])
                    if c['name'].upper() in {'ENA_RUN', 'SRA_RUN'})
        known = {r.get('run') for r in samples.get(sid, {}).get('sra_run', [])}
        runs &= known
        if len(runs) != 1:
            continue
        run = next(iter(runs))
        for node in path['steps']:
            if not node['kind'].endswith('_file'):
                continue
            link = node.get('link', {})
            uri = structural_uri(link.get('value'))
            matched = inventory.get((sid, run, uri), [])
            if not matched:
                continue
            current = defaultdict(set)
            def add(key, value):
                if value not in (None, ''):
                    current[key].add(str(field_value(key, value)))
            add('NAME', node.get('name'))
            add('FORMAT', link.get('type'))
            for key, target in fields.items():
                if key not in ('format', 'checksum'):
                    add(target, link.get(key))
            add('MD5', link.get('checksum'))
            for comment in node.get('comments', []):
                key = str(comment.get('name', '')).upper()
                add(aliases.get(key, key), comment.get('value'))
            combined = {k: set(v) for k, v in current.items()}
            for record in matched:
                for key, value in record.items():
                    if value not in (None, ''):
                        combined.setdefault(key, set()).add(str(field_value(key, value)))
            # Compare MD5 aliases as one assertion before considering completion.
            for facts in [current, combined]:
                if facts.get('CHECKSUM_METHOD') == {'MD5'}:
                    facts.setdefault('MD5', set()).update(facts.get('CHECKSUM', set()))
            if any(len(combined.get(key, set())) > 1 for key in ('NAME', *fields.values())):
                continue
            for key in ('MD5', 'BYTES', 'READ_ROLE', 'READ_INDEX'):
                values = combined.get(key, set())
                if current.get(key) or len(values) != 1:
                    continue
                value = next(iter(values))
                if key == 'MD5' and not re.fullmatch(r'[a-f0-9]{32}', value):
                    continue
                if key == 'BYTES' and not value.isdecimal():
                    continue
                if key in {'READ_ROLE', 'READ_INDEX'} and current.get('READ_TYPE'):
                    continue
                empty = next((c for c in node.get('comments', []) if not c.get('value')
                              and aliases.get(c['name'].upper(), c['name'].upper()) == key), None)
                if empty is not None:
                    empty['value'] = value
                else:
                    node.setdefault('comments', []).append({'name': key, 'value': value})
