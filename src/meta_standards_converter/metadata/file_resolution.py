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
