# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Explicit administrative attribute projections for native archive packages."""
from collections import Counter

_RULES = {
    'ena first public': ('ENA', 'release_date'),
    'ena last update': ('ENA', 'last_update_date'),
    'ena-first-public': ('ENA', 'release_date'),
    'ena-last-update': ('ENA', 'last_update_date'),
    'insdc first public': ('INSDC', 'release_date'),
    'insdc last update': ('INSDC', 'last_update_date'),
}


def _date_key(name, value):
    rule = _RULES.get(str(name).strip().lower())
    return (*rule, value) if rule and isinstance(value, str) and value else None


def normalize_archive_dates(data):
    """Move administrative characteristics to sample status on a mutable export copy.

    Repeated source occurrences survive. Channel/path copies use maximum occurrence
    counts so multiple assay paths do not multiply the sample's timestamps.
    Unbound paths are left intact rather than guessing their sample identity.
    """
    if data.get('source', {}).get('format') not in {'ENA', 'SRA'}:
        return
    samples = {s['iid']: s for s in data.get('sample', [])}
    wanted = {iid: Counter() for iid in samples}

    def clean(owner):
        counts = Counter()
        retained = []
        for item in owner.get('characteristics', []):
            key = _date_key(item.get('name'), item.get('value'))
            if key:
                counts[key] += 1
                # Preserve qualified/unknown annotations alongside their status.
                extra = {k: v for k, v in item.items() if k not in {'name', 'value'}}
                if extra:
                    database, field, value = key
                    extra.update(attribute_name=item['name'])
                    annotations.setdefault(key, []).append(extra)
            else:
                retained.append(item)
        if 'characteristics' in owner:
            owner['characteristics'] = retained
        return counts

    for iid, sample in samples.items():
        annotations = {}
        for channel in sample.get('channel', []):
            wanted[iid] += clean(channel)
        for path in data.get('series', {}).get('assay_paths', []):
            refs = {step.get('sample_ref') for step in path.get('steps', [])} - {None}
            for step in path.get('steps', []):
                ref = step.get('sample_ref') or (next(iter(refs)) if len(refs) == 1 else None)
                if ref == iid:
                    wanted[iid] |= clean(step)
        statuses = sample.setdefault('status', [])
        existing = Counter((s.get('database'), field, s[field])
                           for s in statuses for field in ('release_date', 'last_update_date') if s.get(field))
        for (database, field, value), count in wanted[iid].items():
            key = (database, field, value)
            for _ in range(max(0, count - existing[key])):
                statuses.append({'database': database, field: value})
            if annotations.get(key):
                status = next(s for s in statuses if s.get('database') == database and s.get(field) == value)
                saved = status.setdefault('attribute_annotations', [])
                for extra in annotations[key]:
                    if extra not in saved:
                        saved.append(extra)
    declarations = {d['iid'] for d in data.get('database', [])}
    for sample in samples.values():
        for status in sample.get('status', []):
            database = status.get('database')
            if database and database not in declarations:
                data.setdefault('database', []).append({'iid': database, 'name': database})
                declarations.add(database)
