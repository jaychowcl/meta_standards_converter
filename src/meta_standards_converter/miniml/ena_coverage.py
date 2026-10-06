# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Study-bound ENA membership and fallback coverage, independent of field fidelity."""
from collections import defaultdict
import re

from ..sources.archive_support import identifier


def scoped_rows(records, kind):
    studies = {records.seed.study, records.seed.primary}
    for row in records.indexed.get(kind, []):
        supplied = {row.get(k) for k in ('study_accession', 'secondary_study_accession') if row.get(k)}
        if supplied and supplied <= studies:
            yield row


def fallback_samples(records):
    groups = defaultdict(list)
    for kind in ('read_experiment', 'read_run'):
        for row in scoped_rows(records, kind):
            acc = row.get('sample_accession') or row.get('secondary_sample_accession')
            if acc and re.fullmatch(r'(?:SAM(?:N|EA|D)|[SED]RS)\d+', acc):
                groups[acc].append(row)
    fields = ('secondary_sample_accession', 'sample_title', 'sample_description',
              'scientific_name', 'tax_id', 'isolation_source', 'host', 'host_tax_id',
              'sex', 'tissue_type', 'cell_type', 'disease', 'developmental_stage')
    for acc, rows in groups.items():
        item = {'sample_accession': acc}
        for field in fields:
            values = {row[field] for row in rows if row.get(field) not in (None, '')}
            if len(values) == 1:
                item[field] = values.pop()
            elif len(values) > 1:
                issue = f'{acc}: ambiguous indexed sample {field}; values retained in source evidence'
                if issue not in records.issues:
                    records.issues.append(issue)
        yield item
    known = set(groups)
    known.update(row.get('secondary_sample_accession') for rows in groups.values() for row in rows)
    for acc in scoped_sample_refs(records):
        if acc not in known:
            known.add(acc)
            yield {'sample_accession': acc}


def scoped_sample_refs(records):
    """Experiment descriptors and explicit pool members establish membership."""
    for root in records.xml:
        for experiment in root.findall('.//EXPERIMENT'):
            if identifier(experiment.find('STUDY_REF')) not in {records.seed.study, records.seed.primary}:
                continue
            descriptor = experiment.find('DESIGN/SAMPLE_DESCRIPTOR')
            if descriptor is None:
                continue
            members = [identifier(n) for n in descriptor.findall('POOL/*') if identifier(n)]
            for acc in members or [identifier(descriptor)]:
                if acc and re.fullmatch(r'(?:SAM(?:N|EA|D)|[SED]RS)\d+', acc):
                    yield acc


def record_coverage(records, full_entities, samples):
    expected = {'sample': set(), 'experiment': set(), 'run': set()}
    for kind in ('read_experiment', 'read_run'):
        for row in scoped_rows(records, kind):
            for entity in expected:
                if row.get(entity+'_accession'):
                    expected[entity].add(row[entity+'_accession'])
    aliases = defaultdict(set)
    for row in fallback_samples(records):
        primary = row['sample_accession']
        aliases[primary].add(primary)
        if row.get('secondary_sample_accession'):
            aliases[row['secondary_sample_accession']].add(primary)
    for sample in samples:
        for accession in sample.get('accession', []):
            aliases[accession['value']].add(sample['iid'])
    for ref in scoped_sample_refs(records):
        owners = aliases.get(ref, set())
        expected['sample'].add(next(iter(owners)) if len(owners)==1 else ref)
    received = {kind: set() for kind in expected}
    for kind, tag in (('sample','SAMPLE'),('experiment','EXPERIMENT'),('run','RUN')):
        for node in full_entities.get(tag, {}).values():
            received[kind].add(identifier(node))
            received[kind].update(n.text for n in node.findall('IDENTIFIERS/*') if n.text)
    model = {'sample': {a['value'] for s in samples for a in s.get('accession', [])},
             'run': {r['run'] for s in samples for r in s.get('sra_run', [])},
             'experiment': {r['experiment'] for s in samples for r in s.get('sra_run', []) if r.get('experiment')}}
    model['sample'].update(s['iid'] for s in samples)
    entities = {kind: {'expected': sorted(ids), 'received': sorted(ids & received[kind]),
                      'fallback': sorted((ids & model[kind])-received[kind]),
                      'unresolved': sorted(ids-model[kind])} for kind, ids in expected.items()}
    complete = all(not(ids-received[kind]) and not(ids-model[kind]) for kind,ids in expected.items())
    records.linked[:] = [r for r in records.linked if not (r.get('provider')=='ena'
        and r.get('kind')=='retrieval_coverage' and r.get('accession')==records.seed.study)]
    if not complete:
        metadata = {'complete': False, 'membership_complete': all(not v['unresolved'] for v in entities.values()),
                    'entities': entities, 'unfinished': [{'stage': 'native_core',
                    'reason': 'full Browser records or model identities unavailable'}]}
        records.linked.append({'provider':'ena', 'kind':'retrieval_coverage',
                               'accession':records.seed.study, 'metadata':metadata})
        issue = f'{records.seed.study}: native identity/metadata coverage incomplete'
        if issue not in records.issues:
            records.issues.append(issue)
