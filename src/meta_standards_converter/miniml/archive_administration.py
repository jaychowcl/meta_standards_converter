# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Scope supplied administrative characteristics without inventing actors."""
from collections import Counter
from copy import deepcopy
import hashlib
import json
import re

_CENTERS = {'insdc center name', 'insdc center alias', 'broker name'}
_IDS = {'external id', 'sra accession', 'insdc secondary accession'}
_COMMENTS = {'ena-checklist', 'biosamplemodel', 'ncbi submission model', 'ncbi submission package',
             'ena-submission-tool', 'submission model', 'submission package'}
_NAMES = {'insdc status', 'submitter id'} | _CENTERS | _IDS | _COMMENTS


def administrative_destination(item, provider):
    """The same explicit field contract is used by mapping and residual matching."""
    name = str(item.get('name') or '').strip().lower()
    value = str(item.get('value') or '')
    annotations = {k:deepcopy(v) for k,v in item.items() if k not in {'name','value'}}
    if name in _IDS and re.fullmatch(r'(?:SAM(?:N|EA|D)\d+|[SED]RS\d+)', value):
        from .insdc_support import database_for
        return 'accession', {'value':value, 'database':database_for(value), 'label':item['name'], **annotations}
    if name == 'submitter id':
        return 'relation', {'type':item['name'], 'target':value, 'namespace':f'{provider} submitter', **annotations}
    if name in _COMMENTS | _IDS:
        return 'comments', deepcopy(item)
    return None


def normalize_administration(data):
    """Project native sample administration, preserving occurrence counts and references."""
    provider = data.get('source', {}).get('format')
    if provider not in {'ENA', 'SRA'}:
        return
    organizations = data.setdefault('organization', [])
    org_ids = {o['iid'] for o in organizations}
    for sample in data.get('sample', []):
        def clean(node):
            wanted = Counter()
            kept = []
            for item in node.get('characteristics', []):
                if item.get('name', '').strip().lower() in _NAMES:
                    wanted[json.dumps(item, sort_keys=True)] += 1
                else:
                    kept.append(item)
            if 'characteristics' in node:
                node['characteristics'] = kept
            return wanted
        wanted = Counter()
        for channel in sample.get('channel', []):
            wanted += clean(channel)
        for path in data.get('series', {}).get('assay_paths', []):
            refs = {s.get('sample_ref') for s in path['steps']} - {None}
            for step in path['steps']:
                if step.get('sample_ref') == sample['iid'] or (not step.get('sample_ref') and refs == {sample['iid']}):
                    wanted |= clean(step)
        for literal, count in wanted.items():
            item = json.loads(literal)
            destination = administrative_destination(item, provider)
            if destination:
                field, record = destination
                values = sample.setdefault(field, [])
                values.extend(deepcopy(record) for _ in range(max(0, count - values.count(record))))
            elif item['name'].strip().lower() == 'insdc status':
                statuses = sample.setdefault('status', [])
                status = {'database': 'INSDC', 'comment': [item]}
                statuses.extend(deepcopy(status) for _ in range(max(0, count - statuses.count(status))))
            else:
                for occurrence in range(count):
                    iid = f"{provider.lower()}:sample:{sample['iid']}:center:{hashlib.sha256(literal.encode()).hexdigest()[:16]}:{occurrence}"
                    if iid not in org_ids:
                        field = 'alias' if item['name'].strip().lower().endswith('alias') else 'name'
                        organization = {'iid': iid, field: item['value'], 'role': item['name']}
                        extras = {k:v for k,v in item.items() if k not in {'name','value'}}
                        if extras:
                            organization['attribute_annotations'] = extras
                        organizations.append(organization)
                        org_ids.add(iid)
                    relation = {'type': 'archive broker' if item['name'].strip().lower() == 'broker name' else 'archive center', 'target': iid}
                    if relation not in sample.setdefault('relation', []):
                        sample['relation'].append(relation)
        accessions = {sample['iid'], *[a['value'] for a in sample.get('accession', [])]}
        for organization in organizations:
            if organization.get('sample_accession') in accessions:
                relation = {'type': 'archive center', 'target': organization['iid']}
                if relation not in sample.setdefault('relation', []):
                    sample['relation'].append(relation)
    declared = {d['iid'] for d in data.get('database', [])}
    for sample in data.get('sample', []):
        for status in sample.get('status', []):
            db = status.get('database')
            if db and db not in declared:
                data.setdefault('database', []).append({'iid': db, 'name': db})
                declared.add(db)
