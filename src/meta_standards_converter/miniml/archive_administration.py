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
_BIOLOGICAL_IDS = {'sample name', 'alias'}
_NAMES = {'insdc status', 'submitter id'} | _CENTERS | _IDS | _COMMENTS | _BIOLOGICAL_IDS


def administrative_destination(item, provider, sample=None):
    """The same explicit field contract is used by mapping and residual matching."""
    name = str(item.get('name') or '').strip().lower()
    value = str(item.get('value') or '')
    annotations = {k:deepcopy(v) for k,v in item.items() if k not in {'name','value'}}
    verified_name = name == 'sample name' and sample is not None and value in {
        sample.get('iid'), *[a.get('value') for a in sample.get('accession', [])]}
    if (name in _IDS or verified_name) and re.fullmatch(r'(?:SAM(?:N|EA|D)\d+|[SED]RS\d+)', value):
        from .insdc_support import database_for
        return 'accession', {'value':value, 'database':database_for(value), 'label':item['name'], **annotations}
    if name == 'alias':
        namespace = re.fullmatch(r'(GSE\d+|E-(?:MTAB|GEOD)-\d+):(.+)', value)
        if namespace:
            return 'relation', {'type':item['name'], 'target':value, 'namespace':namespace[1], **annotations}
        return 'comments', deepcopy(item)
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
        def clean(node, *, biological=True):
            wanted = Counter()
            kept = []
            for item in node.get('characteristics', []):
                name = item.get('name', '').strip().lower()
                recognized = (name in _CENTERS or name == 'insdc status'
                              or administrative_destination(item, provider, sample) is not None)
                if recognized and (biological or name not in _BIOLOGICAL_IDS):
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
                    wanted |= clean(step, biological=step.get('kind') in {'source', 'sample'})
        for literal, count in wanted.items():
            item = json.loads(literal)
            destination = administrative_destination(item, provider, sample)
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
            occurrences = [o for o in organization.get('source_occurrences') or [organization]
                           if o.get('sample_accession') in accessions]
            if not occurrences:
                continue
            types = set()
            for occurrence in occurrences:
                roles = [occurrence.get('role'), *occurrence.get('roles', [])]
                if not any(roles):
                    # Early coalesced packages omitted roles on occurrences.
                    # Only a sole aggregate role is safe to recover at this scope.
                    aggregate = [r for r in [organization.get('role'), *organization.get('roles', [])] if r]
                    if len(aggregate) == 1:
                        roles = aggregate
                for role in roles:
                    literal = role.get('value') if isinstance(role, dict) else role
                    literal = str(literal or '').strip().lower()
                    if literal == 'owner':
                        types.add('sample owner')
                    elif literal in {'center', 'centre', 'center_name', 'centre_name', 'insdc center name', 'insdc center alias'}:
                        types.add('archive center')
                    elif literal in {'broker', 'broker name'}:
                        types.add('archive broker')
            relations = sample.setdefault('relation', [])
            old = {'type': 'archive center', 'target': organization['iid']}
            # Repair the exact legacy projection, without discarding independently
            # annotated relationships or applying another sample's centre role.
            if 'archive center' not in types:
                relations[:] = [r for r in relations if r != old]
            for kind in sorted(types or {'organization'}):
                relation = {'type': kind, 'target': organization['iid']}
                if relation not in relations:
                    relations.append(relation)
    declared = {d['iid'] for d in data.get('database', [])}
    for sample in data.get('sample', []):
        for status in sample.get('status', []):
            db = status.get('database')
            if db and db not in declared:
                data.setdefault('database', []).append({'iid': db, 'name': db})
                declared.add(db)
