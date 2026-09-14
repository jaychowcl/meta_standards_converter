# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Source-bound archive actors and ontology declarations (no retrieval)."""
from copy import deepcopy
import logging
import re


def _address(node):
    value = {}
    lines = []
    for child in node:
        text = ''.join(child.itertext()).strip()
        if not text: continue
        if child.tag in ('City','Country'): value[child.tag.lower()] = text
        else: lines.append(text)
    if lines: value['line'] = lines
    if node.get('postal_code'): value['postal_code'] = node.get('postal_code')
    return value


def actors(records, provider):
    from .insdc_support import text
    organizations, contributors = [], []
    seen = set()
    from ..sources.archive_support import identifier
    def organization_nodes(node, owner, path=()):
        accession = identifier(node) or (node.get('uid') if node.tag == 'DocumentSummary' else None)
        if node.tag == 'EXPERIMENT_PACKAGE': accession = identifier(node.find('EXPERIMENT'))
        if accession:
            owner, path = (node.tag, accession), ()
        if node.tag == 'Organization' or (node.tag == 'Owner' and owner[0] == 'BioSample'):
            yield node, owner, path
        for index, child in enumerate(node):
            yield from organization_nodes(child, owner, (*path, index))
    for ri, root in enumerate(records.xml):
        for org, owner, path in organization_nodes(root, (root.tag, str(ri))):
            identity = (owner, path)
            if identity in seen:
                continue
            seen.add(identity)
            iid = f'{provider}:{owner[0]}:{owner[1]}:organization-' + '-'.join(map(str, path))
            value = {'iid': iid, 'name': text(org, 'Name')}
            sample_accession = owner[1] if owner[0] in ('BioSample', 'SAMPLE') else None
            if sample_accession: value['sample_accession'] = sample_accession
            if org.tag == 'Owner':
                value['role'] = 'owner'
                if text(org, 'Name') and org.find('Name').get('url'):
                    value['web_link'] = org.find('Name').get('url')
            for source, target in [('url', 'web_link'), ('role', 'role'), ('type', 'type')]:
                if org.get(source): value[target] = org.get(source)
            address = org.find('Address')
            if address is not None: value['address'] = _address(address)
            organizations.append(value)
            contacts = org.findall('Contacts/Contact') if org.tag == 'Owner' else org.findall('Contact')
            for ci, contact in enumerate(contacts):
                person = {k: text(contact, 'Name/' + v) for k, v in [('first','First'), ('middle','Middle'), ('last','Last')]}
                item = {'iid': f'{iid}:contact-{ci}', 'organization_ref': {'ref': iid}}
                if any(person.values()): item['person'] = {k:v for k,v in person.items() if v}
                email = contact.get('email')
                if email and re.fullmatch(r'[^\s@]+@[^\s@]+\.[^\s@]+', email): item['email'] = email
                secondary = contact.get('sec_email')
                if secondary and re.fullmatch(r'[^\s@]+@[^\s@]+\.[^\s@]+', secondary): item.setdefault('extensions', {})['secondary_email'] = secondary
                for source, target in [('phone','phone'), ('fax','fax'), ('url','web_link')]:
                    if contact.get(source): item[target] = contact.get(source)
                if contact.get('role'): item['roles'] = [{'value': contact.get('role')}]
                address = contact.find('Address')
                if address is not None: item['address'] = _address(address)
                if set(item) - {'iid', 'organization_ref', 'sample_accession'}: contributors.append(item)
    centers = set()
    for root in records.xml:
        for node in root.iter():
            if node.tag not in ('STUDY','PROJECT','SAMPLE') or not node.get('center_name'):
                continue
            iid = f'{provider}:{node.tag}:{identifier(node)}:center_name'
            if iid not in centers:
                centers.add(iid)
                organizations.append({'iid':iid, 'name':node.get('center_name'), 'role':'center_name',
                                      **({'sample_accession':identifier(node)} if node.tag == 'SAMPLE' else {})})
    return organizations, contributors


def declare_ontologies(data, issues=None):
    """Resolve only namespaces explicit in identifiers; preserve disagreements."""
    declarations = {v['iid'] for v in data.get('database', [])}
    residuals = []
    def visit(value, path=()):
        if isinstance(value, list):
            for i, v in enumerate(value): visit(v, (*path, i))
        elif isinstance(value, dict):
            term = value.get('term_accession_number', '')
            match = re.fullmatch(r'https?://purl.obolibrary.org/obo/([A-Za-z][A-Za-z0-9]*)_\d+', str(term))
            if match and value.get('term_source_ref') != match[1]:
                if value.get('term_source_ref'):
                    residuals.append({'provider': data.get('source', {}).get('format'), 'kind': 'annotation',
                        'accession': data['series']['iid'], 'metadata': {'path': list(path), 'annotation': deepcopy(value)}})
                    message = f'{data["series"]["iid"]}: corrected explicit ontology namespace {value["term_source_ref"]} to {match[1]}'
                    logging.getLogger(__name__).warning(message)
                value['term_source_ref'] = match[1]
            publication_source = value.get('status_term_source_ref')
            if publication_source and publication_source not in declarations:
                declarations.add(publication_source)
                data.setdefault('database', []).append({'iid': publication_source, 'name': publication_source})
            if value.get('term_source_ref') and value['term_source_ref'] not in declarations:
                name = value['term_source_ref']; declarations.add(name)
                data.setdefault('database', []).append({'iid': name, 'name': name})
            for key, v in list(value.items()):
                if key not in ('extensions', 'database'): visit(v, (*path, key))
    visit(data)
    if residuals:
        data.setdefault('extensions', {}).setdefault('insdc', {'version':'1.0','records':[]})['records'].extend(residuals)


def local_platforms(data):
    """Declare local sequencing platforms from exact explicit instrument facts."""
    if data.get('source', {}).get('format') not in ('ENA', 'SRA'):
        return
    import hashlib
    import json
    from ..metadata.archive_enrichment import informative
    provider = data['source']['format'].lower()
    platforms = data.setdefault('platform', [])
    known = {p['iid'] for p in platforms}
    by_experiment = {}
    for sample in data.get('sample', []):
        refs = []
        for run in sample.get('sra_run', []):
            model = run.get('instrument_model')
            if not informative(model):
                refs.append(None)
                continue
            facts = {'instrument_model': model}
            if informative(run.get('instrument_platform')):
                facts['instrument_platform'] = run['instrument_platform']
            iid = provider + ':platform:' + hashlib.sha256(json.dumps(facts, sort_keys=True).encode()).hexdigest()[:16]
            if iid not in known:
                platforms.append({'iid': iid, 'title': model, 'technology': 'high-throughput sequencing', **facts})
                known.add(iid)
            run['platform_ref'] = {'ref': iid}
            refs.append(iid)
            by_experiment.setdefault((sample['iid'], run.get('experiment')), set()).add(iid)
        old = sample.get('platform_ref', {}).get('ref', '')
        if not old or old.startswith(provider + ':platform:'):
            if refs and None not in refs and len(set(refs)) == 1:
                sample['platform_ref'] = {'ref': refs[0]}
            else:
                sample.pop('platform_ref', None)
    for path in data.get('series', {}).get('assay_paths', []):
        for step in path['steps']:
            if step['kind'] == 'assay':
                refs = by_experiment.get((step.get('sample_ref'), step.get('name')), set())
                if len(refs) == 1:
                    platform = next(p for p in platforms if p['iid'] == next(iter(refs)))
                    if not any(c.get('name') == 'Platform_title' for c in step.get('comments', [])):
                        step.setdefault('comments', []).append({'name': 'Platform_title', 'value': platform['title']})


def coalesce_organizations(data):
    """Share identical organization facts, retaining independent source occurrences."""
    if data.get('source', {}).get('format') not in ('ENA', 'SRA'):
        return
    import json
    groups, mappings = {}, {}
    for organization in data.get('organization', []):
        facts = {k: v for k, v in organization.items() if k not in ('iid', 'sample_accession', 'source_occurrences')}
        # Unknown unnamed entities have no sufficient factual identity to share.
        identity = json.dumps(facts, sort_keys=True) if facts.get('name') or facts.get('alias') else organization['iid']
        occurrences = deepcopy(organization.get('source_occurrences') or [
            {k: v for k, v in organization.items() if k in ('iid', 'sample_accession')}])
        if identity not in groups:
            groups[identity] = {**deepcopy(organization), 'source_occurrences': []}
        target = groups[identity]
        for occurrence in occurrences:
            mappings[occurrence['iid']] = target['iid']
            if occurrence not in target['source_occurrences']:
                target['source_occurrences'].append(occurrence)
        mappings[organization['iid']] = target['iid']
    data['organization'] = list(groups.values())
    for organization in data['organization']:
        occurrences = organization['source_occurrences']
        if len(occurrences) == 1:
            organization.pop('source_occurrences')
        elif len({v.get('sample_accession') for v in occurrences}) > 1:
            organization.pop('sample_accession', None)
    def rewrite(value):
        if isinstance(value, list):
            for child in value: rewrite(child)
        elif isinstance(value, dict):
            for key, child in value.items():
                if key in ('extensions', 'source_occurrences'): continue
                if key in ('ref', 'target') and isinstance(child, str): value[key] = mappings.get(child, child)
                else: rewrite(child)
    rewrite(data)


def repository_databases(data):
    """Separate participating repositories from source vocabulary declarations."""
    if data.get('source', {}).get('format') not in ('ENA', 'SRA'):
        return []
    repositories = {'ENA', 'SRA', 'DRA', 'DDBJ', 'GEO', 'ArrayExpress', 'BioProject', 'BioSample', 'PubMed', 'PMC', 'DOI'}
    used = set()
    def visit(value):
        if isinstance(value, list):
            for child in value: visit(child)
        elif isinstance(value, dict):
            for accession in value.get('accession', []) if isinstance(value.get('accession'), list) else []:
                if isinstance(accession, dict) and accession.get('database'): used.add(accession['database'])
            for key, child in list(value.items()):
                if (key in ('term_source_ref', 'term_accession_number')
                        or key.endswith(('_term_source_ref', '_term_accession_number'))):
                    if isinstance(child, str) and not child.strip():
                        value.pop(key)
                        continue
                if key not in ('extensions', 'database'): visit(child)
    visit(data)
    kept, records = [], []
    for declaration in data.get('database', []):
        iid = declaration['iid']
        if iid != 'INSDC' and (iid in repositories or iid in used):
            kept.append(declaration)
        elif set(declaration) - {'iid', 'name'} or declaration.get('name', iid) != iid:
            records.append({'provider': data['source']['format'], 'kind': 'term_source_declaration',
                            'accession': iid, 'metadata': deepcopy(declaration)})
    data['database'] = kept
    return records
