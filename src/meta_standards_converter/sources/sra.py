# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Native Entrez discovery and retrieval; never implicitly calls ENA."""
import xml.etree.ElementTree as ET
from .archive_support import ArchiveHTTP, Resolution, StudyRecords, StudySeed, accession_kind, attempt, chunks, identifier


class SRASource:
    base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
    page_size = 10000

    def __init__(self, http=None, requester=None, resource_profile='standard', evidence_dir=None):
        self.http = http or ArchiveHTTP('ncbi_eutils', requester, resource_profile, evidence_dir)
        self._projects = {}

    def search(self, term, db='sra'):
        ids, issues, start, count = [], [], 0, None
        history = {}
        while count is None or start < count:
            try:
                params = {'db': db, 'term': term, 'retmode': 'json', 'usehistory': 'y',
                          'retstart': start, 'retmax': self.page_size, **history}
                result = self.http.get(self.base + 'esearch.fcgi', params, 'json')['esearchresult']
                if not history and result.get('webenv') and result.get('querykey'):
                    history = {'WebEnv': result['webenv'], 'term': '#' + str(result['querykey'])}
                if result.get('errorlist', {}).get('fieldsnotfound'):
                    raise ValueError('Unsupported Entrez search field')
                count = int(result['count'])
                page = result.get('idlist', [])
                if not page and start < count:
                    raise ValueError('Entrez inventory ended early')
                ids.extend(page)
                start += len(page)
            except Exception as error:
                issues.append(f'{db} inventory: {type(error).__name__}')
                break
        unique = list(dict.fromkeys(ids))
        if len(unique) != len(ids):
            issues.append(f'{db} inventory contains repeated identifiers')
        return unique, issues

    def xml(self, db, ids):
        return self.http.get(self.base + 'efetch.fcgi', {'db': db, 'id': ','.join(ids), 'retmode': 'xml'})

    def project_xml(self, accession, result):
        """Resolve the accession namespace before EFetch, then verify the record."""
        if accession in self._projects:
            return self._projects[accession]
        ids = [accession]
        if not accession.isdigit():
            ids, issues = self.search(accession + '[PRJA]', db='bioproject')
            result.issues.extend(issues)
        if len(ids) != 1:
            result.issues.append(f'BioProject {accession}: missing or ambiguous UID')
            return None
        root = attempt(result, f'BioProject {accession}', lambda: self.xml('bioproject', ids))
        if root is None:
            return None
        nodes = root.findall('.//ProjectID/ArchiveID')
        valid = [n for n in nodes if n.get('id') == ids[0] and n.get('accession')
                 and (accession.isdigit() or n.get('accession') == accession)]
        if len(nodes) != 1 or len(valid) != 1 or any(n.tag.lower() == 'error' for n in root.iter()):
            result.issues.append(f'BioProject {accession}: missing, error or mismatched identity')
            return None
        self._projects[accession] = root
        self._projects[ids[0]] = root
        self._projects[valid[0].get('accession')] = root
        return root

    def linked_xml(self, db, ids, result):
        """Keep only requested, identifiable linked entities from a batch."""
        root = attempt(result, db, lambda: self.xml(db, ids))
        if root is None:
            return None
        paths = {'biosample': './/BioSample', 'taxonomy': './/Taxon', 'pubmed': './/PubmedArticle'}
        clean = ET.Element(root.tag)
        found = set()
        for node in root.findall(paths[db]):
            if any(n.tag.lower() == 'error' for n in node.iter()):
                continue
            aliases = ({node.get('accession'), node.get('id')} if db == 'biosample' else
                       {node.findtext('TaxId')} if db == 'taxonomy' else
                       {node.findtext('MedlineCitation/PMID')})
            matched = set(ids) & aliases
            if matched:
                clean.append(node)
                found.update(matched)
        for missing in sorted(set(ids) - found):
            result.issues.append(f'{db} {missing}: missing, error or mismatched identity')
        return clean if len(clean) else None

    def links(self, dbfrom, db, ids, name):
        root = self.http.get(self.base + 'elink.fcgi', {'dbfrom': dbfrom, 'db': db,
            'id': ','.join(ids), 'linkname': name, 'cmd': 'neighbor'})
        return list(dict.fromkeys(n.text for n in root.findall('.//LinkSetDb/Link/Id') if n.text))

    def resolve(self, accession):
        accession, kind = accession_kind(accession)
        result = Resolution()
        if kind == 'study':
            result.studies.append(StudySeed(accession, accession))
            return result
        queries, visited = [], set()
        def project(acc):
            if acc in visited:
                return
            visited.add(acc)
            root = self.project_xml(acc, result)
            if root is None:
                queries.append(f'{acc}[GPRJ]')
                return
            if root.find('.//ProjectTypeTopAdmin') is None:
                queries.append(f'{acc}[GPRJ]')
                return
            node = root.find('.//ProjectID/ArchiveID')
            uid = node.get('id') if node is not None else acc
            children = attempt(result, f'BioProject children {acc}', lambda: self.links('bioproject', 'bioproject', [uid], 'bioproject_bioproject_u2d')) or []
            for child in children:
                doc = self.project_xml(child, result)
                if doc is not None:
                    for item in doc.findall('.//ProjectID/ArchiveID'):
                        if item.get('accession'):
                            project(item.get('accession'))
        if kind == 'project':
            project(accession)
        else:
            queries.append(f'{accession}[{"BSPL" if kind == "biosample" else "ACCN"}]')
        for query in queries:
            ids, issues = self.search(query)
            result.issues.extend(issues)
            for batch in chunks(ids):
                root = attempt(result, 'SRA accession resolution', lambda: self.xml('sra', batch))
                if root is None:
                    continue
                for node in root.findall('.//EXPERIMENT/STUDY_REF'):
                    acc = identifier(node)
                    if acc and StudySeed(acc, acc) not in result.studies:
                        result.studies.append(StudySeed(acc, acc))
        if not result.studies:
            result.issues.append(f'{accession}: no public read study resolved')
        return result

    def fetch(self, seed):
        records = StudyRecords(seed)
        ids, records.issues = self.search(f'{seed.study}[ACCN]')
        for batch in chunks(ids):
            root = attempt(records, f'{seed.study} experiments', lambda: self.xml('sra', batch))
            if root is not None:
                records.xml.append(root)
                if len(root.findall('.//EXPERIMENT_PACKAGE')) != len(batch):
                    records.issues.append(f'{seed.study}: incomplete experiment package batch')
        observed = {identifier(n) for root in records.xml for n in root.findall('.//EXPERIMENT') if identifier(n)}
        if len(observed) != len(ids):
            records.issues.append(f'{seed.study}: experiment inventory mismatch')
        if not ids:
            records.issues.append(f'{seed.study}: no public experiment records')
        self._pool_samples(records)
        linked = {'biosample': set(), 'bioproject': set(), 'pubmed': set(), 'taxonomy': set()}
        for root in records.xml:
            for node in root.findall('.//EXTERNAL_ID'):
                key = (node.get('namespace') or '').lower()
                if key in linked and node.text:
                    linked[key].add(node.text.strip())
            for node in root.findall('.//XREF_LINK'):
                key, value = node.findtext('DB', '').lower(), node.findtext('ID')
                if key in linked and value:
                    linked[key].add(value)
            linked['taxonomy'].update(n.text for n in root.findall('.//SAMPLE_NAME/TAXON_ID') if n.text)
        accepted_projects = set()
        for db in ('biosample', 'bioproject'):
            for batch in chunks(sorted(linked[db]), size=1 if db == 'bioproject' else 100):
                root = self.project_xml(batch[0], records) if db == 'bioproject' else self.linked_xml(db, batch, records)
                if root is not None and db == 'bioproject':
                    uid = root.find('.//ProjectID/ArchiveID').get('id')
                    if uid in accepted_projects:
                        continue
                    accepted_projects.add(uid)
                if root is not None:
                    records.xml.append(root)
                    if db == 'bioproject':
                        linked['pubmed'].update(n.get('id') for n in root.findall('.//Publication') if n.get('id') and n.findtext('DbType') == 'ePubmed')
        for db in ('taxonomy', 'pubmed'):
            for batch in chunks(sorted(linked[db])):
                root = self.linked_xml(db, batch, records)
                if root is not None:
                    records.xml.append(root)
        assembly_ids = []
        for batch in chunks(ids):
            assembly_ids.extend(attempt(records, 'linked assemblies', lambda: self.links('sra', 'assembly', batch, 'sra_assembly')) or [])
        project_uids = {n.get('id') for root in records.xml for n in root.findall('.//ProjectID/ArchiveID') if n.get('id')}
        for batch in chunks(sorted(project_uids)):
            assembly_ids.extend(attempt(records, 'project assemblies', lambda: self.links('bioproject', 'assembly', batch, 'bioproject_assembly_all')) or [])
        for batch in chunks(list(dict.fromkeys(assembly_ids))):
            data = attempt(records, 'assembly metadata', lambda: self.http.get(self.base + 'esummary.fcgi',
                {'db': 'assembly', 'id': ','.join(batch), 'retmode': 'json', 'report': 'full'}, 'json'))
            if data is not None:
                for uid in batch:
                    item = data.get('result', {}).get(uid)
                    if item and not item.get('error') and item.get('assemblyaccession') and str(item.get('uid', uid)) == uid:
                        records.linked.append({'provider': 'sra', 'kind': 'assembly', 'accession': item.get('assemblyaccession'), 'metadata': item})
                    else:
                        records.issues.append(f'assembly {uid}: missing summary')
        return records

    def _pool_samples(self, records):
        present = {identifier(n) for root in records.xml for n in root.findall('.//SAMPLE')}
        required = {identifier(n) for root in records.xml for path in ('.//Pool/Member', './/SAMPLE_DESCRIPTOR/POOL/*') for n in root.findall(path) if identifier(n)}
        for accession in sorted(required - present):
            ids, issues = self.search(accession + '[ACCN]')
            records.issues.extend(issues)
            found = False
            for batch in chunks(ids):
                root = attempt(records, f'pool sample {accession}', lambda: self.xml('sra', batch))
                if root is not None:
                    # Only explicitly referenced sample records enter this dataset.
                    sample_set = ET.Element('SAMPLE_SET')
                    for node in root.findall('.//SAMPLE'):
                        aliases = {identifier(node)} | {n.text for n in node.findall('IDENTIFIERS/*')}
                        if accession in aliases and identifier(node) not in present:
                            sample_set.append(node)
                            present.add(identifier(node))
                            found = True
                    if len(sample_set):
                        records.xml.append(sample_set)
            if not found:
                records.issues.append(f'{accession}: unresolved pool sample')
