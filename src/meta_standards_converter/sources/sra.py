"""Native Entrez discovery and retrieval; never implicitly calls ENA."""
from .archive_support import ArchiveHTTP, Resolution, StudyRecords, StudySeed, accession_kind, attempt, chunks, identifier


class SRASource:
    base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
    page_size = 10000

    def __init__(self, http=None, requester=None, resource_profile='standard', evidence_dir=None):
        self.http = http or ArchiveHTTP('ncbi_eutils', requester, resource_profile, evidence_dir)

    def search(self, term, db='sra'):
        ids, issues, start, count = [], [], 0, None
        while count is None or start < count:
            try:
                result = self.http.get(self.base + 'esearch.fcgi', {'db': db, 'term': term,
                    'retmode': 'json', 'usehistory': 'y', 'retstart': start, 'retmax': self.page_size}, 'json')['esearchresult']
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
            root = attempt(result, f'BioProject {acc}', lambda: self.xml('bioproject', [acc]))
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
                doc = attempt(result, f'BioProject child {child}', lambda: self.xml('bioproject', [child]))
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
        if not ids:
            records.issues.append(f'{seed.study}: no public experiment records')
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
        for db in ('biosample', 'bioproject'):
            for batch in chunks(sorted(linked[db])):
                root = attempt(records, db, lambda: self.xml(db, batch))
                if root is not None:
                    records.xml.append(root)
                    if db == 'bioproject':
                        linked['pubmed'].update(n.get('id') for n in root.findall('.//Publication') if n.get('id') and n.findtext('DbType') == 'ePubmed')
        for db in ('taxonomy', 'pubmed'):
            for batch in chunks(sorted(linked[db])):
                root = attempt(records, db, lambda: self.xml(db, batch))
                if root is not None:
                    records.xml.append(root)
        assembly_ids = []
        for batch in chunks(ids):
            assembly_ids.extend(attempt(records, 'linked assemblies', lambda: self.links('sra', 'assembly', batch, 'sra_assembly')) or [])
        for batch in chunks(list(dict.fromkeys(assembly_ids))):
            data = attempt(records, 'assembly metadata', lambda: self.http.get(self.base + 'esummary.fcgi',
                {'db': 'assembly', 'id': ','.join(batch), 'retmode': 'json', 'report': 'full'}, 'json'))
            if data is not None:
                for uid in batch:
                    item = data.get('result', {}).get(uid)
                    if item:
                        records.linked.append({'provider': 'sra', 'kind': 'assembly', 'accession': item.get('assemblyaccession'), 'metadata': item})
                    else:
                        records.issues.append(f'assembly {uid}: missing summary')
        return records
