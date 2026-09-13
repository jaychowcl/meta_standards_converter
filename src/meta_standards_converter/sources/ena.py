"""ENA Portal membership discovery and full Browser record retrieval."""
import csv
import io
from .archive_support import ArchiveHTTP, Resolution, StudyRecords, StudySeed, accession_kind, attempt, chunks, identifier


class ENASource:
    portal = 'https://www.ebi.ac.uk/ena/portal/api/'
    browser = 'https://www.ebi.ac.uk/ena/browser/api/xml/'

    def __init__(self, http=None, requester=None, resource_profile='standard', evidence_dir=None):
        self.http = http or ArchiveHTTP('ena_portal', requester, resource_profile, evidence_dir)
        self._catalogues = {}

    def search(self, result, query, fields='all'):
        rows = self.http.get(self.portal + 'search', {'result': result, 'query': query,
            'fields': fields, 'format': 'json', 'limit': 0, 'includeMetagenomes': 'true'}, 'json')
        if not isinstance(rows, list):
            raise ValueError('ENA search must return an array')
        return rows

    def xml(self, accessions):
        return self.http.get(self.browser + ','.join(accessions), {'includeLinks': 'true'})

    def resolve(self, accession):
        accession, kind = accession_kind(accession)
        result, visited, queries = Resolution(), set(), []
        def project(acc):
            if acc in visited:
                return
            visited.add(acc)
            root = attempt(result, acc, lambda: self.xml([acc]))
            if root is not None and root.find('.//UMBRELLA_PROJECT') is not None:
                for child in root.findall('.//CHILD_PROJECT'):
                    if child.get('accession'):
                        project(child.get('accession'))
            else:
                queries.append(f'study_accession="{acc}"')
        if kind == 'project':
            project(accession)
        else:
            field = {'study': 'secondary_study_accession', 'sample': 'secondary_sample_accession',
                     'biosample': 'sample_accession', 'experiment': 'experiment_accession', 'run': 'run_accession'}[kind]
            queries.append(f'{field}="{accession}"')
        for query in queries:
            rows = attempt(result, 'ENA study resolution', lambda: self.search('read_experiment', query, 'study_accession,secondary_study_accession')) or []
            for row in rows:
                study = row.get('secondary_study_accession') or row.get('study_accession')
                primary = row.get('study_accession') or study
                if study and StudySeed(study, primary) not in result.studies:
                    result.studies.append(StudySeed(study, primary))
        if not result.studies:
            result.issues.append(f'{accession}: no public read study resolved')
        return result

    def fetch(self, seed):
        records = StudyRecords(seed)
        query = f'secondary_study_accession="{seed.study}"' if seed.study.startswith(('SRP', 'ERP', 'DRP')) else f'study_accession="{seed.primary}"'
        for accession in dict.fromkeys((seed.study, seed.primary)):
            root = attempt(records, accession, lambda: self.xml([accession]))
            if root is not None:
                records.xml.append(root)
        for kind in ('study', 'sample', 'read_experiment', 'read_run', 'analysis', 'assembly'):
            if kind not in self._catalogues:
                catalogue = attempt(records, f'{kind} catalogue', lambda: self.http.get(self.portal + 'returnFields', {'result': kind}, 'text'))
                if catalogue:
                    self._catalogues[kind] = list(csv.DictReader(io.StringIO(catalogue), delimiter='\t'))
            scoped = f'study_accession="{seed.primary}"' if kind in ('study', 'sample', 'assembly') else query
            rows = attempt(records, f'{kind} inventory', lambda: self.search(kind, scoped))
            if rows is not None:
                records.indexed[kind] = rows
        # File report is a separate contract; do not concatenate the same Portal
        # run projection with itself and accidentally manufacture file duplicates.
        files = attempt(records, 'ENA file report', lambda: self.http.get(self.portal + 'filereport',
            {'accession': seed.study, 'result': 'read_run', 'fields': 'all', 'format': 'json', 'limit': 0}, 'json'))
        if files is not None:
            by_run = {r.get('run_accession'): r for r in records.indexed.get('read_run', [])}
            for row in files:
                previous = by_run.get(row.get('run_accession'), {})
                by_run[row.get('run_accession')] = {**previous, **row}
            records.indexed['read_run'] = list(by_run.values())
        inventory = {}
        for kind, field in [('sample', 'sample_accession'), ('read_experiment', 'experiment_accession'),
                            ('read_run', 'run_accession'), ('analysis', 'analysis_accession'), ('assembly', 'assembly_accession')]:
            accessions = list(dict.fromkeys(r[field] for r in records.indexed.get(kind, []) if r.get(field)))
            inventory[kind] = accessions
            for batch in chunks(accessions):
                root = attempt(records, f'{kind} XML', lambda: self.xml(batch))
                if root is not None:
                    records.xml.append(root)
                    if len(root) < len(batch):
                        records.issues.append(f'{kind}: incomplete XML batch')
        count = attempt(records, 'ENA run count', lambda: self.http.get(self.portal + 'count',
            {'result': 'read_run', 'query': query, 'format': 'json', 'includeMetagenomes': 'true'}, 'json'))
        if count is not None:
            try:
                expected = int(count[0]['count'] if isinstance(count, list) else count.get('count', count))
                if expected != len(inventory['read_run']):
                    records.issues.append('ENA run inventory/count mismatch')
            except (ValueError, TypeError, KeyError, IndexError, AttributeError):
                records.issues.append('ENA run count response unrecognized')
        if not inventory['read_run']:
            records.issues.append(f'{seed.study}: no public run records')
        taxa, pmids = set(), set()
        for root in records.xml:
            taxa.update(n.text for n in root.findall('.//SAMPLE_NAME/TAXON_ID') if n.text)
            pmids.update(n.findtext('ID') for n in root.findall('.//XREF_LINK') if n.findtext('DB', '').lower() == 'pubmed' and n.findtext('ID'))
        for acc in inventory['sample']:
            if acc.startswith('SAM'):
                obj = attempt(records, f'BioSamples {acc}', lambda: self.http.get('https://www.ebi.ac.uk/biosamples/samples/' + acc, fmt='json'))
                if obj is not None:
                    records.linked.append({'provider': 'biosamples', 'kind': 'sample', 'accession': acc, 'metadata': obj})
        for taxid in sorted(taxa):
            obj = attempt(records, f'taxonomy {taxid}', lambda: self.http.get('https://www.ebi.ac.uk/ena/taxonomy/rest/tax-id/' + taxid, fmt='json'))
            if obj is not None:
                records.linked.append({'provider': 'ena', 'kind': 'taxonomy', 'accession': taxid, 'metadata': obj})
        for batch in chunks(sorted(pmids)):
            root = attempt(records, 'PubMed', lambda: self.http.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi',
                {'db': 'pubmed', 'id': ','.join(batch), 'retmode': 'xml'}))
            if root is not None:
                records.xml.append(root)
        return records
