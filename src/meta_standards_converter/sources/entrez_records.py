# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Supporting Entrez inventories and identity-checked linked records."""
import xml.etree.ElementTree as ET
from .archive_support import attempt, chunks


class EntrezRecords:
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

    def __init__(self, http, page_size=10000):
        self.http, self.page_size = http, page_size

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

    def linked_xml(self, db, ids, result):
        """Keep only requested, identifiable linked entities from a batch."""
        fetch_ids = list(ids)
        if db == 'biosample':
            accessions = [a for a in ids if not a.isdigit()]
            fetch_ids = [a for a in ids if a.isdigit()]
            if accessions:
                resolved, issues = self.search(' OR '.join(a + '[Accession]' for a in accessions), db=db)
                result.issues.extend(issues)
                fetch_ids.extend(resolved)
            fetch_ids = list(dict.fromkeys(fetch_ids))
            if not fetch_ids:
                result.issues.append('biosample: no requested accession resolved to a UID')
                return None
        root = ET.Element({'biosample': 'BioSampleSet', 'taxonomy': 'TaxaSet', 'pubmed': 'PubmedArticleSet'}[db])
        for batch in chunks(fetch_ids):
            response = attempt(result, db, lambda: self.xml(db, batch))
            if response is not None:
                root.extend(response)
        paths = {'biosample': 'BioSample', 'taxonomy': 'Taxon', 'pubmed': 'PubmedArticle'}
        clean = ET.Element(root.tag)
        found = set()
        for node in root.findall(paths[db]):
            if any(n.tag.lower() == 'error' for n in node.iter()):
                continue
            if db == 'biosample' and node.get('id') not in fetch_ids:
                continue
            aliases = ({node.get('accession'), node.get('id')} |
                       {n.text for n in node.findall('Ids/Id') if n.get('db') == 'BioSample'} if db == 'biosample' else
                       {node.findtext('TaxId')} if db == 'taxonomy' else
                       {node.findtext('MedlineCitation/PMID')})
            matched = set(ids) & aliases
            if matched:
                clean.append(node)
                found.update(matched)
        for missing in sorted(set(ids) - found):
            result.issues.append(f'{db} {missing}: missing, error or mismatched identity')
        return clean if len(clean) else None

