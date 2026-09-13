"""Explicit, entity-bound citation identifiers; no literature discovery."""
import re
from urllib.parse import urlsplit, unquote
from .archive_support import identifier, attempt

EUTILS = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'


def citation_identifier(namespace, value):
    namespace = str(namespace or '').strip().casefold()
    value = str(value or '').strip()
    if namespace in ('pubmed', 'epubmed', 'pmid') and re.fullmatch(r'[1-9]\d*', value):
        return {'pubmed_id': value}
    if namespace in ('pmc', 'pmcid', 'epmc') and re.fullmatch(r'PMC\d+(?:\.\d+)?', value, re.I):
        return {'pmcid': value.upper()}
    if namespace in ('doi', 'edoi') and re.fullmatch(r'10\.\d{4,9}/\S+', value):
        return {'doi': value}
    try:
        parsed = urlsplit(value)
    except ValueError:
        return {}
    host, path = parsed.hostname, unquote(parsed.path).strip('/')
    if host in ('doi.org', 'dx.doi.org'): return citation_identifier('doi', path)
    if host == 'pubmed.ncbi.nlm.nih.gov': return citation_identifier('pubmed', path)
    if host in ('www.ncbi.nlm.nih.gov', 'ncbi.nlm.nih.gov') and path.startswith('pubmed/'):
        return citation_identifier('pubmed', path[7:])
    if host in ('europepmc.org', 'www.europepmc.org', 'pmc.ncbi.nlm.nih.gov', 'www.ncbi.nlm.nih.gov'):
        parts = path.split('/')
        if len(parts) >= 3 and parts[:2] == ['abstract', 'MED']: return citation_identifier('pubmed', parts[2])
        for part in parts:
            if re.fullmatch(r'PMC\d+(?:\.\d+)?', part, re.I): return citation_identifier('pmc', part)
    return {}


def references(records):
    """Normalized references retain their source entity, not their discovery path."""
    result = []
    def add(acc, value):
        if value and any(value.get(k) for k in ('pubmed_id','pmcid','doi','title')):
            item = {'accession': acc or records.seed.study, **value}
            if item not in result: result.append(item)
    entity_tags = {'STUDY','PROJECT','Project','SAMPLE','BioSample','EXPERIMENT','RUN','ANALYSIS','ASSEMBLY'}
    def visit(node, acc):
        if node.tag in ('PubmedArticle','PubmedBookArticle'): return
        if node.tag in entity_tags:
            candidate = identifier(node) or node.findtext('ProjectID/ArchiveID')
            archive = node.find('ProjectID/ArchiveID')
            if archive is not None: candidate = archive.get('accession') or candidate
            acc = candidate or acc
        if node.tag == 'XREF_LINK': add(acc, citation_identifier(node.findtext('DB'), node.findtext('ID')))
        if node.tag == 'EXTERNAL_ID': add(acc, citation_identifier(node.get('namespace'), node.text))
        if node.tag == 'URL_LINK': add(acc, citation_identifier('url', node.findtext('URL')))
        if node.tag == 'Link': add(acc, citation_identifier(node.get('type') or node.get('target'), node.text))
        if node.tag == 'Publication':
            value = citation_identifier(node.findtext('DbType'), node.get('id'))
            for source, target in [('Title','title'),('AuthorList','author_list'),('DOI','doi')]:
                if node.findtext(source): value[target] = node.findtext(source)
            add(acc, value)
        for child in node: visit(child, acc)
    for root in records.xml: visit(root, records.seed.study)
    for kind, rows in records.indexed.items():
        for row in rows:
            acc = row.get(kind.removeprefix('read_')+'_accession') or row.get('assembly_set_accession') or records.seed.study
            if row.get('pubmed_id'): add(acc, citation_identifier('pubmed',row['pubmed_id']))
    for record in records.linked:
        acc = record.get('accession')
        if record['kind'] == 'publication_reference': add(acc, record['metadata'])
        elif record['kind'] == 'cross_references':
            for row in record['metadata']:
                value = {}
                source = row.get('Source','').casefold()
                if source in ('europepmc','pubmed','citation'):
                    value.update(citation_identifier('pubmed', row.get('Source Secondary Accession')))
                    value.update(citation_identifier('pmc' if source=='europepmc' else source, row.get('Source Primary Accession')))
                for field in ('Source URL','Source Secondary URL','url'):
                    value.update(citation_identifier('url',row.get(field)))
                add(acc,value)
        elif record['provider'] == 'biosamples' and record['kind'] == 'sample':
            for link in record['metadata'].get('externalReferences',[]):
                add(acc,citation_identifier('url',link.get('url')))
    return result


def study_accessions(records):
    values = {records.seed.study, records.seed.primary}
    for root in records.xml:
        for node in root.iter():
            if node.tag in ('STUDY','PROJECT'): values.add(identifier(node))
            if node.tag == 'ArchiveID' and node.get('accession'): values.add(node.get('accession'))
    return values - {None}


def resolve_identifiers(records, http, provider):
    """Resolve supplied PMCIDs/DOIs, verifying exact identifiers before acceptance."""
    cache = {}
    for ref in references(records):
        if ref.get('pubmed_id'): continue
        kind = 'pmcid' if ref.get('pmcid') else 'doi' if ref.get('doi') else None
        if not kind: continue
        value = ref[kind]
        def resolve():
            if kind == 'pmcid':
                payload = http.get('https://pmc.ncbi.nlm.nih.gov/tools/idconv/api/v1/articles/',
                                   {'ids':value,'idtype':'pmcid','format':'json'}, 'json')
                matches = [r for r in payload.get('records',[]) if r.get('requested-id')==value
                           and r.get('pmcid')==value and citation_identifier('pubmed',r.get('pmid'))]
                if len(matches)!=1: raise ValueError('missing or mismatched PMCID mapping')
                return matches[0]['pmid']
            payload = http.get(EUTILS+'esearch.fcgi', {'db':'pubmed','term':'"'+value+'"[AID]', 'retmode':'json','retmax':100},'json')
            search = payload['esearchresult']; ids = search.get('idlist',[])
            if int(search.get('count',len(ids))) > len(ids): raise ValueError('incomplete DOI lookup')
            if not ids: return None
            root = http.get(EUTILS+'efetch.fcgi', {'db':'pubmed','id':','.join(ids),'retmode':'xml'})
            matches = {a.findtext('MedlineCitation/PMID') for a in root.findall('PubmedArticle')
                       if a.findtext('MedlineCitation/PMID') in ids and any(
                           n.get('IdType')=='doi' and str(n.text).casefold()==value.casefold()
                           for n in a.findall('PubmedData/ArticleIdList/ArticleId'))}
            if len(matches)!=1: raise ValueError('ambiguous or mismatched DOI lookup')
            return matches.pop()
        key = (kind,value)
        if key not in cache: cache[key] = attempt(records, f'{kind} {value}', resolve)
        if cache[key]:
            metadata = {k:v for k,v in ref.items() if k != 'accession'}
            metadata['pubmed_id'] = cache[key]
            records.linked.append({'provider':provider,'kind':'publication_reference','accession':ref['accession'],'metadata':metadata})
