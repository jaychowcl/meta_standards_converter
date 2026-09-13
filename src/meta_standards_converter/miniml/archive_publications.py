"""Pure citation projection, retaining the scope of each explicit source link."""
from copy import deepcopy
from ..sources.archive_publications import references, study_accessions


def project_publications(records, series, samples):
    from .insdc_support import text, database_for
    from ..metadata.ontology_mappings import Harmonizer
    from ..sources.archive_support import identifier
    aliases = study_accessions(records)
    entities = {a['value']:s for s in samples for a in s.get('accession',[])}
    entities.update({r['run']:r for s in samples for r in s.get('sra_run',[])})
    scoped = {identifier(n):n.tag.lower() for root in records.xml for n in root.iter()
              if n.tag in ('EXPERIMENT','ASSEMBLY','ANALYSIS') and identifier(n)}
    scoped.update({r['accession']:r['kind'] for r in records.linked if r['kind'] in ('analysis','assembly')})
    articles = {}
    for root in records.xml:
        for article in root.findall('.//PubmedArticle'):
            pmid = text(article,'MedlineCitation/PMID')
            if not pmid: continue
            authors = [text(a,'CollectiveName') or ' '.join(filter(None,[text(a,'ForeName'),text(a,'LastName')]))
                       for a in article.findall('MedlineCitation/Article/AuthorList/Author')]
            value = {'pubmed_id':pmid, 'title':text(article,'MedlineCitation/Article/ArticleTitle'), 'author_list':', '.join(authors)}
            value.update(zip(('status','status_term_source_ref','status_term_accession_number'),
                             Harmonizer().pubstatus2efo(text(article,'PubmedData/PublicationStatus'))))
            for aid in article.findall('PubmedData/ArticleIdList/ArticleId'):
                if aid.get('IdType')=='doi': value['doi']=aid.text
                if aid.get('IdType')=='pmc': value['pmcid']=aid.text
            articles[pmid] = {k:v for k,v in value.items() if v}
    refs = references(records)
    # A resolved identifier enriches its explicit DOI/PMCID reference as one group.
    for ref in refs:
        if ref.get('pubmed_id'): continue
        matches = [r for r in refs if r['accession']==ref['accession'] and r.get('pubmed_id') and any(
            ref.get(k) and ref[k]==r.get(k) for k in ('doi','pmcid'))]
        if len({r['pubmed_id'] for r in matches})==1: ref['pubmed_id']=matches[0]['pubmed_id']
    for ref in refs:
        acc = ref['accession']
        target = series if acc in aliases else entities.get(acc)
        pub = {k:v for k,v in ref.items() if k!='accession'}
        pub.setdefault('pubmed_id','')
        for key,value in articles.get(pub['pubmed_id'],{}).items():
            if not pub.get(key): pub[key]=value
        if target is not None:
            publications = target.setdefault('pubmed_publication',[])
            existing = next((p for p in publications if any(pub.get(k) and p.get(k)==pub[k] for k in ('pubmed_id','doi','pmcid'))),None)
            if existing is None: publications.append(deepcopy(pub))
            else:
                for key,value in pub.items():
                    if not existing.get(key): existing[key]=value
            if pub['pubmed_id']:
                ids = target.setdefault('pubmed_id',[])
                if pub['pubmed_id'] not in ids: ids.append(pub['pubmed_id'])
            for key, db in [('pubmed_id','PubMed'),('pmcid','PMC'),('doi','DOI')]:
                if pub.get(key):
                    relation = {'type':db,'target':pub[key]}
                    if relation not in target.setdefault('relation',[]):target['relation'].append(relation)
        elif acc in scoped:
            value = pub.get('pubmed_id') or pub.get('doi') or pub.get('pmcid')
            if value:
                relation = {'type':scoped[acc]+' publication','target':value,
                            scoped[acc]+'_ref':acc,'publication':pub}
                if relation not in series.setdefault('relation',[]):series['relation'].append(relation)
    # Cross-reference identifiers enter core before residual pruning/link discovery.
    for record in records.linked:
        if record['kind']!='cross_references' or record['accession'] not in aliases: continue
        for row in record['metadata']:
            source, acc = row.get('Source'), row.get('Source Primary Accession','')
            import re
            if (source=='ArrayExpress' and re.fullmatch(r'E-[A-Z]+-\d+',acc)) or (source=='GEO' and re.fullmatch(r'GSE\d+',acc)):
                value = {'value':acc,'database':database_for(acc)}
                if value not in series.setdefault('accession',[]): series['accession'].append(value)
