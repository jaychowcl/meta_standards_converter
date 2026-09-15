"""Local publication identifier cleanup; no lookup or citation reassignment."""
from copy import deepcopy
from ..sources.archive_publications import citation_identifier


def valid_pubmed_ids(values):
    result = []
    for value in values:
        pmid = citation_identifier('pubmed', value).get('pubmed_id')
        if pmid and pmid not in result:
            result.append(pmid)
    return result


def clean_publication_identifiers(data):
    records = []
    entities = [data.get('series', {})]
    entities.extend(v for s in data.get('sample', []) for v in [s, *s.get('sra_run', [])])
    provider = data.get('source', {}).get('format')
    for entity in entities:
        identity = entity.get('iid') or entity.get('run') or data.get('series', {}).get('iid')
        def retain(value):
            if value in (None, ''):
                return
            record = {'provider':provider, 'kind':'invalid_publication_identifier', 'accession':identity,
                      'metadata':{'pubmed_id':deepcopy(value)}}
            if record not in records:
                records.append(record)
            if provider not in {'ENA','SRA'}:
                extra = entity.setdefault('unresolved_pubmed_id', [])
                if value not in extra:
                    extra.append(deepcopy(value))
        if 'pubmed_id' in entity:
            values = entity['pubmed_id']
            values = values if isinstance(values, list) else [values]
            for value in values:
                if not valid_pubmed_ids([value]): retain(value)
            entity['pubmed_id'] = valid_pubmed_ids(values)
        pubs = [*entity.get('pubmed_publication', []),
                *[r['publication'] for r in entity.get('relation', []) if isinstance(r.get('publication'), dict)]]
        for publication in pubs:
            value = publication.get('pubmed_id')
            ids = valid_pubmed_ids([value])
            if not ids: retain(value)
            publication['pubmed_id'] = ids[0] if ids else ''
        for relation in entity.get('relation', []):
            if str(relation.get('type', '')).casefold() == 'pubmed' and not valid_pubmed_ids([relation.get('target')]):
                retain(relation.get('target'))
                citation = relation.get('publication')
                if isinstance(citation, dict) and citation not in entity.setdefault('pubmed_publication', []):
                    entity['pubmed_publication'].append(deepcopy(citation))
                siblings = {k:deepcopy(v) for k,v in relation.items() if k != 'publication'}
                records.append({'provider':provider, 'kind':'invalid_publication_identifier', 'accession':identity,
                                'metadata':{'relation':siblings}})
        if 'relation' in entity:
            entity['relation'] = [r for r in entity['relation'] if str(r.get('type', '')).casefold() != 'pubmed'
                                  or valid_pubmed_ids([r.get('target')])]
    return records
