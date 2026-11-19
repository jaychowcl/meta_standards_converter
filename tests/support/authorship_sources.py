# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Independent raw-source evidence extraction for the authorship audit.

This module does not import MSC. It does not claim that a literal occurring
somewhere in a response proves ownership or correct scientific interpretation.
"""
from collections import Counter, defaultdict
import hashlib
import io
import json
from pathlib import Path
import re
import tarfile
import xml.etree.ElementTree as ET
from urllib.parse import parse_qs, urlsplit

from tests.support.authorship_fidelity import EQUIVALENT, LOSS, UNRESOLVED, row_identity


def documents(case, mode):
    """Hash-check recorded bodies and retain XML member identity and source URL."""
    case = Path(case)
    selected = case / 'http' / (mode + '-json-selected.manifest.json')
    manifest = json.loads((selected if selected.exists() else case / 'http' / (mode + '.manifest.json')).read_text())
    seen, result = set(), []
    for request in manifest['requests']:
        for attempt in request['attempts']:
            if attempt.get('status') != 200 or 'response' not in attempt:
                continue
            file_id = attempt['response']
            identity = (request['identity']['url'], file_id)
            if identity in seen: continue
            seen.add(identity)
            item = manifest['files'][file_id]
            raw = (case / 'http' / item['path']).read_bytes()
            assert hashlib.sha256(raw).hexdigest() == item['sha256'], 'corrupt frozen body'
            evidence = {'url': identity[0], 'sha256': item['sha256'], 'body': 'http/' + item['path'],
                        'reference_assisted': '/biostudies/' in identity[0] or '/arrayexpress/' in identity[0]}
            blobs = [(None, raw)]
            if raw.startswith(b'\x1f\x8b'):
                with tarfile.open(fileobj=io.BytesIO(raw), mode='r:gz') as archive:
                    blobs = [(member.name, archive.extractfile(member).read()) for member in archive
                             if member.isfile() and member.name.endswith('.xml')]
            for member, content in blobs:
                entry = dict(evidence, member=member)
                try:
                    xml = ET.fromstring(content)
                    for node in xml.iter(): node.tag = node.tag.rsplit('}', 1)[-1]
                    entry.update(kind='xml', xml=xml)
                except (ET.ParseError, UnicodeDecodeError):
                    try: entry.update(kind='json', data=json.loads(content))
                    except (ValueError, UnicodeDecodeError): entry.update(kind='text', text=content.decode('utf-8', errors='replace'))
                result.append(entry)
    snapshots = case / 'source-snapshots.json'
    if snapshots.exists():
        for item in json.loads(snapshots.read_text()):
            raw = (case / item['body']).read_bytes()
            assert hashlib.sha256(raw).hexdigest() == item['sha256'], 'corrupt native snapshot'
            snapshot = json.loads(raw)
            for index, content in enumerate(snapshot['xml']):
                xml = ET.fromstring(content)
                for node in xml.iter(): node.tag = node.tag.rsplit('}', 1)[-1]
                result.append(dict(item, kind='xml', xml=xml, member=f'xml/{index}',
                                   reference_assisted=False))
    return result


def provenance(doc, path):
    return {key: doc[key] for key in ('url', 'sha256', 'body', 'member', 'reference_assisted')} | {'path': path}


def leaf_inventory(docs):
    facts = []
    for doc in docs:
        if doc['kind'] == 'xml':
            def walk(node, path):
                current = path + '/' + node.tag
                identifiers = {key: value for key, value in node.attrib.items() if key in ('iid', 'accession', 'Name', 'tag', 'position')}
                if identifiers: current += '[' + json.dumps(identifiers, sort_keys=True) + ']'
                if node.text and node.text.strip():
                    facts.append({'value': node.text, 'evidence': provenance(doc, current)})
                for name, value in node.attrib.items():
                    facts.append({'value': value, 'evidence': provenance(doc, current + '/@' + name)})
                for child in node: walk(child, current)
            walk(doc['xml'], '')
        elif doc['kind'] == 'json':
            def walk_json(value, path):
                if isinstance(value, dict):
                    for key, child in value.items(): walk_json(child, path + '/' + str(key))
                elif isinstance(value, list):
                    for index, child in enumerate(value): walk_json(child, path + '/' + str(index))
                elif value is not None and value != '': facts.append({'value': value, 'evidence': provenance(doc, path)})
            walk_json(doc['data'], '')
    return facts


def citation_facts(docs):
    """PubMed's own identifiers, ordered authors and supplied name components."""
    result = defaultdict(list)
    for doc in docs:
        if doc['kind'] != 'xml': continue
        root = doc['xml']
        for article in root.iter('PubmedArticle'):
            pmid = article.findtext('MedlineCitation/PMID')
            if not pmid: continue
            authors, structured = [], []
            for author in article.findall('MedlineCitation/Article/AuthorList/Author'):
                fields = {child.tag: ''.join(child.itertext()) for child in author
                          if child.tag in {'ForeName', 'Initials', 'LastName', 'Suffix', 'CollectiveName'}}
                structured.append(fields)
                authors.append(fields.get('CollectiveName') or ' '.join(filter(None, (
                    fields.get('ForeName') or fields.get('Initials'), fields.get('LastName'), fields.get('Suffix')))))
            title = article.find('MedlineCitation/Article/ArticleTitle')
            value = {'pubmed_id': pmid, 'author_list': ', '.join(authors), 'authors': structured,
                     'title': ''.join(title.itertext()) if title is not None else '',
                     'format': 'efetch', 'evidence': provenance(doc, f'PubmedArticle[PMID={pmid}]')}
            for node in article.findall('PubmedData/ArticleIdList/ArticleId'):
                if node.get('IdType') == 'doi': value['doi'] = node.text
            result[pmid].append(value)
        for article in root.iter('DocSum'):
            pmid = article.findtext('Id')
            authors = article.findall("Item[@Name='AuthorList']/Item[@Name='Author']")
            if not pmid or not authors: continue
            value = {'pubmed_id': pmid, 'author_list': ', '.join(n.text or '' for n in authors),
                     'title': article.findtext("Item[@Name='Title']") or '', 'format': 'esummary',
                     'evidence': provenance(doc, f'DocSum[Id={pmid}]')}
            doi = article.find(".//Item[@Name='doi']")
            if doi is not None: value['doi'] = doi.text
            result[pmid].append(value)
    return dict(result)


def comparison(scope, field, expected, actual, evidence):
    return {'scope': scope, 'field': field, 'source': expected, 'generated': actual, 'reference': None,
            'classification': EQUIVALENT if expected == actual else LOSS, 'evidence': evidence,
            'reason': 'Independent explicit source mapping; values and entity identity compared.'}


def geo_checks(docs, package):
    checks = []
    roots = [d for d in docs if d.get('kind') == 'xml' and d['xml'].tag == 'MINiML']
    if not roots: return checks
    doc = roots[0]
    root = doc['xml']
    for node in root.findall('Series'):
        if node.get('iid') != package['series']['iid']: continue
        for tag, field in [('Title', 'title'), ('Summary', 'summary'), ('Overall-Design', 'overall_design')]:
            raw = node.findtext(tag)
            if raw and raw.strip():
                checks.append(comparison('GEO series ' + node.get('iid'), field, raw,
                    package['series'].get(field), provenance(doc, 'Series/' + tag)))
        pmids = [n.text for n in node.findall('Pubmed-ID') if n.text]
        checks.append(comparison('GEO series ' + node.get('iid'), 'pubmed_id', pmids,
                                 package['series'].get('pubmed_id', []), provenance(doc, 'Series/Pubmed-ID')))
    samples = {s['iid']: s for s in package.get('sample', [])}
    for node in root.findall('Sample'):
        iid = node.get('iid')
        actual = samples.get(iid, {})
        for tag, field in [('Title', 'title'), ('Description', 'description'), ('Type', 'type'), ('Channel-Count', 'channel_count'),
                           ('Hybridization-Protocol', 'hybridization_protocol'), ('Scan-Protocol', 'scan_protocol'),
                           ('Data-Processing', 'data_processing')]:
            raw = node.findtext(tag)
            if raw and raw.strip(): checks.append(comparison(iid, field, raw, actual.get(field), provenance(doc, iid + '/' + tag)))
        for index, channel in enumerate(node.findall('Channel')):
            target = (actual.get('channel', []) + [{}] * (index + 1))[index]
            for tag, field in [('Source', 'source'), ('Molecule', 'molecule'), ('Label', 'label')]:
                raw = channel.findtext(tag)
                if raw and raw.strip(): checks.append(comparison(f'{iid}/channel/{index + 1}', field, raw,
                    (target.get(field) or {}).get('value'), provenance(doc, f'{iid}/Channel[{index + 1}]/{tag}')))
            for tag, field in [('Growth-Protocol', 'growth_protocol'), ('Treatment-Protocol', 'treatment_protocol'),
                               ('Extract-Protocol', 'extract_protocol'), ('Label-Protocol', 'label_protocol')]:
                raw = channel.findtext(tag)
                if raw and raw.strip(): checks.append(comparison(f'{iid}/channel/{index + 1}', field, raw,
                    target.get(field), provenance(doc, f'{iid}/Channel[{index + 1}]/{tag}')))
            for tag, field in [('Organism', 'organism'), ('Characteristics', 'characteristics')]:
                expected = []
                for child in channel.findall(tag):
                    if not child.text or not child.text.strip(): continue
                    value = {'value': child.text}
                    if tag == 'Characteristics': value['name'] = child.get('tag') or 'unspecified'
                    if child.get('taxid'): value['taxid'] = child.get('taxid')
                    expected.append(value)
                if expected:
                    observed = [{k: v for k, v in value.items() if k in {'value', 'name', 'taxid'}}
                                for value in target.get(field, [])]
                    checks.append(comparison(f'{iid}/channel/{index + 1}', field, expected, observed,
                                             provenance(doc, f'{iid}/Channel[{index + 1}]/{tag}')))
    contributors = {c['iid']: c for c in package.get('contributor', []) if 'iid' in c}
    for node in root.findall('Contributor'):
        iid = node.get('iid')
        actual = contributors.get(iid, {})
        for tag, field in [('First', 'first'), ('Middle', 'middle'), ('Last', 'last')]:
            raw = node.findtext('Person/' + tag)
            if raw and raw.strip(): checks.append(comparison('contact/' + iid, 'person/' + field, raw,
                actual.get('person', {}).get(field), provenance(doc, f'Contributor[{iid}]/Person/{tag}')))
        for tag, field in [('Email', 'email'), ('Phone', 'phone'), ('Fax', 'fax'), ('Laboratory', 'laboratory'),
                           ('Department', 'department'), ('Organization', 'organization')]:
            raw = node.findtext(tag)
            if raw and raw.strip(): checks.append(comparison('contact/' + iid, field, raw, actual.get(field),
                                                            provenance(doc, f'Contributor[{iid}]/{tag}')))
    return checks


def geo_sdrf_checks(docs, document):
    """Source-derived sample/channel characteristics, with row identity checks."""
    checks = []
    roots = [d for d in docs if d.get('kind') == 'xml' and d['xml'].tag == 'MINiML']
    if not roots: return checks
    doc = roots[0]
    sdrf = next(row[1] for row in document if row[0] == 'SDRF File')
    headers = sdrf[0]
    for sample in doc['xml'].findall('Sample'):
        iid = sample.get('iid')
        channels = sample.findall('Channel')
        for channel_number, channel in enumerate(channels, 1):
            rows = [(index, row) for index, row in enumerate(sdrf[1:], 1)
                    if iid in row_identity(headers, row)['sample']]
            if len(channels) > 1:
                source = channel.findtext('Source')
                index = headers.index('Comment[Sample_source_name]') if 'Comment[Sample_source_name]' in headers else None
                rows = [(i, row) for i, row in rows if index is not None and row[index] == source]
            if not rows:
                checks.append({'scope': iid, 'field': f'channel/{channel_number}/SDRF row', 'source': channel.findtext('Source'),
                    'generated': None, 'reference': None, 'classification': UNRESOLVED,
                    'reason': 'No unique source/channel row evidence; not guessed from physical row order.',
                    'evidence': provenance(doc, f'{iid}/Channel[{channel_number}]')})
                continue
            expected = defaultdict(list)
            for annotation in channel.findall('Characteristics'):
                if annotation.text and annotation.text.strip():
                    expected['Characteristics[' + (annotation.get('tag') or 'unspecified') + ']'].append(annotation.text)
            for annotation in channel.findall('Organism'):
                if annotation.text: expected['Characteristics[organism]'].append(annotation.text)
            for index, row in rows:
                for header, values in expected.items():
                    actual = [row[j] for j, name in enumerate(headers) if name == header and row[j] != '']
                    checks.append(comparison(f'{iid}/channel/{channel_number}/SDRF-row/{index}', header, values, actual,
                        provenance(doc, f'{iid}/Channel[{channel_number}]/{header}')))
    return checks


def geo_contact_idf_checks(docs, document):
    roots = [d for d in docs if d.get('kind') == 'xml' and d['xml'].tag == 'MINiML']
    if not roots: return []
    doc = roots[0]
    root = doc['xml']
    eligible = {node.get('ref') for entity in [*root.findall('Series'), *root.findall('Sample')]
                for node in entity.iter() if node.tag in {'Contact-Ref', 'Contributor-Ref'}}
    rows = {row[0]: row[1:] for row in document if row and row[0].startswith('Person ')}
    count = max(map(len, rows.values()), default=0)
    contacts = [{key: values[i] if i < len(values) else None for key, values in rows.items()} for i in range(count)]
    checks = []
    for person in root.findall('Contributor'):
        if person.get('iid') not in eligible or person.find('Person') is None: continue
        expected = {key: person.findtext(path) for key, path in [
            ('Person First Name', 'Person/First'), ('Person Mid Initials', 'Person/Middle'),
            ('Person Last Name', 'Person/Last'), ('Person Email', 'Email'), ('Person Phone', 'Phone'), ('Person Fax', 'Fax')]}
        parts = [person.findtext(tag) for tag in ('Laboratory', 'Department', 'Organization') if person.findtext(tag)]
        expected['Person Affiliation'] = ', '.join(dict.fromkeys(parts)) or None
        address = person.find('Address')
        if address is not None:
            expected['Person Address'] = ', '.join(value for value in [expected['Person Affiliation'], *[child.text for child in address]] if value)
        candidates = [contact for contact in contacts if all((contact.get(key) or '') == (expected.get(key) or '')
                      for key in ('Person First Name', 'Person Mid Initials', 'Person Last Name', 'Person Email'))]
        for key, value in expected.items():
            if value:
                actual = candidates[0].get(key) if len(candidates) == 1 else None
                checks.append(comparison('GEO contact ' + person.get('iid') + '/IDF', key, value, actual,
                    provenance(doc, 'Contributor[' + person.get('iid') + ']')))
    return checks


def native_contact_checks(docs, package, document):
    """Bind BioSample owner contacts to their own sample, then IDF presentation."""
    checks, organizations_only = [], []
    contributors = {item.get('iid'): item for item in package.get('contributor', [])}
    organizations = {item.get('iid'): item for item in package.get('organization', [])}
    rows = {row[0]: row[1:] for row in document if row and row[0].startswith('Person ')}
    idf_contacts = [{key: values[i] if i < len(values) else None for key, values in rows.items()}
                    for i in range(max(map(len, rows.values()), default=0))]
    seen = set()
    for doc in docs:
        if doc['kind'] != 'xml' or doc['reference_assisted']: continue
        for sample in doc['xml'].iter('BioSample'):
            accession = sample.get('accession')
            if accession in seen: continue
            seen.add(accession)
            owner = sample.find('Owner')
            if owner is None: continue
            people = owner.findall('Contacts/Contact')
            if not people:
                organizations_only.append({'sample': accession, 'organization': owner.findtext('Name'),
                    'evidence': provenance(doc, f'BioSample[{accession}]/Owner'),
                    'assessment': 'Owner organization is not a named publication author or personal contact.'})
                continue
            candidates = [item for item in package.get('sample', []) if accession in {a.get('value') for a in item.get('accession', [])}]
            if len(candidates) != 1: continue
            refs = {r.get('ref') for r in candidates[0].get('contact_ref', [])}
            eligible = [contributors[ref] for ref in refs if ref in contributors]
            for person in people:
                expected = {field: person.findtext('Name/' + tag) for field, tag in [('first', 'First'), ('middle', 'Middle'), ('last', 'Last')]}
                found = [item for item in eligible if all((item.get('person', {}).get(key) or '') == (value or '') for key, value in expected.items())]
                actual = found[0] if len(found) == 1 else {}
                for key, value in expected.items():
                    if value:
                        checks.append(comparison(f'BioSample/{accession}/contact/JSON', key, value,
                            actual.get('person', {}).get(key), provenance(doc, f'BioSample[{accession}]/Owner/Contacts/Contact/Name')))
                if person.get('email'):
                    checks.append(comparison(f'BioSample/{accession}/contact/JSON', 'email', person.get('email'), actual.get('email'),
                        provenance(doc, f'BioSample[{accession}]/Owner/Contacts/Contact/@email')))
                displayed = [item for item in idf_contacts if all((item.get('Person ' + label) or '') == (expected.get(field) or '')
                    for field, label in [('first', 'First Name'), ('middle', 'Mid Initials'), ('last', 'Last Name')])]
                if len(displayed) == 1:
                    for field, label in [('first', 'First Name'), ('middle', 'Mid Initials'), ('last', 'Last Name')]:
                        if expected[field]: checks.append(comparison(f'BioSample/{accession}/contact/IDF', 'Person ' + label,
                            expected[field], displayed[0].get('Person ' + label), provenance(doc, f'BioSample[{accession}]/Owner/Contacts')))
                    if owner.findtext('Name'):
                        checks.append(comparison(f'BioSample/{accession}/contact/IDF', 'Person Affiliation', owner.findtext('Name'),
                            displayed[0].get('Person Affiliation'), provenance(doc, f'BioSample[{accession}]/Owner/Name')))
    return {'checks': checks, 'organization_only_evidence': organizations_only}


def publication_checks(docs, package, document, provider):
    facts = citation_facts(docs)
    checks = []
    series = package.get('series', {})
    rows = {row[0]: row[1:] for row in document if row}
    idf_pmids = rows.get('PubMed ID', [])
    for index, publication in enumerate(series.get('pubmed_publication') or []):
        pmid = publication.get('pubmed_id')
        options = facts.get(pmid, [])
        if not options: continue  # absent hydration is visible in evidence coverage
        selected = next((x for x in options if x['format'] == ('esummary' if provider == 'geo' else 'efetch')), options[0])
        for field, label in [('author_list', 'Publication Author List'), ('doi', 'Publication DOI'), ('title', 'Publication Title')]:
            if not selected.get(field): continue
            checks.append(comparison(f'publication/{pmid}/JSON', field, selected[field], publication.get(field), selected['evidence']))
            # Find by the paper identifier, never by the position of another list.
            positions = [i for i, value in enumerate(idf_pmids) if value == pmid]
            if len(positions) == 1:
                values = rows.get(label, [])
                actual = values[positions[0]] if positions[0] < len(values) else None
                checks.append(comparison(f'publication/{pmid}/IDF', label, selected[field], actual, selected['evidence']))
        for check in checks:
            if not check['scope'].startswith(f'publication/{pmid}/') or check['field'] not in {'author_list', 'Publication Author List'}:
                continue
            if author_representation_equivalent(selected, check['generated']):
                check.update(classification=EQUIVALENT,
                    reason='Ordered author identities agree with PubMed; full names/explicit initials or final list punctuation differ only in representation.')
    return checks


def author_representation_equivalent(source, actual):
    if not isinstance(actual, str): return False
    expected = source['author_list']
    if actual.rstrip('.') == expected.rstrip('.'): return True
    authors = source.get('authors')
    if not authors: return False
    # Only PubMed's explicit Initials establish abbreviated given names. Never
    # infer initials from a surname or waive a supplied suffix/consortium name.
    abbreviated = []
    for author in authors:
        if author.get('CollectiveName'):
            abbreviated.append(author['CollectiveName'])
        elif author.get('LastName') and author.get('Initials'):
            abbreviated.append(' '.join(filter(None, [author['LastName'], author['Initials'], author.get('Suffix')])))
        else:
            return False
    return actual.rstrip('.') == ', '.join(abbreviated).rstrip('.')


def native_coverage(docs, package):
    xml_runs, listed_runs, expected_samples = set(), set(), set()
    evidence = []
    for doc in docs:
        if doc['reference_assisted']: continue
        if doc['kind'] == 'xml':
            for node in doc['xml'].iter('RUN'):
                value = node.get('accession')
                if value and re.fullmatch(r'[SED]RR\d+', value): xml_runs.add(value)
            for node in doc['xml'].iter('SAMPLE'):
                value = node.get('accession')
                if value and re.fullmatch(r'[SED]RS\d+', value): expected_samples.add(value)
        elif (doc['kind'] == 'json' and isinstance(doc['data'], list)
              and parse_qs(urlsplit(doc['url']).query).get('result', [''])[0] in {'read_run', 'read_experiment'}):
            for row in doc['data']:
                if not isinstance(row, dict): continue
                if row.get('run_accession'):
                    listed_runs.update(value for value in row['run_accession'].split(';') if re.fullmatch(r'[SED]RR\d+', value))
                if row.get('secondary_sample_accession'): expected_samples.add(row['secondary_sample_accession'])
        evidence.append(provenance(doc, 'RUN/@accession; SAMPLE/@accession; read_run accession fields'))
    actual_runs = {r.get('run') for sample in package.get('sample', []) for r in sample.get('sra_run', []) if r.get('run')}
    actual_samples = {a['value'] for sample in package.get('sample', []) for a in sample.get('accession', [])
                      if re.fullmatch(r'[SED]RS\d+', a.get('value', ''))}
    expected_runs = listed_runs or xml_runs
    return {'source_runs': sorted(expected_runs), 'generated_runs': sorted(actual_runs),
            'source_samples': sorted(expected_samples), 'generated_samples': sorted(actual_samples),
            'runs_equal': expected_runs == actual_runs, 'samples_equal': expected_samples == actual_samples,
            'run_enumerator': 'study-scoped ENA read_run/read_experiment listings' if listed_runs else 'native RUN XML',
            'additional_retrieved_xml_runs_not_in_study_listing': sorted(xml_runs - listed_runs) if listed_runs else [],
            'evidence': evidence}
