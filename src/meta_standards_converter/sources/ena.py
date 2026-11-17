# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""ENA Portal membership discovery and full Browser record retrieval."""
from meta_standards_converter.provider_issues import provider_boundary
import re
import xml.etree.ElementTree as ET
from .archive_support import ArchiveHTTP, Resolution, StudyRecords, StudySeed, accession_kind, attempt, chunks, identifier


class ENASource:
    xref_page_size = 1000
    portal = 'https://www.ebi.ac.uk/ena/portal/api/'
    browser = 'https://www.ebi.ac.uk/ena/browser/api/xml/'

    def __init__(self, http=None, requester=None, resource_profile='standard', evidence_dir=None):
        self.http = http or ArchiveHTTP('ena_portal', requester, resource_profile, evidence_dir)
        self._xrefs = {}

    def cross_references(self, accession, records):
        """Fetch bounded JSON pages; TSV has no header and loses its first row."""
        def family_scope(row):
            return (isinstance(row, dict) and row.get('Target')=='assembly'
                    and re.fullmatch(r'GC[AF]_\d+\.\d+', accession)
                    and row.get('Target Primary Accession')==accession.rsplit('.',1)[0]
                    and row.get('Target Secondary Accession') in ('',None,accession.rsplit('.',1)[0]))
        if accession not in self._xrefs:
            rows, issues, failures, offset = [], [], [], 0
            seen_pages = set()
            while True:
                try:
                    page = self.http.get('https://www.ebi.ac.uk/ena/xref/rest/json/search',
                        {'accession':accession,'limit':self.xref_page_size,'offset':offset}, 'json')
                    if not isinstance(page, list): raise ValueError('expected cross-reference array')
                    signature = repr(page)
                    if page and signature in seen_pages: raise ValueError('repeated cross-reference page')
                    seen_pages.add(signature)
                    for row in page:
                        if not isinstance(row, dict) or (accession not in {row.get('Target Primary Accession'),row.get('Target Secondary Accession')} and not family_scope(row)):
                            issues.append(f'ENA cross references {accession}: mismatched target')
                            continue
                        rows.append(row)
                    if len(page) < self.xref_page_size: break
                    offset += len(page)
                except Exception as error:
                    from ..metadata.retrieval_budget import _BudgetExhausted
                    if isinstance(error, _BudgetExhausted): raise
                    issues.append(f'ENA cross references {accession}: {type(error).__name__}')
                    failures.append(error)
                    break
            self._xrefs[accession] = rows, issues, failures
        rows, issues, failures = self._xrefs[accession]
        records.issues.extend(i for i in issues if i not in records.issues)
        records.failures.extend(failures)
        grouped = {}
        for row in rows:
            owner = row['Target Primary Accession'] if family_scope(row) else accession
            grouped.setdefault(owner, []).append(row)
        for owner, values in grouped.items():
            record = {'provider':'ena','kind':'cross_references','accession':owner,'metadata':values}
            if record not in records.linked: records.linked.append(record)

    def search(self, result, query, fields='all'):
        rows = self.http.get(self.portal + 'search', {'result': result, 'query': query,
            'fields': fields, 'format': 'json', 'limit': 0, 'includeMetagenomes': 'true'}, 'json')
        if not isinstance(rows, list):
            raise ValueError('ENA search must return an array')
        return rows

    def xml(self, accessions):
        return self.http.get(self.browser + ','.join(accessions), {'includeLinks': 'true'})

    def verified_xml(self, accessions, result, kind='record'):
        from collections import deque
        import requests
        from ..helpers.metadata_stream import (AttemptAllowance, attempt_allowance,
            split_payload, BatchPayloadFailure, BatchRequestLimit, _UnchangedRequestLimit,
            _BatchHeaderTimeout)
        from ..helpers.request_helper import RateLimitedRequester
        from ..metadata.retrieval_budget import _BudgetExhausted, _check_cancelled
        from ..conversion_logging import emit
        allowance = AttemptAllowance()
        token = attempt_allowance.set(allowance)
        queue = deque([(list(dict.fromkeys(accessions)), 0, None)])
        clean, found, failures = None, set(), []
        headers_available = False
        accounted = isinstance(self.http, ArchiveHTTP) and isinstance(self.http.requester, RateLimitedRequester)
        try:
            while queue and allowance.used < allowance.maximum:
                _check_cancelled()
                batch, repeated, probe_followup = queue.popleft()
                batch = [acc for acc in batch if acc not in found]
                if not batch:
                    continue
                split = split_payload.set(len(batch) > 25)
                payload_failure = False
                header_stall = False
                try:
                    if not accounted:
                        allowance.charge(tuple(batch))
                    root = self.xml(batch)
                    if clean is None:
                        clean = ET.Element(root.tag)
                    wanted = set(batch)
                    for node in root:
                        aliases = {identifier(node)} | {n.text for n in node.findall('IDENTIFIERS/PRIMARY_ID')} | {n.text for n in node.findall('IDENTIFIERS/SECONDARY_ID')}
                        # Explicit archive accessions are aliases; arbitrary
                        # submitter names and unrelated external labels are not.
                        aliases.update(n.text for n in node.findall('IDENTIFIERS/EXTERNAL_ID')
                                       if re.fullmatch(r'(?:SAM(?:N|EA|D)\d+|PRJ(?:NA|EB|DB|DA)\d+|[SED]R[PSRXZ]\d+|GC[AF]_\d+(?:\.\d+)?)',n.text or ''))
                        matches = wanted & aliases
                        matches.update(acc for acc in batch if re.fullmatch(r'GC[AF]_\d+', acc)
                            and any(re.fullmatch(re.escape(acc) + r'\.\d+', a or '') for a in aliases))
                        allowed = {'sample':{'SAMPLE'},'experiment':{'EXPERIMENT'},'run':{'RUN'},
                                   'study':{'STUDY','PROJECT'},'project':{'PROJECT','STUDY'},
                                   'analysis':{'ANALYSIS'},'assembly':{'ASSEMBLY'}}.get(kind.removeprefix('read_'))
                        if allowed is None:
                            allowed=set()
                            for acc in matches:
                                for pattern,tags in [(r'SAM(?:N|EA|D)\d+|[SED]RS\d+',{'SAMPLE'}),
                                    (r'[SED]RX\d+',{'EXPERIMENT'}),(r'[SED]RR\d+',{'RUN'}),
                                    (r'[SED]RP\d+',{'STUDY','PROJECT'}),(r'PRJ(?:NA|EB|DB|DA)\d+',{'PROJECT','STUDY'}),
                                    (r'[SED]RZ\d+',{'ANALYSIS'}),(r'GC[AF]_\d+(?:\.\d+)?',{'ASSEMBLY'})]:
                                    if re.fullmatch(pattern,acc):allowed.update(tags)
                        if matches and node.tag in allowed and not any(n.tag.lower() == 'error' for n in node.iter()):
                            # Preserve duplicate occurrences within a complete
                            # response; subsequent retries request only missing IDs.
                            clean.append(node)
                            found.update(matches)
                except _BudgetExhausted:
                    raise
                except _UnchangedRequestLimit:
                    continue
                except BatchRequestLimit:
                    break
                except Exception as error:
                    failures.append(error)
                    # Native streaming already spent its unchanged-request retry
                    # budget. Large body/XML failures permit repair. A header
                    # stall requires verified single-record availability first.
                    payload_failure = isinstance(error, (BatchPayloadFailure, ET.ParseError)) or (
                        not accounted and isinstance(error, requests.exceptions.ChunkedEncodingError))
                    header_stall = isinstance(error, _BatchHeaderTimeout)
                    emit('xml_batch_recovery', level='WARNING', kind=kind,
                         requested=len(batch), physical_attempts=allowance.used,
                         error_type=type(error).__name__, split=payload_failure and len(batch)>25)
                    if not payload_failure:
                        continue
                finally:
                    split_payload.reset(split)
                missing = [acc for acc in batch if acc not in found]
                if probe_followup is not None:
                    if missing:
                        continue  # No verified small response: do not fan out.
                    headers_available = True
                    missing = [acc for acc in probe_followup if acc not in found]
                elif header_stall and not headers_available and missing:
                    queue.appendleft((missing[:1], 0, missing[1:]))
                    continue
                if not missing:
                    continue
                if len(missing) > 25:
                    middle = len(missing)//2
                    queue.appendleft((missing[middle:], 0, None))
                    queue.appendleft((missing[:middle], 0, None))
                elif repeated < (3 if payload_failure else 1):
                    queue.appendleft((missing, repeated+1, None))
        finally:
            attempt_allowance.reset(token)
        missing_ids = sorted(set(accessions) - found)
        if missing_ids and failures:
            result.issues.append(f'{kind} XML: {type(failures[-1]).__name__}')
            if hasattr(result, 'failures'):
                result.failures.extend(failures)
        for missing in missing_ids:
            result.issues.append(f'{kind}: missing XML record {missing}')
        return clean if clean is not None and len(clean) else None

    def linked_json(self, url, accession, id_field, result):
        obj = attempt(result, f'linked {accession}', lambda: self.http.get(url, fmt='json'))
        if obj is not None and (not isinstance(obj, dict) or str(obj.get(id_field)) != accession):
            result.issues.append(f'linked {accession}: missing or mismatched identity')
            return None
        return obj

    def fetch_biosamples(self, accessions, records):
        """Retrieve owner/contact evidence for the explicitly bound BioSamples."""
        from .entrez_records import EntrezRecords
        client = EntrezRecords(self.http)
        for batch in chunks(sorted(set(accessions))):
            root = client.linked_xml('biosample', batch, records)
            if root is not None:
                records.xml.append(root)

    @provider_boundary(provider='pubmed', operation='publication_retrieval')
    def publications(self, ids, result):
        root = attempt(result, 'PubMed', lambda: self.http.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi',
                {'db': 'pubmed', 'id': ','.join(ids), 'retmode': 'xml'}))
        if root is None:
            return None
        clean, found = ET.Element('PubmedArticleSet'), set()
        for article in root.findall('PubmedArticle'):
            pmid = article.findtext('MedlineCitation/PMID')
            if pmid in ids and not any(n.tag.lower() == 'error' for n in article.iter()):
                clean.append(article)
                found.add(pmid)
        for missing in sorted(set(ids) - found):
            result.issues.append(f'PubMed {missing}: missing or mismatched identity')
        return clean if len(clean) else None

    @provider_boundary(provider='ena', operation='study_resolution')
    def resolve(self, accession):
        accession, kind = accession_kind(accession)
        result, visited, queries = Resolution(), set(), []
        def project(acc):
            if acc in visited:
                return
            visited.add(acc)
            root = self.verified_xml([acc], result, 'project')
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

    def fetch_reference_ranges(self, records):
        """Validate bounded ENA accession ranges; keep their records outside membership."""
        kinds = {'ENA-RUN': ('RUN', r'[SED]RR\d+'),
                 'ENA-EXPERIMENT': ('EXPERIMENT', r'[SED]RX\d+'),
                 'ENA-STUDY': ('STUDY', r'[SED]RP\d+')}
        verified = {tag: set() for tag, _ in kinds.values()}
        patterns = dict(kinds.values())

        def remember(root):
            # Only full records establish identity, never Portal rows or *_REF
            # elements. Provider primary/secondary identifiers are exact aliases;
            # submitter labels, external namespaces and version guesses are not.
            if root.tag.rsplit('}', 1)[-1].lower() == 'error':
                return
            nodes = [root] if root.tag in verified else list(root)
            for node in nodes:
                if node.tag not in verified or any(n.tag.rsplit('}', 1)[-1].lower() == 'error' for n in node.iter()):
                    continue
                aliases = {identifier(node)} | {n.text for n in node.findall('IDENTIFIERS/PRIMARY_ID')}
                aliases.update(n.text for n in node.findall('IDENTIFIERS/SECONDARY_ID'))
                verified[node.tag].update(a for a in aliases if a and re.fullmatch(patterns[node.tag], a))

        for root in records.xml:
            remember(root)
        for root in list(records.xml):
            for sample in root.findall('.//SAMPLE'):
                for link in sample.findall('.//XREF_LINK'):
                    db, literal = link.findtext('DB',''), link.findtext('ID','')
                    if db not in ('ENA-RUN','ENA-EXPERIMENT','ENA-STUDY') or '-' not in literal:
                        continue
                    candidates=[]
                    for part in literal.split(','):
                        match=re.fullmatch(r'([SED]R[PRX])(\d+)-\1(\d+)', part.strip())
                        if match and 0 <= int(match[3])-int(match[2]) < 1000:
                            candidates.extend(match[1]+str(i).zfill(len(match[2])) for i in range(int(match[2]),int(match[3])+1))
                        elif re.fullmatch(r'[SED]R[PRX]\d+', part.strip()): candidates.append(part.strip())
                        else:
                            candidates=[]
                            break
                    if not candidates or len(candidates)>1000: continue
                    tag, _ = kinds[db]
                    missing = list(dict.fromkeys(a for a in candidates if a not in verified[tag]))
                    for batch in chunks(missing):
                        response=self.verified_xml(batch,records,'accession range')
                        if response is not None:
                            remember(response)
                        for accession in batch:
                            issue = f'accession range: missing XML record {accession}'
                            if accession not in verified[tag] and issue not in records.issues:
                                records.issues.append(issue)
                    if set(candidates)<=verified[tag]:
                        records.linked.append({'provider':'ena','kind':'accession_range','accession':identifier(sample),
                            'metadata':{'database':db,'literal':literal,'accessions':candidates}})

    def fetch_associated_analyses(self, records):
        """Fetch one hop of explicit analysis links without importing their read studies."""
        members, requests, existing = set(), set(), set()
        for root in records.xml:
            for node in root.iter():
                if node.tag in ('SAMPLE', 'RUN'):
                    members.add(identifier(node))
                    members.update(n.text for n in node.findall('IDENTIFIERS/*'))
                    for link in node.findall('.//XREF_LINK'):
                        if link.findtext('DB') == 'ENA-ANALYSIS':
                            requests.update(a.strip() for a in (link.findtext('ID') or '').split(',') if re.fullmatch(r'[SED]RZ\d+', a.strip()))
                if node.tag == 'ANALYSIS': existing.add(identifier(node))
        for batch in chunks(sorted(requests - existing)):
            root = self.verified_xml(batch, records, 'associated analysis')
            if root is None: continue
            accepted = ET.Element(root.tag)
            for node in root:
                refs = {identifier(n) for tag in ('SAMPLE_REF','RUN_REF') for n in node.findall(tag)}
                if refs & members:
                    accepted.append(node)
                else:
                    records.issues.append(f'{identifier(node)}: associated analysis has no verified dataset binding')
            if len(accepted):
                records.xml.append(accepted)
                for node in accepted:
                    acc = identifier(node)
                    rows = attempt(records, f'{acc} indexed analysis', lambda: self.search('analysis', f'analysis_accession="{acc}"'))
                    for row in rows or []:
                        refs = {v for key in ('sample_accession','run_accession') for v in str(row.get(key) or '').split(';') if v}
                        if row.get('analysis_accession') == acc and (not refs or refs & members):
                            records.indexed.setdefault('analysis', []).append(row)
                        else:
                            records.issues.append(f'{acc}: mismatched indexed analysis identity or membership')

    @provider_boundary(provider='ena', operation='metadata_retrieval')
    def fetch(self, seed, *, _on_core=None, publications=True):
        records = StudyRecords(seed)
        query = f'secondary_study_accession="{seed.study}"' if seed.study.startswith(('SRP', 'ERP', 'DRP')) else f'study_accession="{seed.primary}"'
        for accession in dict.fromkeys((seed.study, seed.primary)):
            root = self.verified_xml([accession], records)
            if root is not None:
                records.xml.append(root)
        kinds = ('study', 'sample', 'read_experiment', 'read_run')
        if _on_core is None: kinds += ('analysis', 'assembly')
        for kind in kinds:
            self._inventory(kind, records, query)
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
            if kind not in kinds: continue
            accessions = list(dict.fromkeys((r.get('assembly_set_accession') or r[field]) if kind == 'assembly' else r[field]
                                           for r in records.indexed.get(kind, []) if r.get(field)))
            if kind == 'sample':
                members = {r.get('sample_accession') for k in ('read_experiment', 'read_run') for r in records.indexed.get(k, []) if r.get('sample_accession')}
                if members:
                    accessions = sorted(members)
            inventory[kind] = accessions
            for batch in chunks(accessions):
                root = self.verified_xml(batch, records, kind)
                if root is not None:
                    records.xml.append(root)
        if _on_core is not None: _on_core(records)
        from ..metadata.retrieval_budget import _budget_stage
        with _budget_stage('native_supplementary'):
            return self._supplement(records, inventory, query, extras=_on_core is not None, publications=publications)

    @provider_boundary(provider='ena', operation=lambda self, kind, *args, **kwargs: {'sample': 'sample_inventory', 'read_run': 'run_inventory', 'read_experiment': 'experiment_inventory'}.get(kind, 'metadata_retrieval'))
    def _inventory(self, kind, records, query):
        scoped = f'study_accession="{records.seed.primary}"' if kind in ('study', 'sample', 'assembly') else query
        rows = attempt(records, f'{kind} inventory', lambda: self.search(kind, scoped))
        if rows is not None: records.indexed[kind] = rows

    def _supplement(self, records, inventory, query, *, extras=False, publications=True):
        seed = records.seed
        if extras:
            for kind, field in [('analysis', 'analysis_accession'), ('assembly', 'assembly_accession')]:
                self._inventory(kind, records, query)
                accessions = list(dict.fromkeys((r.get('assembly_set_accession') or r[field]) if kind == 'assembly' else r[field]
                                               for r in records.indexed.get(kind, []) if r.get(field)))
                inventory[kind] = accessions
                for batch in chunks(accessions):
                    root = self.verified_xml(batch, records, kind)
                    if root is not None: records.xml.append(root)
        self.fetch_reference_ranges(records)
        self.fetch_associated_analyses(records)
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
        references = {seed.study, seed.primary}
        for root in records.xml:
            references.update(identifier(n) for n in root.iter() if n.tag in
                              ('STUDY','PROJECT','SAMPLE','EXPERIMENT','RUN','ANALYSIS','ASSEMBLY') and identifier(n))
        self._fetch_independent_links(records, sorted(references), inventory['sample'])
        from .archive_support import publication_ids
        taxa = set()
        for root in records.xml:
            taxa.update(n.text for n in root.findall('.//SAMPLE_NAME/TAXON_ID') if n.text)
        self.fetch_biosamples([acc for acc in inventory['sample'] if acc.startswith('SAM')], records)
        for taxid in sorted(taxa):
            obj = self.linked_json('https://www.ebi.ac.uk/ena/taxonomy/rest/tax-id/' + taxid, taxid, 'taxId', records)
            if obj is not None:
                records.linked.append({'provider': 'ena', 'kind': 'taxonomy', 'accession': taxid, 'metadata': obj})
        if publications:
            from .archive_publications import resolve_identifiers
            resolve_identifiers(records, self.http, 'ena')
            for batch in chunks(sorted(publication_ids(records))):
                root = self.publications(batch, records)
                if root is not None:
                    records.xml.append(root)
        return records

    def _fetch_independent_links(self, records, references, samples):
        """Two bounded lanes; only the caller merges evidence into StudyRecords."""
        from ..metadata.retrieval_budget import _cancelled, _check_cancelled, _RetrievalCancelled
        def xrefs(target):
            for accession in references:
                _check_cancelled()
                self.cross_references(accession, target)
        def biosamples(target):
            for accession in samples:
                _check_cancelled()
                if accession.startswith('SAM'):
                    obj = self.linked_json('https://www.ebi.ac.uk/biosamples/samples/' + accession,
                                           accession, 'accession', target)
                    if obj is not None:
                        target.linked.append({'provider': 'biosamples', 'kind': 'sample',
                                              'accession': accession, 'metadata': obj})
        workers = getattr(getattr(self.http, 'profile', None), 'network_workers', 1)
        if workers < 2 or not references or not any(s.startswith('SAM') for s in samples):
            xrefs(records)
            biosamples(records)
            return
        from concurrent.futures import ThreadPoolExecutor, wait, FIRST_EXCEPTION
        from threading import Event
        from ..conversion_logging import context_callable
        stopped = Event()
        token = _cancelled.set(stopped)
        parts = [StudyRecords(records.seed), StudyRecords(records.seed)]
        executor = ThreadPoolExecutor(max_workers=2, thread_name_prefix='msc-ena-links')
        futures = []
        def run(fn, part):
            try:
                fn(part)
            except BaseException:
                stopped.set()
                raise
        try:
            futures = [executor.submit(context_callable(run), fn, part)
                       for fn, part in zip((xrefs, biosamples), parts)]
            wait(futures, return_when=FIRST_EXCEPTION)
            errors = [f.exception() for f in futures if f.done() and f.exception() is not None]
            if errors:
                raise next((e for e in errors if not isinstance(e, _RetrievalCancelled)), errors[0])
        finally:
            stopped.set()
            executor.shutdown(wait=True, cancel_futures=True)
            _cancelled.reset(token)
            records.linked.extend(row for row in parts[0].linked if row not in records.linked)
            records.issues.extend(issue for issue in parts[0].issues if issue not in records.issues)
            records.linked.extend(parts[1].linked)
            records.issues.extend(parts[1].issues)
            records.failures.extend(parts[0].failures)
            records.failures.extend(parts[1].failures)
