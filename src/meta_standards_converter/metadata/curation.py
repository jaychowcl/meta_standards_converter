# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Network-free, advisory assessment of final metadata for human curation.

This is a checklist, not schema validation. It never fills scientific fields or
interprets missing biology as a negative assertion. Locations are JSON pointers
within the prepared package; context identifies its dataset and package ordinal.
"""
from collections import defaultdict
from collections.abc import Mapping
import re

from ..converters.unified.contracts import CurationContext, Diagnostic
from .fields import BIOLOGICAL_ALIASES
from .preparation import preparation_operation
from ..miniml.publication_identifiers import publication_identity


_UNKNOWN = {'unknown', 'not provided', 'not available', 'not collected',
            'not reported', 'unspecified', 'missing'}
_BIOLOGY = ('organism', 'organism_part', 'developmental_stage', 'disease', 'genotype')
_LIBRARY = ('library_strategy', 'library_source', 'library_selection', 'library_layout')
_PROTOCOL_TYPES = {
    'collection': r'collect|sampling', 'extraction': r'extract',
    'library_preparation': r'library', 'sequencing': r'sequenc',
    'hybridization': r'hybridi', 'labeling': r'label', 'scan': r'scan',
    'processing': r'processing|normaliz|transformation',
}
_RAW_NODES = {'array_data_file', 'array_data_matrix_file', 'image_file'}
_PROCESSED_NODES = {'derived_array_data_file', 'derived_array_data_matrix_file'}


def _items(value):
    return list(value) if isinstance(value, (list, tuple)) else [] if value is None else [value]


def _slug(value):
    return re.sub(r'[\s_-]+', '_', str(value).strip().casefold())


def missing_reason(value):
    """Only explicit unknown markers count; 'none', zero and N/A remain facts."""
    if value is None:
        return 'missing'
    if isinstance(value, str):
        text = ' '.join(value.split()).casefold()
        return 'missing' if not text else 'reported_unknown' if text in _UNKNOWN else None
    if isinstance(value, Mapping):
        return missing_reason(value['value']) if 'value' in value else missing_reason(list(value.values()))
    if isinstance(value, (list, tuple)):
        reasons = [missing_reason(item) for item in value]
        if any(reason is None for reason in reasons):
            return None
        return 'reported_unknown' if 'reported_unknown' in reasons else 'missing'
    return None


def _pointer(key):
    return str(key).replace('~', '~0').replace('/', '~1')


def _value(value):
    return value.get('value') if isinstance(value, Mapping) else value


def _biological_field(name):
    return next((field for field in _BIOLOGY if name in BIOLOGICAL_ALIASES.get(field, {field})), None)


def _lookup_key(kind, identifier):
    return kind, str(identifier).casefold() if kind == 'doi' else str(identifier)


class _Assessment:
    def __init__(self, data, dataset_id, provider, package_index, preparation):
        self.data, self.series = data, data.get('series', {})
        self.dataset_id, self.provider, self.package_index = dataset_id, provider, package_index
        self.result, self.seen = [], set()
        from ..magetab.technology import series_identity
        identity = series_identity(data) or dataset_id
        self.preparation = [r for r in preparation if r.get('dataset') in {identity, dataset_id}
                            and r.get('dataset_id', dataset_id) == dataset_id
                            and r.get('package_index', package_index) == package_index]
        self.lookups = {_lookup_key(r['identifier_kind'], r['identifier']): r for r in self.preparation
                        if r.get('operation') == 'publication_lookup'}
        self.publications_disabled = any(r.get('operation') == 'standard' and r.get('reason') == 'disabled'
                                         for r in self.preparation)
        self.runs_disabled = self.publications_disabled or any(
            r.get('operation') == 'run_metadata' and r.get('reason') == 'disabled' for r in self.preparation)
        self.protocols = {p['name']: p for p in self.series.get('protocols', [])}
        self.paths = defaultdict(list)
        self.step_locations = {}
        for i, path in enumerate(self.series.get('assay_paths', [])):
            self.step_locations.update({id(s): f'/series/assay_paths/{i}/steps/{j}'
                                        for j, s in enumerate(path['steps'])})
            owners = {s['sample_ref'] for s in path['steps'] if s.get('sample_ref')}
            if len(owners) == 1:
                self.paths[next(iter(owners))].append(path)
        self.platforms = {p.get('iid'): p for p in data.get('platform', [])}

    def add(self, field, reason, path, kind, identity, *, action=None, lookup_status=None, lookup=None):
        key = (field, reason, path, kind, identity)
        if key in self.seen:
            return
        self.seen.add(key)
        label = field.split('.', 1)[-1].replace('_', ' ')
        message = f'{label} {reason.replace("_", " ")}'
        action = action or f'Supply verified {label} metadata or confirm it is not applicable.'
        if lookup_status == 'disabled':
            action += ' Relevant enrichment was disabled; enable it to attempt retrieval.'
        lookup = lookup or {}
        context = CurationContext(self.dataset_id, kind, str(identity), field, reason,
            self.package_index, lookup.get('status', lookup_status), lookup.get('provider'),
            lookup.get('identifier'), lookup.get('identifier_kind'))
        self.result.append(Diagnostic('curation_' + reason, message, 'curation', 'warning',
            action=action, provider=lookup.get('provider', self.provider), operation='curation_assessment',
            location=path, curation=context))

    def check(self, value, field, path, kind, identity, **kwargs):
        reason = missing_reason(value)
        if reason:
            self.add(field, reason, path, kind, identity, **kwargs)

    def companions(self, value, path, kind, identity, *, skip=()):
        """Inspect supplied ontology/unit groups, excluding opaque residual data."""
        if isinstance(value, Mapping):
            for prefix in ('', 'status_'):
                source, accession = prefix + 'term_source_ref', prefix + 'term_accession_number'
                if value.get(source) is not None or value.get(accession) is not None:
                    for field in (source, accession):
                        self.check(value.get(field), ('publication.' if prefix else 'ontology.') + field,
                                   path + '/' + field, kind, identity)
            if isinstance(value.get('unit'), Mapping):
                self.check(value['unit'].get('value'), 'unit.value', path + '/unit/value', kind, identity)
            for key, item in value.items():
                if key not in {'extensions', 'records', *skip} and not str(key).startswith('hz_'):
                    self.companions(item, path + '/' + _pointer(key), kind, identity)
        elif isinstance(value, (list, tuple)):
            for i, item in enumerate(value):
                self.companions(item, path + '/' + str(i), kind, identity)

    def publications(self, entity, path, kind, identity, *, required=False):
        from ..sources.archive_publications import citation_identifier
        pubs = [(p, f'{path}/pubmed_publication/{i}', True, identity)
                for i, p in enumerate(entity.get('pubmed_publication', [])) if missing_reason(p) != 'missing']
        known = {tuple(publication_identity(p).items()) for p, *_ in pubs}
        pmids = {publication_identity(p).get('pubmed_id') for p, *_ in pubs}
        # Keep exact locations for unhydrated IDs and unstructured citations.
        for i, value in enumerate(_items(entity.get('pubmed_id'))):
            pub = citation_identifier('pubmed', value) or {'pubmed_id': value}
            if pub.get('pubmed_id') not in pmids and not missing_reason(value):
                pubs.append((pub, f'{path}/pubmed_id/{i}', False, identity))
                pmids.add(pub.get('pubmed_id'))
        if not pubs:
            for i, value in enumerate(_items(entity.get('citation'))):
                if missing_reason(value) != 'missing':
                    pub = value if isinstance(value, Mapping) else citation_identifier('url', value)
                    pubs.append((pub, f'{path}/citation/{i}', isinstance(value, Mapping), identity))
        has_direct_citation = bool(pubs)
        for i, relation in enumerate(entity.get('relation', [])):
            pub = relation.get('publication')
            if not isinstance(pub, Mapping):
                pub = citation_identifier(relation.get('type'), relation.get('target'))
            if not pub:
                continue
            ids = publication_identity(pub)
            scope = next((relation[k] for k in ('experiment_ref', 'analysis_ref', 'assembly_ref') if relation.get(k)), identity)
            if scope == identity and ids and (tuple(ids.items()) in known or
                    ids.get('pubmed_id') and ids['pubmed_id'] in pmids):
                continue
            supplied = isinstance(relation.get('publication'), Mapping)
            pubs.append((pub, f'{path}/relation/{i}/' + ('publication' if supplied else 'target'), supplied, scope))
            if scope == identity:
                has_direct_citation = True
        if required and not has_direct_citation:
            self.check(None, 'publication.citation', path + '/pubmed_publication', kind, identity,
                action='Check the source for a linked paper or confirm that no publication is available.',
                lookup_status='disabled' if self.publications_disabled else None)
        for i, (pub, location, supplied, scope) in enumerate(pubs):
            ids = publication_identity(pub)
            publication_id = str(pub.get('pubmed_id') or pub.get('doi') or pub.get('pmcid') or i)
            pub_identity = f'{scope}/publication[{publication_id}]'
            evidence = [self.lookups[k, v] for k, v in ids.items() if (k, v) in self.lookups]
            lookup = next((r for status in ('failed', 'empty', 'not_found', 'incomplete', 'complete')
                           for r in evidence if r['status'] == status), None)
            context = dict(lookup=lookup, lookup_status='disabled' if self.publications_disabled else None)
            self.check(ids, 'publication.identifier', location, 'publication', pub_identity, **context)
            if lookup and lookup['status'] == 'not_found':
                self.add('publication.identifier', 'lookup_no_match', location, 'publication', pub_identity,
                    action='Verify the supplied DOI/PMCID and curate the citation manually if it is not indexed in PubMed.',
                    lookup=lookup)
            if lookup and lookup['status'] == 'incomplete':
                self.add('lookup.metadata', 'lookup_incomplete', location, 'publication', pub_identity,
                    action='Check the returned citation and supply any fields still absent from the final metadata.',
                    lookup=lookup)
            for field, key in (('title', 'title'), ('authors', 'author_list'), ('doi', 'doi'), ('status', 'status')):
                self.check(pub.get(key), 'publication.' + field, location + ('/' + key if supplied else ''),
                           'publication', pub_identity, **context)
            if supplied:
                self.companions(pub, location, 'publication', pub_identity)

    def study(self):
        series, identity = self.series, self.dataset_id
        for field, keys in (
            ('title', ('title',)), ('description', ('summary', 'description')),
            ('overall_design', ('overall_design',)), ('experimental_type', ('type',))):
            self.check([series.get(k) for k in keys], 'study.' + field, '/series/' + keys[0], 'study', identity)
        for singular, field in (('contributor', 'contributors'), ('contact', 'contacts')):
            people = list(series.get(singular, []))
            definitions = {p.get('iid'): p for key in ('contributor', 'contact') for p in self.data.get(key, [])}
            people += [definitions.get(ref.get('ref'), {}) for ref in series.get(singular + '_ref', [])]
            information = [{k: v for k, v in p.items() if k not in {'iid', 'extensions', 'status'}} for p in people]
            self.check(information, 'study.' + field, '/series/' + singular, 'study', identity)
        self.publications(series, '/series', 'study', identity, required=True)
        self.companions(series, '/series', 'study', identity, skip=('pubmed_publication', 'relation', 'protocols', 'assay_paths'))
        for i, protocol in enumerate(series.get('protocols', [])):
            for key in ('type', 'description'):
                self.check(protocol.get(key), 'protocol.' + key, f'/series/protocols/{i}/{key}',
                           'protocol', protocol['name'])
            self.companions(protocol, f'/series/protocols/{i}', 'protocol', protocol['name'])
        for path in series.get('assay_paths', []):
            owners = {s['sample_ref'] for s in path['steps'] if s.get('sample_ref')}
            for step in path['steps']:
                sp = self.step_locations[id(step)]
                kind, owner = ('sample', next(iter(owners))) if len(owners) == 1 else ('assay', step.get('name', identity))
                self.companions(step, sp, kind, owner)
                self.declared_values(step.get('characteristics', []), sp + '/characteristics', kind, owner)
                self.declared_values(step.get('factor_values', []), sp + '/factor_values', kind, owner, factor=True)
        for i, platform in enumerate(self.data.get('platform', [])):
            platform_id = platform.get('iid') or str(i)
            self.publications(platform, f'/platform/{i}', 'platform', platform_id)
            self.companions(platform, f'/platform/{i}', 'platform', platform_id, skip=('pubmed_publication', 'relation'))
        if self.processed(series, ()):
            descriptions = [p.get('description') for p in self.protocols.values()
                            if re.search(_PROTOCOL_TYPES['processing'], str(_value(p.get('type'))), re.I)]
            self.check(descriptions, 'protocol.processing', '/series/protocols', 'study', identity)

    def declared_values(self, values, path, kind, identity, *, factor=False, factor_names=()):
        for i, value in enumerate(values):
            name = _slug(value.get('name') or value.get('tag'))
            if name.startswith('hz_'):
                continue
            biological = _biological_field(name)
            location = f'{path}/{i}/value'
            if factor or name in factor_names:
                self.check(value.get('value'), 'factor.' + name, location, kind, identity)
            if not factor:
                field = 'sample.' + biological if biological else 'characteristic.' + name
                if biological or name not in factor_names:
                    self.check(value.get('value'), field, location, kind, identity)

    def scoped_steps(self, sample, run=None):
        paths = self.paths.get(sample.get('iid'), ())
        if run is not None:
            paths = [p for p in paths if any(s.get('kind') == 'scan' and s.get('name') == run.get('run')
                                             for s in p['steps']) and
                     (not run.get('experiment') or any(s.get('kind') == 'assay' and s.get('name') == run['experiment']
                                                       for s in p['steps']))]
        return [s for p in paths for s in p['steps']]

    def protocol_text(self, steps, category):
        definitions = [self.protocols[s['protocol_ref']] for s in steps
                       if s.get('protocol_ref') in self.protocols]
        return [p.get('description') for p in definitions
                if re.search(_PROTOCOL_TYPES[category], str(_value(p.get('type')) or ''), re.I)]

    def sample(self, sample, index):
        identity = sample.get('iid') or str(index)
        base = f'/sample/{index}'
        steps = self.scoped_steps(sample)
        self.check(sample.get('title'), 'sample.title', base + '/title', 'sample', identity)
        self.publications(sample, base, 'sample', identity)
        self.companions(sample, base, 'sample', identity, skip=('pubmed_publication', 'relation', 'sra_run'))
        factors = [v.get('name') or v.get('factor') for v in self.series.get('variable', [])
                   if not v.get('sample_ref') or identity in {r.get('ref') for r in v['sample_ref']}]
        factor_names = {_slug(f) for f in factors if f}
        biological = [s for s in steps if s.get('kind') in {'source', 'sample'}]
        for ci, channel in enumerate(sample.get('channel') or [{}]):
            cp = base + f'/channel/{ci}' if sample.get('channel') else base + '/channel'
            entity = str(identity)
            self.check([channel.get('source'), *[s.get('name') for s in biological if s['kind'] == 'source']],
                       'sample.source', cp + '/source', 'sample', entity)
            chars = [(v, cp + f'/characteristics/{j}') for j, v in enumerate(channel.get('characteristics', []))]
            for step in biological:
                # Semantic paths are additional authored evidence, never other samples.
                chars.extend((v, self.step_locations[id(step)] + f'/characteristics/{j}')
                             for j, v in enumerate(step.get('characteristics', [])))
            self.declared_values(channel.get('characteristics', []), cp + '/characteristics',
                                 'sample', entity, factor_names=factor_names)
            by_name = defaultdict(list)
            for value, path in chars:
                name = value.get('name') or value.get('tag')
                if name:
                    by_name[_slug(name)].append((value.get('value'), path + '/value'))
            for field in _BIOLOGY:
                names = sorted(BIOLOGICAL_ALIASES.get(field, {field}))
                candidates = [v for n in names for v in by_name[n]]
                if field == 'organism':
                    candidates += [(v, cp + '/organism') for v in _items(channel.get('organism'))]
                # Explicit empty occurrences were assessed at their own paths.
                if not candidates or field == 'organism' and not by_name[field]:
                    self.check([v for v, _ in candidates], 'sample.' + field,
                               candidates[0][1] if candidates else cp + '/characteristics', 'sample', entity)
            self.check([channel.get('molecule'), *[s.get('material_type') for s in steps if s.get('kind') == 'extract']],
                       'sample.molecule', cp + '/molecule', 'sample', entity)
            factor_values = defaultdict(list)
            for step in steps:
                for v in step.get('factor_values', []):
                    factor_values[_slug(v['name'])].append(v.get('value'))
            for factor in filter(None, factors):
                values = [v for v, _ in by_name[_slug(factor)]] + factor_values[_slug(factor)]
                if not by_name[_slug(factor)] and not factor_values[_slug(factor)]:
                    self.check(values, 'factor.' + _slug(factor), cp + '/characteristics', 'sample', entity)
            for category, values, location in (
                ('collection', self.protocol_text(steps, 'collection') + [channel.get('sample_collection_protocol')], cp),
                ('extraction', self.protocol_text(steps, 'extraction') + [channel.get('extract_protocol')], cp + '/extract_protocol')):
                self.check(values, 'protocol.' + category, location, 'sample', entity)
            self.assays(sample, channel, base, cp, entity, steps)

    @staticmethod
    def processed(entity, steps):
        from ..magetab.sdrf.handlers.base import classify_file
        return (any(s.get('kind') in _PROCESSED_NODES for s in steps)
                or not missing_reason(entity.get('data_table'))
                or any(classify_file(str(_value(v))) == 'matrix_or_derived'
                       for v in entity.get('supplementary_data', [])))

    def assays(self, sample, channel, base, cp, identity, steps):
        from ..magetab.technology import resolve_technology
        from ..magetab.chemistry import resolve_chemistry
        from ..magetab.sdrf.handlers.base import classify_file
        runs = sample.get('sra_run') or []
        platform = self.platforms.get((sample.get('platform_ref') or {}).get('ref'), {})
        scopes = [(run, base + f'/sra_run/{i}', f'{identity}/run[{run.get("run") or i}]',
                   'run', self.scoped_steps(sample, run), False) for i, run in enumerate(runs)]
        covered = {id(s) for _, _, _, _, group, _ in scopes for s in group}
        # Imported MAGE-TAB can declare several assays on the same sample without
        # SRA records. Respect those explicit nodes before applying legacy routing.
        for path in self.paths.get(sample.get('iid'), ()):
            declared = [s for s in path['steps'] if s.get('technology_type') or s.get('array_design_ref')]
            if declared and not any(id(s) in covered for s in path['steps']):
                assay = declared[-1]
                scopes.append((None, self.step_locations[id(assay)], f'{identity}/assay[{assay["name"]}]',
                               'assay', path['steps'], True))
        if not scopes:
            scopes = [(None, base, identity, 'sample', steps, False)]
        for run, rp, entity, kind, scoped, native in scopes:
            decision = resolve_technology(sample, channel, run, data=self.data)
            handler = decision.handler
            technologies = {_slug(_value(s['technology_type'])) for s in scoped if s.get('technology_type')}
            if technologies == {'array_assay'} or (not technologies and any(s.get('array_design_ref') for s in scoped)):
                handler = 'array'
            elif (technologies == {'sequencing_assay'} and 'sequencing' not in handler) or (run is not None and handler == 'generic'):
                handler = 'sequencing'
            comments = defaultdict(list)
            for step in scoped:
                for c in step.get('comments', []):
                    comments[_slug(c.get('name'))].append(c.get('value'))
            if run is not None:
                self.check(run.get('run'), 'run.identifier', rp + '/run', kind, entity)
                self.publications(run, rp, kind, entity)
                self.companions(run, rp, kind, entity, skip=('pubmed_publication', 'relation'))
            if handler == 'array':
                self.check([sample.get('platform_ref'), *[s.get('array_design_ref') for s in scoped]],
                           'assay.array_design', rp + '/array_design_ref' if native else base + '/platform_ref', kind, entity)
                for category, value, path in (
                    ('hybridization', sample.get('hybridization_protocol'), base + '/hybridization_protocol'),
                    ('scan', sample.get('scan_protocol'), base + '/scan_protocol')):
                    self.check([value, *self.protocol_text(scoped, category)], 'protocol.' + category, path, kind, entity)
                if not missing_reason(channel.get('label')):
                    self.check([channel.get('label_protocol'), *self.protocol_text(scoped, 'labeling')],
                               'protocol.labeling', cp + '/label_protocol', kind, entity)
            elif 'sequencing' in handler:
                lookup_context = {'lookup_status': 'disabled' if self.runs_disabled else None}
                if run is None:
                    named_runs = [s['name'] for s in scoped if s.get('kind') == 'scan'
                                  and re.fullmatch(r'[SED]RR\d+', s['name'])]
                    self.check([*comments['ena_run'], *comments['run_accession'], *named_runs], 'run.identifier',
                               rp + '/comments' if native else base + '/sra_run', kind, entity, **lookup_context)
                for field in _LIBRARY:
                    self.check([sample.get(field), (run or {}).get(field), *comments[field]],
                               'assay.' + field, rp + ('/comments' if native else '/' + field), kind, entity, **lookup_context)
                self.check([sample.get('instrument_model'), platform.get('title'), (run or {}).get('instrument_model'),
                            *comments['instrument_model']], 'assay.instrument_model',
                           rp + ('/comments' if native else '/instrument_model'), kind, entity, **lookup_context)
                from .preparation import bound_protocols
                library_text = [text for _, text in bound_protocols(sample, run, self.series)
                                if re.search(r'\blibrar(?:y|ies)\b', text, re.I)]
                extraction = channel.get('extract_protocol')
                if extraction and re.search(r'\blibrar(?:y|ies)\b', extraction, re.I):
                    library_text.append(extraction)
                self.check([*self.protocol_text(scoped, 'library_preparation'), *library_text,
                            (run or {}).get('library_protocol'), (run or {}).get('library_construction_protocol'),
                            ], 'protocol.library_preparation', rp, kind, entity)
                self.check(self.protocol_text(scoped, 'sequencing'), 'protocol.sequencing', rp, kind, entity)
                if 'single_cell' in handler or 'spatial' in handler:
                    chemistry = resolve_chemistry(sample, channel=channel, run=run, series=self.series)
                    if chemistry.manufacturer == '10x Genomics':
                        self.check(chemistry.family, 'assay.chemistry', rp, kind, entity)
            raw = _items((run or {}).get('fastq_files')) + _items((run or {}).get('files'))
            raw += _items((run or {}).get('submitted_file_name'))
            if len(scopes) == 1:
                raw += _items(sample.get('raw_data'))
            raw += [s.get('link') or s.get('name') for s in scoped if s.get('kind') in _RAW_NODES]
            for link in sample.get('supplementary_data', []) if len(scopes) == 1 else []:
                value = link.get('value') if isinstance(link, Mapping) else link
                if value and classify_file(str(value)) in {'sequencing_raw', 'array_raw'}:
                    raw.append(value)
            references = [v.get('uri') or v.get('filename') or v.get('value') if isinstance(v, Mapping) else v for v in raw]
            if handler != 'generic' or raw:
                self.check(references, 'data.raw_file', rp if native else rp + ('/fastq_files' if run is not None else '/raw_data'), kind, entity)
            if self.processed(sample, scoped):
                self.check([sample.get('data_processing'), *self.protocol_text(scoped, 'processing')],
                           'protocol.processing', base + '/data_processing', kind, entity)

    def lookup_findings(self):
        for record in self.preparation:
            if record.get('operation') != 'run_lookup':
                continue
            status = record['status']
            # Failures already have their correlated operational diagnostic.
            # Publication outcomes are attached to their explicit citations above.
            if status not in {'incomplete', 'not_found'}:
                continue
            entity = record.get('entity_id') or self.dataset_id
            location = record.get('location', '/series/pubmed_publication')
            self.add('lookup.metadata', 'lookup_no_match' if status == 'not_found' else 'lookup_incomplete',
                     location, 'lookup', entity + '/' + record['identifier'], lookup=record,
                     action='Verify the requested identifier and supply any metadata still absent from the final record.')


@preparation_operation
def assess_curation(package, *, dataset_id, provider=None, package_index=0, preparation=()):
    """Return advisories for one validated prepared package; perform no retrieval."""
    data = package.to_mapping() if hasattr(package, 'to_mapping') else package
    audit = _Assessment(data, dataset_id, provider, package_index, preparation)
    audit.study()
    for index, sample in enumerate(data.get('sample', [])):
        audit.sample(sample, index)
    audit.lookup_findings()
    return audit.result
