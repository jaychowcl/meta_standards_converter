# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Occurrence-scoped residual projection for the INSDC 2.0 extension.

Source trees remain private while an import is enriched. Serialization contains
only unmatched fields. Rules bind source paths to entity fields, never to a
package-wide collection of equal strings.
"""
from copy import deepcopy
import json
import re
from urllib.parse import urlsplit, unquote


def children(node, tag):
    return [c for c in node.get('children', []) if c.get('tag') == tag]


def child_text(node, path):
    current = node
    for part in path.split('/'):
        current = next(iter(children(current, part)), {})
    return str(current.get('text') or '').strip()


def all_text(node):
    return str(node.get('text') or '') + ''.join(all_text(c) + str(c.get('tail') or '') for c in node.get('children', []))


def accession(node):
    return node.get('attributes', {}).get('accession') or child_text(node, 'IDENTIFIERS/PRIMARY_ID')


def contains(expected, actual):
    """A coupled source value is represented by a compatible richer value."""
    if isinstance(expected, dict):
        return isinstance(actual, dict) and all(k in actual and
            (isinstance(actual[k], list) and len(v) == len(actual[k])
             and all(contains(a, b) for a, b in zip(v, actual[k]))
             if k == 'steps' and isinstance(v, list) else contains(v, actual[k]))
            for k,v in expected.items() if v not in (None, '', [], {}))
    if isinstance(expected, list):
        if not isinstance(actual, list): return False
        available = list(actual)
        for value in expected:
            match = next((i for i, item in enumerate(available) if contains(value, item)), None)
            if match is None: return False
            available.pop(match)
        return True
    return expected == actual


def _path_signature(path):
    if not isinstance(path, dict) or not isinstance(path.get('steps'), list):
        return None
    return json.dumps([{k: s[k] for k in ('kind', 'name', 'sample_ref', 'protocol_ref') if k in s}
                       | ({'uri': s['link']['value']} if isinstance(s.get('link'), dict) and 'value' in s['link'] else {})
                       for s in path['steps']], sort_keys=True)


def _path_diff(source, target):
    # Equal ordered identities establish occurrence position, including repeated
    # applications. Keep the skeleton only when an unmapped sibling remains.
    steps = []
    changed = False
    for item, other in zip(source['steps'], target['steps']):
        delta = diff(item, other)
        changed |= delta is not None
        steps.append(delta if delta is not None else
                     {k: deepcopy(item[k]) for k in ('kind', 'name', 'sample_ref', 'protocol_ref') if k in item})
    result = {k: delta for k, value in source.items() if k != 'steps'
              and (delta := diff(value, target.get(k), field=k)) is not None}
    if changed or result:
        result['steps'] = steps
    return result or None


def _mapped_date_view(metadata, target, workflows=False):
    """Prune original occurrences proved by same-sample administrative projections."""
    from .archive_dates import _date_key
    from .archive_administration import _NAMES
    projection = Projection(target)
    view = deepcopy(metadata)
    samples = {s['iid']: s for s in target.get('sample', [])}

    def clean(owner, sample, used):
        retained = []
        for item in owner.get('characteristics', []):
            if str(item.get('name') or '').strip().lower() in _NAMES:
                unit = item.get('unit', {}).get('value') if isinstance(item.get('unit'), dict) else item.get('unit')
                terms = {k:v for k,v in item.items() if k not in {'name','value','unit'}}
                if projection.character(sample, item.get('name'), item.get('value'), unit, terms):
                    continue
            key = _date_key(item.get('name'), item.get('value'))
            match = None
            if key:
                database, field, value = key
                wanted = {'database': database, field: value}
                extra = {k: v for k, v in item.items() if k not in ('name', 'value')}
                if extra:
                    wanted['attribute_annotations'] = [{**extra, 'attribute_name': item['name']}]
                match = next(((i, field) for i, status in enumerate(sample.get('status', []))
                              if (i, field) not in used and contains(wanted, status)), None)
            if match is None:
                retained.append(item)
            else:
                used.add(match)
        if 'characteristics' in owner:
            owner['characteristics'] = retained

    if not workflows:
        for sample in view.get('sample', []):
            used = set()
            projection._mapped_characters.clear()
            for channel in sample.get('channel', []):
                clean(channel, samples.get(sample.get('iid'), {}), used)
    series = view if workflows else view.get('series', {})
    for path in series.get('assay_paths', []):
        bound = {s.get('sample_ref') for s in path.get('steps', [])} - {None}
        for step in path.get('steps', []):
            ref = step.get('sample_ref') or (next(iter(bound)) if len(bound) == 1 else None)
            if ref in samples:
                projection._mapped_characters.clear()
                clean(step, samples[ref], set())
    return view


def diff(source, target, *, field=None):
    """Residual MINiML fields, using scoped identities and coherent value groups."""
    if contains(source, target): return None
    if isinstance(source, dict):
        if any(k in source for k in ('value', 'protocol_ref', 'link', 'unit')):
            return deepcopy(source)
        result = {k: delta for k,v in source.items() if k not in ('miniml_schema_version', 'extensions')
                  and not (k == 'source' and field is None)
                  and (delta := diff(v, target.get(k) if isinstance(target,dict) else None, field=k)) is not None}
        if result:
            for k in ('iid', 'kind', 'name', 'sample_ref'):
                if k in source: result.setdefault(k, source[k])
        return result or None
    if isinstance(source, list):
        result = []
        available = list(target) if isinstance(target,list) else []
        source_paths = [_path_signature(p) for p in source] if field == 'assay_paths' else []
        target_paths = [_path_signature(p) for p in available] if field == 'assay_paths' else []
        for item in source:
            match = next((i for i,x in enumerate(available) if contains(item,x)), None)
            if match is not None:
                available.pop(match)
                continue
            if field == 'assay_paths':
                signature = _path_signature(item)
                if signature is not None and source_paths.count(signature) == target_paths.count(signature) == 1:
                    match = next((i for i, p in enumerate(available) if _path_signature(p) == signature), None)
                    if match is not None:
                        delta = _path_diff(item, available.pop(match))
                        if delta is not None: result.append(delta)
                        continue
            candidate = None
            # Native channel projection is unambiguous only for the single
            # channel of an already matched sample, never by biological title.
            if field == 'channel' and len(source) == len(available) == 1:
                candidate = available[0]
            if isinstance(item,dict):
                identity = next((k for k in ('iid','run','name') if item.get(k)), None)
                if identity:
                    options = [x for x in available if isinstance(x,dict) and x.get(identity)==item[identity]]
                    if len(options)==1: candidate=options[0]
            delta = diff(item, candidate, field=field)
            if delta is not None: result.append(delta)
            if candidate is not None:
                available.remove(candidate)
        return result or None
    return deepcopy(source) if source not in (None, '', [], {}) else None


class Projection:
    def __init__(self, data, records=()):
        self.ranges = [r for r in records if r["kind"] == "accession_range"]
        self.data = data
        self.series = data['series']
        self.samples = {a['value']: s for s in data.get('sample', []) for a in s.get('accession', [])}
        self.samples.update({s['iid']: s for s in data.get('sample', [])})
        self.runs = {}
        self.experiments = {}
        for s in data.get('sample', []):
            for run in s.get('sra_run', []):
                self.runs[run['run']] = run
                self.experiments.setdefault(run.get('experiment'), []).append(run)
        self.paths = self.series.get('assay_paths', [])
        from .archive_factors import ExperimentBindings
        bindings = ExperimentBindings(data)
        self.factor_nodes = {}
        for assay_path in self.paths:
            experiment = bindings.experiment(assay_path)
            if experiment:
                self.factor_nodes.setdefault(experiment, []).extend(
                    node['factor_values'] for node in assay_path['steps'] if node.get('factor_values'))
        self.actors = {v['iid']: v for k in ('organization','contributor') for v in data.get(k, [])}
        self.organizations = {v['iid']:v for v in data.get('organization', [])}
        for organization in data.get('organization', []):
            for occurrence in organization.get('source_occurrences', []):
                self.actors[occurrence['iid']] = {**organization, **occurrence}
        self.publications = {str(p['pubmed_id']): p for p in self.series.get('pubmed_publication', [])}
        self.ids = {a['value'] for a in self.series.get('accession', [])}
        self._mapped_characters = set()
        self.all_publications = [p for entity in [self.series, *data.get('sample', []), *self.runs.values()]
                                 for p in entity.get('pubmed_publication', [])]
        self.all_publications.extend(r['publication'] for entity in [self.series, *data.get('sample', []), *self.runs.values()]
                                     for r in entity.get('relation', []) if r.get('publication'))
        self.publications = {str(p['pubmed_id']):p for p in self.all_publications if p.get('pubmed_id')}

    def citations(self, acc):
        target = self.series if acc in self.ids else self.samples.get(acc, self.runs.get(acc, {}))
        values = list(target.get('pubmed_publication', []))
        values.extend(r['publication'] for r in target.get('relation', []) if r.get('publication'))
        values.extend(r['publication'] for r in self.series.get('relation', []) if r.get('publication')
                      and acc in (r.get('experiment_ref'),r.get('assembly_ref'),r.get('analysis_ref')))
        return values

    def entity(self, kind, acc):
        if kind in ('SAMPLE','BioSample','sample'): return self.samples.get(acc, {})
        if kind in ('RUN','read_run'): return self.runs.get(acc, {})
        if kind in ('STUDY','PROJECT','DocumentSummary','study') and acc in self.ids: return self.series
        return {}

    def organisms(self):
        return [o for s in self.data.get('sample', []) for c in s.get('channel', []) for o in c.get('organism', [])]

    def character(self, sample, name, value, unit=None, terms=None):
        from .archive_administration import administrative_destination, _CENTERS
        wanted = {'name':name, 'value':value}
        if unit: wanted['unit'] = {'value':unit}
        if terms: wanted.update(terms)
        destination = administrative_destination(wanted, self.data.get('source', {}).get('format'))
        if destination:
            field, record = destination
            for index, actual in enumerate(sample.get(field, [])):
                occurrence = (sample.get('iid'), 'administration', field, index)
                if contains(record, actual) and occurrence not in self._mapped_characters:
                    self._mapped_characters.add(occurrence)
                    return True
        from .archive_dates import _date_key
        date_key = _date_key(name, value)
        if date_key and not unit and not terms:
            database, field, literal = date_key
            for index, status in enumerate(sample.get('status', [])):
                occurrence = (sample.get('iid'), 'status', index, field)
                if status.get('database') == database and status.get(field) == literal and occurrence not in self._mapped_characters:
                    self._mapped_characters.add(occurrence)
                    return True
        if name and name.strip().lower() == 'insdc status' and not unit and not terms:
            for index, status in enumerate(sample.get('status', [])):
                occurrence = (sample.get('iid'), 'administration', index)
                if status.get('database') == 'INSDC' and {'name':name, 'value':value} in status.get('comment', []) and occurrence not in self._mapped_characters:
                    self._mapped_characters.add(occurrence)
                    return True
        if name and name.strip().lower() in _CENTERS and not unit and not terms:
            for relation in sample.get('relation', []):
                organization = self.organizations.get(relation.get('target'), {})
                for index, original in enumerate(organization.get('source_occurrences') or [organization]):
                    actor = {**organization, **original}
                    scope = actor.get('sample_accession')
                    if scope and scope not in {sample.get('iid'), *[a['value'] for a in sample.get('accession', [])]}:
                        continue
                    occurrence = ('administration', organization.get('iid'), index)
                    field = 'alias' if name.strip().lower().endswith('alias') else 'name'
                    if name in [original.get('role'), *original.get('roles', [])] and actor.get(field) == value and occurrence not in self._mapped_characters:
                        self._mapped_characters.add(occurrence)
                        return True
        wanted = {'name': name, 'value': value}
        if unit: wanted['unit'] = {'value': unit}
        if terms: wanted.update(terms)
        for channel_index, channel in enumerate(sample.get('channel', [])):
            for item_index, item in enumerate(channel.get('characteristics', [])):
                occurrence = (sample.get('iid'), channel_index, item_index)
                if occurrence in self._mapped_characters: continue
                represented = contains(wanted, item)
                # A supplied value may have been replaced by an explicit linked
                # value/unit pair (e.g. Age "1 days" -> 1 + days).
                if item.get('name') == name and not unit and item.get('unit', {}).get('value'):
                    represented |= value == str(item.get('value')) + ' ' + item['unit']['value']
                if represented:
                    self._mapped_characters.add(occurrence)
                    return True
        return False

    def factor(self, experiment, wanted):
        def matches(expected, actual):
            if contains(expected, actual):
                return True
            return (expected.get('name') == actual.get('name') and not expected.get('unit')
                    and bool(actual.get('unit', {}).get('value'))
                    and expected.get('value') == str(actual.get('value')) + ' ' + actual['unit']['value'])

        def facts(value):
            # Definition annotations can enrich a projection; source factors here
            # supply only their name, value and optional literal unit.
            return {k: value[k] for k in ('name', 'value') if k in value} | (
                {'unit': {'value': value['unit']['value']}} if value.get('unit', {}).get('value') else {})

        def available(node_index, values):
            return [(i, value) for i, value in enumerate(values)
                    if (experiment, 'factor', node_index, i) not in self._mapped_characters]

        nodes = self.factor_nodes.get(experiment, [])
        candidates = [(not contains(wanted, value), facts(value))
                      for ni, values in enumerate(nodes) for _, value in available(ni, values)
                      if matches(wanted, value)]
        if not candidates:
            return False
        # Prefer exact facts before a literal/value+unit projection. Reserve one
        # equivalent occurrence on each node, so copies never multiply capacity.
        witness = min(candidates, key=lambda item: item[0])[1]
        for ni, values in enumerate(nodes):
            equivalent = [(facts(value) != witness, i) for i, value in available(ni, values)
                          if matches(witness, facts(value)) or matches(facts(value), witness)]
            if equivalent:
                _, index = min(equivalent)
                self._mapped_characters.add((experiment, 'factor', ni, index))
        return True

    def relations(self, entity):
        return entity.get('relation', [])

    def reference(self, value, kind, acc):
        entity = self.entity(kind, acc)
        if value in {a['value'] for a in entity.get('accession', [])}: return True
        if value in [r.get('target') for r in entity.get('relation', [])]: return True
        runs = [entity] if kind in ('RUN','read_run') else self.experiments.get(acc, [])
        if any(value == r.get(k) for r in runs for k in ('run','study','sample','biosample','experiment','geo_sample')): return True
        return False

    def xml(self, node, kind, acc, provider, path=(), owner=None, actor=None, file=None, citation=None):
        tag = node['tag']; attrs=node.get('attributes', {})
        text=str(node.get('text') or '').strip()
        entity = self.entity(kind, acc)
        if owner is None: owner=(kind, attrs.get('uid') or acc)
        node_acc=accession(node) or (attrs.get('uid') if tag=='DocumentSummary' else None)
        if tag=='EXPERIMENT_PACKAGE': node_acc=accession(next(iter(children(node,'EXPERIMENT')),{})) or acc
        if node_acc: owner, path=(tag,node_acc), ()
        actor_id = f'{provider}:{owner[0]}:{owner[1]}:organization-' + '-'.join(map(str,path))
        if tag=='Organization' or (tag=='Owner' and owner[0]=='BioSample'): actor=self.actors.get(actor_id, {})
        if tag=='Contact' and actor is not None:
            # Contacts are numbered within their organization's Contact list.
            actor=self.actors.get(actor.get('iid','') + ':contact-' + str(getattr(self,'contact_index',0)), {})
        runs = self.experiments.get(acc, []) if kind=='EXPERIMENT' else [entity] if kind in ('RUN','read_run') else []
        if tag=='RUN' and entity.get('scan_name')==attrs.get('alias'): mapped_alias=True
        else: mapped_alias=False
        if tag=='PubmedArticle': self.publication=self.publications.get(child_text(node,'MedlineCitation/PMID'),{})
        if tag=='Statistics' and any(r.get('statistics') == node for r in runs): return None
        mapped_text=False; mapped_attrs={'alias'} if mapped_alias else set(); drop=False
        from ..sources.archive_publications import citation_identifier
        if tag=='Publication':
            wanted = citation_identifier(child_text(node,'DbType'),attrs.get('id'))
            citation = next((p for p in (self.citations(owner[1]) or self.citations(acc)) if wanted and contains(wanted,p)), {})
            if citation:
                mapped_attrs.add('id')
        if citation:
            field = {'Title':'title','AuthorList':'author_list','DOI':'doi'}.get(tag)
            if field and citation.get(field)==text: mapped_text=True
            if tag=='DbType': mapped_text=True
        if tag=='Link':
            wanted = citation_identifier(attrs.get('type') or attrs.get('target'),text)
            if wanted and any(contains(wanted,p) for p in self.citations(owner[1])):
                mapped_text=True; mapped_attrs.update(('type','target'))
        if kind=='PubmedArticle':
            pub=getattr(self,'publication',{})
            if tag == 'PublicationStatus':
                from ..metadata.ontology_mappings import Harmonizer
                expected = Harmonizer().pubstatus2efo(text)
                if expected[0] and all(pub.get(k) == v for k, v in zip(
                        ('status', 'status_term_source_ref', 'status_term_accession_number'), expected)):
                    mapped_text = True
            if tag=='PMID' and pub.get('pubmed_id')==text:mapped_text=True;mapped_attrs.update(('Version',))
            if tag=='ArticleTitle' and pub.get('title')==all_text(node).strip():return None
            if tag=='ArticleId' and attrs.get('IdType')=='doi' and pub.get('doi')==text:mapped_text=True;mapped_attrs.add('IdType')
            if tag=='AuthorList':
                names=[]
                for author in children(node,'Author'):
                    names.append(child_text(author,'CollectiveName') or ' '.join(filter(None,[child_text(author,'ForeName'),child_text(author,'LastName')])))
                if pub.get('author_list')==', '.join(names):
                    node=deepcopy(node)
                    for author in node.get('children',[]):
                        author['children']=[c for c in author.get('children',[]) if c['tag'] not in ('ForeName','LastName','CollectiveName')]
        if kind=='TaxaSet' and tag=='Taxon':
            taxid=child_text(node,'TaxId');name=child_text(node,'ScientificName')
            if any(o.get('taxid')==taxid and o.get('value')==name for o in self.organisms()):
                node=deepcopy(node);node['children']=[c for c in node.get('children',[]) if c['tag']!='ScientificName']
        center = self.actors.get(f'{provider}:{tag}:{node_acc}:center_name', {})
        if attrs.get('center_name') and center.get('name')==attrs['center_name']:mapped_attrs.add('center_name')
        if tag in ('SAMPLE_ATTRIBUTE','Attribute'):
            name = child_text(node,'TAG') if tag=='SAMPLE_ATTRIBUTE' else attrs.get('attribute_name') or attrs.get('display_name')
            value = child_text(node,'VALUE') if tag=='SAMPLE_ATTRIBUTE' else text
            unit = child_text(node,'UNITS') if tag=='SAMPLE_ATTRIBUTE' else attrs.get('unit')
            if self.character(entity,name,value,unit):
                if tag=='SAMPLE_ATTRIBUTE':
                    node=deepcopy(node);node['children']=[c for c in node.get('children',[]) if c['tag'] not in ('TAG','VALUE','UNITS')]
                    drop=not node['children'] and not attrs
                else: mapped_text=True;mapped_attrs.update(('attribute_name','display_name','unit'))
        if tag == 'EXPERIMENT_ATTRIBUTE' and owner[0] == 'EXPERIMENT':
            from .archive_factors import declared_factor
            wanted = declared_factor(child_text(node, 'TAG'), child_text(node, 'VALUE'), child_text(node, 'UNITS'))
            if wanted and self.factor(owner[1], wanted):
                node = deepcopy(node)
                left = []
                for child in node.get('children', []):
                    if child['tag'] in ('TAG', 'VALUE', 'UNITS'):
                        child.pop('text', None)
                        if not child.get('attributes') and not child.get('children') and not child.get('tail', '').strip():
                            continue
                    left.append(child)
                # Keep the factor name to interpret unmapped occurrence siblings.
                if left or attrs:
                    identity = next((c for c in left if c['tag'] == 'TAG'), None)
                    if identity is None:
                        identity = {'tag': 'TAG'}
                        left.insert(0, identity)
                    identity['text'] = 'Experimental Factor: ' + wanted['name']
                node['children'] = left
                drop = not node['children'] and not attrs
        if tag == 'STUDY_ATTRIBUTE':
            name, value = child_text(node, 'TAG'), child_text(node, 'VALUE')
            date_field = {'ENA-FIRST-PUBLIC': 'release_date', 'ENA-LAST-UPDATE': 'last_update_date'}.get(name)
            represented = bool(date_field and any(s.get('database') == 'ENA' and s.get(date_field) == value
                                                  for s in entity.get('status', [])))
            represented |= {'type': name, 'target': value} in self.relations(entity)
            if represented:
                node = deepcopy(node)
                node['children'] = [c for c in node.get('children', []) if c['tag'] not in ('TAG', 'VALUE')]
                drop = not node['children'] and not attrs
        if tag in ('XREF_LINK','URL_LINK'):
            from .reference_targets import parse_reference_targets
            db=child_text(node,'DB') if tag=='XREF_LINK' else child_text(node,'LABEL') or 'external'
            value=child_text(node,'ID') if tag=='XREF_LINK' else child_text(node,'URL')
            ranges = [r['metadata'] for r in self.ranges if r['accession'] == acc]
            parts = parse_reference_targets(db, value, ranges)
            unmapped = [p for p in parts if not any(r.get('type') == db and r.get('target') == p for r in self.relations(entity))]
            if len(unmapped) < len(parts):
                node = deepcopy(node)
                kept = []
                for child in node.get('children', []):
                    if child['tag'] in (('DB','ID') if tag == 'XREF_LINK' else ('LABEL','URL')):
                        if child['tag'] in ('ID','URL') and unmapped:
                            child['text'] = ','.join(unmapped)
                        elif not unmapped:
                            child.pop('text', None)
                        if not child.get('text') and not child.get('attributes') and not child.get('children') and not child.get('tail', '').strip():
                            continue
                    kept.append(child)
                node['children'] = kept
                drop = not kept and not attrs
        if tag in ('PRIMARY_ID','SECONDARY_ID','EXTERNAL_ID') and self.reference(text,kind,acc): mapped_text=True;mapped_attrs.update(('namespace','label'))
        field={'STUDY_TITLE':'title','STUDY_ABSTRACT':'summary','STUDY_DESCRIPTION':'summary','TITLE':'title','DESCRIPTION':'description'}.get(tag)
        if kind in ('PROJECT','STUDY','DocumentSummary','study') and tag in ('DESCRIPTION','Description'): field='summary'
        if tag=='Title' and kind=='DocumentSummary': field='title'
        if field and text and ' '.join(text.split()) == ' '.join(str(entity.get(field) or '').split()): mapped_text=True
        channel=entity.get('channel',[{}])[0]
        if kind=='BioSample' and tag=='Title' and entity.get('title')==text: mapped_text=True
        if kind=='BioSample' and tag=='OrganismName' and any(o.get('value')==text for o in channel.get('organism',[])): mapped_text=True
        if kind=='BioSample' and tag=='Id' and self.reference(text,kind,acc): mapped_text=True;mapped_attrs.update(('db','is_primary'))
        if tag in ('SCIENTIFIC_NAME','TAXON_ID') and kind=='SAMPLE':
            key='value' if tag=='SCIENTIFIC_NAME' else 'taxid'
            if text and any(text==o.get(key) for o in channel.get('organism',[])): mapped_text=True
        if kind=='BioSample' and tag=='Organism':
            if any(o.get('taxid')==attrs.get('taxonomy_id') and o.get('value')==attrs.get('taxonomy_name') for o in channel.get('organism',[])):
                mapped_attrs.update(('taxonomy_id','taxonomy_name'))
        if kind=='BioSample':
            if tag == 'BioSample' and any(s.get('database') == 'BioSample' and {'name':'access','value':attrs.get('access')} in s.get('comment', []) for s in entity.get('status', [])): mapped_attrs.add('access')
            if tag == 'Status' and any(s.get('database') == 'BioSample' and {'name':'record status','value':attrs.get('status')} in s.get('comment', []) for s in entity.get('status', [])): mapped_attrs.add('status')
            for attr,field in [('publication_date','release_date'),('submission_date','submission_date'),('last_update','last_update_date')]:
                if attrs.get(attr) and any(s.get(field)==attrs[attr] and s.get('database')=='BioSample' for s in entity.get('status',[])):mapped_attrs.add(attr)
        if tag=='STUDY_TYPE' and any(t.get('value')==attrs.get('existing_study_type') for t in entity.get('type',[])):mapped_attrs.add('existing_study_type')
        lib={'LIBRARY_STRATEGY':'library_strategy','LIBRARY_SOURCE':'library_source','LIBRARY_SELECTION':'library_selection','INSTRUMENT_MODEL':'instrument_model'}
        if tag in lib and runs and all(r.get(lib[tag])==text for r in runs):mapped_text=True
        if tag in ('PAIRED','SINGLE') and runs and all(r.get('library_layout')==tag for r in runs):drop=not attrs
        if tag=='DESIGN_DESCRIPTION':
            if any(s.get('kind')=='assay' and s.get('name')==acc and s.get('description')==text for p in self.paths for s in p['steps']):mapped_text=True
        if tag in ('LIBRARY_CONSTRUCTION_PROTOCOL','PROTOCOL','METHOD','ASSEMBLY_METHOD'):
            if any(p.get('description')==text and p['name'].split(':')[0]==acc for p in self.series.get('protocols',[])):mapped_text=True
            if tag == 'LIBRARY_CONSTRUCTION_PROTOCOL':
                from .archive_protocols import library_description
                designs = {s.get('description') for p in self.paths for s in p['steps']
                           if s.get('kind') == 'assay' and s.get('name') == acc}
                combined = {library_description(text, design) for design in designs}
                if any(p.get('description') in combined and p['name'] == acc + ':library'
                       for p in self.series.get('protocols', [])):
                    mapped_text = True
        if tag=='SRAFile':
            file=[f for r in runs for f in r.get('files',[]) if f.get('filename')==attrs.get('filename')]
            for source,target in [('filename','filename'),('size','bytes'),('md5','md5'),('semantic_name','format'),('supertype','role'),('url','uri')]:
                if attrs.get(source) and any(f.get(target)==attrs[source] for f in file):mapped_attrs.add(source)
        if tag=='Alternatives' and file:
            for key,value in attrs.items():
                if any(f.get('uri')==attrs.get('url') and f.get(key)==value for f in file):mapped_attrs.add(key)
        if tag=='FILE':
            files=[f for r in runs for f in r.get('files',[])]
            targets={s['iid'] for s in self.data.get('sample',[]) if any(r.get('target')==acc for r in s.get('relation',[]))}
            links=[s.get('link',{}) for p in self.paths if any(n.get('sample_ref') in targets for n in p['steps']) for s in p['steps']]
            if any(r.get('target')==acc for r in self.series.get('relation',[])): links.extend(self.series.get('supplementary_data',[]))
            for source,target in [('filename','filename'),('filetype','format'),('checksum','checksum'),('checksum_method','checksum_method')]:
                if attrs.get(source) and any(f.get('filename')==attrs.get('filename') and f.get(target)==attrs[source] for f in files):mapped_attrs.add(source)
            matching_links = [l for l in links if l.get('value')==attrs.get('filename') or
                              (attrs.get('filename') and unquote(urlsplit(l.get('value','')).path).endswith('/'+attrs['filename']))]
            if len({l.get('value') for l in matching_links}) == 1:
                mapped_attrs.add('filename')
                if any(str(l.get('type','')).casefold()==str(attrs.get('filetype','')).casefold() for l in matching_links): mapped_attrs.add('filetype')
                if attrs.get('checksum_method','').upper()=='MD5' and any(l.get('checksum')==attrs.get('checksum') for l in matching_links):mapped_attrs.update(('checksum','checksum_method'))
        if tag in ('STUDY_REF','SAMPLE_REF','EXPERIMENT_REF','RUN_REF','SAMPLE_DESCRIPTOR'):
            value=accession(node)
            if self.reference(value,kind,acc):mapped_attrs.add('accession')
            if kind in ('ASSEMBLY','ANALYSIS') and tag=='SAMPLE_REF' and any(r.get('target')==acc for r in self.samples.get(value,{}).get('relation',[])):
                mapped_attrs.add('accession')
                node=deepcopy(node);node['children']=[c for c in node.get('children',[]) if c['tag']!='IDENTIFIERS']
        if actor is not None:
            if tag=='Organization':
                for source,target in [('url','web_link'),('role','role'),('type','type')]:
                    if attrs.get(source) and actor.get(target)==attrs[source]:mapped_attrs.add(source)
            if tag=='Contact':
                for source,target in [('email','email'),('phone','phone'),('fax','fax'),('url','web_link')]:
                    if attrs.get(source) and actor.get(target)==attrs[source]:mapped_attrs.add(source)
                if attrs.get('sec_email') and actor.get('extensions', {}).get('secondary_email')==attrs['sec_email']:mapped_attrs.add('sec_email')
                if attrs.get('role') and {'value':attrs['role']} in actor.get('roles',[]):mapped_attrs.add('role')
            if tag=='Name' and actor.get('name')==text:mapped_text=True
            if tag=='Name' and attrs.get('url') and actor.get('web_link')==attrs['url']:mapped_attrs.add('url')
            if tag in ('First','Middle','Last') and actor.get('person',{}).get(tag.lower())==text:mapped_text=True
            address = actor.get('address')
            if isinstance(address, dict):
                if text and (text in address.get('line',address.get('lines',[])) or (tag in ('City','Country') and address.get(tag.lower())==text)): mapped_text=True
                if tag=='Address' and attrs.get('postal_code') and address.get('postal_code')==attrs['postal_code']: mapped_attrs.add('postal_code')
        # Identity is carried by the record envelope, not counted as residual data.
        if not path: mapped_attrs.update(('accession','uid'))
        out={'tag':tag}; left={k:v for k,v in attrs.items() if k not in mapped_attrs}
        if left:out['attributes']=left
        remaining=[]; ci=0
        for i,c in enumerate(node.get('children',[])):
            if tag=='EXPERIMENT_PACKAGE' and c['tag'] in ('STUDY','SAMPLE','EXPERIMENT','SUBMISSION','RUN_SET'): continue
            if c['tag']=='Contact':self.contact_index=ci;ci+=1
            v=self.xml(c,kind,acc,provider,(*path,i),owner,actor,file,citation)
            if tag == 'PLATFORM' and c['tag'] != 'INDEXED' and not (
                    runs and all(r.get('instrument_platform') == c['tag'] for r in runs)):
                # The XML element name itself is a supplied family value.
                v = v or {'tag': c['tag']}
            if v:remaining.append(v)
        if remaining:out['children']=remaining
        if text and not mapped_text:out['text']=node['text']
        if node.get('tail','').strip():out['tail']=node['tail']
        if len(out)==1:
            return None if drop or mapped_text or mapped_attrs or node.get('children') else out
        for key in ('accession','filename','attribute_name'):
            if attrs.get(key): out.setdefault('attributes',{})[key]=attrs[key]
        return out

    def indexed(self, metadata, kind, acc):
        from .insdc_support import files_from_ena
        entity=self.entity(kind,acc); result=deepcopy(metadata)
        runs=self.experiments.get(acc,[]) if kind=='read_experiment' else [entity] if kind=='read_run' else []
        fields={'study_title':'title','study_description':'summary','sample_title':'title','sample_description':'description',
                'library_strategy':'library_strategy','library_source':'library_source','library_selection':'library_selection',
                'library_layout':'library_layout','instrument_model':'instrument_model',
                'instrument_platform':'instrument_platform'}
        for source,target in fields.items():
            if metadata.get(source) and (entity.get(target)==metadata[source] or (runs and all(r.get(target)==metadata[source] for r in runs))):result.pop(source,None)
        for field in ('study_accession','secondary_study_accession','sample_accession','secondary_sample_accession','experiment_accession','run_accession'):
            if metadata.get(field) and self.reference(metadata[field],kind,acc):result.pop(field,None)
        for field,key in [('scientific_name','value'),('tax_id','taxid')]:
            if any(str(o.get(key))==str(metadata.get(field)) for c in entity.get('channel',[]) for o in c.get('organism',[])):result.pop(field,None)
        for field,key in [('first_public','release_date'),('last_updated','last_update_date')]:
            if metadata.get(field) and any(s.get(key)==metadata[field] for s in entity.get('status',[])):result.pop(field,None)
        for field in ('read_count','base_count'):
            if entity.get('indexed_statistics',{}).get(field)==metadata.get(field):result.pop(field,None)
        files=files_from_ena(metadata)
        actual=[f for r in runs for f in r.get('files',[])]
        if files and all(any(contains(f,a) for a in actual) for f in files):
            for family in ('fastq','submitted','sra','bam'):
                for suffix in ('ftp','md5','bytes','file_role','format','aspera','galaxy'):result.pop(family+'_'+suffix,None)
        return {k:v for k,v in result.items() if v not in ('',None,[],{})}


def source_records(package):
    return list(deepcopy(getattr(package, '_archive_source_records', package.to_mapping().get('extensions',{}).get('insdc',{}).get('records',[]))))


def finalize(data, records=None):
    """Finalize after each mapping/merge; private source occurrences survive only in memory."""
    from .codec import MINiMLCodec
    from .archive_dates import normalize_archive_dates
    normalize_archive_dates(data)
    from .archive_administration import normalize_administration
    normalize_administration(data)
    records=deepcopy(records if records is not None else data.get('extensions',{}).get('insdc',{}).get('records',[]))
    from .archive_factors import restore_native_factors
    restore_native_factors(data, records)
    from .publication_identifiers import clean_publication_identifiers
    records.extend(r for r in clean_publication_identifiers(data) if r not in records)
    from .archive_results import normalize_result_bundles, comparison_view
    records.extend(r for r in normalize_result_bundles(data) if r not in records)
    from .archive_entities import local_platforms, coalesce_organizations, repository_databases
    local_platforms(data)
    coalesce_organizations(data)
    records.extend(r for r in repository_databases(data) if r not in records)
    from ..metadata.archive_enrichment import linked_accessions
    for acc in linked_accessions(data):
        if acc not in {a['value'] for a in data['series'].get('accession', [])} and not any(r.get('target')==acc for r in data['series'].get('relation', [])):
            data['series'].setdefault('relation', []).append({'type': 'GEO' if acc.startswith('GSE') else 'ArrayExpress', 'target': acc})
    projection=Projection(data, records);residual=[];seen=set()
    comparison = comparison_view(data)
    for record in records:
        # One destination occurrence can represent at most one occurrence within
        # this source record. Independent records keep their own matching scope.
        projection._mapped_characters.clear()
        kind,acc,metadata=record['kind'],record.get('accession'),record['metadata']
        if kind in ('MINiML','MINiML_workflows'):
            # Incoming identifiers have been remapped by the merger before this
            # snapshot is made. Saved packages are compared by those exact IDs.
            view = _mapped_date_view(metadata, comparison, workflows=kind=='MINiML_workflows')
            left=diff(view,comparison if kind=='MINiML' else comparison['series'])
            if kind == 'MINiML':
                from .archive_workflow_residuals import prune_bound_workflows
                left = prune_bound_workflows(view, data, left, record['provider'], acc)
        elif isinstance(metadata,dict) and 'tag' in metadata:
            if kind=='PubmedArticle' and not acc:
                acc=child_text(metadata,'MedlineCitation/PMID') or None
            if kind=='DocumentSummary' and not acc:
                project=next(iter(children(metadata,'Project')),{}); acc=next(iter(children(next(iter(children(project,'ProjectID')),{}),'ArchiveID')),{}).get('attributes',{}).get('accession')
            left=projection.xml(metadata,kind,acc,record['provider'])
            if left and set(left)=={'tag'}:left=None
        elif kind=='publication_reference':
            candidates=projection.citations(acc)
            candidate=next((p for p in candidates if any(metadata.get(k) and metadata[k]==p.get(k) for k in ('pubmed_id','doi','pmcid'))),{})
            left=diff(metadata,candidate)
        elif kind=='cross_references':
            from ..sources.archive_publications import citation_identifier
            left=[]
            for row in metadata:
                item=deepcopy(row)
                pubs=projection.citations(acc)
                mapped=False
                for key in ('Source Primary Accession','Source Secondary Accession','Source URL','Source Secondary URL','url'):
                    value=row.get(key)
                    if not value: item.pop(key,None); continue
                    wanted=citation_identifier('pmc' if key=='Source Primary Accession' else 'pubmed' if key=='Source Secondary Accession' else 'url',value)
                    represented=bool(wanted and any(contains(wanted,p) for p in pubs))
                    represented |= acc in projection.ids and (value in projection.ids or
                        (row.get('Source') in ('ArrayExpress','GEO') and row.get('Source Primary Accession') in projection.ids and key=='Source URL'))
                    if represented: item.pop(key,None); mapped=True
                if mapped:
                    for key in ('Source','Target','Target Primary Accession','Target Secondary Accession','Target URL'):
                        item.pop(key,None)
                item={k:v for k,v in item.items() if v not in ('',None,[],{})}
                if item and mapped:
                    for key in ('Source','Source Primary Accession','Source Secondary Accession','Target','Target Primary Accession','Target Secondary Accession'):
                        if row.get(key): item[key]=row[key]
                if item:left.append(item)
        elif kind in ('study','sample','read_run','read_experiment') and record['provider']!='biosamples':
            left=projection.indexed(metadata,kind,acc)
        elif record['provider']=='biosamples' and kind=='sample':
            left=deepcopy(metadata);sample=projection.samples.get(acc,{})
            attrs={}
            for name,values in metadata.get('characteristics',{}).items():
                from .insdc_support import biosample_characteristic
                rest=[]
                for value in values:
                    mapped = biosample_characteristic(name, value)
                    terms = {k:v for k,v in mapped.items() if k not in ('name','value','unit')}
                    # Alias recognition is limited to the import's explicit
                    # underscore/space name rule, within this resolved sample.
                    names = {a['name'] for c in sample.get('channel', []) for a in c.get('characteristics', [])
                             if a['name'].replace('_',' ').casefold() == name.replace('_',' ').casefold()}
                    resolved_name = next(iter(names)) if len(names) == 1 else name
                    if projection.character(sample,resolved_name,mapped['value'],value.get('unit'),terms or None):
                        remaining=deepcopy(value)
                        for key in ('text','unit'):
                            remaining.pop(key,None)
                        if 'term_accession_number' in terms:
                            remaining.pop('ontologyTerms',None)
                        if remaining:
                            # This literal identifies which occurrence owns the
                            # residual siblings (including multiple ontology terms).
                            if 'text' in value:
                                remaining['text'] = value['text']
                            rest.append(remaining)
                    else:
                        rest.append(deepcopy(value))
                if name=='organism':
                    remaining=[]
                    for value in rest:
                        candidates=[(ci,oi,o) for ci,c in enumerate(sample.get('channel',[]))
                                    for oi,o in enumerate(c.get('organism',[]))
                                    if o.get('value')==value.get('text') and
                                    (sample.get('iid'),'organism',ci,oi) not in projection._mapped_characters]
                        if candidates and not value.get('unit'):
                            ci,oi,organism=candidates[0]
                            projection._mapped_characters.add((sample.get('iid'),'organism',ci,oi))
                            matches=[organism]
                            item=deepcopy(value)
                            item.pop('text',None)
                            terms=value.get('ontologyTerms',[])
                            identifiers={str(o[k]) for o in matches for k in ('term_accession_number','taxid') if o.get(k)}
                            for o in matches:
                                if o.get('taxid'):
                                    identifiers.update(f'{prefix}{o["taxid"]}' for prefix in (
                                        'http://purl.obolibrary.org/obo/NCBITaxon_',
                                        'https://purl.obolibrary.org/obo/NCBITaxon_', 'NCBITaxon:'))
                            if len(terms)==1 and terms[0] in identifiers:
                                item.pop('ontologyTerms',None)
                            if item:
                                item['text']=value.get('text','')
                                remaining.append(item)
                        else:
                            remaining.append(value)
                    rest=remaining
                if rest:attrs[name]=rest
            if attrs:
                left['characteristics']=attrs
            else:
                left.pop('characteristics',None)
            if any(s.get('database') == 'BioSamples' and {'name':'status','value':metadata.get('status')} in s.get('comment', []) for s in sample.get('status', [])):
                left.pop('status', None)
            for key,field in [('name','title'),('description','description')]:
                if sample.get(field)==metadata.get(key):left.pop(key,None)
            if set(left) <= {'accession'} and left.get('accession',acc) == acc:
                left = {}
        elif kind in ('analysis','assembly') and record['provider']=='ena':
            from .insdc_support import files_from_ena
            left=deepcopy(metadata)
            for key in ('analysis_accession','assembly_set_accession','assembly_accession','accession'):
                if any(r.get('target')==metadata.get(key) for r in projection.series.get('relation',[])):left.pop(key,None)
            targets=[s for s in data.get('sample',[]) if any(r.get('target')==acc for r in s.get('relation',[]))]
            links=[l for s in targets for l in s.get('supplementary_data',[])]
            files=files_from_ena(metadata)
            if files and all(any(l.get('value')==f.get('uri') for l in links) for f in files):
                for family in ('fastq','submitted','sra','bam'):
                    for suffix in ('ftp','md5','bytes','file_role','format','aspera','galaxy'):
                        left.pop(family+'_'+suffix,None)
            if metadata.get('sample_accession') in projection.samples and targets: left.pop('sample_accession',None)
            left={k:v for k,v in left.items() if v not in ('',None,[],{})}
        elif kind=='accession_range':
            entity=projection.samples.get(acc,{})
            left=None if all({'type':metadata['database'],'target':a} in entity.get('relation',[]) for a in metadata['accessions']) else deepcopy(metadata)
        elif kind=='taxonomy':
            left=deepcopy(metadata)
            for o in projection.organisms():
                if str(o.get('taxid'))==str(metadata.get('taxId')):
                    left.pop('taxId',None)
                    if o.get('value')==metadata.get('scientificName'):left.pop('scientificName',None)
        elif kind=='assembly' and record['provider']=='sra':
            left=deepcopy(metadata)
            for key in ('assemblyaccession','lastmajorreleaseaccession'):
                if any(r.get('target')==metadata.get(key) for r in projection.series.get('relation',[])):left.pop(key,None)
            sample=projection.samples.get(metadata.get('biosampleaccn'),{})
            if any(r.get('target')==acc for r in sample.get('relation',[])):left.pop('biosampleaccn',None)
            for key in ('ftppath_genbank','ftppath_refseq'):
                if any(r.get('target')==metadata.get(key) for r in projection.series.get('relation',[])):left.pop(key,None)
            for key in ('ftppath_stats_rpt','ftppath_regions_rpt','ftppath_assembly_rpt'):
                if any(r.get('value')==metadata.get(key) for r in (sample or projection.series).get('supplementary_data',[])):left.pop(key,None)
            if 'synonym' in left:
                for k in ('genbank','refseq'):
                    if any(r.get('target')==left['synonym'].get(k) for r in projection.series.get('relation',[])):left['synonym'].pop(k,None)
        else:left=deepcopy(metadata)
        if not left:continue
        item={**record,'accession':acc,'metadata':left}
        if kind == 'term_source_declaration':
            compatible = next((r for r in residual if r['kind'] == kind and r['metadata'].get('iid') == left.get('iid')
                              and all(not r['metadata'].get(k) or not v or r['metadata'][k] == v for k,v in left.items())), None)
            if compatible is not None:
                compatible['metadata'].update({k:v for k,v in left.items() if v and not compatible['metadata'].get(k)})
                continue
        key=json.dumps(item,sort_keys=True)
        if key not in seen:seen.add(key);residual.append(item)
    data.setdefault('extensions',{})['insdc']={'version':'2.0','records':residual}
    package=MINiMLCodec().decode(data).package
    object.__setattr__(package,'_archive_source_records',tuple(records))
    return package
