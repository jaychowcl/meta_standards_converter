# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Complete structural protocol coverage on the detached MAGE-TAB export copy.

Authored methods remain bound to their supplied workflows. A missing application
uses a wholly blank definition, never another sample's described procedure.
"""
from .protocols import SEQUENCING_TECHNOLOGIES
from ..metadata.ontology_mappings import Harmonizer


_COLLECTION = 'sample collection protocol'
_SEQUENCING = 'nucleic acid sequencing protocol'
_LABELS = {_COLLECTION: 'Sample-Collection-Protocol',
           _SEQUENCING: 'Nucleic-Acid-Sequencing-Protocol'}
_ACCESSIONS = {'EFO_0005518': _COLLECTION, 'EFO_0004170': _SEQUENCING}


def _type(protocol):
    value = protocol.get('type') or {}
    if not isinstance(value, dict):
        value = {'value': value}
    accession = str(value.get('term_accession_number') or '').rsplit('/', 1)[-1].replace(':', '_').upper()
    label = ' '.join(str(value.get('value') or '').replace('-', ' ').replace('_', ' ').split()).casefold()
    return _ACCESSIONS.get(accession, label)


def _blank(protocol):
    return all(not (value.strip() if isinstance(value, str) else value)
               for key, value in protocol.items() if key not in {'name', 'type'})


class _Definitions:
    def __init__(self, handler):
        self.registry = handler.protocol_registry
        self.protocols = handler.data['series'].setdefault('protocols', [])
        self.used = {p['name'] for p in self.protocols} | {p['ref'] for p in self.registry.records()}
        self.by_name = {p['name']: p for p in self.protocols}
        self.blanks = {}
        for protocol in self.protocols:
            if _blank(protocol):
                self.blanks.setdefault(_type(protocol), protocol['name'])
        self.next_index = 1

    def add(self, protocol):
        while f'P-{self.registry.series_accession}-{self.next_index}' in self.used:
            self.next_index += 1
        name = f'P-{self.registry.series_accession}-{self.next_index}'
        protocol = dict(protocol, name=name)
        self.used.add(name)
        self.by_name[name] = protocol
        self.protocols.append(protocol)
        return name

    def placeholder(self, kind):
        if kind not in self.blanks:
            label, source, accession = Harmonizer().geoprotocols2efo(_LABELS[kind])
            self.blanks[kind] = self.add({'type': {'value': label, 'term_source_ref': source,
                                                  'term_accession_number': accession}, 'description': ''})
        return self.blanks[kind]


def _sequencing(handler, path, by_sample, decisions):
    """Reuse routing decisions, narrowing by the path's sample/channel/run joins."""
    if not by_sample:  # Explicit platform handler (or empty input).
        return getattr(handler, 'technology_type', None) in SEQUENCING_TECHNOLOGIES
    if len(decisions) == 1:
        return decisions == {True}
    steps = path.get('steps', [])
    sample_ids = {str(s.get('sample_ref') or s.get('name')) for s in steps
                  if s.get('kind') in {'source', 'sample'}}
    matches = [item for sid in sample_ids for item in by_sample.get(sid, ())]
    # Unbound paths can use a unanimous operation decision, never one neighbour.
    candidates = matches
    from .harmonized import _binding_from_samples
    binding = _binding_from_samples(handler.sample_by_id, path) if matches else None
    if binding:
        candidates = [(s, d) for s, d in candidates if s.channel is binding[1]]
    run_ids = {str(s.get('name')) for s in steps if s.get('kind') == 'scan'}
    experiment_ids = {str(s.get('name')) for s in steps if s.get('kind') == 'assay'}
    for step in steps:
        for comment in step.get('comments', []):
            if comment.get('name') in {'ENA_RUN', 'SRA_RUN'}:
                run_ids.add(str(comment.get('value')))
            if comment.get('name') in {'ENA_EXPERIMENT', 'SRA_EXPERIMENT'}:
                experiment_ids.add(str(comment.get('value')))
    for identities, fields in ((run_ids, ('accession', 'run')), (experiment_ids, ('experiment',))):
        narrowed = [(s, d) for s, d in candidates if s.run and
                    any(str(s.run.get(key)) in identities for key in fields)]
        if narrowed:
            candidates = narrowed
    decisions = {d.handler in SEQUENCING_TECHNOLOGIES for _, d in candidates}
    if len(decisions) != 1:
        _warn(handler, 'nucleic acid sequencing protocol: ambiguous workflow technology; Protocol REF not added.')
    return decisions == {True}


def _boundary(steps, kind):
    if kind == _COLLECTION:
        for node in ('source', 'sample'):
            index = next((i for i, step in enumerate(steps) if step.get('kind') == node), None)
            if index is not None:
                return index + 1
    else:
        index = next((i for i, step in enumerate(steps)
                      if step.get('kind') in {'scan', 'array_data_file', 'array_data_matrix_file'}), None)
        if index is not None:
            return index
        return next((i + 1 for i in range(len(steps) - 1, -1, -1)
                     if steps[i].get('kind') in {'assay', 'extract'}), None)
    return None


def _warn(handler, message):
    if message not in handler.audit.warnings:
        handler.audit.warnings.append(message)


def complete_canonical_protocols(handler, paths):
    """Add only absent required applications; retain authored step/attribute order."""
    definitions = _Definitions(handler)
    by_sample, decisions = {}, set()
    for scope, decision in getattr(handler, 'technology_scopes', ()):
        by_sample.setdefault(str(scope.sample.get('iid')), []).append((scope, decision))
        decisions.add(decision.handler in SEQUENCING_TECHNOLOGIES)
    for path in paths:
        steps = path['steps']
        applied = {_type(definitions.by_name[step['protocol_ref']]) for step in steps
                   if step.get('protocol_ref') in definitions.by_name}
        required = [_COLLECTION, _SEQUENCING] if _sequencing(handler, path, by_sample, decisions) else [_COLLECTION]
        for kind in required:
            if kind in applied:
                continue
            reference = definitions.placeholder(kind)
            position = _boundary(steps, kind)
            if position is None:
                identity = path.get('id') or next((s.get('name') for s in steps if s.get('name')), 'unnamed')
                _warn(handler, f'{kind}: no usable workflow boundary for {identity}; definition retained, Protocol REF not added.')
                continue
            steps.insert(position, {'kind': 'protocol_application', 'protocol_ref': reference})


def reconcile_generated_protocols(handler, paths):
    """Keep generated references declared when supplied IDF definitions take priority.

    Without explicit paths, handler protocols describe the generated workflow.
    Match an authored definition only by exact supplied type/text; structural
    blanks never borrow authored details. Retain unused supplied definitions too.
    """
    from .sdrf.model import SDRFEdge
    definitions = _Definitions(handler)
    bindings = {}
    records = {r['ref']: r for r in handler.protocol_registry.records()}
    for path in paths:
        for edge in path.parts:
            if not isinstance(edge, SDRFEdge) or edge.protocol_ref not in records:
                continue
            old = edge.protocol_ref
            if old not in bindings:
                record = records[old]
                label, source, accession = Harmonizer().geoprotocols2efo(record['label'])
                if label in _LABELS and not record['text']:
                    reference = definitions.placeholder(label)
                else:
                    matches = [p['name'] for p in definitions.protocols if _type(p) == label
                               and p.get('description') == record['text']]
                    reference = matches[0] if len(matches) == 1 else definitions.add({
                        'type': {k: v for k, v in {'value': label, 'term_source_ref': source,
                                                   'term_accession_number': accession}.items() if v},
                        'description': record['text']})
                bindings[old] = reference
            edge.protocol_ref = bindings[old]
