# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Deterministic, network-free preparation of the canonical investigation graph.

Provider adapters still supply their technology-specific construction rules. This
service converts their in-memory fields to the existing typed model, before JSON
serialization or MAGE-TAB projection. It neither renders a table nor retrieves
metadata. Supplied paths retain their order, applications and annotation owners.
"""
from ..miniml import MINiMLCodec
from ..miniml.model import Series
from dataclasses import replace


def _prepared(package, data):
    """Rebuild only changed graph records; reuse validated immutable entities."""
    result = replace(package, series=Series.from_mapping(data['series']),
                     extensions=data.get('extensions', {}))
    result._raise_structural_errors()
    return result


def prepare_workflows(package, *, platform_handler=None):
    package = MINiMLCodec().decode(package).package
    data = package.to_mapping()
    from ..miniml.archive_paths import complete_native_paths
    complete_native_paths(data)
    if data.get('series', {}).get('assay_paths') or not data.get('sample'):
        return _prepared(package, data)

    from ..magetab.sdrf.constructor import SDRFConstructor
    from ..magetab.sdrf.model import SDRFEdge
    from ..magetab.protocols import ProtocolRegistry
    from ..magetab.technology import series_identity
    from ..magetab.semantics import _assay_path
    from ..miniml.migration import MINiMLV1Migrator
    from .ontology_mappings import Harmonizer
    constructor = SDRFConstructor()
    registry = ProtocolRegistry(series_identity(data))
    if platform_handler is None:
        handler, _ = constructor.create_operation_handler(data, registry)
    else:
        handler = constructor.create_handler(data, registry, platform_handler)
    paths = handler.prepare_paths()
    if handler.audit.warnings:
        data.setdefault('extensions', {})['msc_workflow_preparation'] = {
            'warnings': list(dict.fromkeys(handler.audit.warnings))}
    if data['series'].get('protocols'):
        from ..magetab.required_protocols import reconcile_generated_protocols
        reconcile_generated_protocols(handler, paths)
    else:
        definitions = []
        for record in registry.records():
            label, source, accession = Harmonizer().geoprotocols2efo(record['label'])
            definitions.append({'name': record['ref'], 'description': record['text'],
                'type': {k: v for k, v in {'value': label, 'term_source_ref': source,
                    'term_accession_number': accession}.items() if v is not None}})
        data['series']['protocols'] = definitions
    identities = {handler.material_identity(sample, channel): sample['iid']
                  for sample in handler.samples for channel in handler._channel_records(sample)}
    canonical = []
    for index, path in enumerate(paths):
        headers, values, required = [], [], []
        def append(label, value, attrs=(), mandatory=False):
            headers.append(label)
            values.append('' if value is None else value)
            required.append(mandatory)
            for attr in attrs:
                append(attr.label, attr.value, attr.attrs, attr.required)
        for part in path.parts:
            append('Protocol REF' if isinstance(part, SDRFEdge) else part.kind,
                   part.protocol_ref if isinstance(part, SDRFEdge) else part.value, part.attrs)
        # Reuse the semantic field decoder, including repeated qualifier owners;
        # no column planning, rendering, TSV serialization or re-import occurs.
        unbound = []
        fields = _assay_path('generated', index, headers, values)['steps']
        for field in fields:
            if required[field['column_index']]:
                field['_required'] = True
        steps = MINiMLV1Migrator._assay_steps(fields, unbound=unbound)
        if unbound:
            raise ValueError('Generated workflow contains unbound annotation fields')
        for step in steps:
            if step['kind'] in {'source', 'sample'} and step.get('name') in identities:
                step['sample_ref'] = identities[step['name']]
        if steps:
            canonical.append({'steps': steps, 'comments': [
                {'name': 'msc_workflow_origin', 'value': 'generated:' + str(data.get('source', {}).get('format', 'unknown'))}]})
    if canonical:
        data['series']['assay_paths'] = canonical
    return _prepared(package, data)
