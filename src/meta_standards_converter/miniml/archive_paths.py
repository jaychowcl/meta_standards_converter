# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Complete native paths from explicit sample fields without inferring methods."""
from copy import deepcopy


def _literal(value):
    return value.get('value') if isinstance(value, dict) else value


def _comment(node, name, value):
    value = _literal(value)
    if value not in (None, ''):
        node['comments'] = [c for c in node.get('comments', []) if c.get('name', '').casefold() != name.casefold()]
        node['comments'].append({'name': name, 'value': str(value)})


def complete_native_paths(data):
    """Project explicit channel materials/protocols and scoped descriptive fields."""
    if data.get('source', {}).get('format') not in ('ENA', 'SRA'):
        return
    from ..metadata.archive_enrichment import informative
    samples = {s['iid']: s for s in data.get('sample', [])}
    platforms = {p['iid']: p for p in data.get('platform', [])}
    series = data['series']
    definitions = series.setdefault('protocols', [])
    types = {'growth_protocol': 'growth protocol', 'treatment_protocol': 'treatment protocol',
             'extract_protocol': 'nucleic acid extraction protocol', 'label_protocol': 'labeling protocol',
             'data_processing': 'normalization data transformation protocol',
             'scan_protocol': 'nucleic acid sequencing protocol', 'hybridization_protocol': 'hybridization protocol'}
    refs = {}
    for sid, sample in samples.items():
        channel = sample.get('channel', [None])[0] if len(sample.get('channel', [])) == 1 else None
        if channel is None:
            continue
        for field, kind in types.items():
            description = sample.get(field) if field in ('data_processing', 'scan_protocol', 'hybridization_protocol') else channel.get(field)
            if not informative(description):
                continue
            base = f"{series['iid']}:{sid}:{field}"
            old = next((p for p in definitions if p['name'].startswith(base) and p.get('description') == description), None)
            if old is None:
                name = base
                count = 1
                while any(p['name'] == name for p in definitions):
                    count += 1
                    name = f'{base}:{count}'
                old = {'name': name, 'description': description, 'type': {'value': kind}}
                definitions.append(old)
            refs[(sid, field)] = old['name']
    for path in series.get('assay_paths', []):
        steps = path['steps']
        bound = {s['sample_ref'] for s in steps if s.get('sample_ref')}
        if len(bound) != 1:
            continue
        sid = next(iter(bound)); sample = samples.get(sid)
        if sample is None or len(sample.get('channel', [])) != 1:
            continue
        channel = sample['channel'][0]
        for step in steps:
            if step['kind'] in ('source', 'sample'):
                _comment(step, 'Sample_title', sample.get('title'))
                _comment(step, 'Sample_source_name', channel.get('source'))
                _comment(step, 'Sample_description', sample.get('description'))
            elif step['kind'] == 'assay':
                platform = platforms.get(sample.get('platform_ref', {}).get('ref'), {})
                _comment(step, 'Platform_title', platform.get('title'))
        # Legacy generated processing applications on raw branches are not an
        # evidenced raw-to-result relationship. Definitions themselves survive.
        has_result = any(s['kind'].startswith('derived_') for s in steps)
        generated_processing = {p['name'] for p in definitions
                                if ':' + sid + ':data_processing' in p['name']}
        steps[:] = [s for s in steps if s.get('protocol_ref') not in generated_processing]
        # Refresh this helper's deterministic generated material identity.
        # Authored AE extracts retain their supplied (namespaced) identities.
        steps[:] = [s for s in steps if not (s['kind'] == 'extract' and s.get('name') == sid + ':extract')]
        # Complete explicit AE material workflows take precedence. GEO supplies
        # a channel molecule/extraction field, which is one intact method.
        if any(s['kind'] == 'assay' for s in steps) and not any(s['kind'] in ('extract', 'labeled_extract') for s in steps):
            fields = ('growth_protocol', 'treatment_protocol', 'extract_protocol', 'label_protocol')
            generated = {p['name'] for p in definitions if any(':' + sid + ':' + field in p['name'] for field in fields)}
            steps[:] = [s for s in steps if s.get('protocol_ref') not in generated]
            position = next((i for i, s in enumerate(steps) if s['kind'] not in ('source', 'sample')), len(steps))
            addition = []
            for field in fields:
                if (sid, field) in refs:
                    addition.append({'kind': 'protocol_application', 'protocol_ref': refs[(sid, field)]})
                if field == 'extract_protocol' and informative(channel.get('molecule')):
                    addition.append({'kind': 'extract', 'name': sid + ':extract', 'sample_ref': sid,
                                     'material_type': deepcopy(channel['molecule'])})
            steps[position:position] = addition
        for field, boundary in (('scan_protocol', 'scan'), ('hybridization_protocol', 'hybridization')):
            if (sid, field) not in refs or not any(s['kind'] == boundary for s in steps):
                continue
            generated = {p['name'] for p in definitions if ':' + sid + ':' + field in p['name']}
            steps[:] = [s for s in steps if s.get('protocol_ref') not in generated]
            position = next(i for i, s in enumerate(steps) if s['kind'] == boundary)
            steps.insert(position, {'kind': 'protocol_application', 'protocol_ref': refs[(sid, field)]})
        if has_result and (sid, 'data_processing') in refs:
            by_name = {p['name']: p for p in definitions}
            processing = [s for s in steps if s.get('protocol_ref') and any(word in str(by_name.get(s['protocol_ref'], {}).get('type', {})).lower()
                          for word in ('processing', 'transformation', 'alignment', 'normalization'))]
            if not processing:
                position = next(i for i, s in enumerate(steps) if s['kind'].startswith('derived_'))
                steps.insert(position, {'kind': 'protocol_application', 'protocol_ref': refs[(sid, 'data_processing')]})
