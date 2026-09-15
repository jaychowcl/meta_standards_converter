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
from .protocol_text import comparable_protocol_text


def _literal(value):
    return value.get('value') if isinstance(value, dict) else value


def _comment(node, name, value):
    value = _literal(value)
    if value not in (None, ''):
        node['comments'] = [c for c in node.get('comments', []) if c.get('name', '').casefold() != name.casefold()]
        node['comments'].append({'name': name, 'value': str(value)})


def _authored_method(steps, by_name, field, description, base):
    families = {
        'growth_protocol': {'grow', 'growth', 'growth protocol'},
        'treatment_protocol': {'treatment', 'treatment protocol'},
        'extract_protocol': {'nucleic acid extraction', 'nucleic acid extraction protocol', 'extraction'},
        'label_protocol': {'labeling', 'labelling', 'labeling protocol'},
        'scan_protocol': {'scanning', 'sequencing', 'nucleic acid sequencing protocol'},
        'hybridization_protocol': {'hybridization', 'hybridization protocol'},
        'data_processing': {'normalization', 'data processing', 'normalization data transformation protocol'},
    }
    for step in steps:
        method = by_name.get(step.get('protocol_ref'), {})
        name = method.get('name', '')
        kind = str(_literal(method.get('type', ''))).lower().replace('_', ' ')
        if (name and name != base and not name.startswith(base + ':')
                and kind in families[field]
                and ' '.join(str(method.get('description', '')).split()) == ' '.join(description.split())):
            return True
    return False


def _compound_application(steps, sample, methods, description):
    """Recognize duplicate scalar projection at an unambiguous native boundary."""
    assays = [s for s in steps if s['kind'] == 'assay']
    scans = [s for s in steps if s['kind'] == 'scan']
    if not description or len(assays) != 1 or len(scans) != 1:
        return None
    experiment = assays[0].get('name')
    matches = [r for r in sample.get('sra_run', [])
               if r.get('experiment') == experiment and r.get('run') == scans[0].get('name')]
    if len(matches) != 1 or not experiment:
        return None
    ref = experiment + ':library'
    applications = [s for s in steps if s.get('protocol_ref') == ref]
    method = methods.get(ref, {})
    kind = str(_literal(method.get('type', ''))).lower().replace('_', ' ')
    if (len(applications) != 1 or set(applications[0]) != {'kind', 'protocol_ref'}
            or comparable_protocol_text(method.get('description') or '') != comparable_protocol_text(description)
            or kind not in {'library construction protocol', 'nucleic acid library construction protocol'}):
        return None
    application = applications[0]
    # Never move an application across an authored procedure. Insert the
    # generated material only at the native library-to-assay boundary.
    if steps.index(application) + 1 != steps.index(assays[0]):
        return None
    return application


def complete_native_paths(data):
    """Project explicit channel materials/protocols and scoped descriptive fields."""
    if data.get('source', {}).get('format') not in ('ENA', 'SRA'):
        return
    from ..metadata.archive_enrichment import informative
    samples = {s['iid']: s for s in data.get('sample', [])}
    platforms = {p['iid']: p for p in data.get('platform', [])}
    series = data['series']
    definitions = series.setdefault('protocols', [])
    methods = {p['name']: p for p in definitions}
    paths_by_sample = {}
    for path in series.get('assay_paths', []):
        bound = {s['sample_ref'] for s in path['steps'] if s.get('sample_ref')}
        if len(bound) == 1:
            paths_by_sample.setdefault(next(iter(bound)), []).append(path)
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
            boundary = {'scan_protocol': 'scan', 'hybridization_protocol': 'hybridization',
                        'data_processing': 'derived_array_data_file'}.get(field, 'assay')
            relevant = [p for p in paths_by_sample.get(sid, [])
                        if any(s['kind'].startswith('derived_') if field == 'data_processing'
                               else s['kind'] == boundary for s in p['steps'])]
            if relevant and all(_authored_method(p['steps'], methods, field, description, base) for p in relevant):
                obsolete = {p['name'] for p in definitions
                            if (p['name'] == base or p['name'].startswith(base + ':'))
                            and set(p) <= {'name', 'type', 'description'}
                            and ' '.join(str(p.get('description', '')).split()) == ' '.join(description.split())}
                definitions[:] = [p for p in definitions if p['name'] not in obsolete]
                for name in obsolete: methods.pop(name, None)
                for path in series.get('assay_paths', []):
                    path['steps'][:] = [s for s in path['steps'] if s.get('protocol_ref') not in obsolete]
                continue
            old = next((p for p in definitions if p['name'].startswith(base) and p.get('description') == description), None)
            if old is None:
                name = base
                count = 1
                while any(p['name'] == name for p in definitions):
                    count += 1
                    name = f'{base}:{count}'
                old = {'name': name, 'description': description, 'type': {'value': kind}}
                definitions.append(old)
                methods[old['name']] = old
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
                for accession in sample.get('accession', []):
                    if accession.get('database') == 'BioSample' and accession.get('value'):
                        comment = {'name': 'BioSD_SAMPLE', 'value': accession['value']}
                        if comment not in step.setdefault('comments', []):
                            step['comments'].append(comment)
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
            compound = _compound_application(steps, sample, methods, channel.get('extract_protocol'))
            position = next((i for i, s in enumerate(steps) if s['kind'] not in ('source', 'sample')), len(steps))
            addition = []
            after_compound = []
            for field in fields:
                destination = after_compound if compound is not None and field in ('extract_protocol', 'label_protocol') else addition
                if (sid, field) in refs and not (compound is not None and field == 'extract_protocol'):
                    destination.append({'kind': 'protocol_application', 'protocol_ref': refs[(sid, field)]})
                if field == 'extract_protocol' and informative(channel.get('molecule')):
                    destination.append({'kind': 'extract', 'name': sid + ':extract', 'sample_ref': sid,
                                     'material_type': deepcopy(channel['molecule'])})
            steps[position:position] = addition
            if compound is not None:
                position = steps.index(compound) + 1
                steps[position:position] = after_compound
        for field, boundary in (('scan_protocol', 'scan'), ('hybridization_protocol', 'hybridization')):
            if (sid, field) not in refs or not any(s['kind'] == boundary for s in steps):
                continue
            generated = {p['name'] for p in definitions if ':' + sid + ':' + field in p['name']}
            steps[:] = [s for s in steps if s.get('protocol_ref') not in generated]
            base = f"{series['iid']}:{sid}:{field}"
            if _authored_method(steps, methods, field, sample[field], base):
                continue
            position = next(i for i, s in enumerate(steps) if s['kind'] == boundary)
            steps.insert(position, {'kind': 'protocol_application', 'protocol_ref': refs[(sid, field)]})
        if has_result and (sid, 'data_processing') in refs:
            by_name = {p['name']: p for p in definitions}
            processing = [s for s in steps if s.get('protocol_ref') and any(word in str(by_name.get(s['protocol_ref'], {}).get('type', {})).lower()
                          for word in ('processing', 'transformation', 'alignment', 'normalization'))]
            if not processing:
                position = next(i for i, s in enumerate(steps) if s['kind'].startswith('derived_'))
                steps.insert(position, {'kind': 'protocol_application', 'protocol_ref': refs[(sid, 'data_processing')]})
