# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Conservative native result-object normalization; no retrieval or text rewriting."""
from copy import deepcopy
import re
from urllib.parse import urlsplit


def companion_node(value):
    """Extras are not schema-validated AssayNodes; interpret only known shapes."""
    return (isinstance(value, dict) and value.get('kind') in ('derived_array_data_file', 'derived_array_data_matrix_file')
            and isinstance(value.get('name'), str) and isinstance(value.get('link', {}), dict)
            and isinstance(value.get('comments', []), list)
            and all(isinstance(c, dict) and isinstance(c.get('name', ''), str) for c in value.get('comments', [])))


def _cell_ranger(definition):
    """Recognize affirmative supplied software/method forms, not a bare mention."""
    text = str(definition.get('description', ''))
    if re.search(r'\b(?:not|never|without|rather than|instead of)\b[^.!?;\n]{0,160}\bcell\s*ranger\b|'
                 r'\bcell\s*ranger\s+(?:was|is)\s+not\s+used\b', text, re.I):
        return False
    if any(re.fullmatch(r'cell\s*ranger(?:\s+v?\d[\w. -]*)?', str(s), re.I) for s in definition.get('software', [])):
        return True
    return bool(re.search(r'\b(?:processed|aligned|generated)\s+(?:with|using)\s+cell\s*ranger\b|'
                          r'\bcell\s*ranger(?:\s+v?[\d.]+)?\s+(?:generates|generated|was used)\b', text, re.I))


def normalize_result_bundles(data):
    """Recognize a source-declared Cell Ranger MEX object within one acquisition.

    Require one bound sample/run, one three-member result tail, common literal
    member identity/directory, and identical full application blocks. General
    repeated protocols or unrelated same-name files are never consolidated.
    """
    if data.get('source', {}).get('format') not in ('ENA', 'SRA'):
        return []
    definitions = {p['name']: p for p in data.get('series', {}).get('protocols', [])}
    layouts = []
    for path in data.get('series', {}).get('assay_paths', []):
        steps = path['steps']
        indices = [i for i, s in enumerate(steps) if s['kind'].startswith('derived_')]
        bound = {s['sample_ref'] for s in steps if s.get('sample_ref')}
        runs = [s['name'] for s in steps if s['kind'] == 'scan']
        if len(indices) != 3 or len(bound) != 1 or len(runs) != 1 or indices[-1] != len(steps)-1:
            continue
        files = [steps[i] for i in indices]
        if any(s.get('link', {}).get('companion_files') for s in files):
            continue
        matches = [re.fullmatch(r'(.*?)(barcodes\.tsv(?:\.gz)?|features\.tsv(?:\.gz)?|matrix\.mtx(?:\.gz)?)', s['name']) for s in files]
        if not all(matches) or [m[2].split('.')[0] for m in matches] != ['barcodes', 'features', 'matrix']:
            continue
        if len({m[1] for m in matches}) != 1:
            continue
        locations = [s.get('link', {}).get('value', '') for s in files]
        if any(not v for v in locations) or len({v.rsplit('/', 1)[0] if '/' in v else '' for v in locations}) != 1:
            continue
        if any(urlsplit(v).path.rsplit('/', 1)[-1] != s['name'] for v, s in zip(locations, files)):
            continue
        start = indices[0]
        while start and steps[start-1]['kind'] == 'protocol_application':
            start -= 1
        blocks = [steps[start:indices[0]], steps[indices[0]+1:indices[1]], steps[indices[1]+1:indices[2]]]
        if not blocks[0] or any(block != blocks[0] for block in blocks[1:]):
            continue
        if any(s['kind'] != 'protocol_application' for s in blocks[0]):
            continue
        if not any(_cell_ranger(definitions.get(s['protocol_ref'], {})) for s in blocks[0]):
            continue
        layout = {'provider': files[-1].get('link', {}).get('repository', data['source']['format']), 'kind': 'result_layout',
                  'accession': files[-1].get('link', {}).get('source_accession', data['series']['iid']),
                  'metadata': {'sample_ref': next(iter(bound)), 'run': runs[0],
                               'sequence': [{k: s[k] for k in ('kind', 'name', 'protocol_ref') if k in s} for s in steps[start:]]}}
        if layout not in layouts:
            layouts.append(layout)
        matrix = deepcopy(files[-1])
        matrix['link']['companion_files'] = [{'role': role, 'node': deepcopy(node)}
                                            for role, node in zip(('barcodes', 'features'), files[:2])]
        steps[start:] = deepcopy(blocks[0]) + [matrix]
    return layouts


def comparison_view(data):
    """Expose represented companion leaves to residual matching, not to export.

    Original serial topology is retained separately as result_layout. Re-expanding
    the known source layout here avoids retaining entire mapped source workflows.
    """
    if not any(s.get('link', {}).get('companion_files') for p in data.get('series', {}).get('assay_paths', []) for s in p['steps']):
        return data
    result = deepcopy(data)
    for path in result['series']['assay_paths']:
        steps = path['steps']
        for i in range(len(steps)-1, -1, -1):
            members = steps[i].get('link', {}).get('companion_files')
            if (not isinstance(members, list) or not all(isinstance(m, dict) and companion_node(m.get('node')) for m in members)
                    or [m.get('role') for m in members] != ['barcodes', 'features']):
                continue
            start = i
            while start and steps[start-1]['kind'] == 'protocol_application': start -= 1
            block = deepcopy(steps[start:i])
            matrix = deepcopy(steps[i]); matrix['link'].pop('companion_files')
            steps[start:i+1] = [s for node in [*(m['node'] for m in members), matrix] for s in [*deepcopy(block), deepcopy(node)]]
    return result
