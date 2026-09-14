# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Exact inverse projection of verified ArrayExpress workflow residuals."""
from copy import deepcopy
import re


def prune_bound_workflows(source, target, residual, provider, accession):
    """Prune mapped siblings without discarding source aliases or event positions."""
    if provider != 'MAGE-TAB' or not isinstance(residual, dict):
        return residual
    from meta_standards_converter.miniml.archive_residuals import contains, _path_signature, _path_diff
    from meta_standards_converter.metadata.archive_workflows import _prepare_path, path_ids, mark_file_origins, is_file
    old_paths = source.get('series', {}).get('assay_paths', [])
    left_paths = residual.get('series', {}).get('assay_paths', [])
    if not old_paths or not left_paths:
        return residual
    samples = {s['iid']: s for s in target.get('sample', [])}
    from meta_standards_converter.miniml.archive_results import comparison_view, normalize_result_bundles
    expanded = comparison_view(target)
    actual = expanded.get('series', {}).get('assay_paths', [])
    prepared = []
    for path in old_paths:
        if any(c.get('name', '').upper() in {'FASTQ_URI', 'FASTQ_FILE_NAME'} and set(c) - {'name', 'value'}
               for node in path['steps'] for c in node.get('comments', [])):
            prepared.append(None); continue
        refs = {s.get('sample_ref') for s in path['steps']} - {None}
        scope = path_ids(path)
        runs = {x for x in scope if re.fullmatch(r'[SED]RR\d+', x)}
        experiments = {x for x in scope if re.fullmatch(r'[SED]RX\d+', x)}
        if len(refs) != 1 or len(runs) != 1 or len(experiments) != 1:
            prepared.append(None); continue
        sid = next(iter(refs)); run = next(iter(runs)); experiment = next(iter(experiments))
        if sid not in samples or not any(r.get('run') == run and r.get('experiment') == experiment for r in samples[sid].get('sra_run', [])):
            prepared.append(None); continue
        if any(sum(s['kind'] == kind for s in path['steps']) != 1 for kind in ('source', 'assay', 'scan')) or any(s['kind'] == 'sample' for s in path['steps']):
            prepared.append(None); continue
        try:
            view = _prepare_path(path, sid, {})
        except ValueError:
            prepared.append(None); continue
        for step in view['steps']:
            if step['kind'] in ('source', 'assay', 'scan'):
                step.update(name={'source': sid, 'assay': experiment, 'scan': run}[step['kind']], sample_ref=sid)
        mark_file_origins({'source': {'format': 'MAGE-TAB'}, 'series': {'iid': accession, 'assay_paths': [view]}})
        files = [s for s in view['steps'] if is_file(s)]
        if not any(s['kind'] == 'array_data_file' and s.get('link', {}).get('value') for s in files):
            prepared.append(None); continue
        signature = _path_signature(view)
        prepared.append((view, signature))
    identities = [p[1] if p else None for p in prepared]
    replacements = []
    for index, entry in enumerate(prepared):
        if entry is None or identities.count(entry[1]) != 1:
            continue
        original = old_paths[index]
        # Only replace the whole path retained by the generic conservative diff.
        if sum(p == original for p in left_paths) != 1:
            continue
        view, identity = entry
        matches = []
        for actual_index, path in enumerate(actual):
            candidate = deepcopy(path)
            missing = [i for i, s in enumerate(view['steps']) if is_file(s) and not s.get('link', {}).get('value')]
            if missing:
                # Only a positively verified final MEX object authorizes the
                # filename-only representation inverse, within this exact path.
                probe = {'source': target['source'], 'series': {**target['series'], 'assay_paths': [deepcopy(path)]}}
                if not normalize_result_bundles(probe) or probe['series']['assay_paths'][0] != target['series']['assay_paths'][actual_index]:
                    continue
                if len(candidate['steps']) != len(view['steps']):
                    continue
                valid = True
                for position in missing:
                    a, b = view['steps'][position], candidate['steps'][position]
                    if not a['kind'].startswith('derived_') or b.get('link', {}).get('value') != a.get('name') or b.get('name') != a.get('name'):
                        valid = False; break
                    b['link'].pop('value')
                if not valid:
                    continue
            if (_path_signature(candidate) == identity
                    and all(contains(a, b) for a, b in zip(view['steps'], candidate['steps']) if is_file(a))):
                matches.append(path)
        if len(matches) != 1:
            continue
        inverse = {k: deepcopy(v) for k, v in matches[0].items() if k != 'steps'}
        inverse['steps'] = []
        cursor = 0
        for node in original['steps']:
            projected = view['steps'][cursor]
            restored = deepcopy(matches[0]['steps'][cursor])
            # _prepare_path removes only contract-defined, completely aligned
            # file groups. Their full projected nodes were proved above.
            remaining = deepcopy(projected.get('comments', []))
            moved = []
            for comment in node.get('comments', []):
                if comment in remaining:
                    remaining.remove(comment)
                else:
                    moved.append(comment)
            saved = restored.setdefault('comments', [])
            for comment in moved:
                needed = moved.count(comment)
                while saved.count(comment) < needed:
                    saved.append(deepcopy(comment))
            if node['kind'] == 'source':
                # The biological organism label is canonicalized by mapping;
                # its value and complete ontology/unit group must still agree.
                for character in node.get('characteristics', []):
                    if character.get('name') != 'Organism':
                        continue
                    wanted = {**character, 'name': 'organism'}
                    candidates = [c for c in restored.get('characteristics', []) if contains(wanted, c)]
                    # An exact occurrence remains identifiable beside an
                    # independently supplied, richer taxonomy assertion.
                    exact = [c for c in candidates if c == wanted]
                    if exact:
                        candidates = exact
                    if len(candidates) == 1 and sum(c == character for c in node['characteristics']) == 1:
                        candidates[0]['name'] = 'Organism'
            inverse['steps'].append(restored)
            cursor += 1 + sum(c.get('name', '').upper() == 'FASTQ_URI' and bool(c.get('value')) for c in node.get('comments', []))
        delta = _path_diff(original, inverse)
        if delta:
            # A scan-comment read URI identifies the original row/mate even
            # after mapped checksum/format payload has been pruned.
            for old, left in zip(original['steps'], delta['steps']):
                for c in old.get('comments', []):
                    if c.get('name', '').upper() == 'FASTQ_URI' and c.get('value') and c not in left.setdefault('comments', []):
                        left['comments'].append(deepcopy(c))
                if not left.get('comments'):left.pop('comments', None)
        replacements.append((original, delta))
    if not replacements:
        return residual
    result = deepcopy(residual)
    updated = []
    for path in left_paths:
        replacement = next((d for p, d in replacements if p == path), path)
        if replacement is not None:updated.append(replacement)
    result['series']['assay_paths'] = updated
    return result
