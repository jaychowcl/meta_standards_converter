# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Pure sample selection with whole-path closure and frozen selection evidence."""
from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from types import MappingProxyType
from typing import Any, Iterable, Mapping
import re

from .codec import MINiMLCodec
from .model import MINiMLModelError, MINiMLPackage, MINiMLValidationIssue
from .patches import miniml_source_fingerprint

SAMPLE_SELECTION_EXTENSION_KEY = 'msc_sample_selection'
SAMPLE_SELECTION_CONTRACT_VERSION = '1.0'


@dataclass(frozen=True)
class MINiMLSampleSelection:
    """Immutable package plus identities and selection-time scalar provenance."""
    package: MINiMLPackage
    analysis_ids: tuple[str, ...]
    context_ids: tuple[str, ...]
    excluded_ids: tuple[str, ...]
    diagnostics: tuple[MINiMLValidationIssue, ...]
    selected_to_source: Mapping[str, str]


def _mapping(package):
    return package.to_mapping() if isinstance(package, MINiMLPackage) else deepcopy(dict(package))


def selection_content_fingerprint(package: MINiMLPackage | Mapping[str, Any]) -> str:
    """Canonical content hash excluding selection evidence and retained patches.

    The selected hash describes the selection-time snapshot, not later annotation
    changes. Excluding this extension avoids a self-referential hash.
    """
    value = _mapping(package)
    extensions = value.get('extensions', {})
    extensions.pop(SAMPLE_SELECTION_EXTENSION_KEY, None)
    if not extensions:
        value.pop('extensions', None)
    return miniml_source_fingerprint(value)


def validate_sample_selection(package: Mapping[str, Any]) -> Mapping[str, Any] | None:
    """Check role coverage without treating a historical hash as a live checksum."""
    marker = package.get('extensions', {}).get(SAMPLE_SELECTION_EXTENSION_KEY)
    if marker is None:
        return None
    if not isinstance(marker, Mapping) or marker.get('contract_version') != SAMPLE_SELECTION_CONTRACT_VERSION:
        raise MINiMLModelError('unsupported sample selection contract')
    groups = []
    for name in ('analysis_ids', 'context_ids', 'excluded_ids'):
        ids = marker.get(name)
        if not isinstance(ids, (list, tuple)) or any(not isinstance(i, str) or not i.strip() for i in ids):
            raise MINiMLModelError(f'invalid sample selection {name}')
        groups.extend(ids)
    if len(groups) != len(set(groups)):
        raise MINiMLModelError('sample selection identities must be disjoint and unique')
    ids = [s.get('iid') for s in package.get('sample', [])]
    if len(ids) != len(set(ids)) or set(ids) != set(marker['analysis_ids']) | set(marker['context_ids']):
        raise MINiMLModelError('sample selection roles must cover retained sample identities')
    for name in ('source_sha256', 'selected_sha256'):
        if not isinstance(marker.get(name), str) or not re.fullmatch('[0-9a-f]{64}', marker[name]):
            raise MINiMLModelError(f'invalid sample selection {name}')
    pointers = marker.get('selected_to_source')
    if not isinstance(pointers, Mapping) or any(
        not isinstance(k, str) or not k.startswith('/') or not isinstance(v, str) or not v.startswith('/')
        for k, v in pointers.items()
    ):
        raise MINiMLModelError('invalid sample selection pointer map')
    return marker


def sample_role(package: MINiMLPackage | Mapping[str, Any], sample_id: str) -> str:
    """Return analysis/context; legacy unmarked packages retain analysis semantics."""
    value = _mapping(package) if isinstance(package, MINiMLPackage) else package
    marker = validate_sample_selection(value)
    if marker is None:
        return 'analysis'
    if sample_id in marker['analysis_ids']:
        return 'analysis'
    if sample_id in marker['context_ids']:
        return 'context'
    raise MINiMLModelError(f'sample {sample_id!r} is not retained by selection')


def select_miniml_samples(
    package: MINiMLPackage | Mapping[str, Any], analysis_ids: Iterable[str],
) -> MINiMLSampleSelection:
    """Retain complete paths touching analysis samples, closing over sample refs.

    Validate before pruning as well as after selection. Unbound paths and all
    declarations survive as study context. A marked selection can be reused with
    the same analysis set, but a different selection requires the frozen source;
    this prevents promotion of context or stale patch-pointer reinterpretation.
    """
    codec = MINiMLCodec()
    decoded = codec.decode(package, strict=True)
    source = codec.encode(decoded.package)
    wanted = set(analysis_ids)
    ids = [s.get('iid') for s in source.get('sample', [])]
    if any(not isinstance(i, str) or not i.strip() for i in ids):
        raise MINiMLModelError('sample selection requires explicit sample iid identities')
    if wanted - set(ids):
        raise MINiMLModelError(f'unknown analysis sample identities: {sorted(wanted - set(ids))}')
    marker = validate_sample_selection(source)
    if marker:
        if wanted != set(marker['analysis_ids']):
            raise MINiMLModelError('reselect from frozen source; context samples cannot be promoted')
        return MINiMLSampleSelection(decoded.package, tuple(marker['analysis_ids']), tuple(marker['context_ids']),
            tuple(marker['excluded_ids']), decoded.diagnostics, MappingProxyType(dict(marker['selected_to_source'])))

    paths = source['series'].get('assay_paths', [])
    refs = [{s['sample_ref'] for s in path['steps'] if s.get('sample_ref')} for path in paths]
    retained = set(wanted)
    while True:
        expanded = retained | {i for group in refs if group & retained for i in group}
        if expanded == retained:
            break
        retained = expanded
    selected = deepcopy(source)
    # Track old list indices while pruning. Other arrays retain occurrence order.
    indices: dict[str, list[int]] = {}
    def keep(items, predicate, pointer):
        positions = [i for i, item in enumerate(items) if predicate(item)]
        indices[pointer] = positions
        return [items[i] for i in positions]
    selected['sample'] = keep(selected.get('sample', []), lambda s: s['iid'] in retained, '/sample')
    series = selected['series']
    if 'assay_paths' in series:
        positions = [i for i, group in enumerate(refs) if not group or group & retained]
        indices['/series/assay_paths'] = positions
        series['assay_paths'] = [series['assay_paths'][i] for i in positions]
    if 'sample_ref' in series:
        series['sample_ref'] = keep(series['sample_ref'], lambda r: r['ref'] in retained, '/series/sample_ref')
    for key in ('variable', 'repeats'):
        if key not in series:
            continue
        subsets, positions = [], []
        for old_index, subset in enumerate(series[key]):
            members = subset.get('sample_ref', [])
            member_positions = [i for i, r in enumerate(members) if r['ref'] in retained]
            if members and not member_positions:
                continue
            new_index = len(subsets)
            positions.append(old_index)
            if 'sample_ref' in subset:
                subset['sample_ref'] = [members[i] for i in member_positions]
                indices[f'/series/{key}/{new_index}/sample_ref'] = member_positions
            subsets.append(subset)
        indices[f'/series/{key}'] = positions
        series[key] = subsets
    output = codec.encode(codec.decode(selected, strict=True).package)
    pointers = {}
    def walk(value, out='', original=''):
        if isinstance(value, dict):
            for key, item in value.items():
                segment = key.replace('~', '~0').replace('/', '~1')
                walk(item, f'{out}/{segment}', f'{original}/{segment}')
        elif isinstance(value, list):
            positions = indices.get(out, list(range(len(value))))
            for i, item in enumerate(value):
                walk(item, f'{out}/{i}', f'{original}/{positions[i]}')
        else:
            pointers[out] = original
    walk(output)
    analysis = tuple(i for i in ids if i in wanted)
    context = tuple(i for i in ids if i in retained - wanted)
    excluded = tuple(i for i in ids if i not in retained)
    output.setdefault('extensions', {})[SAMPLE_SELECTION_EXTENSION_KEY] = {
        'contract_version': SAMPLE_SELECTION_CONTRACT_VERSION,
        'analysis_ids': list(analysis), 'context_ids': list(context), 'excluded_ids': list(excluded),
        'source_sha256': selection_content_fingerprint(source),
        'selected_sha256': selection_content_fingerprint(output),
        'selected_to_source': pointers,
    }
    result = codec.decode(output, strict=True)
    return MINiMLSampleSelection(result.package, analysis, context, excluded,
        result.diagnostics, MappingProxyType(pointers))
