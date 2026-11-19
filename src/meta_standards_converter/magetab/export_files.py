# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Exclude processed files from private MAGE-TAB plans, retaining row facts."""
import json

from .sdrf.model import SDRFNode
DERIVED_LABELS = frozenset({'Derived Array Data File', 'Derived Array Data Matrix File'})
INDEPENDENT = ('Factor Value[', 'Characteristics[', 'Parameter Value[')


def _source_file_annotations(step):
    """Retain typed row facts from native source-only file comment bundles."""
    from ..miniml.model import NamedValue
    kept = []
    for comment in step.get('comments', []):
        name = comment.get('name', '')
        if not name.startswith('SAMPLE_FILE_'):
            kept.append(comment)
        elif name == 'SAMPLE_FILE_NODE':
            try:
                node = json.loads(comment['value'])
                fields = {field: [NamedValue.from_mapping(v).to_mapping() for v in node.get(field, [])]
                          for field in ('characteristics', 'factor_values')}
            except (ValueError, TypeError, AttributeError, KeyError):
                # Uninterpretable file-owned metadata remains in the input JSON.
                continue
            for field, values in fields.items():
                if values:
                    step.setdefault(field, []).extend(values)
    if 'comments' in step:
        step['comments'] = kept


def prepare_export_files(data):
    """Prevent generated supplementary-file comments from bypassing projection.

    Only the constructor's detached export view is passed here. Canonical file
    nodes are filtered after native file consolidation, preserving its existing
    acquisition occurrences, and before handler routing and column planning.
    Prose and provenance are not searched or rewritten.
    """
    from .sdrf.handlers.base import classify_file
    series = data.get('series', {})
    from ..metadata.file_resolution import structural_uri
    samples = data.get('sample', [])
    sample_objects = {id(sample) for sample in samples}
    for entity in [series, *samples, *data.get('platform', [])]:
        if 'supplementary_data' in entity:
            entity['supplementary_data'] = [link for link in entity['supplementary_data']
                if link.get('type') != 'assembly report'
                and classify_file(link) != 'matrix_or_derived']
        if id(entity) not in sample_objects:
            for link in entity.get('supplementary_data', []):
                if classify_file(link) == 'supplementary':
                    comment = {'name': 'Unclassified File' + (f" [{entity['iid']}]" if entity is not series and entity.get('iid') else ''),
                               'value': structural_uri(link.get('value'))}
                    if comment not in series.setdefault('comments', []):
                        series['comments'].append(comment)
    # Native projection may represent source-only results as file comment
    # bundles. Strip those alongside nodes, without changing row occurrences.
    for path in series.get('assay_paths', []):
        for step in path['steps']:
            _source_file_annotations(step)
        removed = [step for step in path['steps'] if step['kind'].startswith('derived_array_data')]
        if not removed:
            continue
        retained = [step for step in path['steps'] if not step['kind'].startswith('derived_array_data')]
        nodes = [step for step in retained if step['kind'] != 'protocol_application']
        target = next((step for step in reversed(nodes) if step['kind'] in {'assay', 'hybridization'}),
                      nodes[0] if nodes else None)
        for step in removed:
            for field in ('characteristics', 'factor_values'):
                if step.get(field):
                    if target is None:
                        target = {'kind': 'source', 'name': ''}
                        retained.insert(0, target)
                    target.setdefault(field, []).extend(step[field])
        path['steps'] = retained


def _independent(attrs):
    for attr in attrs:
        if attr.label.startswith(INDEPENDENT):
            yield attr
        else:
            yield from _independent(attr.attrs)


def exclude_derived_paths(paths):
    """Filter typed nodes and attributes before planning columns, without dedup.

    Scientific annotations carry their complete qualifier trees to a retained
    assay (or source) node. File-specific comments/checksums leave with the file.
    """
    for path in paths:
        retained, annotations = [], []
        for part in path.parts:
            if isinstance(part, SDRFNode) and part.kind in DERIVED_LABELS:
                annotations.extend(_independent(part.attrs))
                continue
            def filtered(attrs):
                result = []
                for attr in attrs:
                    if attr.label in DERIVED_LABELS:
                        annotations.extend(_independent(attr.attrs))
                    else:
                        attr.attrs = filtered(attr.attrs)
                        result.append(attr)
                return result
            part.attrs = filtered(part.attrs)
            retained.append(part)
        if annotations:
            nodes = [part for part in retained if isinstance(part, SDRFNode)]
            target = next((part for part in reversed(nodes) if part.kind in {'Assay Name', 'Hybridization Name'}),
                          nodes[0] if nodes else None)
            if target is None:
                # A blank carrier preserves annotations without inventing identity.
                target = SDRFNode('Source Name', 'source:unidentified', '')
                retained.insert(0, target)
            target.attrs.extend(annotations)
        path.parts = retained
    return paths
