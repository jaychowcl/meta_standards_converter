# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Occurrence-local harmonized groups shared by SDRF construction and parsing."""
from dataclasses import replace
from collections.abc import Mapping
import re

from meta_standards_converter.miniml.harmonization import HarmonizedValue, iter_harmonized_values, parse_harmonized_key


def columns(value, prefix="Characteristics", used=None):
    """Render complete value groups; indexes are scoped to their output occurrence."""
    used = set() if used is None else used
    result = []
    for item in _ordered_values(value):
        while (item.field, item.index) in used:
            item = replace(item, index=item.index + 1)
        used.add((item.field, item.index))
        label = f"hz_{item.field}{item.suffix}"
        result.append((f"{prefix}[{label}]", item.value))
        if item.term_source_ref is not None:
            result.append(("Term Source REF", item.term_source_ref))
        if item.term_accession_number is not None:
            result.append(("Term Accession Number", item.term_accession_number))
        if item.hierarchy_depth is not None:
            result.append((f"Comment[hz_{item.field}_hierarchy_depth{item.suffix}]", item.hierarchy_depth))
    return result


def _ordered_values(container):
    values = list(iter_harmonized_values(container))
    keys = list(container) if isinstance(container, Mapping) else [row.get("name", "") for row in container]
    order = {}
    for key in keys:
        if str(key).startswith("hz_"):
            field, role, index = parse_harmonized_key(key)
            if role == "value":
                order.setdefault((field, index), len(order))
    return sorted(values, key=lambda item: order.get((item.field, item.index), len(order)))


def _path_binding(package, path):
    """Identify a biological channel using explicit sample, source and label evidence."""
    samples = {str(s.get("iid")): s for s in package.get("sample", [])}
    steps = path.get("steps", [])
    bindings = [(node, samples.get(str(node.get("sample_ref") or node.get("name"))))
                for node in steps if node.get("kind") in {"sample", "source"}]
    bindings = [(node, sample) for node, sample in bindings if sample is not None]
    if not bindings:
        return None
    identities = {str(sample.get("iid")) for _, sample in bindings}
    if len(identities) != 1:
        import logging
        logging.getLogger(__name__).warning("Ambiguous sample association in explicit assay path; local evidence retained")
        return None
    # Source is the biological root when both source and sample nodes are present.
    node, sample = bindings[0]
    channels = sample.get("channel", [])
    markers = {str(c.get("value")) for step in steps for c in step.get("comments", [])
               if c.get("name") == "msc_channel"}
    if len(channels) > 1 and markers:
        channels = [c for i, c in enumerate(channels)
                    if str(c.get("extensions", {}).get("msc_channel", i)) in markers]
    if len(channels) > 1:
        def text(value):
            return str(value.get("value", "")) if isinstance(value, Mapping) else str(value or "")
        labels = {text(s.get("label")) for s in steps if s.get("label")}
        sources = {str(s.get("name")) for s in steps if s.get("kind") == "source"}
        candidates = [c for c in channels if (labels and text(c.get("label")) in labels)
                      or (sources and text(c.get("source")) in sources)]
        channels = candidates
    if len(channels) != 1:
        import logging
        logging.getLogger(__name__).warning("Ambiguous channel association for sample %s; explicit assay-local evidence retained", sample.get("iid"))
        return None
    return node, channels[0]


def channel_groups(channel):
    """Yield supported containers without recursively losing biological scope."""
    yield channel
    for key in ("source", "molecule", "material_type"):
        value = channel.get(key)
        if isinstance(value, Mapping):
            yield value
    for organism in channel.get("organism", []) or []:
        if isinstance(organism, Mapping):
            yield organism
    rows = channel.get("characteristics", []) or []
    yield rows
    for row in rows:
        if isinstance(row, Mapping):
            yield row


def channel_columns(channel):
    used, result = set(), []
    for container in channel_groups(channel):
        result.extend(columns(container, used=used))
    return result


def sdrf_attrs(pairs):
    from .sdrf.model import SDRFAttr
    result = []
    for label, value in pairs:
        if result and (label in {"Term Source REF", "Term Accession Number"} or "_hierarchy_depth" in label):
            result[-1].attrs.append(SDRFAttr(label, value))
        else:
            result.append(SDRFAttr(label, value))
    return result


def read_group(header, row, index):
    """Read a harmonized value plus its adjacent standard ontology companions."""
    match = re.fullmatch(r"(Characteristics|Factor Value|Parameter Value|Comment)\[(hz_[^]]+)\]", header[index], re.I)
    if not match:
        return None
    field, role, ordinal = parse_harmonized_key(match[2])
    if role != "value":
        return None
    end = index + 1
    source = accession = depth = None
    while end < len(header):
        label = header[end]
        value = row[end] if end < len(row) else ""
        if label == "Term Source REF":
            source = value or None
        elif label == "Term Accession Number":
            accession = value or None
        elif label == f"Comment[hz_{field}_hierarchy_depth{'('+str(ordinal)+')' if ordinal else ''}]":
            depth = int(value) if value != "" else None
        else:
            break
        end += 1
    value = row[index] if index < len(row) else ""
    return match[1].lower(), (None if value == "" else HarmonizedValue(field, value, source, accession, depth, ordinal)), end


def bind_sample_groups(package, paths):
    """Attach sample evidence to explicit paths only through unambiguous identities."""
    import copy
    import logging
    from meta_standards_converter.miniml.harmonization import named_harmonized_rows
    result = copy.deepcopy(paths or [])
    for path in result:
        binding = _path_binding(package, path)
        if binding is None:
            continue
        node, channel = binding
        rows = node.setdefault("characteristics", [])
        existing = list(iter_harmonized_values(rows))
        existing.extend(v for row in rows for v in iter_harmonized_values(row))
        used = {(v.field, v.index) for v in existing}
        for container in channel_groups(channel):
            for value in iter_harmonized_values(container):
                if any(replace(v, index=0) == replace(value, index=0) for v in existing):
                    continue
                while (value.field, value.index) in used:
                    value = replace(value, index=value.index+1)
                used.add((value.field, value.index))
                existing.append(value)
                rows.extend(named_harmonized_rows([value]))
    return result
