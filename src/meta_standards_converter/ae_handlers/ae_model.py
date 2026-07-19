# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Typed, editable MAGE-TAB extension for MINiML-compatible JSON packages."""

from __future__ import annotations

import copy
import re


PROTOCOL_FIELDS = {
    "Protocol Name": "name",
    "Protocol Type": "type",
    "Protocol Type Term Source REF": "type_term_source_ref",
    "Protocol Type Term Accession Number": "type_term_accession_number",
    "Protocol Description": "description",
    "Protocol Hardware": "hardware",
    "Protocol Software": "software",
    "Protocol Parameters": "parameters",
    "Protocol Contact": "contact",
    "Protocol Performer": "performer",
}

PROTOCOL_FIELD_ALIASES = {
    "type_term_source_ref": (
        "Protocol Type Term Source REF",
        "Protocol Term Source REF",
    ),
    "type_term_accession_number": (
        "Protocol Type Term Accession Number",
        "Protocol Term Accession Number",
    ),
}

DECLARATION_FIELDS = {
    "quality_control": (
        "Quality Control Type",
        "Quality Control Term Source REF",
        "Quality Control Term Accession Number",
    ),
    "replicate": (
        "Replicate Type",
        "Replicate Term Source REF",
        "Replicate Term Accession Number",
    ),
    "normalization": (
        "Normalization Type",
        "Normalization Term Source REF",
        "Normalization Term Accession Number",
    ),
}

NODE_HEADERS = {
    "Source Name", "Sample Name", "Extract Name", "Labeled Extract Name",
    "Hybridization Name", "Assay Name", "Scan Name", "Normalization Name",
}

FILE_HEADERS = {
    "Array Data File", "Array Data Matrix File", "Derived Array Data File",
    "Derived Array Data Matrix File", "Image File",
}


def build_model(idf_rows: list[list], sdrfs: list[tuple[str, list[list]]]) -> dict:
    """Build a typed model without projecting unsupported values into MINiML fields."""
    values = {_normalized(row[0]): list(row[1:]) for row in idf_rows if row}
    names = values.get(_normalized("Protocol Name"), [])
    inferred_protocol_count = max(
        [
            len(names),
            *[
                len(_protocol_values(values, label, key)[0])
                for label, key in PROTOCOL_FIELDS.items()
            ],
        ],
        default=0,
    )
    protocol_count = len(names) if names else inferred_protocol_count
    protocol_labels = {}
    protocol_widths = {}
    protocols = []
    for position in range(protocol_count):
        record = {"id": f"protocol:{position + 1}", "position": position}
        for label, key in PROTOCOL_FIELDS.items():
            row_values, source_label = _protocol_values(values, label, key)
            if source_label:
                protocol_labels[key] = source_label
                protocol_widths[key] = len(row_values)
            record[key] = row_values[position] if position < len(row_values) else ""
        if not record.get("name"):
            record["name"] = record["id"]
        protocols.append(record)

    declarations = {
        kind: _declarations(values, labels)
        for kind, labels in DECLARATION_FIELDS.items()
    }
    declaration_widths = {
        kind: {
            field: len(values.get(_normalized(label), []))
            for label, field in zip(
                labels,
                ("value", "term_source_ref", "term_accession_number"),
            )
        }
        for kind, labels in DECLARATION_FIELDS.items()
    }
    typed_labels = {
        *(_normalized(label) for label in PROTOCOL_FIELDS),
        *(
            _normalized(label)
            for labels in PROTOCOL_FIELD_ALIASES.values()
            for label in labels
        ),
        *(
            _normalized(label)
            for labels in DECLARATION_FIELDS.values()
            for label in labels
        ),
    }
    investigation_fields = [
        {
            "row_index": index,
            "label": row[0],
            "values": copy.deepcopy(row[1:]),
        }
        for index, row in enumerate(idf_rows)
        if row and _normalized(row[0]) not in typed_labels
    ]

    model_sdrfs = []
    assay_paths = []
    for sdrf_name, table in sdrfs:
        if not table:
            continue
        header = table[0]
        columns = _columns(header)
        model_sdrfs.append({"name": sdrf_name, "columns": columns})
        for row_index, row in enumerate(table[1:], start=1):
            assay_paths.append(_assay_path(sdrf_name, row_index, header, row))

    return {
        "schema_version": 1,
        "idf_layout": [
            {"row_index": index, "label": row[0]}
            for index, row in enumerate(idf_rows)
            if row
        ],
        "protocols": protocols,
        "protocol_field_labels": protocol_labels,
        "protocol_field_widths": protocol_widths,
        "declarations": declarations,
        "declaration_field_widths": declaration_widths,
        "assay_paths": assay_paths,
        "sdrfs": model_sdrfs,
        "investigation_fields": investigation_fields,
    }


def render_model(model: dict) -> list | None:
    """Render a version-1 typed model into the constructor's in-memory MAGE-TAB form."""
    if not isinstance(model, dict) or model.get("schema_version") != 1:
        return None
    sdrfs = model.get("sdrfs") or []
    if not sdrfs:
        return None

    typed_rows = _typed_idf_rows(model)
    investigation = {
        item.get("row_index"): [item.get("label"), *(item.get("values") or [])]
        for item in model.get("investigation_fields", [])
        if isinstance(item, dict) and item.get("label") is not None
    }
    rows = []
    for layout in sorted(model.get("idf_layout", []), key=lambda item: item.get("row_index", 0)):
        label = layout.get("label")
        row = typed_rows.get(_normalized(label)) or investigation.get(layout.get("row_index"))
        if row is not None:
            rows.append(copy.deepcopy(row))

    known = {_normalized(row[0]) for row in rows if row}
    for row in [*typed_rows.values(), *investigation.values()]:
        if row and _normalized(row[0]) not in known:
            rows.append(copy.deepcopy(row))
            known.add(_normalized(row[0]))

    rendered_sdrfs = [_render_sdrf(model, descriptor) for descriptor in sdrfs]
    sdrf = _consolidate_sdrfs(rendered_sdrfs)
    for index, row in enumerate(rows):
        if row and _normalized(row[0]) == _normalized("SDRF File"):
            rows[index] = ["SDRF File", sdrf, *row[2:]]
            break
    else:
        rows.append(["SDRF File", sdrf])
    return rows


def overlay_core(model_rows: list, core_rows: list) -> list:
    """Apply exact MINiML projections while retaining model-only graph information."""
    result = copy.deepcopy(model_rows)
    replace_labels = {
        _normalized(label)
        for label in (
            "Investigation Title", "Investigation Accession",
            "Investigation Accession Term Source REF", "Comment[SecondaryAccession]",
            "Comment[SecondaryAccessionTermSourceRef]", "Experimental Design",
            "Experimental Design Term Source REF", "Experimental Design Term Accession Number",
            "Experimental Factor Name", "Experimental Factor Type",
            "Experimental Factor Term Source REF", "Experimental Factor Term Accession Number",
            "Person Last Name", "Person First Name", "Person Mid Initials", "Person Email",
            "Person Phone", "Person Fax", "Person Address", "Person Affiliation",
            "Date of Experiment", "Public Release Date", "Comment[GEOReleaseDate]",
            "Comment[GEOLastUpdateDate]", "PubMed ID", "Publication DOI",
            "Publication Author List", "Publication Title", "Publication Status",
            "Status Term Source Ref", "Status Term Accession Number",
            "Experiment Description", "Term Source Name", "Term Source File", "Term Source Version",
        )
    }
    core_by_label = {_normalized(row[0]): row for row in core_rows if row}
    _overlay_protocol_rows(result, core_rows)
    for index, row in enumerate(result):
        label = _normalized(row[0]) if row else ""
        if label in replace_labels and label in core_by_label:
            result[index] = copy.deepcopy(core_by_label[label])

    model_sdrf_index = _sdrf_index(result)
    core_sdrf_index = _sdrf_index(core_rows)
    if model_sdrf_index is not None and core_sdrf_index is not None:
        model_sdrf = result[model_sdrf_index][1]
        core_sdrf = core_rows[core_sdrf_index][1]
        if _table(model_sdrf) and _table(core_sdrf):
            _overlay_sdrf(model_sdrf, core_sdrf)
    return result


def _declarations(values: dict, labels: tuple[str, str, str]) -> list[dict]:
    terms = values.get(_normalized(labels[0]), [])
    sources = values.get(_normalized(labels[1]), [])
    accessions = values.get(_normalized(labels[2]), [])
    present = any(_normalized(label) in values for label in labels)
    count = max(len(terms), len(sources), len(accessions), 1 if present else 0)
    return [
        {
            "position": index,
            "value": terms[index] if index < len(terms) else "",
            "term_source_ref": sources[index] if index < len(sources) else "",
            "term_accession_number": accessions[index] if index < len(accessions) else "",
        }
        for index in range(count)
    ]


def _columns(header: list[str]) -> list[dict]:
    occurrences = {}
    result = []
    for index, label in enumerate(header):
        key = _normalized(label)
        occurrences[key] = occurrences.get(key, 0) + 1
        result.append({"index": index, "header": label, "occurrence": occurrences[key]})
    return result


def _assay_path(sdrf_name: str, row_index: int, header: list[str], row: list[str]) -> dict:
    binding = {
        _snake(label): row[index]
        for label in ("Source Name", "Sample Name", "Assay Name", "Comment[ENA_RUN]")
        for index, candidate in enumerate(header)
        if _normalized(candidate) == _normalized(label) and index < len(row) and row[index]
    }
    steps = []
    occurrences = {}
    consumed = set()
    for index, label in enumerate(header):
        if index in consumed:
            continue
        value = row[index] if index < len(row) else ""
        key = _normalized(label)
        occurrences[key] = occurrences.get(key, 0) + 1
        base = {
            "column_index": index,
            "header": label,
            "occurrence": occurrences[key],
            "value": value,
        }
        annotation = re.fullmatch(r"\s*(Characteristics|Factor\s+Value|Parameter\s+Value)\s*\[(.*)]\s*", label, re.I)
        if annotation:
            base.update({
                "kind": "attribute",
                "attribute_type": " ".join(annotation.group(1).split()).lower(),
                "name": annotation.group(2).strip(),
            })
            companions = {}
            for companion_index in range(index + 1, len(header)):
                companion_label = header[companion_index]
                normalized = _normalized(companion_label)
                if normalized not in {
                    _normalized("Unit"), _normalized("Term Source REF"),
                    _normalized("Term Accession Number"),
                }:
                    break
                companion_value = row[companion_index] if companion_index < len(row) else ""
                field = {
                    _normalized("Unit"): "unit",
                    _normalized("Term Source REF"): "term_source_ref",
                    _normalized("Term Accession Number"): "term_accession_number",
                }[normalized]
                base[field] = companion_value
                companions[field] = companion_index
                consumed.add(companion_index)
            if companions:
                base["companion_columns"] = companions
        elif label in NODE_HEADERS:
            base.update({"kind": "node", "node_type": label})
        elif _normalized(label) == _normalized("Protocol REF"):
            base.update({"kind": "protocol_ref", "name": value})
        elif re.fullmatch(r"\s*Comment\[(.*)]\s*", label, re.I):
            base.update({"kind": "comment", "name": re.fullmatch(r"\s*Comment\[(.*)]\s*", label, re.I).group(1)})
        elif label in FILE_HEADERS:
            base.update({"kind": "file", "file_type": label})
        else:
            base.update({"kind": "field"})
        steps.append(base)
    return {
        "id": f"{sdrf_name}:row:{row_index}",
        "sdrf": sdrf_name,
        "row_index": row_index,
        "binding": binding,
        "steps": steps,
    }


def _typed_idf_rows(model: dict) -> dict[str, list]:
    protocols = sorted(model.get("protocols", []), key=lambda item: item.get("position", 0))
    layout_labels = {
        _normalized(item.get("label"))
        for item in model.get("idf_layout", [])
        if isinstance(item, dict)
    }
    rows = {}
    if protocols:
        protocol_labels = model.get("protocol_field_labels") or {}
        protocol_widths = model.get("protocol_field_widths") or {}
        for canonical_label, key in PROTOCOL_FIELDS.items():
            label = protocol_labels.get(key) or canonical_label
            values = [item.get(key, "") for item in protocols]
            if _normalized(label) in layout_labels or any(value not in (None, "") for value in values):
                width = _edited_width(values, protocol_widths.get(key, 0))
                rows[_normalized(label)] = [label, *values[:width]]
    for kind, labels in DECLARATION_FIELDS.items():
        records = sorted(
            model.get("declarations", {}).get(kind, []),
            key=lambda item: item.get("position", 0),
        )
        fields = ("value", "term_source_ref", "term_accession_number")
        declaration_widths = model.get("declaration_field_widths", {}).get(kind, {})
        for label, field in zip(labels, fields):
            values = [item.get(field, "") for item in records]
            if records and (
                _normalized(label) in layout_labels
                or any(value not in (None, "") for value in values)
            ):
                width = _edited_width(values, declaration_widths.get(field, 0))
                rows[_normalized(label)] = [label, *values[:width]]
    return rows


def _render_sdrf(model: dict, descriptor: dict) -> list[list]:
    columns = sorted(descriptor.get("columns", []), key=lambda item: item.get("index", 0))
    width = max((item.get("index", 0) for item in columns), default=-1) + 1
    header = [""] * width
    for column in columns:
        header[column["index"]] = column.get("header", "")
    paths = sorted(
        (path for path in model.get("assay_paths", []) if path.get("sdrf") == descriptor.get("name")),
        key=lambda item: item.get("row_index", 0),
    )
    rows = [header]
    for path in paths:
        row = [""] * width
        for step in path.get("steps", []):
            index = step.get("column_index")
            if isinstance(index, int) and index < width:
                row[index] = step.get("value", "")
            for field, companion_index in (step.get("companion_columns") or {}).items():
                if isinstance(companion_index, int) and companion_index < width:
                    row[companion_index] = step.get(field, "")
        rows.append(row)
    return rows


def _consolidate_sdrfs(tables: list[list[list]]) -> list[list]:
    if len(tables) == 1:
        return tables[0]
    union = []
    seen = set()
    table_keys = []
    for table in tables:
        occurrences = {}
        keys = []
        for header in table[0]:
            label = _normalized(header)
            occurrences[label] = occurrences.get(label, 0) + 1
            key = (label, occurrences[label])
            keys.append(key)
            if key not in seen:
                seen.add(key)
                union.append((key, header))
        table_keys.append(keys)
    union_indexes = {key: index for index, (key, _) in enumerate(union)}
    result = [[header for _, header in union]]
    for table, keys in zip(tables, table_keys):
        for source_row in table[1:]:
            row = [""] * len(union)
            for index, value in enumerate(source_row):
                if index < len(keys):
                    row[union_indexes[keys[index]]] = value
            result.append(row)
    return result


def _overlay_protocol_rows(model_rows: list, core_rows: list) -> None:
    labels = {
        key: (_protocol_row(model_rows, label, key), _protocol_row(core_rows, label, key))
        for label, key in PROTOCOL_FIELDS.items()
    }
    model_types, core_types = labels["type"]
    if not model_types or not core_types:
        return
    used = set()
    supported_fragments = (
        "manufactur", "treatment", "growth", "extract", "label", "hybrid",
        "scan", "processing", "normalization",
    )
    for core_position, core_type in enumerate(core_types[1:], start=1):
        normalized_type = _normalized(core_type)
        if not any(fragment in normalized_type for fragment in supported_fragments):
            continue
        model_position = next(
            (
                index
                for index, model_type in enumerate(model_types[1:], start=1)
                if index not in used and _normalized(model_type) == normalized_type
            ),
            None,
        )
        if model_position is None:
            continue
        used.add(model_position)
        for key in (
            "type", "type_term_source_ref", "type_term_accession_number", "description",
        ):
            model_row, core_row = labels[key]
            if model_row and core_row and core_position < len(core_row):
                while len(model_row) <= model_position:
                    model_row.append("")
                model_row[model_position] = core_row[core_position]


def _overlay_sdrf(model: list[list], core: list[list]) -> None:
    model_header, core_header = model[0], core[0]
    core_indexes = {_normalized(label): index for index, label in enumerate(core_header)}
    protected = {_normalized(label) for label in NODE_HEADERS | {"Protocol REF"}}
    core_rows_by_identity = {}
    for row in core[1:]:
        for identity in _identities(core_header, row):
            core_rows_by_identity.setdefault(identity, row)
    for row in model[1:]:
        core_row = next(
            (core_rows_by_identity[value] for value in _identities(model_header, row) if value in core_rows_by_identity),
            None,
        )
        if core_row is None:
            continue
        for model_index, label in enumerate(model_header):
            key = _normalized(label)
            core_index = core_indexes.get(key)
            if key in protected or core_index is None or core_index >= len(core_row):
                continue
            row[model_index] = core_row[core_index]


def _identities(header: list, row: list) -> list[str]:
    values = []
    for label in ("Sample Name", "Source Name", "Comment[ENA_RUN]"):
        for index, candidate in enumerate(header):
            if _normalized(candidate) == _normalized(label) and index < len(row) and row[index]:
                values.append(row[index])
    return values


def _sdrf_index(rows: list) -> int | None:
    return next(
        (index for index, row in enumerate(rows) if row and _normalized(row[0]) == _normalized("SDRF File")),
        None,
    )


def _row(rows: list, label: str) -> list | None:
    return next(
        (row for row in rows if row and _normalized(row[0]) == _normalized(label)),
        None,
    )


def _protocol_row(rows: list, canonical_label: str, key: str) -> list | None:
    for label in PROTOCOL_FIELD_ALIASES.get(key, (canonical_label,)):
        value = _row(rows, label)
        if value is not None:
            return value
    return None


def _protocol_values(values: dict, canonical_label: str, key: str) -> tuple[list, str | None]:
    for label in PROTOCOL_FIELD_ALIASES.get(key, (canonical_label,)):
        normalized = _normalized(label)
        if normalized in values:
            return values[normalized], label
    return [], None


def _edited_width(values: list, original_width: int) -> int:
    last_nonblank = max(
        (index + 1 for index, value in enumerate(values) if value not in (None, "")),
        default=0,
    )
    return max(original_width, last_nonblank)


def _table(value) -> bool:
    return isinstance(value, list) and bool(value) and all(isinstance(row, list) for row in value)


def _normalized(value) -> str:
    return "".join(str(value).split()).casefold()


def _snake(value) -> str:
    return re.sub(r"[^a-z0-9]+", "_", value.casefold()).strip("_")
