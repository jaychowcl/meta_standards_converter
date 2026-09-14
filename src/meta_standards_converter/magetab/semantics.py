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
from meta_standards_converter.metadata.provenance import (
    patch_provenance_columns,
)


class MAGETabModelError(ValueError):
    """Raised when the enriched MAGE-TAB model contract is invalid."""


MODEL_COLLECTIONS = (
    "idf_layout", "protocols", "assay_paths", "sdrfs", "investigation_fields",
)
HZ_FIELDS = (
    "hz_value", "hz_value_id", "hz_value_onto", "hz_field",
    "hz_unit", "hz_unit_id", "hz_unit_onto",
    "hz_value_hierarchy_depth", "hz_unit_hierarchy_depth",
)


def validate_model(model: dict) -> dict:
    """Validate and return an enriched MAGE-TAB model version 1."""
    if not isinstance(model, dict):
        raise MAGETabModelError("mage_tab.model must be an object")
    if model.get("schema_version") != 1:
        raise MAGETabModelError("mage_tab.model schema_version must be 1")
    for field in MODEL_COLLECTIONS:
        if not isinstance(model.get(field), list):
            raise MAGETabModelError(f"mage_tab.model {field} must be a list")
    if not isinstance(model.get("declarations"), dict):
        raise MAGETabModelError("mage_tab.model declarations must be an object")
    sdrf_names = {
        item.get("name") for item in model["sdrfs"]
        if isinstance(item, dict) and isinstance(item.get("name"), str)
    }
    if len(sdrf_names) != len(model["sdrfs"]):
        raise MAGETabModelError("mage_tab.model sdrfs require unique names")
    assay_ids: set[str] = set()
    for assay in model["assay_paths"]:
        if not isinstance(assay, dict) or not isinstance(assay.get("steps"), list):
            raise MAGETabModelError("mage_tab.model assay_paths require step lists")
        assay_id = assay.get("id")
        if not isinstance(assay_id, str) or not assay_id or assay_id in assay_ids:
            raise MAGETabModelError("mage_tab.model assay_paths require unique ids")
        assay_ids.add(assay_id)
        if assay.get("sdrf") not in sdrf_names:
            raise MAGETabModelError("mage_tab.model assay path references an unknown SDRF")
        for step in assay["steps"]:
            if not isinstance(step, dict) or not isinstance(step.get("column_index"), int):
                raise MAGETabModelError("mage_tab.model assay steps require column_index")
            if step.get("kind") == "attribute":
                if step.get("attribute_type") not in {
                    "characteristics", "factor value", "parameter value"
                } or not isinstance(step.get("name"), str):
                    raise MAGETabModelError("mage_tab.model contains an invalid attribute step")
                for field in HZ_FIELDS:
                    if field in step and not isinstance(step[field], (str, int, float)):
                        raise MAGETabModelError(
                            f"mage_tab.model attribute {field} must be scalar"
                        )
    return model


PROTOCOL_FIELDS = {
    "Protocol Name": "name",
    "Protocol Type": "type",
    "Protocol Term Source REF": "type_term_source_ref",
    "Protocol Term Accession Number": "type_term_accession_number",
    "Protocol Description": "description",
    "Protocol Hardware": "hardware",
    "Protocol Software": "software",
    "Protocol Parameters": "parameters",
    "Protocol Contact": "contact",
    "Protocol Performer": "performer",
}

PROTOCOL_FIELD_ALIASES = {
    "type_term_source_ref": (
        "Protocol Term Source REF",
        "Protocol Type Term Source REF",
    ),
    "type_term_accession_number": (
        "Protocol Term Accession Number",
        "Protocol Type Term Accession Number",
    ),
}

IDF_LABEL_ALIASES = {
    "Publication Status Term Source REF": (
        "Publication Status Term Source REF",
        "Status Term Source Ref",
    ),
    "Publication Status Term Accession Number": (
        "Publication Status Term Accession Number",
        "Status Term Accession Number",
    ),
    "Protocol Term Source REF": PROTOCOL_FIELD_ALIASES["type_term_source_ref"],
    "Protocol Term Accession Number": PROTOCOL_FIELD_ALIASES["type_term_accession_number"],
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
    protocol_count = inferred_protocol_count
    protocol_labels = {}
    protocol_widths = {}
    protocols = []
    for position in range(protocol_count):
        record = {"id": f"protocol:{position + 1}", "position": position}
        for label, key in PROTOCOL_FIELDS.items():
            row_values, source_label = _protocol_values(values, label, key)
            if source_label:
                protocol_labels[key] = label
                protocol_widths[key] = len(row_values)
            record[key] = row_values[position] if position < len(row_values) else ""
        if not any(str(record.get(key) or "").strip() for key in PROTOCOL_FIELDS.values()):
            continue
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
            "label": _canonical_idf_label(row[0]),
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
            {"row_index": index, "label": _canonical_idf_label(row[0])}
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
    if not isinstance(model, dict):
        return None
    validate_model(model)
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


def overlay_miniml_semantics(package: dict, core_rows: list) -> list:
    """Render ordered MSC MINiML semantics over constructor-generated IDF/SDRF rows."""
    rows = copy.deepcopy(core_rows)
    series = package.get("series") if isinstance(package, dict) else None
    if not isinstance(series, dict):
        return rows
    _replace_row(rows, "Investigation Accession", [series.get("iid")])
    protocols = [item for item in series.get("protocols", []) if isinstance(item, dict)]
    if protocols or package.get("source", {}).get("format") in {"SRA", "ENA"}:
        fields = (
            ("Protocol Name", lambda item: item.get("name")),
            ("Protocol Type", lambda item: _ontology_text(item.get("type"))),
            ("Protocol Term Source REF", lambda item: _ontology_field(item.get("type"), "term_source_ref")),
            ("Protocol Term Accession Number", lambda item: _ontology_field(item.get("type"), "term_accession_number")),
            ("Protocol Description", lambda item: item.get("description")),
            ("Protocol Hardware", lambda item: " | ".join(str(value) for value in item.get("hardware", []))),
            ("Protocol Software", lambda item: " | ".join(str(value) for value in item.get("software", []))),
            ("Protocol Parameters", lambda item: "; ".join(str(value) for value in item.get("parameters", []))),
            ("Protocol Contact", lambda item: "; ".join(str(value) for value in item.get("contacts", []))),
            ("Protocol Performer", lambda item: "; ".join(str(value) for value in item.get("performers", []))),
        )
        labels = {label for label, _ in fields}
        position = next((i for i, row in enumerate(rows) if row and row[0] in labels), len(rows))
        rows[:] = [row for row in rows if not row or row[0] not in labels]
        rows[position:position] = [[label, *[accessor(item) or "" for item in protocols]]
                                   for label, accessor in fields]
    for field, labels in (
        ("quality_controls", DECLARATION_FIELDS["quality_control"]),
        ("replicate_types", DECLARATION_FIELDS["replicate"]),
        ("normalization_types", DECLARATION_FIELDS["normalization"]),
    ):
        values = [item for item in series.get(field, []) if isinstance(item, dict)]
        if values:
            _replace_row(rows, labels[0], [_ontology_text(item) for item in values])
            _replace_row(rows, labels[1], [_ontology_field(item, "term_source_ref") for item in values])
            _replace_row(rows, labels[2], [_ontology_field(item, "term_accession_number") for item in values])
    for comment in series.get("comments", []) or []:
        if isinstance(comment, dict) and comment.get("name"):
            _replace_row(rows, f"Comment[{comment['name']}]", [comment.get("value", "")])
    if package.get('source', {}).get('format') in {'ENA', 'SRA'}:
        files = [v for v in series.get('supplementary_data', []) if v.get('value')]
        if files:
            _replace_row(rows, 'Comment[Study supplementary file]', [v['value'] for v in files])
            _replace_row(rows, 'Comment[Study supplementary file type]', [v.get('type', '') for v in files])
    from .harmonized import bind_sample_groups
    assay_table = _render_miniml_assay_paths(bind_sample_groups(package, series.get("assay_paths")),
        preserve_order=package.get('source', {}).get('format') in {'ENA', 'SRA'})
    if assay_table:
        _replace_row(rows, "SDRF File", [assay_table])
    _insert_retained_patch_comments(package, rows)
    return rows


def _insert_retained_patch_comments(package: dict, rows: list) -> None:
    """Add applied-patch evidence comments without manufacturing attributes."""

    if not package.get("extensions", {}).get("msc_harmonization"):
        return
    sdrf_row = next(
        (
            row
            for row in rows
            if row and _normalized(row[0]) == _normalized("SDRF File")
        ),
        None,
    )
    if not sdrf_row or len(sdrf_row) < 2 or not isinstance(sdrf_row[1], list):
        return
    table = sdrf_row[1]
    if not table or not isinstance(table[0], list):
        return
    samples = [item for item in package.get("sample", []) if isinstance(item, dict)]
    from meta_standards_converter.miniml import iter_harmonization_operations
    grouped = {}
    for operation in iter_harmonization_operations(package):
        match = re.match(r"/sample/(\d+)/", operation['path'])
        if match:
            grouped.setdefault(int(match.group(1)), []).append(operation)
    row_columns = [patch_provenance_columns(package, sample, operations=tuple(grouped.get(i, [])))
                   for i, sample in enumerate(samples)]
    by_identity = {s['iid']: columns for s, columns in zip(samples, row_columns)}
    for path in package.get('series', {}).get('assay_paths', []):
        for step in path.get('steps', []):
            if step.get('kind') in ('source', 'sample') and step.get('sample_ref') in by_identity:
                by_identity.setdefault(step.get('name'), by_identity[step['sample_ref']])
    union: list[str] = []
    for columns in row_columns:
        for key in columns:
            if key not in union:
                union.append(key)
    if not union:
        return

    header = table[0]
    characteristic_anchors = [
        index
        for index, label in enumerate(header)
        if str(label).startswith("Characteristics[")
    ]
    anchor = characteristic_anchors[-1] if characteristic_anchors else len(header) - 1
    while anchor + 1 < len(header) and header[anchor + 1] in {
        "Term Source REF",
        "Term Accession Number",
    }:
        anchor += 1
    offset = anchor + 1
    labels = [
        "Comment[msc_harmonization_"
        + key.removeprefix("msc.harmonization.").replace(".", "_")
        + "]"
        for key in union
    ]
    header[offset:offset] = labels
    identity_indexes = [i for i, label in enumerate(header) if label in ('Sample Name', 'Source Name')]
    for row_index, row in enumerate(table[1:]):
        columns = next((by_identity[row[i]] for i in identity_indexes if i < len(row) and row[i] in by_identity), {})
        row[offset:offset] = [columns.get(key, "") or "" for key in union]


def _replace_row(rows: list, label: str, values: list) -> None:
    normalized = _normalized(label)
    for index, row in enumerate(rows):
        if row and _normalized(row[0]) == normalized:
            rows[index] = [label, *values]
            return
    rows.append([label, *values])


def _ontology_text(value):
    return value.get("value", "") if isinstance(value, dict) else (value or "")


def _ontology_field(value, field):
    return value.get(field, "") if isinstance(value, dict) else ""


def _render_miniml_assay_paths(paths, *, preserve_order=False) -> list | None:
    documents = render_miniml_assay_documents(paths, preserve_order=preserve_order)
    if not documents:
        return None
    tables = list(documents.values())
    if len(tables) == 1:
        return tables[0]
    header = tables[0][0]
    if any(table[0] != header for table in tables[1:]):
        raise ValueError("Multiple SDRF documents with different graph layouts cannot be consolidated safely.")
    return [header, *(row for table in tables for row in table[1:])]


def render_miniml_assay_documents(paths, *, preserve_order=False) -> dict[str, list[list]]:
    grouped = {}
    for index, path in enumerate(paths if isinstance(paths, list) else []):
        if isinstance(path, dict):
            grouped.setdefault(str(path.get("document") or "study.sdrf.txt"), []).append(path)
    return {name: _render_assay_path_group(values, preserve_order=preserve_order) for name, values in grouped.items()}


def _render_assay_path_group(paths, *, preserve_order=False) -> list[list]:
    from .sdrf.model import SDRFPath, SDRFNode, SDRFEdge, SDRFAttr
    from .sdrf.renderer import SDRFRenderer
    rendered = []
    for path in paths:
        parts = []
        for step in path.get("steps", []):
            pairs = _miniml_path_columns([step])
            if not pairs:
                continue
            label, value = pairs[0]
            part = SDRFEdge(value) if label == "Protocol REF" else SDRFNode(label, label, value)
            last = named = unit = None
            for label, value in pairs[1:]:
                attr = SDRFAttr(label, value, required=True)
                if last is not None and (label in {"Term Source REF", "Term Accession Number"} or "_hierarchy_depth" in label):
                    last.attrs.append(attr)
                    continue
                if label.startswith("Unit") and named is not None:
                    named.attrs.append(attr)
                    unit = attr
                elif label.startswith("Comment[hz_unit") and unit is not None:
                    unit.attrs.append(attr)
                elif label.startswith(("Parameter Value[hz_", "Factor Value[hz_")) and named is not None:
                    named.attrs.append(attr)
                else:
                    part.attrs.append(attr)
                    if label.startswith(("Characteristics[", "Parameter Value[", "Factor Value[")):
                        named = attr
                        unit = None
                last = attr
            parts.append(part)
        rendered.append(SDRFPath(parts))
    renderer = SDRFRenderer(preserve_order=preserve_order)
    return renderer.render_paths(renderer.plan_columns(rendered), rendered)


def _miniml_path_columns(steps) -> list[tuple[str, object]]:
    result = []
    headers = {
        "source": "Source Name", "sample": "Sample Name", "extract": "Extract Name",
        "labeled_extract": "Labeled Extract Name", "hybridization": "Hybridization Name",
        "assay": "Assay Name", "scan": "Scan Name", "normalization": "Normalization Name",
        "array_data_file": "Array Data File", "array_data_matrix_file": "Array Data Matrix File",
        "derived_array_data_file": "Derived Array Data File",
        "derived_array_data_matrix_file": "Derived Array Data Matrix File", "image_file": "Image File",
    }
    for step in steps if isinstance(steps, list) else []:
        if not isinstance(step, dict):
            continue
        if step.get("kind") == "protocol_application":
            result.append(("Protocol REF", step.get("protocol_ref", "")))
            if step.get("performer") not in (None, ""):
                result.append(("Performer", step["performer"]))
            if step.get("date") not in (None, ""):
                result.append(("Date", step["date"]))
            result.extend(_named_values_columns("Parameter Value", step.get("parameter_values", [])))
            for comment in step.get("comments", []) or []:
                result.append((f"Comment[{comment.get('name', '')}]", comment.get("value", "")))
            continue
        header = headers.get(step.get("kind"))
        if not header:
            continue
        result.append((header, step.get("name", "")))
        if step.get('link', {}).get('value'):
            result.append(('Comment[File URI]', step['link']['value']))
        members = step.get('link', {}).get('companion_files', [])
        from ..miniml.archive_results import companion_node
        for member in members if isinstance(members, list) else []:
            if isinstance(member, dict) and member.get('role') in ('barcodes', 'features') and companion_node(member.get('node')):
                from .native_files import _record, _comments
                result.extend((f"Comment[{c['name']}]", c['value']) for c in
                              _comments(_record(member['node']), 'MATRIX_' + member['role'].upper() + '_'))
        result.extend(_named_values_columns("Characteristics", step.get("characteristics", [])))
        result.extend(_named_values_columns("Factor Value", step.get("factor_values", [])))
        for field, field_header in (
            ("provider", "Provider"), ("material_type", "Material Type"),
            ("description", "Description"), ("label", "Label"),
            ("technology_type", "Technology Type"),
        ):
            if step.get(field) not in (None, ""):
                result.extend(_ontology_columns(field_header, step[field]))
        reference = step.get("array_design_ref")
        if isinstance(reference, dict) and reference.get("ref"):
            result.append(("Array Design REF", reference["ref"]))
        for comment in step.get("comments", []) or []:
            result.append((f"Comment[{comment.get('name', '')}]", comment.get("value", "")))
    return result


def _named_values_columns(prefix, values):
    from .harmonized import columns
    used, result = set(), []
    for value in values or []:
        if str(value.get("name", "")).startswith("hz_"):
            continue
        result.extend(_named_value_columns(prefix, value))
        result.extend(columns(value, prefix, used))
    result.extend(columns(values or [], prefix, used))
    return result


def _named_value_columns(prefix: str, value) -> list[tuple[str, object]]:
    if not isinstance(value, dict):
        return []
    qualifier = f" ({value['qualifier']})" if value.get("qualifier") else ""
    result = [(f"{prefix}[{value.get('name', '')}]{qualifier}", value.get("value", ""))]
    if value.get("term_source_ref"):
        result.append(("Term Source REF", value["term_source_ref"]))
    if value.get("term_accession_number"):
        result.append(("Term Accession Number", value["term_accession_number"]))
    if value.get("unit") is not None:
        unit_header = f"Unit[{value['unit_type']}]" if value.get("unit_type") else "Unit"
        result.extend(_ontology_columns(unit_header, value["unit"]))
    for comment in value.get("comments", []) or []:
        result.append((f"Comment[{comment.get('name', '')}]", comment.get("value", "")))
    return result


def _ontology_columns(header: str, value) -> list[tuple[str, object]]:
    result = [(header, _ontology_text(value))]
    source = _ontology_field(value, "term_source_ref")
    accession = _ontology_field(value, "term_accession_number")
    if source:
        result.append(("Term Source REF", source))
    if accession:
        result.append(("Term Accession Number", accession))
    from .harmonized import columns
    if isinstance(value, dict):
        result.extend(columns(value, "Comment"))
    return result


def overlay_core(model_rows: list, core_rows: list) -> list:
    """Union MINiML projections into model tables while preserving model structure."""
    result = copy.deepcopy(model_rows)
    for row in result:
        if row:
            row[0] = _canonical_idf_label(row[0])
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
            "Publication Status Term Source REF",
            "Publication Status Term Accession Number",
            "Experiment Description", "Term Source Name", "Term Source File", "Term Source Version",
        )
    }
    core_by_label = {
        _normalized(_canonical_idf_label(row[0])): [
            _canonical_idf_label(row[0]), *row[1:]
        ]
        for row in core_rows
        if row
    }
    _overlay_protocol_rows(result, core_rows)
    for index, row in enumerate(result):
        label = _normalized(row[0]) if row else ""
        if label in replace_labels and label in core_by_label:
            result[index] = copy.deepcopy(core_by_label[label])
    _insert_missing_idf_rows(result, core_rows, replace_labels)

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
        from .harmonized import read_group
        group = read_group(header, row, index)
        if group is not None:
            prefix, harmonized, end = group
            consumed.update(range(index + 1, end))
            if harmonized is not None:
                base.update(kind="harmonized", prefix=prefix, harmonized=harmonized.to_annotation_mapping())
                steps.append(base)
            continue
        annotation = re.fullmatch(r"\s*(Characteristics|Factor\s+Value|Parameter\s+Value)\s*\[([^]]*)]\s*(?:\(([^)]*)\))?\s*", label, re.I)
        if annotation:
            base.update({
                "kind": "attribute",
                "attribute_type": " ".join(annotation.group(1).split()).lower(),
                "name": annotation.group(2).strip(),
                **({"qualifier": annotation.group(3).strip()} if annotation.group(3) else {}),
            })
            companions = {}
            unit_seen = False
            for companion_index in range(index + 1, len(header)):
                companion_label = header[companion_index]
                normalized = _normalized(companion_label)
                unit_match = re.fullmatch(r"\s*Unit(?:\[([^]]*)])?\s*", companion_label, re.I)
                if unit_match:
                    field = "unit"
                    unit_seen = True
                    if unit_match.group(1):
                        base["unit_type"] = unit_match.group(1).strip()
                elif normalized == _normalized("Term Source REF"):
                    field = "unit_term_source_ref" if unit_seen else "term_source_ref"
                elif normalized == _normalized("Term Accession Number"):
                    field = "unit_term_accession_number" if unit_seen else "term_accession_number"
                else:
                    break
                companion_value = row[companion_index] if companion_index < len(row) else ""
                base[field] = companion_value
                companions[field] = companion_index
                consumed.add(companion_index)
            if companions:
                base["companion_columns"] = companions
            annotation_index = max([index, *companions.values()]) + 1
            while annotation_index < len(header):
                candidate = read_group(header, row, annotation_index)
                if candidate is not None:
                    # Old flat Comment companions are accepted only when actually
                    # present; new ontology/depth groups are parsed independently.
                    next_label = header[annotation_index + 1] if annotation_index + 1 < len(header) else ""
                    if not re.fullmatch(r"Comment\[hz_(?:value|unit)_(?:id|onto)\]", next_label):
                        break
                harmonized = re.fullmatch(
                    r"\s*Comment\[(hz_(?:value|unit)(?:_id|_onto|_hierarchy_depth)?|hz_field)]\s*",
                    header[annotation_index],
                    re.I,
                )
                if harmonized is None:
                    break
                field = harmonized.group(1).casefold()
                base[field] = row[annotation_index] if annotation_index < len(row) else ""
                consumed.add(annotation_index)
                annotation_index += 1
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
    return _insert_harmonization_columns(rows, paths)


def _insert_harmonization_columns(
    rows: list[list], paths: list[dict]
) -> list[list]:
    insertions: dict[int, list[str]] = {}
    for path in paths:
        for step in path.get("steps", []):
            if step.get("kind") != "attribute":
                continue
            fields = [field for field in HZ_FIELDS if step.get(field) not in (None, "")]
            if not fields:
                continue
            after = max([
                step.get("column_index", 0),
                *(step.get("companion_columns") or {}).values(),
            ])
            existing = insertions.setdefault(after, [])
            for field in fields:
                if field not in existing:
                    existing.append(field)
    for after in sorted(insertions, reverse=True):
        fields = insertions[after]
        offset = after + 1
        rows[0][offset:offset] = [f"Comment[{field}]" for field in fields]
        for row, path in zip(rows[1:], paths):
            step = next(
                (
                    item for item in path.get("steps", [])
                    if item.get("kind") == "attribute"
                    and max([
                        item.get("column_index", 0),
                        *(item.get("companion_columns") or {}).values(),
                    ]) == after
                ),
                {},
            )
            row[offset:offset] = [step.get(field, "") for field in fields]
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
    _insert_missing_sdrf_columns(model, core)
    model_header, core_header = model[0], core[0]
    model_keys = _header_keys(model_header)
    core_indexes = {key: index for index, key in enumerate(_header_keys(core_header))}
    protected = {_normalized(label) for label in NODE_HEADERS | {"Protocol REF"}}
    core_rows_by_identity = _rows_by_identity(core_header, core[1:])
    for row in model[1:]:
        core_rows = _matching_core_rows(
            model_header,
            row,
            core[1:],
            core_rows_by_identity,
        )
        if not core_rows:
            continue
        for model_index, key in enumerate(model_keys):
            core_index = core_indexes.get(key)
            if key[0] in protected or core_index is None:
                continue
            values = {
                core_row[core_index]
                for core_row in core_rows
                if core_index < len(core_row)
            }
            if len(values) == 1:
                row[model_index] = values.pop()


def _insert_missing_idf_rows(
    model_rows: list[list],
    core_rows: list[list],
    allowed_labels: set[str],
) -> None:
    core_keys = [_normalized(row[0]) if row else "" for row in core_rows]
    model_keys = [_normalized(row[0]) if row else "" for row in model_rows]
    for core_index, row in enumerate(core_rows):
        key = core_keys[core_index]
        if not row or key not in allowed_labels or key in model_keys:
            continue
        insert_at = _anchored_insert_index(model_keys, core_keys, core_index)
        model_rows.insert(insert_at, copy.deepcopy(row))
        model_keys.insert(insert_at, key)


def _insert_missing_sdrf_columns(model: list[list], core: list[list]) -> None:
    model_header, core_header = model[0], core[0]
    model_keys = _header_keys(model_header)
    core_keys = _header_keys(core_header)
    protected = {_normalized(label) for label in NODE_HEADERS | {"Protocol REF"}}
    for core_index, key in enumerate(core_keys):
        if not key[0] or key[0] in protected or key in model_keys:
            continue
        insert_at = _anchored_insert_index(model_keys, core_keys, core_index)
        model_header.insert(insert_at, core_header[core_index])
        for row in model[1:]:
            row.insert(insert_at, "")
        model_keys.insert(insert_at, key)


def _anchored_insert_index(current_keys: list, source_keys: list, source_index: int) -> int:
    for key in source_keys[source_index + 1:]:
        if key in current_keys:
            return current_keys.index(key)
    for key in reversed(source_keys[:source_index]):
        if key in current_keys:
            return len(current_keys) - current_keys[::-1].index(key)
    return len(current_keys)


def _header_keys(header: list) -> list[tuple[str, int]]:
    occurrences = {}
    keys = []
    for label in header:
        normalized = _normalized(label)
        occurrences[normalized] = occurrences.get(normalized, 0) + 1
        keys.append((normalized, occurrences[normalized]))
    return keys


def _rows_by_identity(header: list, rows: list[list]) -> dict[str, set[int]]:
    result = {}
    for index, row in enumerate(rows):
        for identity in set(_identities(header, row)):
            result.setdefault(identity, set()).add(index)
    return result


def _matching_core_rows(
    model_header: list,
    model_row: list,
    core_rows: list[list],
    rows_by_identity: dict[str, set[int]],
) -> list[list]:
    candidates = None
    for identity in dict.fromkeys(_identities(model_header, model_row)):
        indexes = rows_by_identity.get(identity)
        if not indexes:
            continue
        candidates = set(indexes) if candidates is None else candidates & indexes
    if not candidates:
        return []
    return [core_rows[index] for index in sorted(candidates)]


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


def _canonical_idf_label(value: str) -> str:
    normalized = _normalized(value)
    for canonical, aliases in IDF_LABEL_ALIASES.items():
        if any(normalized == _normalized(alias) for alias in aliases):
            return canonical
    return value


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
