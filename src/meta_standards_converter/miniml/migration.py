# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Explicit one-way migration from legacy MINiML JSON to MSC MINiML 2.0."""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from typing import Any, Mapping

from .model import MINiMLModelError, MINiMLPackage, MINiMLValidationIssue


@dataclass(frozen=True)
class MINiMLMigrationResult:
    package: MINiMLPackage
    diagnostics: tuple[MINiMLValidationIssue, ...] = ()


class MINiMLV1Migrator:
    """Translate legacy/unversioned packages without weakening the v2 codec."""

    def migrate(self, value: Mapping[str, Any]) -> MINiMLMigrationResult:
        if not isinstance(value, Mapping):
            raise MINiMLModelError("legacy MINiML package must be an object")
        version = value.get("miniml_schema_version")
        if version not in (None, "1.0"):
            raise MINiMLModelError(f"cannot migrate MINiML schema version {version!r}")
        legacy = deepcopy(dict(value))
        diagnostics: list[MINiMLValidationIssue] = []
        mage_tab = legacy.pop("mage_tab", None)
        source_format = "MAGE-TAB" if isinstance(mage_tab, Mapping) else "GEO MINiML"
        source_version = legacy.pop("version", None)
        schema_location = legacy.pop("schema_location", None)
        source_documents = self._documents(mage_tab)
        legacy["miniml_schema_version"] = "2.0"
        legacy["source"] = {
            "format": source_format,
            **({"version": str(source_version)} if source_version is not None else {}),
            **({"schema_location": schema_location} if schema_location is not None else {}),
            **({"documents": source_documents} if source_documents else {}),
        }
        self._migrate_core(legacy)
        if isinstance(mage_tab, Mapping):
            self._fold_model(legacy, mage_tab.get("model"))
            if mage_tab.get("roundtrip") is not None:
                diagnostics.append(
                    MINiMLValidationIssue(
                        "/mage_tab/roundtrip",
                        "source_layout_dropped",
                        "Raw IDF/SDRF source layout is not part of MSC MINiML 2.0.",
                    )
                )
        package = MINiMLPackage.from_mapping(legacy)
        diagnostics.extend(package.validate())
        return MINiMLMigrationResult(package, tuple(diagnostics))

    @staticmethod
    def _documents(mage_tab: Any) -> list[dict[str, str]]:
        if not isinstance(mage_tab, Mapping):
            return []
        names: list[tuple[str, str]] = []
        source = mage_tab.get("source")
        if isinstance(source, Mapping):
            if source.get("idf"):
                names.append(("idf", str(source["idf"])))
            for name in source.get("sdrf", []) or []:
                if name:
                    names.append(("sdrf", str(name)))
        model = mage_tab.get("model")
        if isinstance(model, Mapping):
            for item in model.get("sdrfs", []) or []:
                if isinstance(item, Mapping) and item.get("name"):
                    names.append(("sdrf", str(item["name"])))
        roundtrip = mage_tab.get("roundtrip")
        if isinstance(roundtrip, Mapping):
            for item in roundtrip.get("sdrfs", []) or []:
                if isinstance(item, Mapping) and item.get("name"):
                    names.append(("sdrf", str(item["name"])))
        result = [{"kind": kind, "name": name} for kind, name in dict.fromkeys(names)]
        return result

    @classmethod
    def _migrate_core(cls, package: dict[str, Any]) -> None:
        for contributor in cls._mappings(package.get("contributor")):
            contributor.pop("position", None)
            roles = contributor.pop("role", contributor.get("roles"))
            if roles is not None:
                contributor["roles"] = [cls._ontology_value(item) for item in cls._items(roles)]
        for sample in cls._mappings(package.get("sample")):
            for ref_name in ("contact_ref",):
                cls._remove_positions(sample.get(ref_name))
            for channel in cls._mappings(sample.get("channel")):
                channel.pop("position", None)
                cls._fold_harmonization(channel)
                if channel.get("source") is not None:
                    channel["source"] = cls._ontology_value(channel["source"])
                if channel.get("molecule") is not None:
                    channel["molecule"] = cls._ontology_value(channel["molecule"])
                channel["characteristics"] = [
                    cls._named_value(item) for item in cls._items(channel.get("characteristics"))
                ]
        series = package.get("series")
        if isinstance(series, dict):
            for ref_name in ("sample_ref", "contributor_ref", "contact_ref"):
                cls._remove_positions(series.get(ref_name))
            if "type" in series:
                series["type"] = [cls._ontology_value(item) for item in cls._items(series["type"])]

    @classmethod
    def _fold_harmonization(cls, channel: dict[str, Any]) -> None:
        """Move legacy ``hz_*`` values into v2 characteristic annotations."""
        rows = [
            deepcopy(item)
            for item in cls._items(channel.get("characteristics"))
            if isinstance(item, Mapping)
        ]
        grouped: dict[str, dict[str, Any]] = {}

        def collect(label: str, value: Any) -> None:
            if not label.startswith("hz_"):
                return
            body = label[3:]
            qualifier = "value"
            for suffix, candidate in (
                ("_hierarchy_depth", "hierarchy_depth"),
                ("_id", "id"),
                ("_onto", "onto"),
            ):
                if body.endswith(suffix):
                    body = body[:-len(suffix)]
                    qualifier = candidate
                    break
            grouped.setdefault(body, {})[qualifier] = value

        retained = []
        for row in rows:
            name = str(row.get("name", row.get("tag", "")))
            if name.startswith("hz_"):
                collect(name, row.get("value"))
            else:
                retained.append(row)
        for key in tuple(channel):
            if not str(key).startswith("hz_"):
                continue
            value = channel.pop(key)
            if isinstance(value, Mapping):
                body = str(key)[3:]
                grouped.setdefault(body, {}).update({
                    "value": value.get("value"),
                    "id": value.get("id", value.get("term_accession_number")),
                    "onto": value.get("onto", value.get("term_source_ref")),
                    "hierarchy_depth": value.get("hierarchy_depth"),
                })
            else:
                collect(str(key), value)
        raw_label = channel.pop("pre_hz_label", None)

        by_name = {
            str(item.get("name", item.get("tag", ""))).casefold(): item
            for item in retained
        }
        aliases = {
            "species": "organism",
            "organism": "organism",
            "tissue": "tissue",
            "disease": "disease",
            "cell": "cell type",
            "exposure": "exposure",
        }
        for field, item in grouped.items():
            value = item.get("value")
            if value in (None, ""):
                continue
            prefix = field.split("_", 1)[0].casefold()
            target_name = aliases.get(prefix, prefix.replace("_", " "))
            target = by_name.get(target_name.casefold())
            if target is None:
                target = {
                    "name": target_name,
                    "value": str(raw_label if target_name == "organism" and raw_label else value),
                }
                retained.append(target)
                by_name[target_name.casefold()] = target
            annotation = {"field": field, "value": str(value)}
            if item.get("onto") not in (None, ""):
                annotation["term_source_ref"] = str(item["onto"])
            if item.get("id") not in (None, ""):
                annotation["term_accession_number"] = str(item["id"])
            depth = item.get("hierarchy_depth")
            if depth not in (None, ""):
                try:
                    annotation["hierarchy_depth"] = int(depth)
                except (TypeError, ValueError):
                    pass
            target.setdefault("annotations", []).append(annotation)
        channel["characteristics"] = retained

    @classmethod
    def _fold_model(cls, package: dict[str, Any], model: Any) -> None:
        if not isinstance(model, Mapping):
            return
        series = package.get("series")
        if not isinstance(series, dict):
            return
        protocols = []
        for item in model.get("protocols", []) or []:
            if not isinstance(item, Mapping) or not item.get("name"):
                continue
            protocol: dict[str, Any] = {"name": str(item["name"])}
            if item.get("type"):
                protocol["type"] = {
                    "value": str(item["type"]),
                    **({"term_source_ref": item["type_term_source_ref"]} if item.get("type_term_source_ref") else {}),
                    **({"term_accession_number": item["type_term_accession_number"]} if item.get("type_term_accession_number") else {}),
                }
            for source, destination in (
                ("description", "description"),
                ("parameters", "parameters"),
                ("hardware", "hardware"),
                ("software", "software"),
                ("contact", "performers"),
                ("performer", "performers"),
            ):
                if item.get(source):
                    protocol.setdefault(destination, []).extend(cls._items(item[source]))
                    if destination == "description":
                        protocol[destination] = str(item[source])
            protocols.append(protocol)
        if protocols:
            series["protocols"] = protocols
        declarations = model.get("declarations")
        if isinstance(declarations, Mapping):
            for old, new in (
                ("quality_control", "quality_controls"),
                ("replicate", "replicate_types"),
                ("normalization", "normalization_types"),
            ):
                values = [cls._declaration(item) for item in declarations.get(old, []) or []]
                if values:
                    series[new] = values
        paths = []
        for item in model.get("assay_paths", []) or []:
            if not isinstance(item, Mapping):
                continue
            steps = cls._assay_steps(item.get("steps"))
            if steps:
                paths.append({
                    **({"document": str(item["sdrf"])} if item.get("sdrf") else {}),
                    "steps": steps,
                })
        if paths:
            series["assay_paths"] = paths

    @classmethod
    def _assay_steps(cls, values: Any) -> list[dict[str, Any]]:
        result: list[dict[str, Any]] = []
        current_application: dict[str, Any] | None = None
        header_kinds = {
            "source name": "source", "sample name": "sample", "extract name": "extract",
            "labeled extract name": "labeled_extract", "hybridization name": "hybridization",
            "assay name": "assay", "scan name": "scan", "normalization name": "normalization",
            "array data file": "array_data_file", "derived array data file": "derived_array_data_file",
            "array data matrix file": "array_data_matrix_file",
            "derived array data matrix file": "derived_array_data_matrix_file", "image file": "image_file",
        }
        for step in values if isinstance(values, list) else []:
            if not isinstance(step, Mapping):
                continue
            kind = step.get("kind")
            if kind in {"protocol", "protocol_ref"}:
                current_application = {
                    "kind": "protocol_application",
                    "protocol_ref": str(step.get("value") or step.get("protocol_ref") or ""),
                }
                result.append(current_application)
                continue
            if kind == "attribute":
                attribute = cls._legacy_attribute(step)
                if step.get("attribute_type") == "parameter value" and current_application is not None:
                    current_application.setdefault("parameter_values", []).append(attribute)
                elif result and result[-1].get("kind") != "protocol_application":
                    destination = "factor_values" if step.get("attribute_type") == "factor value" else "characteristics"
                    result[-1].setdefault(destination, []).append(attribute)
                continue
            if kind == "comment" and result and result[-1].get("kind") != "protocol_application":
                result[-1].setdefault("comments", []).append({
                    "name": str(step.get("name") or step.get("header") or "comment"),
                    "value": str(step.get("value", "")),
                })
                continue
            if kind == "field" and result and result[-1].get("kind") != "protocol_application":
                cls._attach_node_field(result[-1], step)
                continue
            if kind in {"node", "file"}:
                label = str(step.get("name") or step.get("header") or "").strip().casefold()
                node_kind = header_kinds.get(label)
                if node_kind and step.get("value") not in (None, ""):
                    node = {"kind": node_kind, "name": str(step["value"])}
                    if node_kind == "sample":
                        node["sample_ref"] = str(step["value"])
                    result.append(node)
                    current_application = None
        return result

    @classmethod
    def _attach_node_field(cls, node: dict[str, Any], item: Mapping[str, Any]) -> None:
        header = str(item.get("header") or item.get("name") or "").strip().casefold()
        value = str(item.get("value", ""))
        ontology = cls._ontology_value(value)
        if item.get("term_source_ref"):
            ontology["term_source_ref"] = item["term_source_ref"]
        if item.get("term_accession_number"):
            ontology["term_accession_number"] = item["term_accession_number"]
        destinations = {
            "provider": "provider",
            "description": "description",
            "material type": "material_type",
            "label": "label",
            "technology type": "technology_type",
        }
        destination = destinations.get(header)
        if destination in {"material_type", "label", "technology_type"}:
            node[destination] = ontology
        elif destination:
            node[destination] = value
        elif header == "array design ref":
            node["array_design_ref"] = {"ref": value}
        else:
            node.setdefault("comments", []).append({
                "name": str(item.get("header") or item.get("name") or "field"),
                "value": value,
            })

    @classmethod
    def _legacy_attribute(cls, item: Mapping[str, Any]) -> dict[str, Any]:
        result = {
            "name": str(item.get("name", "")),
            "value": str(item.get("value", "")),
        }
        if item.get("term_source_ref"):
            result["term_source_ref"] = item["term_source_ref"]
        if item.get("term_accession_number"):
            result["term_accession_number"] = item["term_accession_number"]
        if item.get("unit"):
            unit: dict[str, Any] = {"value": str(item["unit"])}
            annotation = cls._legacy_annotation(item, "unit")
            if annotation:
                unit["annotations"] = [annotation]
            result["unit"] = unit
        annotation = cls._legacy_annotation(item, "value")
        if annotation:
            result["annotations"] = [annotation]
        return result

    @staticmethod
    def _legacy_annotation(item: Mapping[str, Any], prefix: str) -> dict[str, Any] | None:
        value = item.get(f"hz_{prefix}")
        if value in (None, ""):
            return None
        result = {"field": str(item.get("hz_field") or prefix), "value": str(value)}
        if item.get(f"hz_{prefix}_onto"):
            result["term_source_ref"] = item[f"hz_{prefix}_onto"]
        if item.get(f"hz_{prefix}_id"):
            result["term_accession_number"] = item[f"hz_{prefix}_id"]
        if item.get(f"hz_{prefix}_hierarchy_depth") is not None:
            result["hierarchy_depth"] = int(item[f"hz_{prefix}_hierarchy_depth"])
        return result

    @classmethod
    def _declaration(cls, value: Any) -> dict[str, Any]:
        if not isinstance(value, Mapping):
            return {"value": str(value)}
        return {
            "value": str(value.get("value", "")),
            **({"term_source_ref": value["term_source_ref"]} if value.get("term_source_ref") else {}),
            **({"term_accession_number": value["term_accession_number"]} if value.get("term_accession_number") else {}),
        }

    @classmethod
    def _named_value(cls, value: Any) -> dict[str, Any]:
        if isinstance(value, Mapping):
            name = value.get("name", value.get("tag"))
            result = {"name": str(name or "unspecified"), "value": str(value.get("value", ""))}
            for key in ("term_source_ref", "term_accession_number", "unit", "annotations", "comments", "qualifier"):
                if value.get(key) is not None:
                    result[key] = deepcopy(value[key])
            return result
        return {"name": "unspecified", "value": str(value)}

    @staticmethod
    def _ontology_value(value: Any) -> dict[str, Any]:
        return deepcopy(dict(value)) if isinstance(value, Mapping) else {"value": str(value)}

    @staticmethod
    def _items(value: Any) -> list[Any]:
        if value is None:
            return []
        return list(value) if isinstance(value, (list, tuple)) else [value]

    @classmethod
    def _mappings(cls, value: Any) -> list[dict[str, Any]]:
        return [item for item in cls._items(value) if isinstance(item, dict)]

    @classmethod
    def _remove_positions(cls, value: Any) -> None:
        for item in cls._mappings(value):
            item.pop("position", None)


__all__ = ["MINiMLMigrationResult", "MINiMLV1Migrator"]
