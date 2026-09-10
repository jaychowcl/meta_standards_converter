# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from __future__ import annotations
import json
import hashlib
import os
import re
from dataclasses import field
from pathlib import Path
from typing import Any, Mapping, Sequence
from urllib.parse import urlparse
from meta_standards_converter.miniml import (
    harmonized_value_mappings,
    iter_harmonization_operations,
    iter_harmonization_patches,
)
from meta_standards_converter.converters.dataset_combination import DatasetCombinationPolicy
from meta_standards_converter.metadata.projection.assay import _parameter_rows, _parameter_summary
from meta_standards_converter.metadata.provenance import patch_provenance_columns
from meta_standards_converter.expression.assets import Asset
from meta_standards_converter.expression.readers import scientific_modules

class AnnDataNormalizer:
    MINIML_SCHEMA_VERSION = "1.0"

    H5AD_METADATA_SCHEMA_VERSION = "2.0"

    PUBLICATION_POLICY = "citation_metadata_only"

    OBS_METADATA_FIELDS = {
        "msc.sample.title": "title",
        "msc.sample.description": "description",
        "msc.sample.channel.organism.value": "organism",
        "msc.sample.channel.organism.taxid": "organism_taxid",
        "msc.sample.channel.organism_part": "organism_part",
        "msc.sample.channel.developmental_stage": "developmental_stage",
        "msc.sample.channel.disease": "disease",
        "msc.sample.channel.genotype": "genotype",
        "msc.sample.channel.source": "source",
        "msc.sample.channel.biomaterial_provider": "biomaterial_provider",
        "msc.sample.channel.material_type": "material_type",
        "msc.sample.channel.molecule": "molecule",
        "msc.platform.accession": "platform_accession",
        "msc.archive.sra_accession": "sra_accession",
        "msc.archive.ena_accession": "ena_accession",
        "msc.archive.biosample_accession": "biosample_accession",
        "msc.archive.sra_run_accessions": "sra_run_accessions",
        "msc.library.strategy": "library_strategy",
        "msc.library.source": "library_source",
        "msc.library.selection": "library_selection",
        "msc.library.layout": "library_layout",
        "msc.instrument.model": "instrument_model",
        "msc.protocol.types": "protocol_types",
        "msc.protocol.term_source_refs": "protocol_term_source_refs",
        "msc.protocol.term_accession_numbers": "protocol_term_accession_numbers",
        "msc.database.identifier": "metadata_source",
        "msc.database.name": "metadata_source_name",
        "msc.database.uri": "metadata_source_uri",
    }

    PUBLICATION_FIELDS = (
        "pubmed_id",
        "doi",
        "title",
        "author_list",
        "status",
        "status_term_source_ref",
        "status_term_accession_number",
    )

    def __init__(self, *, metadata_service, planner, localize, package_version, combination_policy=None):
        self.metadata_service = metadata_service
        self.planner = planner
        self._local_path = localize
        self._package_version = package_version
        self.combination_policy = combination_policy or DatasetCombinationPolicy(
            scientific_modules=scientific_modules,
            attach_sample_values=self._attach_sample_values,
            package_version=package_version,
            metadata_schema_version=self.H5AD_METADATA_SCHEMA_VERSION)

    def normalize(
        self,
        adata,
        sample: dict,
        package: dict,
        study_accession: str,
        asset: Asset,
        characteristic_columns: list[str],
        artifact_parent: Path,
        harmonization_resolution=None,
    ) -> dict:
        sample_id = self.planner.sample_accession(sample)
        metadata_values = self._sample_metadata_values(sample, package)
        metadata = self._render_sample_metadata(metadata_values)
        modality = self._sample_modality(sample)
        original_names = [str(value) for value in adata.obs_names]
        adata.obs["msc.observation.original_id"] = original_names
        candidates = [
            value
            if self._is_sample_qualified(value, sample_id)
            else f"{value}-{sample_id}"
            for value in original_names
        ]
        adata.obs_names = self._deduplicate_observation_ids(candidates)
        source_uri, source_uri_scope = self._portable_location(asset.path, artifact_parent)
        canonical_values = {
            "msc.sample.accession": (sample_id,),
            "msc.series.accession": (study_accession,),
            **{
                column: metadata_values.get(key, ())
                for column, key in self.OBS_METADATA_FIELDS.items()
            },
            "msc.asset.tier": (asset.kind,),
            "msc.asset.uri": (source_uri,) if source_uri else (),
            "msc.asset.uri_scope": (source_uri_scope,) if source_uri_scope else (),
            "msc.expression.modality": (modality,),
        }
        canonical_values.update(_parameter_summary(package, sample))
        for column in characteristic_columns:
            canonical_values[f"msc.characteristics.{column}"] = metadata_values[
                "characteristics"
            ].get(column, ())
        for item in getattr(harmonization_resolution, "selections", ()):
            if item.sample_accession != sample_id:
                continue
            prefix = f"msc.harmonization.{item.destination}"
            canonical_values.update({
                f"{prefix}.value": (item.value,),
                f"{prefix}.id": (item.identifier,) if item.identifier else (),
                f"{prefix}.ontology": (item.ontology,) if item.ontology else (),
                f"{prefix}.source_field": (item.source_field,),
                f"{prefix}.hierarchy_depth": (
                    (item.hierarchy_depth,) if item.hierarchy_depth is not None else ()
                ),
            })
        for key, value in patch_provenance_columns(
            package,
            sample,
            occupied=set(canonical_values),
        ).items():
            canonical_values[key] = () if value is None else (value,)
        for key, values in canonical_values.items():
            adata.obs[key] = self._join_values(values)
        self._attach_sample_values(
            adata,
            {sample_id: canonical_values},
        )
        provenance = {
            "study_accession": study_accession,
            "sample_accession": sample_id,
            "source_tier": asset.kind,
            "source_uri": source_uri,
            "source_uri_scope": source_uri_scope,
            "path_base": "artifact_parent",
            "source_origin": asset.source,
            "source_sha256": self._sha256(asset.path, md5=asset.md5),
            "converter_version": self._package_version(),
            "metadata_schema_version": self.H5AD_METADATA_SCHEMA_VERSION,
            "modality": modality,
        }
        declared_reference = asset.reference or self._declared_reference(adata)
        if declared_reference:
            provenance["reference"] = declared_reference
        for key in (
            "annotation_source",
            "annotation_format",
            "annotation_sha256",
            "effective_annotation",
        ):
            value = getattr(asset, key)
            if value:
                if key in {"annotation_source", "effective_annotation"}:
                    portable, scope = self._portable_location(value, artifact_parent)
                    provenance[key] = portable
                    provenance[f"{key}_scope"] = scope
                else:
                    provenance[key] = value
        existing = adata.uns.get("meta_standards_converter")
        if isinstance(existing, dict):
            provenance = {**existing, **provenance}
        adata.uns["meta_standards_converter"] = provenance
        return metadata

    def _sample_metadata(self, sample: dict, package: dict) -> dict:
        return self.metadata_service.sample_metadata(sample, package)

    def _sample_metadata_values(self, sample: dict, package: dict) -> dict:
        return self.metadata_service.sample_metadata_values(sample, package)

    def _render_sample_metadata(self, values: Mapping[str, Any]) -> dict:
        return self.metadata_service.render_sample_metadata(values)

    def _attach_sample_values(
        self,
        adata,
        samples: Mapping[str, Mapping[str, Sequence[Any]]],
    ) -> None:
        _anndata, _numpy, pandas, _sparse = scientific_modules()
        rows = []
        for sample_accession, fields in samples.items():
            for field_name, values in fields.items():
                for ordinal, value in enumerate(values):
                    if value is None or str(value) == "":
                        continue
                    rows.append(
                        (
                            str(sample_accession),
                            str(field_name),
                            ordinal,
                            str(value),
                            self._metadata_value_type(value),
                        )
                    )
        frame = pandas.DataFrame(
            rows,
            columns=(
                "sample_accession",
                "field",
                "ordinal",
                "value",
                "value_type",
            ),
        )
        frame.index = [f"value_{index:06d}" for index in range(len(frame))]
        adata.uns["msc_metadata"] = {
            "schema_version": self.H5AD_METADATA_SCHEMA_VERSION,
            "sample_values": frame,
        }

    @staticmethod
    def _metadata_value_type(value: Any) -> str:
        if isinstance(value, bool):
            return "boolean"
        if isinstance(value, int):
            return "integer"
        if isinstance(value, float):
            return "number"
        return "string"

    @staticmethod
    def _is_sample_qualified(observation_id: str, sample_id: str) -> bool:
        token = re.escape(str(sample_id))
        return re.search(
            rf"(?:^|[-_.:]){token}(?:$|[-_.:])",
            str(observation_id),
            flags=re.IGNORECASE,
        ) is not None

    @staticmethod
    def _deduplicate_observation_ids(values: Sequence[str]) -> list[str]:
        totals = {}
        for value in values:
            totals[value] = totals.get(value, 0) + 1
        seen = {}
        rendered = []
        for value in values:
            if totals[value] == 1:
                rendered.append(value)
                continue
            seen[value] = seen.get(value, 0) + 1
            rendered.append(f"{value}-{seen[value]}")
        return rendered

    def _ensure_global_observation_ids(self, adatas: Mapping[str, Any]) -> None:
        locations: dict[str, list[tuple[str, int]]] = {}
        for sample_id, adata in adatas.items():
            for position, value in enumerate(adata.obs_names.astype(str)):
                locations.setdefault(value, []).append((sample_id, position))
        collisions = {
            value: entries
            for value, entries in locations.items()
            if len(entries) > 1
        }
        if not collisions:
            return
        for value, entries in collisions.items():
            for sample_id, position in entries:
                names = list(adatas[sample_id].obs_names.astype(str))
                names[position] = f"{value}-{sample_id}"
                adatas[sample_id].obs_names = names
        all_values = [
            value
            for adata in adatas.values()
            for value in adata.obs_names.astype(str)
        ]
        if len(all_values) != len(set(all_values)):
            raise ValueError("Observation identifiers remain non-unique after sample qualification.")

    def ensure_observation_ids(
        self,
        adata,
        sample_id: str,
        used: set[str],
    ) -> None:
        """Make one sample globally unique without retaining earlier matrices."""

        rendered: list[str] = []
        local: set[str] = set()
        for original in adata.obs_names.astype(str):
            candidate = original
            ordinal = 1
            while candidate in used or candidate in local:
                suffix = f"-{sample_id}" if ordinal == 1 else f"-{sample_id}-{ordinal}"
                candidate = f"{original}{suffix}"
                ordinal += 1
            rendered.append(candidate)
            local.add(candidate)
        adata.obs_names = rendered
        used.update(local)

    def characteristic_columns(self, packages: list[dict]) -> list[str]:
        columns = []
        for package in packages:
            for sample in self.planner._as_list(package.get("sample")):
                if not isinstance(sample, dict):
                    continue
                for channel in self.planner._as_list(sample.get("channel")):
                    if not isinstance(channel, dict):
                        continue
                    for annotation in harmonized_value_mappings(channel):
                        if not isinstance(annotation, dict) or not annotation.get("field"):
                            continue
                        annotation_slug = "hz_" + self._metadata_slug(annotation["field"])
                        for candidate in (annotation_slug, f"{annotation_slug}_id", f"{annotation_slug}_onto"):
                            if candidate not in columns:
                                columns.append(candidate)
                    characteristic_rows = self.planner._as_list(
                        channel.get("characteristics")
                    )
                    for item in characteristic_rows:
                        if not isinstance(item, dict):
                            continue
                        slug = self._metadata_slug(item.get("name", item.get("tag")))
                        if str(item.get("name", item.get("tag", ""))).startswith("hz_"):
                            continue
                        if slug and slug not in columns:
                            columns.append(slug)
                    for annotation in harmonized_value_mappings(characteristic_rows):
                        annotation_slug = "hz_" + self._metadata_slug(
                            annotation["field"]
                        )
                        for candidate in (
                            annotation_slug,
                            f"{annotation_slug}_id",
                            f"{annotation_slug}_onto",
                        ):
                            if candidate not in columns:
                                columns.append(candidate)
        return columns

    def _metadata_slug(self, value) -> str:
        return self.metadata_service.metadata_slug(value)

    def _values(self, values) -> list[str]:
        return self.metadata_service.values(values)

    def _join_values(self, values) -> str:
        return self.metadata_service.join_values(values)

    def _metadata_database(self, package: dict) -> Mapping[str, Any]:
        return self.metadata_service.metadata_database(package)

    def _platform_accessions(self, sample: dict, package: dict) -> str:
        return self._join_values(self._platform_accession_values(sample, package))

    def _platform_accession_values(self, sample: dict, package: dict) -> list[str]:
        return self.metadata_service.platform_accession_values(sample, package)

    def attach_miniml(
        self,
        adata,
        packages: list[dict],
        source_json: str,
        source_json_sha256: str | None,
        sample_id: str | None = None,
        artifact_parent: Path | None = None,
    ) -> None:
        _anndata, _numpy, pandas, _sparse = scientific_modules()
        rows = []
        transported_packages = []
        for package_index, package in enumerate(packages):
            if sample_id is not None and not self._package_has_sample(package, sample_id):
                continue
            to_mapping = getattr(package, "to_mapping", None)
            transported_packages.append(
                to_mapping() if callable(to_mapping) else package
            )
            entities = self._metadata_entities(package, sample_id=sample_id)
            for entity_type, entity_id, entity in entities:
                safe_entity = self._publication_safe(entity)
                self._flatten_metadata(
                    safe_entity,
                    rows=rows,
                    package_index=package_index,
                    entity_type=entity_type,
                    entity_id=entity_id,
                )
        fields = pandas.DataFrame(
            rows,
            columns=(
                "package_index",
                "entity_type",
                "entity_id",
                "path",
                "value",
                "value_type",
            ),
        )
        fields.index = [f"field_{index:06d}" for index in range(len(fields))]
        database = next(
            (
                self._metadata_database(package)
                for package in packages
                if sample_id is None or self._package_has_sample(package, sample_id)
            ),
            {},
        )
        portable_source, source_scope = self._portable_location(
            source_json, artifact_parent or Path.cwd()
        )
        adata.uns["msc_miniml"] = {
            "schema_version": self.MINIML_SCHEMA_VERSION,
            "packages_json": json.dumps(
                transported_packages, sort_keys=True, ensure_ascii=False
            ),
            "source_json": portable_source,
            "source_json_scope": source_scope,
            "path_base": "artifact_parent",
            "source_sha256": source_json_sha256 or "",
            "publication_policy": self.PUBLICATION_POLICY,
            "metadata_source": self._join_values(
                database.get("public_id") or database.get("iid") or database.get("name")
            ),
            "metadata_source_name": self._join_values(database.get("name")),
            "metadata_source_uri": self._join_values(database.get("web_link")),
            "fields": fields,
        }
        parameter_occurrences = []
        for package in packages:
            if sample_id is None:
                parameter_occurrences.extend(_parameter_rows(package))
                continue
            sample = next(
                (
                    item for item in self.planner._as_list(package.get("sample"))
                    if isinstance(item, dict)
                    and self.planner.sample_accession(item) == sample_id
                ),
                None,
            )
            if sample is not None:
                parameter_occurrences.extend(_parameter_rows(package, sample=sample))
        if parameter_occurrences:
            parameters = pandas.DataFrame(parameter_occurrences)
            parameters.index = [
                f"parameter_{index:06d}" for index in range(len(parameters))
            ]
            adata.uns["msc_assay"] = {
                "schema_version": "3.0",
                "parameters": parameters,
            }

    def attach_harmonization(
        self,
        adata,
        resolution,
        *,
        packages: list[dict] | tuple[dict, ...] = (),
        sample_id: str | None = None,
    ) -> None:
        fragments = []
        operations = []
        for package in packages:
            if sample_id is not None and not self._package_has_sample(package, sample_id):
                continue
            fragments.extend(iter_harmonization_patches(package))
            operations.extend(
                iter_harmonization_operations(package, sample=sample_id)
            )
        if not fragments and (resolution is None or not resolution.enabled):
            return
        _anndata, _numpy, pandas, _sparse = scientific_modules()
        rows = [
            {
                "sample_accession": item.sample_accession,
                "destination": item.destination,
                "value": item.value,
                "id": item.identifier or "",
                "ontology": item.ontology or "",
                "source_field": item.source_field,
                "hierarchy_depth": (
                    item.hierarchy_depth if item.hierarchy_depth is not None else -1
                ),
                "status": item.status,
            }
            for item in getattr(resolution, "selections", ())
        ]
        selections = pandas.DataFrame(rows, columns=(
            "sample_accession", "destination", "value", "id", "ontology",
            "source_field", "hierarchy_depth", "status",
        ))
        selections.index = [f"selection_{index:06d}" for index in range(len(selections))]
        operation_rows = []
        for item in operations:
            value = item["harmonized_value"]
            evidence = item.get("source_evidence") or {}
            operation_rows.append({
                "path": item["path"],
                "field": value["field"],
                "value": value["value"],
                "id": value.get("term_accession_number") or "",
                "ontology": value.get("term_source_ref") or "",
                "hierarchy_depth": (
                    value.get("hierarchy_depth")
                    if value.get("hierarchy_depth") is not None
                    else -1
                ),
                "source_field": evidence.get("source_field") or "",
                "source_label": evidence.get("source_label") or "",
                "source_path": evidence.get("source_value_path") or "",
                "match_kind": evidence.get("match_kind") or "",
                "patch_id": item["patch_id"],
            })
        operation_table = pandas.DataFrame(operation_rows, columns=(
            "path", "field", "value", "id", "ontology", "hierarchy_depth",
            "source_field", "source_label", "source_path", "match_kind",
            "patch_id",
        ))
        operation_table.index = [
            f"operation_{index:06d}" for index in range(len(operation_table))
        ]
        adata.uns["msc_harmonization"] = {
            "schema_version": "1.0",
            "enabled": bool(resolution is not None and resolution.enabled),
            "applied": bool(getattr(resolution, "applied", False)),
            "profile": dict(getattr(resolution, "profile", None) or {}),
            "selections": selections,
            "warnings": list(getattr(resolution, "warnings", ())),
            "patches_json": json.dumps(
                fragments, sort_keys=True, ensure_ascii=False
            ),
            "operations": operation_table,
        }

    def _package_has_sample(self, package: dict, sample_id: str) -> bool:
        return any(
            self.planner.sample_accession(sample) == sample_id
            for sample in self.planner._as_list(package.get("sample"))
            if isinstance(sample, dict)
        )

    def _metadata_entities(self, package: dict, sample_id: str | None):
        entity_groups = {
            "database": [
                item for item in self.planner._as_list(package.get("database")) if isinstance(item, dict)
            ],
            "contributor": [
                item for item in self.planner._as_list(package.get("contributor")) if isinstance(item, dict)
            ],
            "platform": [
                item for item in self.planner._as_list(package.get("platform")) if isinstance(item, dict)
            ],
            "sample": [
                item for item in self.planner._as_list(package.get("sample")) if isinstance(item, dict)
            ],
        }
        series = package.get("series") if isinstance(package.get("series"), dict) else None
        package_scalars = {
            key: value
            for key, value in package.items()
            if key not in {"database", "contributor", "platform", "sample", "series"}
        }
        entities = [("package", "package", package_scalars)]
        if series is not None:
            entities.append(("series", self._entity_id("series", series, 0), series))

        if sample_id is None:
            for entity_type in ("database", "contributor", "platform", "sample"):
                for index, entity in enumerate(entity_groups[entity_type]):
                    entities.append((entity_type, self._entity_id(entity_type, entity, index), entity))
            return entities

        selected_sample = next(
            (
                sample
                for sample in entity_groups["sample"]
                if self.planner.sample_accession(sample) == sample_id
            ),
            None,
        )
        if selected_sample is None:
            return entities
        selected = {"sample": [selected_sample], "platform": [], "contributor": [], "database": []}
        references = self._metadata_references(selected_sample)
        if series is not None:
            references.update(self._metadata_references(series))
        references.discard("")
        changed = True
        while changed:
            changed = False
            for entity_type in ("platform", "contributor", "database"):
                for index, entity in enumerate(entity_groups[entity_type]):
                    if entity in selected[entity_type]:
                        continue
                    identifiers = self._entity_identifiers(entity_type, entity, index)
                    if not references.intersection(identifiers):
                        continue
                    selected[entity_type].append(entity)
                    references.update(self._metadata_references(entity))
                    changed = True
        database = self._metadata_database(package)
        if database and database not in selected["database"]:
            selected["database"].append(database)
        for entity_type in ("database", "contributor", "platform", "sample"):
            for index, entity in enumerate(selected[entity_type]):
                entities.append((entity_type, self._entity_id(entity_type, entity, index), entity))
        return entities

    def _metadata_references(self, value) -> set[str]:
        references = set()

        def visit(current, key=None):
            if isinstance(current, dict):
                for child_key, child in current.items():
                    if child_key == "sample_ref":
                        continue
                    if child_key.endswith("_ref") or child_key == "database":
                        references.update(self._values(child))
                    else:
                        visit(child, child_key)
            elif isinstance(current, list):
                for child in current:
                    visit(child, key)

        visit(value)
        return references

    def _entity_identifiers(self, entity_type: str, entity: dict, index: int) -> set[str]:
        return {
            value
            for value in (
                self._entity_id(entity_type, entity, index),
                self._join_values(entity.get("iid")),
                *self._values(entity.get("accession")),
                self._join_values(entity.get("public_id")),
                self._join_values(entity.get("name")),
            )
            if value
        }

    def _entity_id(self, entity_type: str, entity: dict, index: int) -> str:
        if entity_type == "sample":
            return self.planner.sample_accession(entity) or self._join_values(entity.get("iid")) or f"sample_{index}"
        prefixes = {"series": "GSE", "platform": "GPL"}
        prefix = prefixes.get(entity_type)
        for value in self._values(entity.get("accession")):
            if not prefix or value.upper().startswith(prefix):
                return value
        return (
            self._join_values(entity.get("public_id"))
            or self._join_values(entity.get("iid"))
            or self._join_values(entity.get("name"))
            or f"{entity_type}_{index}"
        )

    def _publication_safe(self, value):
        if isinstance(value, list):
            return [self._publication_safe(item) for item in value]
        if not isinstance(value, dict):
            return value
        safe = {}
        for key, child in value.items():
            if key == "pubmed_publication":
                publications = []
                for publication in self.planner._as_list(child):
                    if not isinstance(publication, dict):
                        continue
                    publications.append(
                        {
                            field: self._publication_safe(publication[field])
                            for field in self.PUBLICATION_FIELDS
                            if field in publication
                        }
                    )
                safe[key] = publications
            else:
                safe[key] = self._publication_safe(child)
        return safe

    def _flatten_metadata(
        self,
        value,
        rows: list,
        package_index: int,
        entity_type: str,
        entity_id: str,
        path: str = "",
    ) -> None:
        if isinstance(value, dict):
            if not value:
                rows.append((package_index, entity_type, entity_id, path, "", "empty_object"))
                return
            for key, child in value.items():
                child_path = f"{path}.{key}" if path else str(key)
                self._flatten_metadata(
                    child,
                    rows,
                    package_index,
                    entity_type,
                    entity_id,
                    child_path,
                )
            return
        if isinstance(value, list):
            if not value:
                rows.append((package_index, entity_type, entity_id, path, "", "empty_list"))
                return
            for index, child in enumerate(value):
                self._flatten_metadata(
                    child,
                    rows,
                    package_index,
                    entity_type,
                    entity_id,
                    f"{path}[{index}]",
                )
            return
        if value is None:
            value_type, serialized = "null", ""
        elif isinstance(value, bool):
            value_type, serialized = "boolean", "true" if value else "false"
        elif isinstance(value, int):
            value_type, serialized = "integer", str(value)
        elif isinstance(value, float):
            value_type, serialized = "number", repr(value)
        else:
            value_type, serialized = "string", str(value)
        rows.append((package_index, entity_type, entity_id, path, serialized, value_type))

    def _sample_modality(self, sample: dict) -> str:
        return self.metadata_service.sample_modality(sample)

    def _combine(
        self,
        adatas: dict[str, object],
        *,
        allow_unverified: bool = False,
    ):
        return self.combination_policy.combine(
            adatas,
            allow_unverified=allow_unverified,
        )

    def _missing_combination_evidence(
        self,
        adatas: Mapping[str, object],
    ) -> dict[str, list[str]]:
        return self.combination_policy.missing_combination_evidence(adatas)

    def _feature_namespace(self, adata) -> str:
        return self.combination_policy.feature_namespace(adata)

    def _declared_reference(self, adata) -> str | None:
        return self.combination_policy.declared_reference(adata)

    def _sha256(self, path: str, md5: str | None = None) -> str | None:
        local = self._local_path(path, md5=md5)
        if not os.path.isfile(local):
            return None
        digest = hashlib.sha256()
        with open(local, "rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        return digest.hexdigest()

    def _text(self, value) -> str | None:
        return self.metadata_service.text(value)

    def _portable_location(self, value: str | None, base: Path) -> tuple[str | None, str | None]:
        if value in (None, ""):
            return value, None
        rendered = str(value)
        parsed = urlparse(rendered)
        if parsed.scheme and parsed.scheme.lower() != "file":
            return rendered, "remote"
        local = parsed.path if parsed.scheme.lower() == "file" else rendered
        absolute = Path(local)
        if not absolute.is_absolute():
            absolute = Path.cwd() / absolute
        relative = os.path.relpath(absolute.resolve(strict=False), Path(base).resolve(strict=False))
        scope = "external" if relative == ".." or relative.startswith(f"..{os.sep}") else "internal"
        return relative, scope
