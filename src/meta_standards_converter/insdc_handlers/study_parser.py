# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Source-faithful projection of SRA/ENA study graphs into MINiML 3.0."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable, Mapping
import json
from pathlib import PurePosixPath
import re
import xml.etree.ElementTree as ET

from meta_standards_converter.miniml import MINiMLCodec, MINiMLPackage
from meta_standards_converter.xml_safety import parse_xml

from .study_models import ProviderDocument, StudyFetchResult


_INSDC_ACCESSION = re.compile(
    r"^(?:PRJ[A-Z]+\d+|SAM[A-Z]+\d+|[SED]RP\d+|[SED]RS\d+|[SED]RX\d+|[SED]RR\d+)$",
    re.I,
)
_GEO_SERIES = re.compile(r"^GSE\d+$", re.I)
_GEO_SAMPLE = re.compile(r"^GSM\d+$", re.I)
_AE_ACCESSION = re.compile(r"^E-[A-Z0-9]+-\d+$", re.I)


def _clean(value) -> str | None:
    if value is None:
        return None
    rendered = str(value).strip()
    return rendered or None


def _strip_namespaces(root: ET.Element) -> ET.Element:
    for node in root.iter():
        node.tag = node.tag.rsplit("}", 1)[-1]
    return root


def _root(document: ProviderDocument) -> ET.Element:
    return _strip_namespaces(
        parse_xml(document.content, max_bytes=max(len(document.content), 1))
    )


def _text(node: ET.Element | None, path: str, default=None) -> str | None:
    if node is None:
        return default
    child = node.find(path)
    return _clean(child.text if child is not None else default)


def _accession(node: ET.Element | None) -> str | None:
    if node is None:
        return None
    value = _clean(node.get("accession")) or _text(node, "./IDENTIFIERS/PRIMARY_ID")
    return value.upper() if value else None


def _append_unique(values: list, value) -> None:
    if value is not None and value not in values:
        values.append(value)


def _all_identifiers(node: ET.Element | None) -> list[tuple[str, str | None]]:
    if node is None:
        return []
    values: list[tuple[str, str | None]] = []
    direct = _accession(node)
    if direct:
        values.append((direct, None))
    for identifier in node.findall("./IDENTIFIERS/*"):
        value = _clean(identifier.text)
        if not value:
            continue
        namespace = _clean(identifier.get("namespace"))
        item = (value, namespace)
        if item not in values:
            values.append(item)
    return values


def _database(accession: str, namespace: str | None = None) -> str | None:
    value = accession.upper()
    if namespace:
        normalized = namespace.casefold()
        if normalized == "geo":
            return "GEO"
        if "biosample" in normalized:
            return "BioSample"
        if "bioproject" in normalized:
            return "BioProject"
        if "arrayexpress" in normalized:
            return "ArrayExpress"
    if value.startswith("PRJ"):
        return "BioProject"
    if value.startswith("SAM"):
        return "BioSample"
    if _GEO_SERIES.fullmatch(value) or _GEO_SAMPLE.fullmatch(value):
        return "GEO"
    if _AE_ACCESSION.fullmatch(value):
        return "ArrayExpress"
    if re.fullmatch(r"[SED]R[PSXR]\d+", value):
        return "INSDC"
    return namespace


def _accession_rows(values: Iterable[tuple[str, str | None]]) -> list[dict]:
    rows: list[dict] = []
    seen: set[tuple[str, str | None]] = set()
    for value, namespace in values:
        normalized = _clean(value)
        if not normalized:
            continue
        item = (normalized, _database(normalized, namespace))
        if item in seen:
            continue
        seen.add(item)
        row = {"value": normalized}
        if item[1]:
            row["database"] = item[1]
        rows.append(row)
    return rows


def _attributes(node: ET.Element | None, prefix: str) -> list[dict]:
    if node is None:
        return []
    rows = []
    for attribute in node.findall(f".//{prefix}_ATTRIBUTE"):
        name = _text(attribute, "./TAG") or _clean(attribute.get("attribute_name"))
        value = _text(attribute, "./VALUE")
        if value is None and attribute.text:
            value = _clean(attribute.text)
        if not name or value is None:
            continue
        row: dict = {"name": name, "value": value}
        unit = _text(attribute, "./UNITS") or _clean(attribute.get("unit"))
        if unit:
            row["unit"] = {"value": unit}
        rows.append(row)
    return rows


def _biosample_attributes(node: ET.Element | None) -> list[dict]:
    if node is None:
        return []
    rows = []
    for attribute in node.findall("./Attributes/Attribute"):
        name = _clean(attribute.get("attribute_name"))
        value = _clean(attribute.text)
        if not name or value is None:
            continue
        row: dict = {"name": name, "value": value}
        unit = _clean(attribute.get("unit"))
        if unit:
            row["unit"] = {"value": unit}
        rows.append(row)
    return rows


def _library(experiment: ET.Element | None) -> dict[str, str | None]:
    descriptor = experiment.find(".//LIBRARY_DESCRIPTOR") if experiment is not None else None
    layout = None
    if descriptor is not None:
        layout_node = descriptor.find("./LIBRARY_LAYOUT")
        if layout_node is not None and list(layout_node):
            layout = list(layout_node)[0].tag.rsplit("}", 1)[-1].upper()
    instrument = None
    if experiment is not None:
        instrument = _text(experiment, ".//PLATFORM/*/INSTRUMENT_MODEL")
    return {
        "library_layout": layout,
        "library_selection": _text(descriptor, "./LIBRARY_SELECTION"),
        "library_source": _text(descriptor, "./LIBRARY_SOURCE"),
        "library_strategy": _text(descriptor, "./LIBRARY_STRATEGY"),
        "instrument_model": instrument,
    }


def _xref_values(node: ET.Element | None) -> list[tuple[str, str]]:
    values: list[tuple[str, str]] = []
    if node is None:
        return values
    for xref in node.findall(".//XREF_LINK"):
        database = _text(xref, "./DB")
        identifier = _text(xref, "./ID") or _text(xref, "./LABEL")
        if database and identifier:
            values.append((database, identifier))
    return values


def _linked_origin_accessions(values: Iterable[str]) -> tuple[list[str], list[str]]:
    geo: list[str] = []
    ae: list[str] = []
    for raw in values:
        for token in re.findall(r"(?:GSE\d+|E-[A-Z0-9]+-\d+)", raw, re.I):
            normalized = token.upper()
            if normalized.startswith("E-GEOD-"):
                normalized = f"GSE{normalized.rsplit('-', 1)[-1]}"
            if _GEO_SERIES.fullmatch(normalized):
                _append_unique(geo, normalized)
            elif _AE_ACCESSION.fullmatch(normalized):
                _append_unique(ae, normalized)
    return geo, ae


def _submission_origin(
    study: str,
    project: str | None,
    brokers: Iterable[str],
    geo: list[str],
    ae: list[str],
) -> str:
    rendered = " ".join(brokers).casefold()
    if geo or "geo" in rendered:
        return "geo_brokered"
    if ae or "arrayexpress" in rendered or "biostudies" in rendered:
        return "arrayexpress_brokered"
    if study.upper().startswith("DRP") or (project or "").upper().startswith("PRJD") or "ddbj" in rendered:
        return "ddbj_origin"
    if study.upper().startswith("ERP") or (project or "").upper().startswith("PRJE"):
        return "ena_direct"
    return "ncbi_direct"


def _status(node: ET.Element | None, prefix: str) -> list[dict]:
    attributes = {
        row["name"].casefold(): row["value"]
        for row in _attributes(node, prefix)
    }
    row = {
        "database": "INSDC",
        "release_date": attributes.get("ena-first-public"),
        "last_update_date": attributes.get("ena-last-update"),
    }
    return [{key: value for key, value in row.items() if value}] if len(row) > 1 else []


def _address(node: ET.Element | None) -> dict | None:
    if node is None:
        return None
    lines = [
        value
        for path in ("./Department", "./Institution", "./Street", "./Line")
        if (value := _text(node, path))
    ]
    row = {
        "line": lines,
        "city": _text(node, "./City"),
        "state": _text(node, "./State"),
        "province": _text(node, "./Sub") or _text(node, "./Province"),
        "postal_code": _clean(node.get("postal_code")) or _text(node, "./PostalCode"),
        "country": _text(node, "./Country"),
    }
    value = {key: item for key, item in row.items() if item not in (None, [])}
    return value or None


def _file_values(row: Mapping, prefix: str) -> list[dict]:
    locations = str(row.get(f"{prefix}_ftp") or "").split(";")
    checksums = str(row.get(f"{prefix}_md5") or "").split(";")
    sizes = str(row.get(f"{prefix}_bytes") or "").split(";")
    formats = str(row.get(f"{prefix}_format") or "").split(";")
    values: list[dict] = []
    for index, raw_location in enumerate(locations):
        location = raw_location.strip()
        if not location:
            continue
        uri = location if "://" in location else f"ftp://{location}"
        item: dict = {
            "uri": uri,
            "filename": PurePosixPath(location).name,
        }
        if index < len(checksums) and checksums[index].strip():
            item["md5"] = checksums[index].strip()
        if index < len(sizes) and sizes[index].strip():
            item["bytes"] = sizes[index].strip()
        if index < len(formats) and formats[index].strip():
            item["format"] = formats[index].strip()
        values.append(item)
    return values


def _unmapped_object_values(kind: str, node: ET.Element) -> list[dict]:
    accession = _accession(node) or _clean(node.get("accession")) or "unknown"
    prefix = f"{kind.upper()}[{accession}]"
    rows: list[dict] = []

    def visit(current: ET.Element, path: str) -> None:
        for name, value in current.attrib.items():
            if name != "accession":
                rows.append({"source_path": f"{path}/@{name}", "value": value})
        children = list(current)
        value = _clean(current.text)
        if value and not children:
            rows.append({"source_path": path, "value": value})
        for child in children:
            visit(child, f"{path}/{child.tag.rsplit('}', 1)[-1]}")

    for child in list(node):
        visit(child, f"{prefix}/{child.tag.rsplit('}', 1)[-1]}")
    return rows


def _pubmed_summaries(documents: Iterable[ProviderDocument]) -> dict[str, dict]:
    summaries: dict[str, dict] = {}
    for document in documents:
        if document.kind != "pubmed_esummary":
            continue
        root = _root(document)
        for docsum in root.findall(".//DocSum"):
            pmid = _text(docsum, "./Id")
            if not pmid:
                continue
            items = {
                item.get("Name"): _clean(item.text)
                for item in docsum.findall("./Item")
                if item.get("Name") and _clean(item.text)
            }
            authors = [
                _clean(item.text)
                for item in docsum.findall('./Item[@Name="AuthorList"]/Item')
                if _clean(item.text)
            ]
            summaries[pmid] = {
                "pubmed_id": pmid,
                "doi": items.get("DOI"),
                "title": items.get("Title"),
                "author_list": ", ".join(authors) or None,
                "status": items.get("RecordStatus"),
            }
    return summaries


class INSDCStudyParser:
    """Normalize one provider study response into a typed MINiML package."""

    def parse(
        self,
        result: StudyFetchResult,
        *,
        enrichment: str = "none",
    ) -> MINiMLPackage:
        if enrichment not in {"none", "geo", "ae"}:
            raise ValueError("enrichment must be none, geo, or ae")
        graph = (
            self._parse_sra(result)
            if result.provider == "sra"
            else self._parse_ena(result)
        )
        return self._to_package(result, graph, enrichment=enrichment)

    def _parse_sra(self, result: StudyFetchResult) -> dict:
        studies: dict[str, ET.Element] = {}
        samples: dict[str, ET.Element] = {}
        experiments: dict[str, ET.Element] = {}
        runs: dict[str, ET.Element] = {}
        submissions: list[ET.Element] = []
        organizations: list[ET.Element] = []
        provider_objects: list[tuple[str, ET.Element]] = []
        biosamples: dict[str, ET.Element] = {}
        bioprojects: list[ET.Element] = []
        for document in result.documents:
            if document.kind == "sra_efetch":
                root = _root(document)
                for package in root.findall(".//EXPERIMENT_PACKAGE"):
                    study = package.find("./STUDY")
                    experiment = package.find("./EXPERIMENT")
                    sample = package.find("./SAMPLE")
                    if study is not None and _accession(study):
                        studies.setdefault(_accession(study), study)
                    if experiment is not None and _accession(experiment):
                        study_ref = _accession(experiment.find("./STUDY_REF"))
                        if study_ref == result.study_accession:
                            experiments.setdefault(_accession(experiment), experiment)
                    if sample is not None and _accession(sample):
                        samples.setdefault(_accession(sample), sample)
                    for run in package.findall(".//RUN_SET/RUN"):
                        if _accession(run):
                            runs.setdefault(_accession(run), run)
                    submissions.extend(package.findall("./SUBMISSION"))
                    organizations.extend(package.findall("./Organization"))
                for study in root.findall("./STUDY"):
                    if _accession(study):
                        studies.setdefault(_accession(study), study)
                provider_objects.extend(("analysis", node) for node in root.findall(".//ANALYSIS"))
                provider_objects.extend(("assembly", node) for node in root.findall(".//ASSEMBLY"))
            elif document.kind == "biosample_efetch":
                for node in _root(document).findall(".//BioSample"):
                    accession = _clean(node.get("accession"))
                    if accession:
                        biosamples[accession] = node
            elif document.kind == "bioproject_efetch":
                bioprojects.extend(_root(document).findall(".//DocumentSummary"))

        experiment_ids = set(experiments)
        runs = {
            key: node
            for key, node in runs.items()
            if _accession(node.find("./EXPERIMENT_REF")) in experiment_ids
        }
        sample_ids = {
            _accession(experiment.find(".//SAMPLE_DESCRIPTOR"))
            for experiment in experiments.values()
        }
        samples = {key: value for key, value in samples.items() if key in sample_ids}
        return {
            "study": studies.get(result.study_accession),
            "samples": samples,
            "experiments": experiments,
            "runs": runs,
            "submissions": submissions,
            "organizations": organizations,
            "biosamples": biosamples,
            "bioprojects": bioprojects,
            "file_rows": [],
            "provider_objects": provider_objects,
        }

    def _parse_ena(self, result: StudyFetchResult) -> dict:
        studies: dict[str, ET.Element] = {}
        samples: dict[str, ET.Element] = {}
        experiments: dict[str, ET.Element] = {}
        runs: dict[str, ET.Element] = {}
        file_rows: list[dict] = []
        xref_values: list[str] = []
        provider_objects: list[tuple[str, ET.Element]] = []
        for document in result.documents:
            if document.kind.startswith("ena_") and document.media_type == "application/xml":
                root = _root(document)
                for node in root.findall(".//STUDY"):
                    if _accession(node):
                        studies.setdefault(_accession(node), node)
                for node in root.findall(".//SAMPLE"):
                    if _accession(node):
                        samples.setdefault(_accession(node), node)
                for node in root.findall(".//EXPERIMENT"):
                    if _accession(node):
                        experiments.setdefault(_accession(node), node)
                for node in root.findall(".//RUN"):
                    if _accession(node):
                        runs.setdefault(_accession(node), node)
                if document.kind == "ena_analysis":
                    provider_objects.extend(("analysis", node) for node in root.findall(".//ANALYSIS"))
                elif document.kind == "ena_assembly":
                    provider_objects.extend(("assembly", node) for node in root.findall(".//ASSEMBLY"))
            elif document.kind == "ena_file_report":
                payload = json.loads(document.content)
                if isinstance(payload, list):
                    file_rows.extend(row for row in payload if isinstance(row, dict))
            elif document.kind == "ena_xref":
                payload = json.loads(document.content)
                for row in payload if isinstance(payload, list) else []:
                    if isinstance(row, Mapping):
                        xref_values.extend(str(value) for value in row.values())

        experiments = {
            key: value
            for key, value in experiments.items()
            if _accession(value.find("./STUDY_REF")) == result.study_accession
        }
        experiment_ids = set(experiments)
        runs = {
            key: value
            for key, value in runs.items()
            if _accession(value.find("./EXPERIMENT_REF")) in experiment_ids
        }
        sample_ids = {
            _accession(experiment.find(".//SAMPLE_DESCRIPTOR"))
            for experiment in experiments.values()
        }
        samples = {key: value for key, value in samples.items() if key in sample_ids}
        return {
            "study": studies.get(result.study_accession),
            "samples": samples,
            "experiments": experiments,
            "runs": runs,
            "submissions": [],
            "organizations": [],
            "biosamples": {},
            "bioprojects": [],
            "file_rows": file_rows,
            "xref_values": xref_values,
            "provider_objects": provider_objects,
        }

    def _to_package(
        self,
        result: StudyFetchResult,
        graph: dict,
        *,
        enrichment: str,
    ) -> MINiMLPackage:
        study = graph["study"]
        if study is None:
            raise ValueError(
                f"{result.provider} response omitted resolved study {result.study_accession}"
            )
        warnings = [dict(item) for item in result.warnings]
        provenance: list[dict] = []
        conflicts: list[dict] = []
        unmapped: list[dict] = []

        mapped_status_attributes = {
            "ENA-STATUS",
            "ENA-FIRST-PUBLIC",
            "ENA-LAST-UPDATE",
        }
        for row in _attributes(study, "STUDY"):
            if row["name"].upper() not in mapped_status_attributes:
                unmapped.append(
                    {
                        "source_path": f"STUDY/STUDY_ATTRIBUTES/{row['name']}",
                        "value": row["value"],
                        **({"unit": row["unit"]["value"]} if row.get("unit") else {}),
                    }
                )
        for experiment_id, experiment in sorted(graph["experiments"].items()):
            for row in _attributes(experiment, "EXPERIMENT"):
                unmapped.append(
                    {
                        "source_path": f"EXPERIMENT[{experiment_id}]/EXPERIMENT_ATTRIBUTES/{row['name']}",
                        "value": row["value"],
                        **({"unit": row["unit"]["value"]} if row.get("unit") else {}),
                    }
                )
        for run_id, run in sorted(graph["runs"].items()):
            for row in _attributes(run, "RUN"):
                if row["name"].upper() not in mapped_status_attributes:
                    unmapped.append(
                        {
                            "source_path": f"RUN[{run_id}]/RUN_ATTRIBUTES/{row['name']}",
                            "value": row["value"],
                            **({"unit": row["unit"]["value"]} if row.get("unit") else {}),
                        }
                    )
        for kind, node in sorted(
            graph.get("provider_objects", []),
            key=lambda item: (item[0], _accession(item[1]) or ""),
        ):
            unmapped.extend(_unmapped_object_values(kind, node))

        study_identifiers = _all_identifiers(study)
        project = next(
            (value for value, _namespace in study_identifiers if value.upper().startswith("PRJ")),
            None,
        )
        link_values = [value for value, _namespace in study_identifiers]
        for database, value in _xref_values(study):
            link_values.extend([database, value])
        for submission in graph["submissions"]:
            link_values.extend(
                value for value, _namespace in _all_identifiers(submission)
            )
        link_values.extend(graph.get("xref_values", []))
        geo_links, ae_links = _linked_origin_accessions(link_values)
        brokers = [
            value
            for node in [study, *graph["submissions"]]
            for value in (node.get("broker_name"), node.get("center_name"))
            if value
        ]
        origin = _submission_origin(
            result.study_accession, project, brokers, geo_links, ae_links
        )
        if enrichment == "none":
            if geo_links:
                warnings.append(
                    {
                        "code": "linked_geo_enrichment_disabled",
                        "message": f"GEO link {geo_links[0]} is available; --enrich-geo is off.",
                    }
                )
            if ae_links:
                warnings.append(
                    {
                        "code": "linked_ae_enrichment_disabled",
                        "message": f"ArrayExpress link {ae_links[0]} is available; --enrich-ae is off.",
                    }
                )

        title = _text(study, "./DESCRIPTOR/STUDY_TITLE") or _text(study, "./TITLE")
        summary = _text(study, "./DESCRIPTOR/STUDY_ABSTRACT") or _text(study, "./DESCRIPTION")
        study_type_node = study.find("./DESCRIPTOR/STUDY_TYPE")
        study_type = (
            _clean(study_type_node.get("existing_study_type"))
            if study_type_node is not None
            else None
        )
        pubmed_ids: list[str] = []
        web_links: list[str] = []
        for database, value in _xref_values(study):
            if database.casefold() in {"pubmed", "pmid"} and value.isdigit():
                _append_unique(pubmed_ids, value)
            elif value.startswith(("http://", "https://", "ftp://")):
                _append_unique(web_links, value)
        for project_node in graph["bioprojects"]:
            for publication in project_node.findall(".//Publication"):
                pmid = _clean(publication.get("id")) or _text(publication, "./Reference")
                if pmid and pmid.isdigit():
                    _append_unique(pubmed_ids, pmid)
        publications = _pubmed_summaries(result.documents)
        publication_rows = [
            {key: value for key, value in publications[pmid].items() if value is not None}
            for pmid in pubmed_ids
            if pmid in publications
        ]

        experiment_by_sample: dict[str, list[tuple[str, ET.Element]]] = defaultdict(list)
        for experiment_id, experiment in graph["experiments"].items():
            sample_id = _accession(experiment.find(".//SAMPLE_DESCRIPTOR"))
            if sample_id:
                experiment_by_sample[sample_id].append((experiment_id, experiment))

        files_by_run: dict[str, dict[str, list[dict]]] = defaultdict(dict)
        for row in graph["file_rows"]:
            run_id = _clean(row.get("run_accession"))
            if not run_id:
                continue
            files_by_run[run_id] = {
                "fastq": _file_values(row, "fastq"),
                "submitted": _file_values(row, "submitted"),
                "sra": _file_values(row, "sra"),
            }

        runs_by_experiment: dict[str, list[tuple[str, ET.Element]]] = defaultdict(list)
        for run_id, run in graph["runs"].items():
            experiment_id = _accession(run.find("./EXPERIMENT_REF"))
            if experiment_id:
                runs_by_experiment[experiment_id].append((run_id, run))

        samples = []
        platforms: list[dict] = []
        platform_ids: dict[str, str] = {}
        protocols: list[dict] = []
        assay_paths: list[dict] = []
        for sample_id in sorted(experiment_by_sample):
            sample = graph["samples"].get(sample_id)
            biosample_id = None
            geo_sample = None
            identifiers = _all_identifiers(sample)
            for value, namespace in identifiers:
                if value.upper().startswith("SAM"):
                    biosample_id = value
                if _GEO_SAMPLE.fullmatch(value):
                    geo_sample = value.upper()
            biosample = graph["biosamples"].get(biosample_id)
            if biosample is not None:
                for value in biosample.findall("./Ids/Id"):
                    identifier = _clean(value.text)
                    if identifier:
                        item = (identifier, _clean(value.get("db")))
                        if item not in identifiers:
                            identifiers.append(item)
                    if (value.get("db") or "").casefold() == "geo" and identifier:
                        geo_sample = identifier.upper()

            taxid = _text(sample, "./SAMPLE_NAME/TAXON_ID")
            organism_name = _text(sample, "./SAMPLE_NAME/SCIENTIFIC_NAME")
            if biosample is not None:
                organism = biosample.find("./Description/Organism")
                taxid = taxid or (organism.get("taxonomy_id") if organism is not None else None)
                organism_name = organism_name or (organism.get("taxonomy_name") if organism is not None else None)
            characteristics = _attributes(sample, "SAMPLE")
            if not characteristics and biosample is not None:
                characteristics = _biosample_attributes(biosample)

            sample_experiments = sorted(experiment_by_sample[sample_id])
            libraries = [_library(experiment) for _, experiment in sample_experiments]
            promoted: dict[str, str] = {}
            for field in (
                "library_layout",
                "library_selection",
                "library_source",
                "library_strategy",
                "instrument_model",
            ):
                values = [value[field] for value in libraries if value[field] is not None]
                unique = list(dict.fromkeys(values))
                if len(unique) == 1 and len(values) == len(libraries):
                    promoted[field] = unique[0]
                elif len(unique) > 1:
                    warnings.append(
                        {
                            "code": "non_unanimous_experiment_value",
                            "message": f"Sample {sample_id} has conflicting {field} values; values remain run-local.",
                            "sample": sample_id,
                            "field": field,
                            "values": unique,
                        }
                    )
                    conflicts.append(
                        {"target": f"sample/{sample_id}/{field}", "values": unique}
                    )

            sample_runs = []
            for (experiment_id, experiment), library in zip(sample_experiments, libraries):
                design_description = _text(experiment, "./DESIGN/DESIGN_DESCRIPTION")
                protocol_name = None
                if design_description:
                    protocol_name = f"{experiment_id} sequencing design"
                    if not any(item["name"] == protocol_name for item in protocols):
                        protocols.append(
                            {
                                "name": protocol_name,
                                "type": {"value": "sequencing protocol"},
                                "description": design_description,
                            }
                        )
                for run_id, run in sorted(runs_by_experiment.get(experiment_id, [])):
                    file_groups = files_by_run.get(run_id, {})
                    run_row: dict = {
                        "run": run_id,
                        "study": result.study_accession,
                        "experiment": experiment_id,
                        "sample": sample_id,
                        "biosample": biosample_id,
                        "geo_sample": geo_sample,
                        "scan_name": _clean(run.get("alias")) or run_id,
                        "fastq_files": file_groups.get("fastq", []),
                        "submitted_files": file_groups.get("submitted", []),
                        "sra_files": file_groups.get("sra", []),
                    }
                    run_row.update({key: value for key, value in library.items() if value})
                    for key, attribute in (
                        ("total_spots", "total_spots"),
                        ("total_bases", "total_bases"),
                        ("bytes", "size"),
                    ):
                        if run.get(attribute):
                            run_row[key] = run.get(attribute)
                    read_lengths = [
                        value
                        for read in run.findall(".//Statistics/Read")
                        if (value := _clean(read.get("average")))
                    ]
                    if read_lengths:
                        run_row["read_lengths"] = read_lengths
                    sra_files = []
                    for source_file in run.findall(".//SRAFile"):
                        item = {
                            "filename": _clean(source_file.get("filename")),
                            "uri": _clean(source_file.get("url")),
                            "md5": _clean(source_file.get("md5")),
                            "bytes": _clean(source_file.get("size")),
                            "semantic_name": _clean(source_file.get("semantic_name")),
                            "supertype": _clean(source_file.get("supertype")),
                        }
                        sra_files.append({k: v for k, v in item.items() if v})
                    if sra_files:
                        run_row["sra_files"] = [*run_row["sra_files"], *sra_files]
                    sample_runs.append({k: v for k, v in run_row.items() if v is not None})

                    steps = [
                        {"kind": "sample", "name": sample_id, "sample_ref": sample_id},
                    ]
                    if protocol_name:
                        steps.append(
                            {"kind": "protocol_application", "protocol_ref": protocol_name}
                        )
                    steps.extend(
                        [
                            {
                                "kind": "assay",
                                "name": experiment_id,
                                "technology_type": {"value": "sequencing assay"},
                            },
                            {"kind": "scan", "name": run_id},
                        ]
                    )
                    assay_paths.append({"steps": steps})

            instrument = promoted.pop("instrument_model", None)
            if instrument:
                platform_id = platform_ids.setdefault(
                    instrument, f"sequencer-{len(platform_ids) + 1}"
                )
                if not any(item["iid"] == platform_id for item in platforms):
                    platforms.append(
                        {
                            "iid": platform_id,
                            "title": instrument,
                            "technology": "high-throughput sequencing",
                            "description": f"INSDC submitted instrument model: {instrument}",
                        }
                    )
            else:
                platform_id = None
            sample_row: dict = {
                "iid": sample_id,
                "accession": _accession_rows(identifiers),
                "status": _status(sample, "SAMPLE"),
                "title": _text(sample, "./TITLE") or _text(biosample, "./Description/Title") or sample_id,
                "type": "SRA",
                "description": _text(sample, "./DESCRIPTION"),
                "channel": [
                    {
                        "organism": (
                            [{"value": organism_name or "", "taxid": taxid}]
                            if organism_name or taxid
                            else []
                        ),
                        "characteristics": characteristics,
                    }
                ],
                "sra_accession": [sample_id],
                "ena_accession": [sample_id] if result.provider == "ena" else [],
                "sra_run": sample_runs,
                **promoted,
            }
            if instrument:
                sample_row["instrument_model"] = {"value": instrument}
            if platform_id:
                sample_row["platform_ref"] = {"ref": platform_id}
            samples.append(sample_row)
            provenance.append(
                {
                    "target": f"sample/{sample_id}",
                    "source": "SAMPLE",
                    "source_accession": sample_id,
                }
            )

        if not graph["experiments"] and not graph["runs"]:
            warnings.append(
                {
                    "code": "metadata_only_study",
                    "message": "The provider study is valid but has no sequencing hierarchy.",
                }
            )

        organization_nodes = list(graph["organizations"])
        for project_node in graph["bioprojects"]:
            organization_nodes.extend(
                project_node.findall(".//Submission/Description/Organization")
            )
        for biosample in graph["biosamples"].values():
            owner = biosample.find("./Owner")
            if owner is not None:
                organization_nodes.append(owner)
        organization_rows = []
        organization_ids: dict[str, str] = {}
        contributor_rows = []
        contact_refs = []
        for organization in organization_nodes:
            name = _text(organization, "./Name")
            if not name:
                continue
            organization_id = organization_ids.get(name)
            if organization_id is None:
                organization_id = f"organization-{len(organization_ids) + 1}"
                organization_ids[name] = organization_id
                organization_row = {"iid": organization_id, "name": name}
                address = _address(organization.find("./Address"))
                if address:
                    organization_row["address"] = address
                organization_rows.append(organization_row)
            for contact in organization.findall(".//Contact"):
                first = _text(contact, "./Name/First")
                middle = _text(contact, "./Name/Middle")
                last = _text(contact, "./Name/Last")
                email = _clean(contact.get("email")) or _text(contact, "./Email")
                contact_address = _address(contact.find("./Address"))
                signature = (first, middle, last, email, organization_id)
                if any(row.get("_signature") == signature for row in contributor_rows):
                    continue
                contact_id = f"contact-{len(contributor_rows) + 1}"
                row = {
                    "iid": contact_id,
                    "person": {
                        key: value
                        for key, value in {"first": first, "middle": middle, "last": last}.items()
                        if value
                    },
                    "email": email,
                    "address": contact_address,
                    "organization_ref": {"ref": organization_id},
                    "_signature": signature,
                }
                contributor_rows.append(row)
                contact_refs.append({"ref": contact_id})
        for row in contributor_rows:
            row.pop("_signature", None)

        source_documents = [
            {
                "kind": document.kind,
                "name": document.name,
                "uri": document.uri,
                "sha256": document.sha256,
                "media_type": document.media_type,
            }
            for document in result.documents
        ]
        series_accessions = _accession_rows(
            [(result.study_accession, "INSDC"), *study_identifiers]
        )
        database_names = sorted(
            {
                row["database"]
                for row in [
                    *series_accessions,
                    *[
                        accession
                        for sample_row in samples
                        for accession in sample_row.get("accession", [])
                    ],
                ]
                if row.get("database")
            }
        )
        series = {
            "iid": result.study_accession,
            "accession": series_accessions,
            "status": _status(study, "STUDY"),
            "title": title or result.study_accession,
            "summary": summary,
            "type": [study_type] if study_type else [],
            "pubmed_id": pubmed_ids,
            "pubmed_publication": publication_rows,
            "web_link": web_links,
            "contact_ref": contact_refs,
            "protocols": protocols,
            "assay_paths": assay_paths,
        }
        provenance.append(
            {
                "target": "series",
                "source": "STUDY",
                "source_accession": result.study_accession,
            }
        )
        extension = {
            "schema_version": "1.0",
            "provider": result.provider,
            "requested_accession": result.requested_accession,
            "resolved_accessions": {
                "study": [result.study_accession],
                "project": [project] if project else [],
            },
            "submission_origin": origin,
            "enrichment": {
                "requested": enrichment,
                "status": "not_requested" if enrichment == "none" else "requested",
                "linked_accessions": {
                    "geo": geo_links,
                    "arrayexpress": ae_links,
                },
            },
            "warnings": warnings,
            "provenance": provenance,
            "conflicts": conflicts,
            "unmapped": unmapped,
        }
        mapping = {
            "miniml_schema_version": "3.0",
            "source": {
                "format": "NCBI SRA" if result.provider == "sra" else "ENA",
                "version": "INSDC SRA 1.5",
                "documents": source_documents,
            },
            "database": [
                {"iid": name, "name": name}
                for name in database_names
            ],
            "organization": organization_rows,
            "contributor": contributor_rows,
            "platform": platforms,
            "sample": samples,
            "series": {key: value for key, value in series.items() if value is not None},
            "extensions": {"insdc": extension},
        }
        return MINiMLCodec().decode(mapping).package


__all__ = ["INSDCStudyParser"]
