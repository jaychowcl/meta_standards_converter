# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Parse MAGE-TAB IDF/SDRF tables into a MINiML-compatible package."""

from __future__ import annotations

import csv
import io
import logging
import os
import re
from urllib.parse import urlparse

from meta_standards_converter.ae_handlers.ae_webfetcher import MAGETabInput


logger = logging.getLogger(__name__)


def normalized_label(value: str) -> str:
    return "".join(str(value).split()).casefold()


class AEParser:
    """Project standard MAGE-TAB metadata into the parsed MINiML shape."""

    KNOWN_IDF_LABELS = {
        normalized_label(value)
        for value in (
            "MAGE-TAB Version", "Investigation Title", "Investigation Accession",
            "Investigation Accession Term Source REF", "Comment[SecondaryAccession]",
            "Comment[SecondaryAccessionTermSourceRef]", "Comment[ArrayExpressAccession]",
            "Comment[RelatedExperiment]", "Experimental Design",
            "Experimental Design Term Source REF", "Experimental Design Term Accession Number",
            "Experimental Factor Name", "Experimental Factor Type",
            "Experimental Factor Term Source REF", "Experimental Factor Term Accession Number",
            "Person Last Name", "Person First Name", "Person Mid Initials", "Person Email",
            "Person Phone", "Person Fax", "Person Address", "Person Affiliation", "Person Roles",
            "Person Roles Term Source Ref", "Person Roles Term Accession Number",
            "Date of Experiment", "Public Release Date", "Comment[GEOReleaseDate]",
            "Comment[GEOLastUpdateDate]", "Comment[ArrayExpressSubmissionDate]",
            "PubMed ID", "Publication DOI", "Publication Author List", "Publication Title",
            "Publication Status", "Status Term Source Ref", "Status Term Accession Number",
            "Publication Status Term Source REF",
            "Publication Status Term Accession Number",
            "Experiment Description", "Protocol Name", "Protocol Type",
            "Protocol Type Term Source REF", "Protocol Type Term Accession Number",
            "Protocol Description", "Protocol Hardware", "Protocol Software", "Protocol Parameters",
            "Protocol Contact", "SDRF File", "Term Source Name", "Term Source File",
            "Term Source Version", "Comment[AEExperiment]", "Comment[AEExperimentType]",
            "Comment[AECurator]", "Comment[SequenceDataURI]", "Comment[AEExpectedClusters]",
            "Comment[AEAdditionalAttributes]", "Comment[AEBatchEffect]",
        )
    }

    EXACT_SDRF_HEADERS = {
        normalized_label(value)
        for value in (
            "Source Name", "Sample Name", "Extract Name", "Labeled Extract Name",
            "Hybridization Name", "Assay Name", "Scan Name", "Normalization Name",
            "Protocol REF", "Provider", "Material Type", "Description", "Label",
            "Array Design REF", "Array Design File", "Technology Type", "Term Source REF",
            "Term Accession Number", "Unit", "Array Data File", "Array Data Matrix File",
            "Derived Array Data File", "Derived Array Data Matrix File",
            "Image File",
        )
    }

    KNOWN_COMMENTS = {
        "sample_title", "sample_description", "sample_source_name", "biosd_sample",
        "library_layout", "library_selection", "library_source", "library_strategy",
        "ena_sample", "ena_experiment", "ena_run", "submitted_file_name", "md5",
        "instrument_model", "fastq_uri", "derived data file", "technical replicate group",
        "library construction", "index", "read_type", "read_index", "cdna read size",
        "read1 file", "read2 file", "read3 file", "read4 file",
    }

    def parse(self, source: MAGETabInput) -> dict:
        self.warnings = []
        idf_rows = self._table(source.idf.text, source.idf.name, rectangular=False)
        if not idf_rows:
            raise ValueError(f"MAGE-TAB IDF {source.idf.name} is empty.")
        idf = self._idf_index(idf_rows)
        protocols = self._protocols(idf)
        series = self._series(idf)
        contributors = self._contributors(idf)
        databases = self._databases(idf)
        samples = {}
        platforms = {}
        unmapped_columns = []

        for resource in source.sdrfs:
            table = self._table(resource.text, resource.name, rectangular=True)
            if len(table) < 2:
                raise ValueError(f"MAGE-TAB SDRF {resource.name} has no data rows.")
            self._map_sdrf(
                resource.name,
                table[0],
                table[1:],
                protocols,
                samples,
                platforms,
                unmapped_columns,
            )

        if not samples:
            raise ValueError("MAGE-TAB SDRF contains no usable Source, Sample, or Assay identity.")
        sample_values = [state["sample"] for state in samples.values()]
        series["sample_ref"] = [{"ref": sample["iid"]} for sample in sample_values]
        unmapped_rows = [
            {"row_index": index, "label": row[0], "values": row[1:]}
            for index, row in enumerate(idf_rows)
            if row and normalized_label(row[0]) not in self.KNOWN_IDF_LABELS
        ]
        for row in unmapped_rows:
            self._warn(f"Unmapped IDF row {row['label']} preserved in mage_tab metadata.")

        package = {
            "version": None,
            "database": databases,
            "organization": [],
            "contributor": contributors,
            "platform": list(platforms.values()),
            "sample": sample_values,
            "series": series,
            "mage_tab": {
                "version": self._first(idf, "MAGE-TAB Version"),
                "source": {
                    "kind": source.source_kind,
                    "value": source.source,
                    "idf": source.idf.name,
                    "sdrf": [resource.name for resource in source.sdrfs],
                },
                "unmapped_idf_rows": unmapped_rows,
                "unmapped_sdrf_columns": unmapped_columns,
                "warnings": self.warnings,
            },
        }
        return package

    def _table(self, text: str, name: str, rectangular: bool) -> list[list[str]]:
        rows = [row for row in csv.reader(io.StringIO(text), delimiter="\t") if row]
        if rectangular and rows:
            width = len(rows[0])
            for index, row in enumerate(rows[1:], start=2):
                if len(row) != width:
                    raise ValueError(
                        f"MAGE-TAB SDRF {name} row {index} has {len(row)} columns; expected {width}."
                    )
        return rows

    def _idf_index(self, rows: list[list[str]]) -> dict[str, list[list[str]]]:
        result = {}
        for row in rows:
            result.setdefault(normalized_label(row[0]), []).append(row[1:])
        return result

    def _values(self, idf: dict, label: str) -> list[str]:
        rows = idf.get(normalized_label(label), [])
        return rows[0] if rows else []

    def _nonblank(self, idf: dict, label: str) -> list[str]:
        return [value.strip() for value in self._values(idf, label) if value.strip()]

    def _first(self, idf: dict, label: str):
        return next(iter(self._nonblank(idf, label)), None)

    def _series(self, idf: dict) -> dict:
        investigation = self._nonblank(idf, "Investigation Accession")
        arrayexpress = self._nonblank(idf, "Comment[ArrayExpressAccession]")
        secondary = self._nonblank(idf, "Comment[SecondaryAccession]")
        secondary_sources = self._values(idf, "Comment[SecondaryAccessionTermSourceRef]")
        source_by_secondary = {
            accession.casefold(): secondary_sources[index].strip()
            for index, accession in enumerate(secondary)
            if index < len(secondary_sources) and secondary_sources[index].strip()
        }
        candidates = [
            *[value for value in secondary if value.upper().startswith("GSE")],
            *investigation,
            *arrayexpress,
            *secondary,
        ]
        accessions = []
        seen = set()
        for value in candidates:
            key = value.casefold()
            if key in seen:
                continue
            seen.add(key)
            accessions.append({
                "value": value,
                "database": source_by_secondary.get(key) or self._accession_database(value),
            })
        if not accessions:
            raise ValueError("MAGE-TAB IDF contains no usable investigation or secondary accession.")

        series = {"accession": accessions}
        title = self._first(idf, "Investigation Title")
        if title:
            series["title"] = title
        description = self._first(idf, "Experiment Description")
        if description:
            series["summary"] = description
        design = self._nonblank(idf, "Experimental Design")
        if design:
            series["type"] = design
        factors = self._values(idf, "Experimental Factor Name")
        factor_types = self._values(idf, "Experimental Factor Type")
        variables = []
        for index, factor in enumerate(factors):
            if factor.strip():
                variables.append({
                    "factor": factor.strip(),
                    "type": factor_types[index].strip() if index < len(factor_types) else factor.strip(),
                })
        if variables:
            series["variable"] = variables
        related = self._nonblank(idf, "Comment[RelatedExperiment]")
        if related:
            series["relation"] = [
                {"type": "related experiment", "target": value} for value in related
            ]
        statuses = self._statuses(idf)
        if statuses:
            series["status"] = statuses
        publications = self._publications(idf)
        if publications:
            series["pubmed_id"] = [item["pubmed_id"] for item in publications if item.get("pubmed_id")]
            series["pubmed_publication"] = publications
        return series

    def _statuses(self, idf: dict) -> list[dict]:
        submissions = self._nonblank(idf, "Date of Experiment")
        releases = self._nonblank(idf, "Comment[GEOReleaseDate]") or self._nonblank(idf, "Public Release Date")
        updates = self._nonblank(idf, "Comment[GEOLastUpdateDate]")
        count = max(len(submissions), len(releases), len(updates), 0)
        return [
            {
                **({"submission_date": submissions[index]} if index < len(submissions) else {}),
                **({"release_date": releases[index]} if index < len(releases) else {}),
                **({"last_update_date": updates[index]} if index < len(updates) else {}),
            }
            for index in range(count)
        ]

    def _publications(self, idf: dict) -> list[dict]:
        fields = {
            "pubmed_id": self._values(idf, "PubMed ID"),
            "doi": self._values(idf, "Publication DOI"),
            "author_list": self._values(idf, "Publication Author List"),
            "title": self._values(idf, "Publication Title"),
            "status": self._values(idf, "Publication Status"),
            "status_term_source_ref": (
                self._values(idf, "Publication Status Term Source REF")
                or self._values(idf, "Status Term Source Ref")
            ),
            "status_term_accession_number": (
                self._values(idf, "Publication Status Term Accession Number")
                or self._values(idf, "Status Term Accession Number")
            ),
        }
        count = max((len(values) for values in fields.values()), default=0)
        return [
            {
                key: values[index].strip() if index < len(values) and values[index].strip() else None
                for key, values in fields.items()
            }
            for index in range(count)
            if any(index < len(values) and values[index].strip() for values in fields.values())
        ]

    def _contributors(self, idf: dict) -> list[dict]:
        fields = {
            "last": self._values(idf, "Person Last Name"),
            "first": self._values(idf, "Person First Name"),
            "middle": self._values(idf, "Person Mid Initials"),
            "email": self._values(idf, "Person Email"),
            "phone": self._values(idf, "Person Phone"),
            "fax": self._values(idf, "Person Fax"),
            "address": self._values(idf, "Person Address"),
            "organization": self._values(idf, "Person Affiliation"),
        }
        count = max((len(values) for values in fields.values()), default=0)
        contributors = []
        for index in range(count):
            values = {
                key: items[index].strip() if index < len(items) and items[index].strip() else None
                for key, items in fields.items()
            }
            if not any(values.values()):
                continue
            contributor = {
                "iid": f"contributor-{index + 1}",
                "person": {
                    key: values[key]
                    for key in ("last", "first", "middle")
                    if values[key]
                },
            }
            for key in ("email", "phone", "fax", "address", "organization"):
                if values[key]:
                    contributor[key] = values[key]
            contributors.append(contributor)
        return contributors

    def _databases(self, idf: dict) -> list[dict]:
        names = self._values(idf, "Term Source Name")
        files = self._values(idf, "Term Source File")
        versions = self._values(idf, "Term Source Version")
        result = []
        for index, name in enumerate(names):
            if not name.strip():
                continue
            item = {"iid": name.strip(), "name": name.strip()}
            if index < len(files) and files[index].strip():
                item["url"] = files[index].strip()
            if index < len(versions) and versions[index].strip():
                item["version"] = versions[index].strip()
            result.append(item)
        return result

    def _protocols(self, idf: dict) -> dict[str, dict]:
        names = self._values(idf, "Protocol Name")
        types = self._values(idf, "Protocol Type")
        descriptions = self._values(idf, "Protocol Description")
        return {
            name.strip(): {
                "type": types[index].strip() if index < len(types) else "",
                "description": descriptions[index].strip() if index < len(descriptions) else "",
            }
            for index, name in enumerate(names)
            if name.strip()
        }

    def _map_sdrf(self, filename, header, rows, protocols, samples, platforms, unmapped):
        normalized = [normalized_label(value) for value in header]
        identity_indexes = [
            index for label in ("Sample Name", "Source Name", "Assay Name")
            for index, value in enumerate(normalized)
            if value == normalized_label(label)
        ]
        if not identity_indexes:
            raise ValueError(f"MAGE-TAB SDRF {filename} has no Source, Sample, or Assay identity column.")
        occurrence = {}
        for index, label in enumerate(header):
            key = normalized[index]
            occurrence[key] = occurrence.get(key, 0) + 1
            if not self._known_sdrf_header(label):
                item = {
                    "file": filename,
                    "column_index": index,
                    "occurrence": occurrence[key],
                    "header": label,
                    "values": [row[index] for row in rows],
                }
                unmapped.append(item)
                self._warn(f"Unmapped SDRF column {label} preserved in mage_tab metadata.")

        for row in rows:
            identity = self._row_identity(header, row)
            if not identity:
                continue
            state = samples.setdefault(identity, self._new_sample(identity))
            sample = state["sample"]
            source_name = self._cell(header, row, "Source Name") or identity
            label = self._cell(header, row, "Label") or ""
            channel_key = (source_name, label)
            channel = state["channels"].setdefault(channel_key, {"characteristics": []})
            if channel not in sample["channel"]:
                sample["channel"].append(channel)
            self._map_sample_scalars(filename, header, row, sample, channel)
            self._map_characteristics(header, row, channel)
            self._map_protocols(header, row, sample, channel, protocols)
            self._map_platform(header, row, sample, platforms)
            self._map_files_and_runs(header, row, sample, state)

    def _new_sample(self, identity: str) -> dict:
        sample = {"iid": identity, "channel": []}
        if re.fullmatch(r"GSM\d+", identity, flags=re.IGNORECASE):
            sample["accession"] = [{"value": identity.upper(), "database": "GEO"}]
        return {"sample": sample, "channels": {}, "runs": {}}

    def _row_identity(self, header, row):
        for label in ("Sample Name", "Source Name", "Assay Name"):
            values = self._cells(header, row, label)
            for value in reversed(values):
                if value.strip():
                    return value.strip()
        return None

    def _map_sample_scalars(self, filename, header, row, sample, channel):
        mappings = (
            ("Comment[Sample_title]", sample, "title"),
            ("Comment[Sample_description]", sample, "description"),
            ("Comment[Sample_source_name]", channel, "source"),
            ("Description", sample, "description"),
            ("Provider", channel, "biomaterial_provider"),
            ("Material Type", channel, "molecule"),
            ("Comment[LIBRARY_LAYOUT]", sample, "library_layout"),
            ("Comment[LIBRARY_SELECTION]", sample, "library_selection"),
            ("Comment[LIBRARY_SOURCE]", sample, "library_source"),
            ("Comment[LIBRARY_STRATEGY]", sample, "library_strategy"),
            ("Comment[INSTRUMENT_MODEL]", sample, "instrument_model"),
        )
        for label, target, key in mappings:
            value = self._cell(header, row, label)
            if value:
                self._set_scalar(target, key, value, f"{filename} sample {sample['iid']}")

    def _map_characteristics(self, header, row, channel):
        existing = {
            (str(item.get("tag", "")).casefold(), str(item.get("value", "")))
            for item in channel["characteristics"]
        }
        for index, label in enumerate(header):
            match = re.fullmatch(r"\s*(Characteristics|Factor\s+Value)\s*\[(.*)]\s*", label, re.I)
            if not match or not row[index].strip():
                continue
            tag = match.group(2).strip()
            value = row[index].strip()
            if tag.casefold() == "organism":
                organisms = channel.setdefault("organism", [])
                record = {"value": value}
                if record not in organisms:
                    organisms.append(record)
                continue
            key = (tag.casefold(), value)
            if key not in existing:
                channel["characteristics"].append({"tag": tag, "value": value})
                existing.add(key)

    def _map_protocols(self, header, row, sample, channel, protocols):
        for reference in self._cells(header, row, "Protocol REF"):
            protocol = protocols.get(reference.strip())
            if not protocol or not protocol["description"]:
                continue
            ptype = normalized_label(protocol["type"])
            if "treatment" in ptype:
                channel.setdefault("treatment_protocol", protocol["description"])
            elif "growth" in ptype:
                channel.setdefault("growth_protocol", protocol["description"])
            elif "label" in ptype:
                channel.setdefault("label_protocol", protocol["description"])
            elif "hybrid" in ptype:
                sample.setdefault("hybridization_protocol", protocol["description"])
            elif "scan" in ptype:
                sample.setdefault("scan_protocol", protocol["description"])
            elif "processing" in ptype or "normalization" in ptype:
                sample.setdefault("data_processing", protocol["description"])
            elif "extract" in ptype or "library" in ptype or "sequenc" in ptype:
                channel.setdefault("extract_protocol", protocol["description"])

    def _map_platform(self, header, row, sample, platforms):
        reference = (
            self._cell(header, row, "Array Design REF")
            or self._cell(header, row, "Array Design File")
        )
        if not reference:
            return
        technology = self._cell(header, row, "Technology Type")
        platform = platforms.setdefault(reference, {
            "iid": reference,
            "accession": [{"value": reference, "database": "ArrayExpress"}],
        })
        if technology:
            platform.setdefault("technology", technology)
        sample.setdefault("platform_ref", {"ref": reference})

    def _map_files_and_runs(self, header, row, sample, state):
        run_id = self._cell(header, row, "Comment[ENA_RUN]")
        fastqs = self._cells(header, row, "Comment[FASTQ_URI]")
        md5s = self._cells(header, row, "Comment[MD5]")
        read_files = [
            row[index].strip()
            for index, label in enumerate(header)
            if re.fullmatch(r"\s*Comment\[read\d+\s+file]\s*", label, re.I) and row[index].strip()
        ]
        if run_id:
            run = state["runs"].get(run_id)
            if run is None:
                run = {"run": run_id, "fastq_files": []}
                state["runs"][run_id] = run
                sample.setdefault("sra_run", []).append(run)
            for label, key in (
                ("Comment[ENA_SAMPLE]", "sample"),
                ("Comment[ENA_EXPERIMENT]", "experiment"),
                ("Comment[SUBMITTED_FILE_NAME]", "submitted_file_name"),
                ("Comment[INSTRUMENT_MODEL]", "instrument_model"),
            ):
                value = self._cell(header, row, label)
                if value:
                    run.setdefault(key, value)
            for index, uri in enumerate(fastqs):
                if not uri.strip() or any(item.get("uri") == uri.strip() for item in run["fastq_files"]):
                    continue
                filename = read_files[index] if index < len(read_files) else os.path.basename(urlparse(uri).path)
                run["fastq_files"].append({
                    "filename": filename or None,
                    "uri": uri.strip(),
                    "md5": md5s[index].strip() if index < len(md5s) and md5s[index].strip() else None,
                })
        elif fastqs:
            raw = sample.setdefault("raw_data", [])
            for uri in fastqs:
                item = {"value": uri.strip()}
                if uri.strip() and item not in raw:
                    raw.append(item)

        for label in ("Array Data File", "Image File"):
            raw = sample.setdefault("raw_data", [])
            for value in self._cells(header, row, label):
                item = {"value": value.strip()}
                if value.strip() and item not in raw:
                    raw.append(item)
        for label in (
            "Array Data Matrix File",
            "Derived Array Data File",
            "Derived Array Data Matrix File",
        ):
            derived = sample.setdefault("supplementary_data", [])
            for value in self._cells(header, row, label):
                item = {"value": value.strip()}
                if value.strip() and item not in derived:
                    derived.append(item)

    def _set_scalar(self, target, key, value, context):
        value = value.strip()
        if key not in target:
            target[key] = value
        elif target[key] != value:
            self._warn(
                f"conflicting {key} values for {context}; keeping {target[key]!r} and dropping {value!r}."
            )

    def _known_sdrf_header(self, label):
        normalized = normalized_label(label)
        if normalized in self.EXACT_SDRF_HEADERS:
            return True
        if re.fullmatch(r"(characteristics|factorvalue|parameterValue)\[.*]", normalized, re.I):
            return True
        match = re.fullmatch(r"comment\[(.*)]", label.strip(), re.I)
        if not match:
            return False
        key = " ".join(match.group(1).split()).casefold()
        return key in self.KNOWN_COMMENTS or re.fullmatch(r"read\d+ file", key) is not None

    def _cells(self, header, row, label):
        target = normalized_label(label)
        return [row[index] for index, value in enumerate(header) if normalized_label(value) == target]

    def _cell(self, header, row, label):
        return next((value.strip() for value in self._cells(header, row, label) if value.strip()), None)

    def _accession_database(self, accession):
        upper = accession.upper()
        if upper.startswith("GSE"):
            return "GEO"
        if upper.startswith("E-"):
            return "ArrayExpress"
        if upper.startswith("ERP"):
            return "ENA"
        if upper.startswith("SRP"):
            return "SRA"
        if upper.startswith("DRP"):
            return "DRA"
        return None

    def _warn(self, message):
        if message not in self.warnings:
            self.warnings.append(message)
            logger.warning(message)
