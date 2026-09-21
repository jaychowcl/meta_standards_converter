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
import hashlib
import io
import logging
import os
import re
from urllib.parse import urlparse
from meta_standards_converter.magetab.semantics import build_model, validate_model
from meta_standards_converter.magetab.accession_rows import secondary_accession_pairs
from meta_standards_converter.sources.magetab import MAGETabInput
from meta_standards_converter.miniml import MINiMLV1Migrator


logger = logging.getLogger(__name__)

MAGETAB_VERSION = "magetabv1.1"
MAGETAB_SCHEMA_LOCATION = (
    "https://www.ebi.ac.uk/biostudies/misc/MAGE-TABv1.1_2011_07_28.pdf"
)
MINIML_MOLECULES = {
    "genomic DNA",
    "polyA RNA",
    "total RNA",
    "cytoplasmic RNA",
    "nuclear RNA",
    "protein",
    "other",
}
from ..miniml.factor_terms import MINIML_VARIABLE_FACTORS, FACTOR_NORMALIZATION, normalized_factor_category


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
            *[f"Comment[{database}{label}]" for database in ("ArrayExpress", "ENA", "SRA", "INSDC")
              for label in ("ReleaseDate", "LastUpdateDate")],
            "PubMed ID", "Publication DOI", "Publication Author List", "Publication Title",
            "Publication Status", "Status Term Source Ref", "Status Term Accession Number",
            "Publication Status Term Source REF",
            "Publication Status Term Accession Number",
            "Experiment Description", "Protocol Name", "Protocol Type",
            "Protocol Term Source REF", "Protocol Term Accession Number",
            "Protocol Type Term Source REF", "Protocol Type Term Accession Number",
            "Protocol Description", "Protocol Hardware", "Protocol Software", "Protocol Parameters",
            "Protocol Contact", "Protocol Performer",
            "Quality Control Type", "Quality Control Term Source REF",
            "Quality Control Term Accession Number", "Replicate Type",
            "Replicate Term Source REF", "Replicate Term Accession Number",
            "Normalization Type", "Normalization Term Source REF",
            "Normalization Term Accession Number",
            "SDRF File", "Term Source Name", "Term Source File",
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
            "Performer", "Date",
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
        "cdna read", "cdna read offset", "cell barcode offset", "cell barcode read",
        "cell barcode size", "end bias", "input molecule", "library_strand", "primer",
        "sample barcode offset", "sample barcode read", "sample barcode size",
        "single cell isolation", "spike in", "umi barcode offset", "umi barcode read",
        "umi barcode size",
    }

    def parse(self, source: MAGETabInput) -> dict:
        self.warnings = []
        self.file_catalogue = source.file_catalogue
        idf_rows = self._table(source.idf.text, source.idf.name, rectangular=False)
        if not idf_rows:
            raise ValueError(f"MAGE-TAB IDF {source.idf.name} is empty.")
        idf = self._idf_index(idf_rows)
        protocols = self._protocols(idf)
        series = self._series(idf, idf_rows)
        idf_comments = self._idf_comments(idf_rows)
        if idf_comments:
            series["comments"] = idf_comments
        contributors = self._contributors(idf)
        databases = self._databases(idf)
        database_ids = {item["iid"] for item in databases}
        ambiguous_sources = {
            item["name"] for item in databases
            if sum(other["name"] == item["name"] for other in databases) > 1
        }
        for accession in series.get("accession", []):
            database = accession.get("database")
            if database and database not in database_ids and database not in ambiguous_sources:
                databases.append({"iid": database, "name": database})
                database_ids.add(database)
        samples = {}
        platforms = {}
        unmapped_columns = []
        source_sdrfs = []

        for resource in source.sdrfs:
            table = self._table(resource.text, resource.name, rectangular=True)
            if len(table) < 2:
                raise ValueError(f"MAGE-TAB SDRF {resource.name} has no data rows.")
            source_sdrfs.append((resource.name, table))
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
        for platform in platforms.values():
            for accession in platform.get("accession", []):
                database = accession.get("database")
                if database and database not in database_ids and database not in ambiguous_sources:
                    databases.append({"iid": database, "name": database})
                    database_ids.add(database)
        self._warn_ambiguous_term_refs(idf_rows, source_sdrfs, ambiguous_sources)
        sample_values = [state["sample"] for state in samples.values()]
        series["sample_ref"] = [{"ref": sample["iid"]} for sample in sample_values]
        unmapped_rows = [
            {"row_index": index, "label": row[0], "values": row[1:]}
            for index, row in enumerate(idf_rows)
            if row and normalized_label(row[0]) not in self.KNOWN_IDF_LABELS
        ]
        for row in unmapped_rows:
            self._warn(f"Unmapped IDF row {row['label']} is outside the semantic MSC MINiML model.")

        package = {
            "version": MAGETAB_VERSION,
            "schema_location": MAGETAB_SCHEMA_LOCATION,
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
                    "documents": [
                        {
                            "kind": "idf",
                            "name": source.idf.name,
                            "uri": source.idf.origin,
                            "sha256": hashlib.sha256(source.idf.text.encode("utf-8")).hexdigest(),
                            "media_type": "text/tab-separated-values",
                        },
                        *[
                            {
                                "kind": "sdrf",
                                "name": resource.name,
                                "uri": resource.origin,
                                "sha256": hashlib.sha256(resource.text.encode("utf-8")).hexdigest(),
                                "media_type": "text/tab-separated-values",
                            }
                            for resource in source.sdrfs
                        ],
                    ],
                },
                "unmapped_idf_rows": unmapped_rows,
                "unmapped_sdrf_columns": unmapped_columns,
                "warnings": self.warnings,
            },
        }
        package["mage_tab"]["model"] = validate_model(
            build_model(idf_rows=idf_rows, sdrfs=source_sdrfs)
        )
        # Carry the exact SDRF row assignment across the legacy model boundary.
        row_samples = {(name, i): self._row_identity(table[0], row)
                       for name, table in source_sdrfs
                       for i, row in enumerate(table[1:], 1)}
        for path in package["mage_tab"]["model"]["assay_paths"]:
            path["sample_ref"] = row_samples.get((path["sdrf"], path["row_index"]))
            for step in path['steps']:
                if step.get('kind') == 'file':
                    uri = self._catalogue_uri(step.get('value'))
                    if uri:
                        step['link'] = {'value': uri}
        retained_hz = {
            (sample["iid"], i): channel.pop("_msc_hz", [])
            for sample in sample_values for i, channel in enumerate(sample.get("channel", []))
        }
        typed = MINiMLV1Migrator().migrate(package).package
        mapping = typed.to_mapping()
        for occurrence in mapping.get('extensions', {}).get('magetab', {}).get('unbound_annotations', []):
            self._warn(f"Unbound annotation {occurrence.get('header')} at {occurrence.get('sdrf')} "
                       f"row {occurrence.get('row_index')}, column {occurrence.get('column_index')}; preserved without assigning it to another node.")
        from meta_standards_converter.miniml import MINiMLCodec
        for sample in mapping.get("sample", []):
            for i, channel in enumerate(sample.get("channel", [])):
                channel.setdefault("characteristics", []).extend(retained_hz.get((sample["iid"], i), []))
        return MINiMLCodec().decode(mapping).package

    def _table(self, text: str, name: str, rectangular: bool) -> list[list[str]]:
        rows = [row for row in csv.reader(io.StringIO(text), delimiter="\t") if row]
        if rectangular and rows:
            def canonical_header(label):
                match = re.fullmatch(r"\s*(Comment|Characteristics|Factor\s*Value|Parameter\s*Value|Unit)\s*\[(.*?)\](.*)", label, re.I)
                if not match:
                    return label
                kind = {'factorvalue': 'Factor Value', 'parametervalue': 'Parameter Value'}.get(normalized_label(match.group(1)), match.group(1).strip().title())
                return f"{kind}[{match.group(2)}]{match.group(3)}"
            rows[0] = [canonical_header(label) for label in rows[0]]
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

    def _series(self, idf: dict, idf_rows: list) -> dict:
        investigation = self._nonblank(idf, "Investigation Accession")
        arrayexpress = self._nonblank(idf, "Comment[ArrayExpressAccession]")
        pairs = secondary_accession_pairs(idf_rows)
        candidates = [
            *[pair for pair in pairs if pair[0].upper().startswith("GSE")],
            *[(value, self._accession_database(value)) for value in [*investigation, *arrayexpress]],
            *pairs,
        ]
        accessions = []
        seen = set()
        for value, source in candidates:
            database = source or self._accession_database(value)
            key = (value.casefold(), database)
            if key in seen:
                continue
            seen.add(key)
            accessions.append({
                "value": value,
                "database": database,
            })
        if not accessions:
            raise ValueError("MAGE-TAB IDF contains no usable investigation or secondary accession.")

        series = {
            "iid": self._series_iid(
                investigation=investigation,
                arrayexpress=arrayexpress,
                accessions=accessions,
            ),
            "accession": accessions,
        }
        title = self._first(idf, "Investigation Title")
        if title:
            series["title"] = title
        description = self._first(idf, "Experiment Description")
        if description:
            series["summary"] = description
        design = self._values(idf, "Experimental Design")
        design_sources = self._values(idf, "Experimental Design Term Source REF")
        design_accessions = self._values(idf, "Experimental Design Term Accession Number")
        if any(value.strip() for value in design):
            series["type"] = [
                {
                    "value": value.strip(),
                    **({"term_source_ref": design_sources[index].strip()} if index < len(design_sources) and design_sources[index].strip() else {}),
                    **({"term_accession_number": design_accessions[index].strip()} if index < len(design_accessions) and design_accessions[index].strip() else {}),
                }
                for index, value in enumerate(design) if value.strip()
            ]
        factors = self._values(idf, "Experimental Factor Name")
        factor_types = self._values(idf, "Experimental Factor Type")
        factor_sources = self._values(idf, "Experimental Factor Term Source REF")
        factor_accessions = self._values(idf, "Experimental Factor Term Accession Number")
        variables = []
        for index, factor in enumerate(factors):
            if factor.strip():
                variables.append({
                    "factor": self._normalized_variable_factor(factor),
                    "name": factor.strip(),
                    "type": {
                        "value": factor_types[index].strip() if index < len(factor_types) else factor.strip(),
                        **({"term_source_ref": factor_sources[index].strip()} if index < len(factor_sources) and factor_sources[index].strip() else {}),
                        **({"term_accession_number": factor_accessions[index].strip()} if index < len(factor_accessions) and factor_accessions[index].strip() else {}),
                    },
                })
        if variables:
            series["variable"] = variables
        experiment_date = self._first(idf, "Date of Experiment")
        if experiment_date:
            series["experiment_date"] = experiment_date
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
        releases = self._nonblank(idf, "Public Release Date")
        result = [{"submission_date": value} for value in submissions]
        result.extend({"release_date": value, "date_source": "Public Release Date"} for value in releases)
        for database in ("GEO", "ArrayExpress", "ENA", "SRA", "INSDC"):
            for field, label in (("release_date", "ReleaseDate"), ("last_update_date", "LastUpdateDate")):
                result.extend({"database": database, field: value} for value in
                              self._nonblank(idf, f"Comment[{database}{label}]"))
        if len({s["release_date"] for s in result if s.get("release_date")}) > 1:
            self._warn("Distinct release-date statements retained with their source field or repository scope.")
        return result

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
            "role": self._values(idf, "Person Roles"),
            "role_source": self._values(idf, "Person Roles Term Source Ref"),
            "role_accession": self._values(idf, "Person Roles Term Accession Number"),
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
            if values.get("role"):
                contributor["role"] = [{
                    "value": values["role"],
                    **({"term_source_ref": values["role_source"]} if values.get("role_source") else {}),
                    **({"term_accession_number": values["role_accession"]} if values.get("role_accession") else {}),
                }]
            contributors.append(contributor)
        return contributors

    @staticmethod
    def _idf_comments(rows: list[list[str]]) -> list[dict[str, str]]:
        comments = []
        reserved = {
            "secondaryaccession", "secondaryaccessiontermsourceref", "arrayexpressaccession",
            "relatedexperiment", "georeleasedate", "geolastupdatedate", "arrayexpresssubmissiondate",
        }
        reserved.update(f"{database}{label}".casefold()
                        for database in ("ArrayExpress", "ENA", "SRA", "INSDC")
                        for label in ("ReleaseDate", "LastUpdateDate"))
        for row in rows:
            match = re.fullmatch(r"\s*Comment\s*\[(.*)]\s*", row[0], re.I) if row else None
            if match and normalized_label(match.group(1)) not in reserved:
                comments.extend({"name": match.group(1), "value": value.strip()} for value in row[1:] if value.strip())
        return comments

    def _databases(self, idf: dict) -> list[dict]:
        names = self._values(idf, "Term Source Name")
        files = self._values(idf, "Term Source File")
        versions = self._values(idf, "Term Source Version")
        declarations = {}
        for index, name in enumerate(names):
            if not name.strip():
                continue
            item = {"iid": name.strip(), "name": name.strip()}
            if index < len(files) and files[index].strip():
                item["url"] = files[index].strip()
            if index < len(versions) and versions[index].strip():
                item["version"] = versions[index].strip()
            identity = (item["name"], item.get("url"), item.get("version"))
            declarations.setdefault(identity, item)
        groups = {}
        for item in declarations.values():
            groups.setdefault(item["name"], []).append(item)
        reserved = set(groups)
        for name, variants in groups.items():
            if len(variants) == 1:
                continue
            suffix = 1
            for item in variants:
                while f"{name}__msc_{suffix}" in reserved:
                    suffix += 1
                item["iid"] = f"{name}__msc_{suffix}"
                reserved.add(item["iid"])
                suffix += 1
            self._warn(
                f"Conflicting term source {name!r}: "
                f"{[(item['iid'], item.get('url'), item.get('version')) for item in variants]!r}; "
                "bare-name references remain unresolved."
            )
        return list(declarations.values())

    def _warn_ambiguous_term_refs(self, idf_rows, source_sdrfs, ambiguous_sources):
        for index, row in enumerate(idf_rows, start=1):
            if row and "termsourceref" in normalized_label(row[0]):
                for column, value in enumerate(row[1:], start=2):
                    if value.strip() in ambiguous_sources:
                        self._warn(f"Unresolved term source reference {value.strip()!r} at IDF row {index}, column {column} ({row[0]}).")
        for name, table in source_sdrfs:
            for column, label in enumerate(table[0]):
                if "termsourceref" not in normalized_label(label):
                    continue
                for index, row in enumerate(table[1:], start=2):
                    if row[column].strip() in ambiguous_sources:
                        self._warn(f"Unresolved term source reference {row[column].strip()!r} at {name} row {index}, column {column + 1} ({label}).")

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
                self._warn(f"Unmapped SDRF column {label} preserved as an assay-node comment.")

        for row in rows:
            identity = self._row_identity(header, row)
            if not identity:
                continue
            # Legacy sample projections must use the same attachment boundaries
            # as canonical migration. Preserve the original row in the model.
            from .semantics import NODE_HEADERS, FILE_HEADERS
            from ..miniml.cells import has_cell_value
            bound_row, active = list(row), False
            for index, label in enumerate(header):
                if label in NODE_HEADERS or label in FILE_HEADERS or normalized_label(label) == normalized_label('Protocol REF'):
                    active = has_cell_value(row[index])
                elif not active:
                    bound_row[index] = ''
            row = bound_row
            state = samples.setdefault(identity, self._new_sample(identity))
            sample = state["sample"]
            source_name = self._cell(header, row, "Source Name") or identity
            label = self._cell(header, row, "Label") or ""
            channel_marker = self._cell(header, row, "Comment[msc_channel]")
            channel_key = (source_name, label, channel_marker)
            channel = state["channels"].setdefault(channel_key, {"characteristics": []})
            if channel not in sample["channel"]:
                sample["channel"].append(channel)
            if channel_marker:
                channel.setdefault("extensions", {})["msc_channel"] = channel_marker
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
            ("Label", channel, "label"),
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
        self._map_material_types(
            filename=filename,
            header=header,
            row=row,
            sample=sample,
            channel=channel,
        )

    def _map_material_types(self, *, filename, header, row, sample, channel):
        biological = True
        from .semantics import NODE_HEADERS
        for index, label in enumerate(header):
            if label in {"Source Name", "Sample Name"}:
                biological = True
            elif label in NODE_HEADERS:
                biological = False
            if normalized_label(label) != normalized_label("Material Type"):
                continue
            raw_value = row[index]
            value = raw_value.strip()
            if not value:
                continue
            molecule = next(
                (
                    candidate
                    for candidate in MINIML_MOLECULES
                    if candidate.casefold() == value.casefold()
                ),
                None,
            )
            if molecule is not None:
                self._set_scalar(
                    channel,
                    "molecule",
                    molecule,
                    f"{filename} sample {sample['iid']}",
                )
                continue
            if not biological:
                continue
            material = {"tag": "material type", "value": value}
            if material not in channel["characteristics"]:
                channel["characteristics"].append(material)

    def _map_characteristics(self, header, row, channel):
        existing = {
            (str(item.get("tag", "")).casefold(), str(item.get("value", "")))
            for item in channel["characteristics"]
        }
        from .harmonized import read_group
        from .semantics import NODE_HEADERS
        from meta_standards_converter.miniml.harmonization import named_harmonized_rows
        biological = True
        for index, label in enumerate(header):
            if label in {"Source Name", "Sample Name"}:
                biological = True
            elif label in NODE_HEADERS or label == "Protocol REF":
                biological = False
            group = read_group(header, row, index)
            if group is not None:
                prefix, value, _ = group
                if biological and prefix == "characteristics" and value is not None:
                    retained = channel.setdefault("_msc_hz", [])
                    for item in named_harmonized_rows([value]):
                        if item not in retained:
                            retained.append(item)
                continue
            match = re.fullmatch(r"\s*(Characteristics|Factor\s+Value)\s*\[(.*)]\s*", label, re.I)
            if not match or not row[index].strip():
                continue
            if not biological or match.group(1).lower().startswith('factor'):
                continue
            tag = match.group(2).strip()
            value = row[index].strip()
            companions = {}
            companion_target = companions
            for companion_index in range(index + 1, len(header)):
                companion = normalized_label(header[companion_index])
                unit_match = re.fullmatch(r"\s*Unit(?:\[([^]]*)])?\s*", header[companion_index], re.I)
                if unit_match:
                    value_unit = row[companion_index].strip()
                    # An absent unit clears the ontology attachment context.
                    # The canonical semantic path retains these unbound cells.
                    companion_target = companions.setdefault("unit", {"value": value_unit}) if value_unit else {}
                    if value_unit and unit_match.group(1):
                        companions["unit_type"] = unit_match.group(1).strip()
                    continue
                if companion == normalized_label("Term Source REF"):
                    if row[companion_index].strip():
                        companion_target["term_source_ref"] = row[companion_index].strip()
                    continue
                if companion == normalized_label("Term Accession Number"):
                    if row[companion_index].strip():
                        companion_target["term_accession_number"] = row[companion_index].strip()
                    continue
                break
            if tag.casefold() == "organism":
                organisms = channel.setdefault("organism", [])
                record = {"value": value, **companions}
                if record not in organisms:
                    organisms.append(record)
                continue
            key = (tag.casefold(), value)
            if key not in existing:
                channel["characteristics"].append({"tag": tag, "value": value, **companions})
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
            elif "extract" in ptype:
                channel.setdefault("extract_protocol", protocol["description"])

    def _map_platform(self, header, row, sample, platforms):
        reference = (
            self._cell(header, row, "Array Design REF")
            or self._cell(header, row, "Array Design File")
        )
        if not reference:
            return
        technology = self._cell(header, row, "Technology Type")
        from meta_standards_converter.metadata.platforms import platform_namespace
        from meta_standards_converter.miniml.model import TECHNOLOGIES
        declared = None
        for index, label in enumerate(header):
            if normalized_label(label) == "arraydesignref" and index + 1 < len(header):
                if normalized_label(header[index + 1]) == "termsourceref":
                    declared = row[index + 1].strip() or None
                break
        namespace = platform_namespace(reference, declared)
        inferred = platform_namespace(reference)
        if declared and inferred and declared != inferred:
            self._warn(f"Array design {reference}: declared namespace {declared} conflicts with identifier namespace {inferred}; source retained.")
        platform = platforms.setdefault(reference, {
            "iid": reference,
            "accession": [{"value": reference, **({"database": namespace} if namespace else {})}],
        })
        if technology:
            platform.setdefault("technology", technology if technology in TECHNOLOGIES else "other")
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

        file_headers = {'Array Data File', 'Image File', 'Array Data Matrix File',
                        'Derived Array Data File', 'Derived Array Data Matrix File'}
        for index, label in enumerate(header):
            if label not in file_headers or not row[index].strip():
                continue
            value = self._catalogue_uri(row[index]) or row[index].strip()
            for j in range(index + 1, len(header)):
                if header[j] in file_headers or header[j] in {'Source Name', 'Sample Name', 'Assay Name', 'Scan Name', 'Protocol REF'}:
                    break
                if normalized_label(header[j]) == normalized_label('Comment[File URI]') and row[j].strip():
                    value = row[j].strip()
                    break
            key = 'raw_data' if label in {'Array Data File', 'Image File'} else 'supplementary_data'
            item = {'value': value}
            if item not in sample.setdefault(key, []):
                sample[key].append(item)

    def _catalogue_uri(self, value):
        name = str(value or '').strip()
        if not name or urlparse(name).scheme:
            return None
        exact = [item for item in self.file_catalogue if item.get('path') == name]
        matches = exact or [item for item in self.file_catalogue if os.path.basename(item.get('path', '')) == name]
        uris = {item['uri'] for item in matches if item.get('uri')}
        return next(iter(uris)) if len(uris) == 1 else None

    def _set_scalar(self, target, key, value, context):
        value = value.strip()
        if key not in target:
            target[key] = value
        elif target[key] != value:
            self._warn(
                f"conflicting {key} values for {context}; keeping {target[key]!r} and dropping {value!r}."
            )

    def _known_sdrf_header(self, label):
        if label == "Comment[msc_channel]" or re.fullmatch(r"Comment\[hz_[^]]+\]", label):
            return True
        normalized = normalized_label(label)
        if normalized in self.EXACT_SDRF_HEADERS:
            return True
        if re.fullmatch(r"(characteristics|factorvalue|parametervalue)\[.*](?:\(.*\))?", normalized, re.I):
            return True
        if re.fullmatch(r"unit\[.*]", normalized, re.I):
            return True
        match = re.fullmatch(r"comment\[(.*)]", label.strip(), re.I)
        if not match:
            return False
        key = " ".join(match.group(1).split()).casefold()
        return key == "file uri" or key in self.KNOWN_COMMENTS or re.fullmatch(r"read\d+ file", key) is not None

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

    @staticmethod
    def _normalized_variable_factor(value: str) -> str:
        return normalized_factor_category(value)

    def _series_iid(self, investigation, arrayexpress, accessions):
        if arrayexpress:
            return arrayexpress[0]
        for accession in investigation:
            if self._accession_database(accession) == "ArrayExpress":
                return accession
        for accession in accessions:
            if accession.get("database") == "ArrayExpress" and accession.get("value"):
                return accession["value"]
        if investigation:
            return investigation[0]
        return accessions[0]["value"]

    def _warn(self, message):
        if not hasattr(self, "warnings"):
            self.warnings = []
        if message not in self.warnings:
            self.warnings.append(message)
            logger.warning(message)
