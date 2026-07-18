# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
'''
Constructor class for ae MAGETAB idf
'''
from meta_standards_converter.harmonizers.harmonizers import Harmonizer
from meta_standards_converter.helpers.json_helper import JSONHandler
from meta_standards_converter.pubmed_handlers.pubmed_webfetcher import PubmedWebFetcher

import logging
import re
from datetime import date
from dateutil import parser as date_parser


logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class IDFConstructor():
    def __init__(self, pubmed_fetcher=None):
        self.pubmed_fetcher = pubmed_fetcher or PubmedWebFetcher()

    def _as_list(self, value):
        if value is None:
            return []
        if isinstance(value, list):
            return value
        return [value]

    def miniml2idf(self, data: dict, protocol_registry=None, technology_type=None) -> list:
        """
        converts miniml json to magetab idf. Walks through sections of idf to extract from miniml
        """
        idf = []

        idf.append(["MAGE-TAB Version", "1.1"])
        idf.extend(self._idf_investigations(data=data))
        idf.extend(self._idf_experimental(data=data))
        idf.extend(self._idf_persons(data=data))
        # idf.extend(self._idf_qc_rep_norm(data=data))
        idf.extend(self._idf_dates(data=data))
        idf.extend(self._idf_publications(data=data))
        idf.extend(self._idf_experiments(data=data))
        idf.extend(self._idf_protocols(
            data=data,
            protocol_registry=protocol_registry,
            technology_type=technology_type,
        ))
        idf.append(["SDRF File"])
        idf.extend(self._idf_term_source(magetab=idf, data=data))
        idf.extend(self._idf_platform_specific(data=data, technology_type=technology_type))

        idf = self._move_experiment_description_after_title(rows=idf)
        return self._move_comment_rows_to_bottom(rows=idf)

    def _move_experiment_description_after_title(self, rows: list) -> list:
        title_index = self._row_index(rows=rows, label="Investigation Title")
        description_index = self._row_index(rows=rows, label="Experiment Description")
        if title_index is None or description_index is None:
            return rows

        description_row = rows.pop(description_index)
        if description_index < title_index:
            title_index -= 1
        rows.insert(title_index + 1, description_row)
        return rows

    def _row_index(self, rows: list, label: str):
        for index, row in enumerate(rows):
            if row and row[0] == label:
                return index
        return None

    def _move_comment_rows_to_bottom(self, rows: list) -> list:
        non_comment_rows = []
        comment_rows = []

        for row in rows:
            if self._is_comment_row(row=row):
                comment_rows.append(row)
            else:
                non_comment_rows.append(row)

        return non_comment_rows + comment_rows

    def _is_comment_row(self, row) -> bool:
        return (
            isinstance(row, (list, tuple))
            and bool(row)
            and isinstance(row[0], str)
            and row[0].strip().lower().startswith("comment[")
        )

    def _idf_platform_specific(self, data: dict, technology_type=None) -> list:
        handler_class = self._platform_idf_handler_class(technology_type=technology_type)
        handler = handler_class(parent=self, data=data, technology_type=technology_type)
        return handler.build()

    def _platform_idf_handler_class(self, technology_type=None):
        return {
            "plate_single_cell_sequencing": _PlateSingleCellSequencingPlatformIDFHandler,
            "droplet_single_cell_sequencing": _DropletSingleCellSequencingPlatformIDFHandler,
            "tenx_v2_droplet_single_cell_sequencing": _DropletSingleCellSequencingPlatformIDFHandler,
            "tenx_v3_droplet_single_cell_sequencing": _DropletSingleCellSequencingPlatformIDFHandler,
            "single_cell_sequencing": _SingleCellSequencingPlatformIDFHandler,
            "spatial_sequencing": _SpatialSequencingPlatformIDFHandler,
            "bulk_sequencing": _BulkSequencingPlatformIDFHandler,
            "sequencing": _SequencingPlatformIDFHandler,
            "array": _ArrayPlatformIDFHandler,
        }.get(technology_type, _GenericPlatformIDFHandler)

    def _idf_investigations(self, data: dict) -> list:
        """
        Extracts investigation title, accession, accession term source ref,
        and ArrayExpress accession from MINiML JSON.
        """
        handler = JSONHandler()

        titles = handler._from_path(data, "series.title")
        series_accession_values = handler._from_path(data, "series.accession.*.value")
        secondary_accessions = self._secondary_accession_pairs(data=data)
        accession_values = [accession for accession, _source in secondary_accessions]
        accession_databases = [source for _accession, source in secondary_accessions]
        related_experiments = self._related_experiments(data=data)

        arrayexpress_accessions = self._to_arrayexpress_accessions(series_accession_values)

        rows = [
            ["Investigation Title", *titles],
            ["Comment[SecondaryAccession]", *accession_values],
            ["Comment[SecondaryAccessionTermSourceRef]", *accession_databases],
            ["Comment[ArrayExpressAccession]", *arrayexpress_accessions],
        ]
        if related_experiments:
            rows.append(["Comment[RelatedExperiment]", *related_experiments])
        return rows

    def _secondary_accession_pairs(self, data: dict) -> list[tuple[str, str | None]]:
        pairs = []
        seen = set()

        for series in self._as_list(data.get("series")):
            if not isinstance(series, dict):
                continue
            for accession in self._as_list(series.get("accession")):
                if not isinstance(accession, dict):
                    self._append_secondary_accession(
                        pairs=pairs,
                        seen=seen,
                        accession=accession,
                    )
                    continue
                self._append_secondary_accession(
                    pairs=pairs,
                    seen=seen,
                    accession=accession.get("value"),
                    declared_source=accession.get("database"),
                )

        for sample in self._as_list(data.get("sample")):
            if not isinstance(sample, dict):
                continue
            for accession in self._as_list(sample.get("ena_accession")):
                self._append_secondary_accession(
                    pairs=pairs,
                    seen=seen,
                    accession=accession,
                )

        return pairs

    def _append_secondary_accession(
        self,
        pairs: list[tuple[str, str | None]],
        seen: set[str],
        accession,
        declared_source=None,
    ) -> None:
        accession = self._clean_secondary_accession_value(accession)
        if not accession:
            return

        dedupe_key = accession.upper()
        if dedupe_key in seen:
            return

        seen.add(dedupe_key)
        pairs.append((
            accession,
            self._secondary_accession_source(
                accession=accession,
                declared_source=declared_source,
            ),
        ))

    def _secondary_accession_source(self, accession: str, declared_source=None):
        prefix_sources = {
            "GSE": "GEO",
            "ERP": "ENA",
            "SRP": "SRA",
            "DRP": "DRA",
        }
        upper_accession = accession.upper()
        for prefix, source in prefix_sources.items():
            if upper_accession.startswith(prefix):
                return source
        return self._clean_secondary_accession_value(declared_source)

    def _clean_secondary_accession_value(self, value):
        if value is None:
            return None
        return " ".join(str(value).replace("\t", " ").replace("\n", " ").split()) or None

    def _related_experiments(self, data: dict) -> list:
        related_gses = []
        seen = set()
        for relation in self._as_list((data.get("series") or {}).get("relation")):
            if not isinstance(relation, dict):
                continue
            if not self._is_related_series_relation(relation=relation):
                continue
            for value in (relation.get("type"), relation.get("target"), relation.get("comment")):
                if not isinstance(value, str):
                    continue
                for match in re.findall(r"GSE\d+", value, re.IGNORECASE):
                    normalized = match.upper()
                    if normalized not in seen:
                        seen.add(normalized)
                        related_gses.append(normalized)
        return related_gses

    def _is_related_series_relation(self, relation: dict) -> bool:
        relation_text = " ".join(
            value
            for value in (relation.get("type"), relation.get("target"), relation.get("comment"))
            if isinstance(value, str)
        ).lower()
        return "superseries" in relation_text or "subseries" in relation_text

    def _to_arrayexpress_accessions(self, accession_values: list) -> list:
        return [
            accession.replace("GSE", "E-GEOD-")
            if isinstance(accession, str)
            else None
            for accession in accession_values
        ]

    def _idf_experimental(self, data: dict) -> list:
        """
        extracts experimental design and factors from miniml json to return list of lists as header as first element in lists
        """
        def clean(value):
            if value is None:
                return None
            cleaned = " ".join(str(value).replace("\t", " ").replace("\n", " ").split())
            return cleaned or None

        factors = {}
        factor_order = []
        declared_designs = []
        declared_factors = []

        for series in self._as_list(data.get("series")):
            if not isinstance(series, dict):
                continue
            for design in self._as_list(series.get("type")):
                value = design.get("value") or design.get("name") if isinstance(design, dict) else design
                value = clean(value)
                if value and value not in declared_designs:
                    declared_designs.append(value)
            for variable in self._as_list(series.get("variable")):
                if not isinstance(variable, dict):
                    continue
                name = clean(variable.get("factor") or variable.get("name") or variable.get("tag"))
                factor_type = clean(variable.get("type")) or name
                if name and not any(item[0].lower() == name.lower() for item in declared_factors):
                    declared_factors.append((name, factor_type))

        for sample in self._as_list(data.get("sample")):
            if not isinstance(sample, dict):
                continue

            for channel in self._as_list(sample.get("channel")):
                if not isinstance(channel, dict):
                    continue

                for characteristic in self._as_list(channel.get("characteristics")):
                    if not isinstance(characteristic, dict):
                        continue

                    tag = clean(characteristic.get("tag"))
                    value = clean(characteristic.get("value"))
                    if not tag or not value:
                        continue

                    tag_key = tag.lower()
                    if tag_key not in factors:
                        factors[tag_key] = {
                            "name": tag,
                            "values": set(),
                        }
                        factor_order.append(tag_key)
                    factors[tag_key]["values"].add(value)

        factor_names = (
            [item[0] for item in declared_factors]
            if declared_factors
            else [
                factors[tag_key]["name"]
                for tag_key in factor_order
                if len(factors[tag_key]["values"]) > 1
            ]
        )
        factor_types = [item[1] for item in declared_factors] if declared_factors else factor_names
        blanks = [None for _ in factor_names]

        return [
            ["Experimental Design", *declared_designs],
            ["Experimental Design Term Source REF", *([None] * len(declared_designs))],
            ["Experimental Design Term Accession Number", *([None] * len(declared_designs))],
            ["Experimental Factor Name", *factor_names],
            ["Experimental Factor Type", *factor_types],
            ["Experimental Factor Term Source REF", *blanks],
            ["Experimental Factor Term Accession Number", *blanks],
        ]

    def _idf_persons(self, data: dict) -> list:
        """
        Extracts person information from MINiML JSON using JSONHandler.
        """
        handler = JSONHandler()

        contributors = handler._from_path(data, "contributor.*")
        contributor_count = len(contributors)

        last_names = handler._from_path(data, "contributor.*.person.last")
        first_names = handler._from_path(data, "contributor.*.person.first")
        mid_initials = handler._from_path(data, "contributor.*.person.middle")
        emails = handler._from_path(data, "contributor.*.email")
        phones = handler._from_path(data, "contributor.*.phone")
        faxes = handler._from_path(data, "contributor.*.fax")

        addresses = []
        for contributor in contributors:
            address = contributor.get("address")
            address_parts = handler._flatten_values(address) if address else []
            if address and not isinstance(address, str):
                address_parts = [contributor.get("organization"), *address_parts]
            addresses.append(", ".join(str(x) for x in address_parts if x) or None)

        affiliations = [
            contributor.get("organization")
            for contributor in contributors
        ]

        return [
            ["Person Last Name", *last_names],
            ["Person First Name", *first_names],
            ["Person Mid Initials", *mid_initials],
            ["Person Email", *emails],
            ["Person Phone", *phones],
            ["Person Fax", *faxes],
            ["Person Address", *addresses],
            ["Person Affiliation", *affiliations],
            ["Person Roles", *([None] * contributor_count)],
            ["Person Roles Term Source Ref", *([None] * contributor_count)],
            ["Person Roles Term Accession Number", *([None] * contributor_count)],
        ]

    def _idf_qc_rep_norm(self, data: dict) -> list:
        """
        extracts qc, rep, norm information from miniml json to return list of lists as header as first element in lists
        """
        return [
            ["Quality Control Type"],
            ["Quality Control Term Source REF"],
            ["Quality Control Term Accession Number"],
            ["Replicate Type"],
            ["Replicate Term Source REF"],
            ["Replicate Term Accession Number"],
            ["Normalization Type"],
            ["Normalization Term Source REF"],
            ["Normalization Term Accession Number"],
        ]

    def _idf_dates(self, data: dict) -> list:
        """
        Extracts date-related fields from MINiML JSON using JSONHandler.
        """
        handler = JSONHandler()

        submission_dates = [
            self._normalized_idf_date(value)
            for value in handler._from_path(data, "series.status.*.submission_date")
        ]
        release_dates = [
            self._normalized_idf_date(value)
            for value in handler._from_path(data, "series.status.*.release_date")
        ]
        last_update_dates = [
            self._normalized_idf_date(value)
            for value in handler._from_path(data, "series.status.*.last_update_date")
        ]
        public_release_date = self._earliest_idf_date(values=release_dates)

        return [
            ["Date of Experiment", *submission_dates],
            ["Public Release Date", public_release_date],
            ["Comment[GEOReleaseDate]", *release_dates],
            ["Comment[GEOLastUpdateDate]", *last_update_dates],
            ["Comment[ArrayExpressSubmissionDate]", self._current_idf_date()],
        ]

    def _current_idf_date(self):
        return date.today().isoformat()

    def _earliest_idf_date(self, values: list):
        dates = [
            value
            for value in values
            if isinstance(value, str) and re.match(r"^\d{4}-\d{2}-\d{2}$", value)
        ]
        return min(dates) if dates else None

    def _normalized_idf_date(self, value):
        if value is None or value == "":
            return value

        text = str(value)
        is_iso_date = re.match(r"^\d{4}-\d{2}-\d{2}$", text) is not None
        try:
            parsed = date_parser.parse(text, dayfirst=not is_iso_date)
        except (ValueError, TypeError, OverflowError):
            return value
        return parsed.date().isoformat()

    def _idf_publications(self, data: dict) -> list:
        """
        Extracts publication information from MINiML JSON using JSONHandler.
        """
        handler = JSONHandler()
        enriched_publications = [
            publication
            for publication in handler._from_path(data, "series.pubmed_publication.*")
            if isinstance(publication, dict)
        ]
        if enriched_publications:
            return [
                ["PubMed ID", *[publication.get("pubmed_id") for publication in enriched_publications]],
                ["Publication DOI", *[publication.get("doi") for publication in enriched_publications]],
                ["Publication Author List", *[publication.get("author_list") for publication in enriched_publications]],
                ["Publication Title", *[publication.get("title") for publication in enriched_publications]],
                ["Publication Status", *[publication.get("status") for publication in enriched_publications]],
                ["Status Term Source Ref", *[publication.get("status_term_source_ref") for publication in enriched_publications]],
                ["Status Term Accession Number", *[publication.get("status_term_accession_number") for publication in enriched_publications]],
            ]

        pubmed_ids = [
            pubmed_id
            for pubmed_id in handler._from_path(data, "series.pubmed_id.*")
            if pubmed_id
        ]
        if not pubmed_ids:
            logger.warning("PubMed lookup skipped: no PubMed ID found in series metadata.")

        publication_details = [
            self._lookup_pubmed_id(pubmed_id)
            for pubmed_id in pubmed_ids
        ]

        doi = [d[0] for d in publication_details]
        author_list = [d[1] for d in publication_details]
        title = [d[2] for d in publication_details]
        status = [d[3] for d in publication_details]
        status_term_source_ref = [d[4] for d in publication_details]
        status_term_accession_number = [d[5] for d in publication_details]

        return [
            ["PubMed ID", *pubmed_ids],
            ["Publication DOI", *doi],
            ["Publication Author List", *author_list],
            ["Publication Title", *title],
            ["Publication Status", *status],
            ["Status Term Source Ref", *status_term_source_ref],
            ["Status Term Accession Number", *status_term_accession_number],
        ]

    def _lookup_pubmed_id(self, pubmed_id: str) -> dict:
        '''
        lookup pubmed id to get doi, authorlist, title
        '''
        return self.pubmed_fetcher.pubmed_summary(pubmed_id=pubmed_id)

    def _idf_experiments(self, data: dict) -> list:
        """
        Extracts experiment summary and overall-design from MINiML JSON using JSONHandler.
        """
        handler = JSONHandler()

        summaries = handler._from_path(data, "series.summary")
        designs = handler._from_path(data, "series.overall_design")

        # Combine summaries and designs first
        combined = []
        for s, d in zip(summaries, designs):
            s_clean = (s or "").replace(chr(9), " ")
            d_clean = (d or "").replace(chr(9), " ")

            if d_clean:
                combined.append(f"{s_clean}. {d_clean}".strip())
            else:
                combined.append(s_clean.strip())

        experiment_description = [
            "Experiment Description",
            *combined,
        ]

        return [experiment_description]

    def _idf_protocols(self, data: dict, protocol_registry=None, technology_type=None) -> list:
        """
        extracts protocol information from miniml json to return list of lists as header as first element in lists
        """
        if protocol_registry is not None:
            return self._idf_protocols_from_registry(
                protocol_registry=protocol_registry,
                technology_type=technology_type,
            )

        handler = JSONHandler()

        p_type = []
        term_source_ref = []
        term_accession_number = []
        description = []

        protocol_type_paths = {
            "Manufacture-Protocol": "platform.*.manufacture_protocol",
            "Treatment-Protocol": "sample.*.channel.*.treatment_protocol",
            "Growth-Protocol": "sample.*.channel.*.growth_protocol",
            "Extract-Protocol": "sample.*.channel.*.extract_protocol",
            "Label-Protocol": "sample.*.channel.*.label_protocol",
            "Hybridization-Protocol": "sample.*.hybridization_protocol",
            "Scan-Protocol": "sample.*.scan_protocol",
            "Data-Processing": "sample.*.data_processing",
        }

        for protocol_type, path in protocol_type_paths.items():
            protocol_descriptions = list({x for x in handler._from_path(data, path) if x is not None})

            for desc in protocol_descriptions:
                if desc:
                    onto_type, onto_source_ref, onto_accession = Harmonizer().geoprotocols2efo(protocol_type=protocol_type)
                    p_type.append(onto_type)
                    term_source_ref.append(onto_source_ref)
                    term_accession_number.append(onto_accession)
                    description.append(desc.replace("\t", " "))
                else:
                    p_type.append(protocol_type)
                    term_source_ref.append(None)
                    term_accession_number.append(None)
                    description.append(None)

        protocol_count = len(p_type)

        series_accessions = handler._from_path(data, "series.accession.*.value")
        series_accession = next((x for x in series_accessions if x), "GEO")
        name = ["P-" + series_accession + "-" + str(i+1) for i in range(protocol_count)]

        parameters = [None for _ in range(protocol_count)] # Blank for internal curation
        hardware = [None for _ in range(protocol_count)]# Blank for internal curation
        software = [None for _ in range(protocol_count)]# Blank for internal curation
        contact = [None for _ in range(protocol_count)]# Blank for internal curation

        rows = [
            ["Protocol Name", *name],
            ["Protocol Type", *p_type],
            ["Protocol Type Term Source REF", *term_source_ref],
            ["Protocol Type Term Accession Number", *term_accession_number],
            ["Protocol Description", *description],
            # ["Protocol Parameters", *parameters],
            ["Protocol Hardware", *hardware],
            ["Protocol Software", *software],
            # ["Protocol Contact", *contact],
        ]
        self._append_required_protocol_rows(
            rows=rows,
            series_accession=series_accession,
            technology_type=technology_type,
        )
        return rows

    def _idf_protocols_from_registry(self, protocol_registry, technology_type=None) -> list:
        records = protocol_registry.records()

        names = []
        p_type = []
        term_source_ref = []
        term_accession_number = []
        description = []

        for record in records:
            onto_type, onto_source_ref, onto_accession = Harmonizer().geoprotocols2efo(
                protocol_type=record["label"]
            )
            names.append(record["ref"])
            p_type.append(onto_type)
            term_source_ref.append(onto_source_ref)
            term_accession_number.append(onto_accession)
            description.append(record["text"])

        protocol_count = len(records)
        blanks = [None for _ in range(protocol_count)]

        rows = [
            ["Protocol Name", *names],
            ["Protocol Type", *p_type],
            ["Protocol Type Term Source REF", *term_source_ref],
            ["Protocol Type Term Accession Number", *term_accession_number],
            ["Protocol Description", *description],
            # ["Protocol Parameters", *blanks],
            ["Protocol Hardware", *blanks],
            ["Protocol Software", *blanks],
            # ["Protocol Contact", *blanks],
        ]
        self._append_required_protocol_rows(
            rows=rows,
            series_accession=protocol_registry.series_accession,
            technology_type=technology_type,
        )
        return rows

    def _append_required_protocol_rows(self, rows: list, series_accession: str, technology_type=None) -> None:
        required_labels = ["Sample-Collection-Protocol"]
        if self._is_sequencing_technology(technology_type=technology_type):
            required_labels.append("Nucleic-Acid-Sequencing-Protocol")

        protocol_name_row = self._row(rows=rows, label="Protocol Name")
        protocol_type_row = self._row(rows=rows, label="Protocol Type")
        existing_types = set(protocol_type_row[1:])

        for required_label in required_labels:
            onto_type, onto_source_ref, onto_accession = Harmonizer().geoprotocols2efo(
                protocol_type=required_label
            )
            if onto_type in existing_types:
                continue

            next_index = len(protocol_name_row)
            protocol_name_row.append(f"P-{series_accession}-{next_index}")
            protocol_type_row.append(onto_type)
            self._row(rows=rows, label="Protocol Type Term Source REF").append(onto_source_ref)
            self._row(rows=rows, label="Protocol Type Term Accession Number").append(onto_accession)
            self._row(rows=rows, label="Protocol Description").append(None)
            self._append_if_row_exists(rows=rows, label="Protocol Parameters", value=None)
            self._append_if_row_exists(rows=rows, label="Protocol Hardware", value=None)
            self._append_if_row_exists(rows=rows, label="Protocol Software", value=None)
            self._append_if_row_exists(rows=rows, label="Protocol Contact", value=None)
            existing_types.add(onto_type)

    def _append_if_row_exists(self, rows: list, label: str, value) -> None:
        row = self._optional_row(rows=rows, label=label)
        if row is not None:
            row.append(value)

    def _optional_row(self, rows: list, label: str):
        for row in rows:
            if row and row[0] == label:
                return row
        return None

    def _row(self, rows: list, label: str) -> list:
        for row in rows:
            if row and row[0] == label:
                return row
        raise ValueError(f"IDF protocol rows missing {label}")

    def _is_sequencing_technology(self, technology_type=None) -> bool:
        return technology_type in {
            "sequencing",
            "bulk_sequencing",
            "plate_single_cell_sequencing",
            "single_cell_sequencing",
            "droplet_single_cell_sequencing",
            "tenx_v2_droplet_single_cell_sequencing",
            "tenx_v3_droplet_single_cell_sequencing",
            "spatial_sequencing",
        }

    def _idf_term_source(self, magetab: list, data: dict | None = None) -> list:
        '''
        Return term source rows, preferring metadata supplied by the input package.
        '''
        sources = set()
        for row in magetab:
            if row and "source ref" in str(row[0]).lower():
                sources.update(x for x in row[1:] if x)

        databases = {}
        raw_databases = (data or {}).get("database", [])
        if isinstance(raw_databases, dict):
            raw_databases = [raw_databases]
        for database in raw_databases:
            if not isinstance(database, dict):
                continue
            declared_name = database.get("iid") or database.get("name")
            if declared_name:
                sources.add(declared_name)
            for key in (database.get("iid"), database.get("name"), declared_name):
                if key:
                    databases[str(key).casefold()] = database

        sources = sorted(sources)
        harmonized = Harmonizer().ontologies
        files = []
        versions = []
        for source in sources:
            supplied = databases.get(str(source).casefold())
            if supplied is not None:
                files.append(supplied.get("url") or supplied.get("web_link"))
                versions.append(supplied.get("version"))
                continue
            fallback = harmonized.get(source, {})
            files.append(fallback.get("Term Source File"))
            versions.append(fallback.get("Term Source Version"))

        return [
            ["Term Source Name", *sources],
            ["Term Source File", *files],
            ["Term Source Version", *versions],
        ]


class _BasePlatformIDFHandler:
    def __init__(self, parent: IDFConstructor, data: dict, technology_type=None):
        self.parent = parent
        self.data = data
        self.technology_type = technology_type

    def build(self) -> list:
        return []

    def _as_list(self, value):
        if value is None:
            return []
        if isinstance(value, list):
            return value
        return [value]

    def clean(self, value):
        if value is None:
            return None
        return " ".join(str(value).replace("\t", " ").replace("\n", " ").split()) or None


class _SequencingPlatformIDFHandler(_BasePlatformIDFHandler):
    ENA_DATA_VIEW_BASE = "http://www.ebi.ac.uk/ena/data/view/"
    EMPTY_COMMENT_ROWS = [
        ["Comment[AEExperiment]"],
        ["Comment[AEExperimentType]"],
        ["Comment[AECurator]"],
    ]

    def build(self) -> list:
        rows = [list(row) for row in self.EMPTY_COMMENT_ROWS]
        urls = self.sequence_data_uris()
        if urls:
            rows.append(["Comment[SequenceDataURI]", *urls])
        else:
            logger.warning("SequenceDataURI row skipped: no valid sample.sra_run run accessions found.")
        return rows

    def _set_ae_experiment_type(self, rows: list, value: str) -> None:
        for row in rows:
            if row and row[0] == "Comment[AEExperimentType]":
                row[:] = ["Comment[AEExperimentType]", value]
                return
        rows.append(["Comment[AEExperimentType]", value])

    def secondary_accessions(self) -> list:
        accessions = []
        seen = set()
        for sample in self._as_list(self.data.get("sample")):
            if not isinstance(sample, dict):
                continue
            for accession in self._as_list(sample.get("ena_accession")):
                accession = self.clean(accession)
                if accession and accession not in seen:
                    seen.add(accession)
                    accessions.append(accession)
        return accessions

    def sequence_data_uris(self) -> list:
        accessions = self.run_accessions()
        grouped = {}
        for accession in accessions:
            parsed = self._parse_run_accession(accession=accession)
            if not parsed:
                logger.warning("SequenceDataURI skipped malformed run accession: %s", accession)
                continue
            prefix, number = parsed
            grouped.setdefault(prefix, []).append((number, accession))

        urls = []
        for prefix, values in grouped.items():
            values = sorted(values, key=lambda item: item[0])
            first = values[0][1]
            last = values[-1][1]
            target = first if first == last else f"{first}-{last}"
            urls.append(f"{self.ENA_DATA_VIEW_BASE}{target}")
        return urls

    def run_accessions(self) -> list:
        accessions = []
        seen = set()
        for sample in self._as_list(self.data.get("sample")):
            if not isinstance(sample, dict):
                continue
            for run in self._as_list(sample.get("sra_run")):
                if not isinstance(run, dict):
                    continue
                accession = self.clean(run.get("run"))
                if not accession:
                    logger.warning("SRA run accession missing in sample.sra_run entry.")
                    continue
                accession = accession.upper()
                if accession not in seen:
                    seen.add(accession)
                    accessions.append(accession)
        return accessions

    def _parse_run_accession(self, accession: str):
        match = re.match(r"^([A-Z]+)(\d+)$", accession)
        if not match:
            return None
        return match.group(1), int(match.group(2))


class _BulkSequencingPlatformIDFHandler(_SequencingPlatformIDFHandler):
    pass


class _SingleCellSequencingPlatformIDFHandler(_SequencingPlatformIDFHandler):
    AE_EXPERIMENT_TYPE = "RNA-seq of coding RNA from single cells"

    def build(self) -> list:
        rows = super().build()
        self._set_ae_experiment_type(rows=rows, value=self.AE_EXPERIMENT_TYPE)
        return rows


class _DropletSingleCellSequencingPlatformIDFHandler(_SingleCellSequencingPlatformIDFHandler):
    DROPLET_EMPTY_COMMENT_ROWS = [
        ["Comment[AEExpectedClusters]"],
        ["Comment[AEAdditionalAttributes]"],
        ["Comment[AEBatchEffect]"],
    ]

    def build(self) -> list:
        rows = super().build()
        rows.extend(list(row) for row in self.DROPLET_EMPTY_COMMENT_ROWS)
        return rows


class _PlateSingleCellSequencingPlatformIDFHandler(_BulkSequencingPlatformIDFHandler):
    AE_EXPERIMENT_TYPE = "RNA-seq of coding RNA from single cells"

    def build(self) -> list:
        rows = super().build()
        self._set_ae_experiment_type(rows=rows, value=self.AE_EXPERIMENT_TYPE)
        return rows


class _SpatialSequencingPlatformIDFHandler(_SingleCellSequencingPlatformIDFHandler):
    def build(self) -> list:
        return _SequencingPlatformIDFHandler.build(self)


class _ArrayPlatformIDFHandler(_BasePlatformIDFHandler):
    pass


class _GenericPlatformIDFHandler(_BasePlatformIDFHandler):
    pass
