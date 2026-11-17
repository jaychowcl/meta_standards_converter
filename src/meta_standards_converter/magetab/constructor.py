# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
'''
Constructor class for ae MAGETAB idf and sdrf
'''
from meta_standards_converter.conversion_logging import measured
from meta_standards_converter.magetab.idf import IDFConstructor
from meta_standards_converter.magetab.semantics import overlay_miniml_semantics
from meta_standards_converter.magetab.protocols import ProtocolRegistry
from meta_standards_converter.magetab.technology import detect_ae_technology, has_array_files, series_identity, resolve_technology


PLATFORM_HANDLER_KEYS = (
    "plate_single_cell_sequencing",
    "droplet_single_cell_sequencing",
    "tenx_v2_droplet_single_cell_sequencing",
    "tenx_v3_droplet_single_cell_sequencing",
    "single_cell_sequencing",
    "spatial_sequencing",
    "bulk_sequencing",
    "sequencing",
    "array",
    "generic",
)


def validate_platform_handler(value: str) -> str:
    if value not in PLATFORM_HANDLER_KEYS:
        raise ValueError(
            f"Unsupported platform handler: {value}. "
            f"Choose one of: {', '.join(PLATFORM_HANDLER_KEYS)}"
        )
    return value
from meta_standards_converter.magetab.sdrf.constructor import SDRFConstructor
from meta_standards_converter.miniml import MINiMLCodec, MINiMLPackage


from .preparation import preparation_operation
from ._performance import phase


from dataclasses import dataclass
from .sdrf.constructor import SDRFExportPlan


@dataclass(frozen=True)
class MAGETabExportPlan:
    data: dict
    protocol_registry: ProtocolRegistry
    technology_type: str
    prefix_rows: list
    sdrf: SDRFExportPlan


class AEConstructor:
    """Build MAGE-TAB exclusively from supplied metadata; no provider clients."""
    def __init__(self, idf_constructor=None, sdrf_constructor=None):
        self.idf_constructor = idf_constructor if idf_constructor is not None else IDFConstructor()
        self.sdrf_constructor = sdrf_constructor if sdrf_constructor is not None else SDRFConstructor()

    @preparation_operation
    @measured("projection")
    def miniml2magetab(self, data: MINiMLPackage, platform_handler: str | None = None) -> list:
        return self.render(self.prepare(data, platform_handler=platform_handler))

    def prepare(self, data: MINiMLPackage, *, platform_handler=None) -> MAGETabExportPlan:
        """Validate once and build a private export plan, retaining audit effects."""
        with phase("magetab_model_reconstruction"):
            # Decode validates external/changed typed inputs as well as mappings.
            # Its canonical package needs one detached view, not another complete
            # encode/decode/validation cycle before private export mutations.
            data = MINiMLCodec().decode(data).package.to_mapping()
        with phase("magetab_normalization"):
            from ..miniml.publication_identifiers import clean_publication_identifiers
            clean_publication_identifiers(data)
            from meta_standards_converter.miniml.file_references import clean_native_file_placeholders
            clean_native_file_placeholders(data)
            from meta_standards_converter.miniml.archive_dates import normalize_archive_dates
            normalize_archive_dates(data)
            from meta_standards_converter.miniml.archive_administration import normalize_administration
            normalize_administration(data)
            from ..miniml.archive_factors import restore_native_factors
            restore_native_factors(data)
            from ..miniml.archive_entities import local_platforms, coalesce_organizations, repository_databases
            if data.get("source", {}).get("format") in {"ENA", "SRA"}:
                from ..miniml.archive_entities import declare_ontologies
                declare_ontologies(data)
            local_platforms(data)
            coalesce_organizations(data)
            declarations = repository_databases(data)
            if declarations:
                data.setdefault('extensions', {}).setdefault('insdc', {'version': '2.0', 'records': []})['records'].extend(declarations)
            from ..miniml.archive_paths import complete_native_paths
            from ..miniml.archive_libraries import synchronize_library_facts
            synchronize_library_facts(data)
            complete_native_paths(data)
            from meta_standards_converter.miniml.archive_results import normalize_result_bundles
            normalize_result_bundles(data)
        from meta_standards_converter.magetab.protocol_export import prepare_protocols
        prepare_protocols(data)
        from meta_standards_converter.magetab.native_files import project_native_files
        project_native_files(data)
        with phase("magetab_routing"):
            protocol_registry = ProtocolRegistry(series_accession=self._series_accession(data=data))
            if platform_handler is not None:
                technology_type = validate_platform_handler(platform_handler)
                handler = self.sdrf_constructor.create_handler(data=data,
                    protocol_registry=protocol_registry, technology_type=technology_type)
            else:
                handler, technology_type = self.sdrf_constructor.create_operation_handler(data, protocol_registry)
        sdrf = self.sdrf_constructor.prepare(handler,
            canonical_paths=data.get("series", {}).get("assay_paths", []),
            preserve_order=data.get("source", {}).get("format") in {"ENA", "SRA"})
        prefix_rows = self.idf_constructor.prefix_rows(data)
        return MAGETabExportPlan(data, protocol_registry, technology_type, prefix_rows, sdrf)

    def render(self, plan: MAGETabExportPlan) -> list:
        data = plan.data
        sdrf = self.sdrf_constructor.render(plan.sdrf)
        with phase("magetab_idf_construction"):
            idf = self.idf_constructor.miniml2idf(
                data=data, prefix_rows=plan.prefix_rows,
                protocol_registry=plan.protocol_registry,
                technology_type=plan.technology_type,
            )
        sdrf_index = self._sdrf_row_index(rows=idf)
        if sdrf_index is None:
            raise ValueError("IDF does not contain an SDRF File row.")
        idf[sdrf_index] = ["SDRF File", sdrf, *idf[sdrf_index][2:]]
        return self._finalize(data, idf)

    def _finalize(self, data, idf):
        from .accession_rows import expand_secondary_accession_rows
        result = IDFConstructor()._move_comment_rows_to_bottom(
            expand_secondary_accession_rows(overlay_miniml_semantics(data, idf)))
        if data.get('source', {}).get('format') in ('ENA', 'SRA'):
            from .native_file_selection import prune_empty_file_columns
            for row in result:
                if row and row[0] == 'SDRF File':
                    row[1] = prune_empty_file_columns(row[1])
        # Handlers and semantic overlays may emit references absent from the
        # original IDF. Declare the final emitted set, preserving supplied metadata.
        declarations = IDFConstructor()._idf_term_source(result, data)
        positions = {row[0]: i for i, row in enumerate(result) if row}
        for declaration in declarations:
            if declaration[0] in positions:
                result[positions[declaration[0]]] = declaration
            else:
                result.append(declaration)
        # An authored Factor Value is an explicit variable even when MINiML
        # omitted its IDF declaration. Preserve its name without guessing a type.
        from .validation import normalized
        import re
        factor_rows = {row[0]: row for row in result if row and str(row[0]).startswith("Experimental Factor")}
        names = factor_rows.setdefault("Experimental Factor Name", ["Experimental Factor Name"])
        for row in result:
            if row and row[0] == "SDRF File":
                for header in row[1][0]:
                    match = re.fullmatch(r"Factor Value\[(.*)\]", str(header), re.I)
                    if match and normalized(match[1]) not in {normalized(v) for v in names[1:]}:
                        names.append(match[1])
        if len(names) > 1 and names not in result:
            result.append(names)
        for label in ("Experimental Factor Type", "Experimental Factor Term Source REF", "Experimental Factor Term Accession Number"):
            row = factor_rows.get(label)
            if row is not None:
                row.extend([None] * (len(names) - len(row)))
        from .validation import validate_magetab
        with phase("magetab_validation"):
            return validate_magetab(result)

    def _detect_ae_technology(self, data: dict) -> str:
        return detect_ae_technology(data)

    def _has_array_files(self, data: dict) -> bool:
        return has_array_files(data)

    def _series_accession(self, data: dict):
        return series_identity(data) or "GEO"


    def _sdrf_row_index(self, rows: list):
        for index, row in enumerate(rows):
            if row and str(row[0]).strip().lower() == "sdrf file":
                return index
        return None
