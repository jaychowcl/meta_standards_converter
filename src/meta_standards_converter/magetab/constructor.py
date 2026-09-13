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
from meta_standards_converter.magetab.idf import IDFConstructor
from meta_standards_converter.metadata.enrichment import MAGETabEvidenceResolver
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


class AEConstructor:
    def __init__(self, idf_constructor=None, sdrf_constructor=None, *, pubmed_client=None, insdc_client=None, evidence_resolver=None):
        self.evidence = evidence_resolver or MAGETabEvidenceResolver(pubmed_client, insdc_client)
        self.idf_constructor = idf_constructor or IDFConstructor()
        self.sdrf_constructor = sdrf_constructor or SDRFConstructor()

    def miniml2magetab(self, data: MINiMLPackage, platform_handler: str | None = None) -> list:
        """
        converts miniml json to magetab idf. Walks through sections of idf to extract from miniml
        """
        data = MINiMLCodec().encode(MINiMLCodec().decode(data).package)
        from meta_standards_converter.miniml.file_references import clean_native_file_placeholders
        clean_native_file_placeholders(data)
        from meta_standards_converter.miniml.archive_dates import normalize_archive_dates
        normalize_archive_dates(data)
        from meta_standards_converter.miniml.archive_administration import normalize_administration
        normalize_administration(data)
        from ..miniml.archive_paths import complete_native_paths
        complete_native_paths(data)
        from meta_standards_converter.magetab.protocol_export import prepare_protocols
        prepare_protocols(data)
        from meta_standards_converter.magetab.native_files import project_native_files
        project_native_files(data)
        forced = platform_handler is not None
        if forced:
            technology_type = validate_platform_handler(platform_handler)
        else:
            technology_type = self._detect_ae_technology(data=data)
        protocol_registry = ProtocolRegistry(series_accession=self._series_accession(data=data))
        handler = self.sdrf_constructor.create_handler(
            data=data, protocol_registry=protocol_registry, technology_type=technology_type,
        )
        evidence_type = technology_type
        if not forced and technology_type == 'generic' and any(
            resolve_technology(sample, data=data).handler not in {'array', 'generic'}
            for sample in handler.ordered_samples()
        ):
            # A mixed IDF summary must not suppress sequencing source evidence.
            evidence_type = 'sequencing'
        handler.run_evidence = self.evidence.sample_runs(handler, evidence_type)
        if not forced:
            handler, technology_type = self.sdrf_constructor.create_operation_handler(
                data, protocol_registry, handler.run_evidence,
            )
        sdrf = self.sdrf_constructor.build(handler)
        # Keep IDF validation before its publication lookup, as in the old flow.
        prefix_rows = self.idf_constructor.prefix_rows(data)
        publication_details = self.evidence.publications(data)
        idf = self.idf_constructor.miniml2idf(
            data=data, prefix_rows=prefix_rows, publication_details=publication_details,
            protocol_registry=protocol_registry,
            technology_type=technology_type,
        )
        sdrf_index = self._sdrf_row_index(rows=idf)
        if sdrf_index is None:
            raise ValueError("IDF does not contain an SDRF File row.")
        idf[sdrf_index] = ["SDRF File", sdrf, *idf[sdrf_index][2:]]
        return IDFConstructor()._move_comment_rows_to_bottom(overlay_miniml_semantics(data, idf))

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
