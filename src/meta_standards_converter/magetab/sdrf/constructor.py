# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from __future__ import annotations
from dataclasses import dataclass
from meta_standards_converter.magetab.technology import resolve_technology
from meta_standards_converter.magetab.sdrf.handlers.base import _BaseSDRFHandler
from meta_standards_converter.magetab.technology import detect_ae_technology, has_array_files
from meta_standards_converter.magetab.sdrf.handlers.sequencing import _SequencingSDRFHandler
from meta_standards_converter.magetab.sdrf.handlers.sequencing import _BulkSequencingSDRFHandler
from meta_standards_converter.magetab.sdrf.handlers.single_cell import _SingleCellSequencingSDRFHandler
from meta_standards_converter.magetab.sdrf.handlers.single_cell import _DropletSingleCellSequencingSDRFHandler
from meta_standards_converter.magetab.sdrf.handlers.single_cell import _TenXV2DropletSingleCellSequencingSDRFHandler
from meta_standards_converter.magetab.sdrf.handlers.single_cell import _TenXV3DropletSingleCellSequencingSDRFHandler
from meta_standards_converter.magetab.sdrf.handlers.single_cell import _PlateSingleCellSequencingSDRFHandler
from meta_standards_converter.magetab.sdrf.handlers.spatial import _SpatialSequencingSDRFHandler
from meta_standards_converter.magetab.sdrf.handlers.array import _ArraySDRFHandler
from meta_standards_converter.magetab.sdrf.handlers.generic import _GenericSDRFHandler

class SDRFConstructor():
    def _add_sdrf_to_idf(self, idf: list, data: dict) -> list:
        """
        Appends the generated SDRF payload to an IDF row list.
        """
        idf.append(["SDRF File", self._miniml2sdrf(data=data)])
        return idf

    def _miniml2sdrf(self, data: dict, protocol_registry=None, technology_type=None):
        """
        converts miniml json to magetab sdrf.
        """
        if technology_type is None:
            handler, _ = self.create_operation_handler(data, protocol_registry)
        else:
            handler = self.create_handler(data, protocol_registry, technology_type)
        return self.build(handler)

    def create_handler(self, data, protocol_registry=None, technology_type=None):
        tech_type = technology_type or self._detect_sdrf_technology(data=data)
        handler_class = {
            "plate_single_cell_sequencing": _PlateSingleCellSequencingSDRFHandler,
            "droplet_single_cell_sequencing": _DropletSingleCellSequencingSDRFHandler,
            "tenx_v2_droplet_single_cell_sequencing": _TenXV2DropletSingleCellSequencingSDRFHandler,
            "tenx_v3_droplet_single_cell_sequencing": _TenXV3DropletSingleCellSequencingSDRFHandler,
            "single_cell_sequencing": _SingleCellSequencingSDRFHandler,
            "spatial_sequencing": _SpatialSequencingSDRFHandler,
            "bulk_sequencing": _BulkSequencingSDRFHandler,
            "sequencing": _SequencingSDRFHandler,
            "array": _ArraySDRFHandler,
        }.get(tech_type, _GenericSDRFHandler)

        handler = handler_class(parent=self, data=data, protocol_registry=protocol_registry)
        return handler

    def create_operation_handler(self, data, protocol_registry=None, run_evidence=None):
        """Plan scoped dispatch, sharing source identities, registry and audit."""
        operation = _RoutedSDRFHandler(self, data, protocol_registry)
        operation.run_evidence = run_evidence or {}
        scopes = []
        for sample in operation.ordered_samples():
            runs = operation.sra_runs(sample) or [None]
            for channel in operation._channel_records(sample):
                for run in runs:
                    decision = resolve_technology(sample, channel, run, data=data)
                    scopes.append((_PathScope(sample, channel, run), decision))
                    for diagnostic in decision.diagnostics:
                        message = f"Sample {operation.sample_accession(sample)} technology {diagnostic.code}: {', '.join(diagnostic.paths)}"
                        if message not in operation.audit.warnings:
                            operation.audit.warnings.append(message)
        keys = {decision.handler for _, decision in scopes}
        if not keys:
            summary = detect_ae_technology(data)
        elif len(keys) == 1:
            summary = next(iter(keys))
        else:
            summary = 'sequencing' if keys.isdisjoint({'array', 'generic'}) else 'generic'
        # Preserve homogeneous registration/rendering, including array handlers.
        if len(keys) <= 1:
            handler = self.create_handler(data, operation.protocol_registry, summary)
            handler.run_evidence = operation.run_evidence
            handler.audit = operation.audit
            return handler, summary
        for scope, decision in scopes:
            child = self.create_handler(data, operation.protocol_registry, decision.handler)
            child.path_scope = scope
            child.run_evidence = operation.run_evidence
            child.audit = operation.audit
            operation.children.append(child)
        return operation, summary

    def prepare(self, handler, *, canonical_paths=None, preserve_order=False):
        """Register protocols/audit once, then select the final graph before rendering."""
        from .renderer import SDRFRenderer
        paths = handler.prepare_paths()
        if canonical_paths:
            from ..harmonized import bind_sample_groups
            from ..semantics import plan_miniml_assay_documents
            documents = plan_miniml_assay_documents(bind_sample_groups(handler.data, canonical_paths))
        else:
            documents = {"study.sdrf.txt": paths}
        from ..export_files import exclude_derived_paths
        documents = {name: exclude_derived_paths(values) for name, values in documents.items()}
        renderer = SDRFRenderer(preserve_order=preserve_order)
        self.last_sdrf_audit = handler.audit
        return SDRFExportPlan(tuple(SDRFDocumentPlan(name, values, renderer.plan_columns(values), preserve_order)
                                   for name, values in documents.items()), handler.audit)

    def render(self, plan):
        from .renderer import SDRFRenderer
        tables = [SDRFRenderer(preserve_order=document.preserve_order).render_paths(document.columns, document.paths)
                  for document in plan.documents]
        if not tables:
            raise ValueError("SDRF plan contains no documents")
        header = tables[0][0]
        if any(table[0] != header for table in tables[1:]):
            raise ValueError("Multiple SDRF documents with different graph layouts cannot be consolidated safely.")
        return tables[0] if len(tables) == 1 else [header, *(row for table in tables for row in table[1:])]

    def build(self, handler):
        return self.render(self.prepare(handler))

    def _detect_sdrf_technology(self, data: dict) -> str:
        return detect_ae_technology(data)

    def _has_array_files(self, data: dict) -> bool:
        return has_array_files(data)


@dataclass(frozen=True)
class SDRFDocumentPlan:
    name: str
    paths: list
    columns: list
    preserve_order: bool = False


@dataclass(frozen=True)
class SDRFExportPlan:
    documents: tuple[SDRFDocumentPlan, ...]
    audit: object


@dataclass(frozen=True)
class _PathScope:
    sample: dict
    channel: dict
    run: dict | None


class _RoutedSDRFHandler(_BaseSDRFHandler):
    """Collect scoped assay paths before one column plan and one table render."""
    def __init__(self, parent, data, protocol_registry=None):
        super().__init__(parent, data, protocol_registry)
        self.children = []

    def preregister_protocols(self):
        for child in self.children:
            child.preregister_protocols()

    def build_paths(self):
        return [path for child in self.children for path in child.build_paths()]
