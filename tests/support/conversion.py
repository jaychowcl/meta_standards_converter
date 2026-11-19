# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Shared conversion builders and deterministic collaborators."""

from __future__ import annotations

from meta_standards_converter.miniml import (
    MINiMLHarmonizationPatch,
    apply_miniml_harmonization_patch,
    miniml_source_fingerprint,
)
from meta_standards_converter.sources.contracts import RequestMetrics
from meta_standards_converter.converters.unified.contracts import ConversionBatchResult, ConversionItemResult
from meta_standards_converter.converters.archive_results import ArchiveImportResult, StudyImportOutcome
from meta_standards_converter.miniml import MINiMLCodec


def converted(converter, source, **request):
    batch = converter.convert(source, **request)
    assert len(batch.items) == 1, batch.to_dict()
    item = batch.items[0]
    assert item.execution == 'succeeded' and item.validation == 'valid', batch.to_dict()
    return item.payload


def package(study: str = "GSE1", sample: str = "GSM1") -> dict:
    return {
        "miniml_schema_version": "3.0",
        "source": {"format": "test"},
        "database": [{"public_id": "GEO", "name": "Gene Expression Omnibus"}],
        "series": {"accession": [{"value": study}]},
        "sample": [
            {
                "iid": sample,
                "title": "sample title",
                "channel": [
                    {
                        "organism": [{"value": "Mus musculus", "taxid": "10090"}],
                        "characteristics": [
                            {"name": "tissue", "value": "brain"},
                        ],
                    }
                ],
            }
        ],
    }


def atlas_v1_dataset(dataset_id: str, metadata: dict) -> dict:
    return {
        "dataset_id": dataset_id,
        "source_repository": "geo",
        "source_ordinal": 0,
        "status": "harmonized",
        "metadata": metadata,
        "publication_ids": [],
        "review": None,
        "harmonization": None,
        "diagnostics": [],
    }


def atlas_v1_payload(datasets: list[dict]) -> dict:
    return {
        "schema_version": "1.0",
        "atlas": {"atlas_id": "atlas-test", "title": "Test", "theme": "test"},
        "run": {
            "run_id": "run-test",
            "created_at": "2026-07-31T00:00:00Z",
            "config": {"queries": [], "metadata_repositories": [], "options": {}},
            "status": "complete",
        },
        "datasets": datasets,
        "publications": [],
        "summary": {
            "dataset_count": len(datasets),
            "publication_count": 0,
            "completed_dataset_count": len(datasets),
            "failed_dataset_count": 0,
        },
    }


def _applied_package() -> dict:
    package = {
        "miniml_schema_version": "3.0",
        "source": {"format": "test"},
        "database": [],
        "organization": [],
        "contributor": [],
        "platform": [],
        "series": {"iid": "GSE1"},
        "sample": [
            {
                "iid": "GSM1",
                "channel": [
                    {
                        "source": {"value": "blood"},
                        "characteristics": [
                            {"name": "disease", "value": "raw case"}
                        ],
                    }
                ],
            }
        ],
    }
    path = "/sample/0/channel/0/characteristics/0"
    patch = MINiMLHarmonizationPatch(
        base_sha256=miniml_source_fingerprint(package),
        adds=(
            {
                "path": path,
                "harmonized_value": {
                    "field": "sample_disease_name",
                    "value": "disease",
                    "term_source_ref": "DOID",
                    "term_accession_number": "DOID:4",
                    "hierarchy_depth": 0,
                },
                "source_evidence": {
                    "target_id": "target-1",
                    "source_field": "disease",
                    "source_label": "raw case",
                    "source_field_path": f"{path}/name",
                    "source_value_path": f"{path}/value",
                    "derivation": "direct",
                    "match_kind": "exact_value",
                },
            },
        ),
    )
    return apply_miniml_harmonization_patch(package, patch)


class FailingSource:
    kind = "geo_accession"

    def __init__(self, error):
        self.error = error
        self.calls = 0

    def probe(self, source, context):
        return str(source).startswith("GSE")

    def validate(self, spec, context):
        pass

    def load(self, spec, context):
        context.service("counted", lambda: self)
        self.calls += 1
        raise self.error

    def metrics(self):
        return RequestMetrics(self.calls, self.calls, self.calls * 0.25)


def complete():
    return ConversionBatchResult([ConversionItemResult('one', 'source', status='complete',
        execution='succeeded', completeness='complete', validation='valid')])


def geo(accession, children=()):
    value = package(accession, 'GSM' + accession[3:])
    value['series']['iid'] = accession
    value['source']['format'] = 'GEO MINiML'
    value['series']['relation'] = [{'type': 'SubSeries', 'target': a} for a in children]
    return value


def native(provider="SRA"):
    value = package("SRP1", "SRS1")
    value["source"]["format"] = provider
    value["series"]["iid"] = "SRP1"
    value["database"][0]["iid"] = "GEO"
    value["series"]["relation"] = [{"type": "GEO", "target": "GSE2"}]
    return value


def services(calls, *, failure=False):
    class Enricher:
        publication_issues = []
        def enrich_selected(self, data, *, publications, run_metadata, session):
            calls.append(("standard", publications, run_metadata))
            if failure:
                raise OSError("unavailable")
            return data
    class Linked:
        def enrich_selected(self, data, **options):
            calls.append(("linked",))
            return data, []
    class Peer:
        def load(self, accession, **options):
            calls.append(("peer", accession))
            p = MINiMLCodec().decode(native("ENA")).package
            return ArchiveImportResult(accession, [StudyImportOutcome("SRP1", "SRP1", "complete", p)])
    return {"enricher": Enricher(), "linked_enricher": Linked(), "ena_source": Peer()}


def events(caplog):
    return [r.msc_event for r in caplog.records if hasattr(r, "msc_event")]
