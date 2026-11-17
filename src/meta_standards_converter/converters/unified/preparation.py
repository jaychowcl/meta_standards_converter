# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Prepare loaded groups once using existing enrichment and merge services."""
from meta_standards_converter.conversion_logging import measured
from copy import deepcopy
from dataclasses import replace
import hashlib
import json
import re
from meta_standards_converter.conversion_logging import emit, count, stage, labels

from meta_standards_converter.miniml import MINiMLCodec
from meta_standards_converter.miniml.codec import MINiMLCompatibilityError
from .contracts import Diagnostic


def provider_of(package, fallback=None):
    source = package.to_mapping().get("source", {})
    name = str(source.get("format", "")).lower()
    provider = {"geo miniml": "geo", "geo": "geo", "sra": "sra", "ena": "ena", "arrayexpress": "biostudies",
                "ae": "biostudies", "mage-tab": "biostudies", "magetab": "biostudies",
                "biostudies": "biostudies"}.get(name)
    # Secondary accessions alone do not identify the source repository.
    return provider or fallback


def study_ids(package):
    from meta_standards_converter.metadata.archive_enrichment import entity_ids
    return sorted(entity_ids(package.to_mapping().get("series", {})))


class MetadataPreparation:
    def __init__(self, context):
        self.context = context
        self.cache = {}
        from .families import StudyFamilies
        self.families = StudyFamilies(context)

    @measured("preparation")
    def prepare(self, loaded):
        if loaded.metadata is None:
            return loaded
        providers = {provider_of(p, loaded.provider) for g in loaded.metadata.groups for p in g.packages}
        if len(providers) == 1:
            loaded.provider = providers.pop()
        from meta_standards_converter.metadata.retrieval_budget import _active, _budget_stage, _coverage, _coverage_issues
        with _budget_stage('study_expansion'):
            loaded = self.families.expand(loaded)
        from meta_standards_converter.metadata.geo_publications import _series_accessions
        self.available = {}
        for group in loaded.metadata.groups:
            for package in group.packages:
                if provider_of(package, loaded.provider) == "geo":
                    data = package.to_mapping()
                    self.available.update({key: data for key in _series_accessions(data.get("series", {}))})
        available_key = hashlib.sha256(json.dumps({k: v["series"] for k, v in self.available.items()}, sort_keys=True).encode()).hexdigest()
        groups = []
        for group in loaded.metadata.groups:
            packages = []
            for package in group.packages:
                provider = provider_of(package, loaded.provider)
                key = (hashlib.sha256(json.dumps(package.to_mapping(), sort_keys=True).encode()).hexdigest(),
                       provider, self.context.preparation_policy, available_key)
                if key not in self.cache:
                    prepared_result = self._package(package, provider)
                    if not (_active.get() and _active.get().exhausted): self.cache[key] = prepared_result
                else:
                    count("preparation_cache_hits")
                    emit("preparation_cache_hit", level="DEBUG", provider=provider)
                    prepared_result = self.cache[key]
                prepared, records, diagnostics = prepared_result
                prepared = _coverage(prepared, _active.get())
                packages.append(prepared)
                loaded.preparation.extend(deepcopy(records))
                loaded.diagnostics.extend(diagnostics)
                loaded.diagnostics.extend(Diagnostic('source_partial', issue, 'preparation', 'warning') for issue in _coverage_issues(prepared))
            groups.append(replace(group, packages=tuple(packages)))
        loaded.metadata = replace(loaded.metadata, groups=tuple(groups))
        loaded.enrichment_applied = True
        return loaded

    def _package(self, package, provider):
        context = self.context
        policy = context.preparation_policy
        records, diagnostics = [], []
        ids = study_ids(package)
        identity = package.to_mapping().get("series", {}).get("iid") or (ids[0] if ids else "input")

        def run(name, applicable, fn):
            nonlocal package
            record = {"operation": name, "dataset": identity, "provider": provider,
                      "status": "skipped", "reason": "not_applicable"}
            records.append(record)
            if not policy.allows(name):
                record["reason"] = "disabled"
                diagnostics.append(Diagnostic("enrichment_skipped", f"{identity}: {name} disabled", "preparation", "info"))
                return
            try:
                if not (applicable() if callable(applicable) else applicable):
                    diagnostics.append(Diagnostic("enrichment_skipped", f"{identity}: {name} not applicable", "preparation", "info"))
                    return
                record.update(status="attempted", reason=None)
                from meta_standards_converter.metadata.retrieval_budget import _active, _budget_stage
                budget = _active.get()
                if budget and budget.exhausted:
                    record.update(status='skipped', reason='metadata_budget_exhausted')
                    diagnostics.append(Diagnostic('source_partial', f'{identity}: {name} unfinished: metadata budget exhausted', 'preparation', 'warning'))
                    return
                # A collaborator gets an isolated, valid candidate. Commit only
                # after it passes the same policy used by strict consumers.
                with _budget_stage(name):
                    candidate, issues = fn(MINiMLCodec().decode(package, strict=True).package)
                package = MINiMLCodec().decode(candidate, strict=True).package
                record["status"] = "partial" if issues else "completed"
                diagnostics.extend(Diagnostic("source_partial", f"{identity}: {name}: {issue}", "preparation", "warning") for issue in issues)
                diagnostics.append(Diagnostic("enrichment_completed", f"{identity}: {name} {record['status']}", "preparation", "info"))
            except MINiMLCompatibilityError as exc:
                record.update(status="failed", reason="invalid_candidate")
                diagnostics.extend(Diagnostic(
                    "source_partial", f"{identity}: {name}: {d.code} at {d.path}: {d.message}; candidate rejected",
                    "preparation", "warning") for d in exc.diagnostics)
            except Exception as exc:
                record["status"] = "failed"
                from meta_standards_converter.metadata.retrieval_budget import _BudgetExhausted
                if isinstance(exc, _BudgetExhausted):
                    record.update(status='partial', reason='metadata_budget_exhausted')
                    message = f'{identity}: {name}: {exc}; valid candidate retained'
                else: message = f"{identity}: {name} failed ({type(exc).__name__})"
                diagnostics.append(Diagnostic("source_partial", message, "preparation", "warning"))

        from meta_standards_converter.metadata.reconciliation import inspect_recorded_conflicts
        def has_recorded_evidence():
            kinds = {'MINiML', 'source_conflict'} if policy.enrichment == 'curators' else {'source_conflict'}
            return any(r.get('kind') in kinds for r in package.to_mapping().get('extensions', {}).get('insdc', {}).get('records', []))
        if policy.enrichment == "off" and has_recorded_evidence():
            # Reporting supplied disagreements is not supplementary reconciliation.
            _, issues = inspect_recorded_conflicts(package, "off")
            records.append({"operation": "recorded_evidence", "dataset": identity,
                            "provider": provider, "status": "partial", "reason": None})
            diagnostics.extend(Diagnostic("source_partial", issue, "preparation", "warning") for issue in issues)
        elif has_recorded_evidence():
            run("recorded_evidence", True,
                lambda p: inspect_recorded_conflicts(p, policy.enrichment))
        read_id = next((v for v in ids if re.fullmatch(r"[SED]RP\d+", v)), None)
        archive_id = read_id or next((v for v in ids if re.fullmatch(r"PRJ(?:NA|EB|DB|DA)\d+", v)), None)
        run("peer_provider", provider in {"sra", "ena"} and read_id is not None,
            lambda p: self._peer(p, provider, read_id))
        from meta_standards_converter.metadata.archive_enrichment import linked_accessions
        from .families import recorded_neighbors
        links = set()
        def linked_applicable():
            nonlocal links
            if not (provider in {"geo", "biostudies"} or (provider in {"sra", "ena"} and archive_id is not None)):
                return False
            links = set(linked_accessions(package)) - set(recorded_neighbors(package, provider))
            # Same-provider experiments are distinct studies, not peer evidence.
            if provider == "geo":
                links = {a for a in links if a.startswith("E-")}
            elif provider == "biostudies":
                links = {a for a in links if a.startswith("GSE")}
            return bool(links)
        run("linked_metadata", linked_applicable, lambda p: self._linked(p, provider, links))
        run("parent_publications", provider == "geo", self._parent_publications)
        run("standard", True, lambda p: self._standard(p, provider))
        from meta_standards_converter.metadata.source_text import _text_issues
        diagnostics.extend(Diagnostic("source_partial", issue, "preparation", "warning") for issue in _text_issues(package))
        for record in records:
            count("preparation_" + record["status"])
            emit("preparation_outcome", level="DEBUG", **record)
        return package, records, diagnostics

    def _parent_publications(self, package):
        from meta_standards_converter.metadata.geo_publications import GEOParentPublications
        service = GEOParentPublications(self.context.geo_source(), self.context.session, self.available)
        result, issues = service.inherit(packages=[package], remove_empty=self.context.input_options.get("remove_empty", True))
        return result[0], issues

    def _peer(self, package, provider, accession):
        from meta_standards_converter.sources.metadata import NativeMetadataSource
        from meta_standards_converter.metadata.archive_enrichment import merge_archive_metadata
        peer = "ena" if provider == "sra" else "sra"
        source = self.context.service(peer + "_source", lambda: NativeMetadataSource(peer, resource_profile=self.context.profile))
        result = self.context.session.load(source, accession,
            publications=self.context.preparation_policy.allows("publications"),
            evidence_dir=self.context.input_options.get("evidence_dir"))
        issues = [issue for outcome in result.studies for issue in outcome.issues]
        if not result.packages:
            issues.append("No peer metadata returned")
        for extra in result.packages:
            if extra == package:
                continue
            package, found = merge_archive_metadata(package, extra, prefer=False,
                mode=self.context.preparation_policy.enrichment)
            issues.extend(found)
        return package, issues

    def _linked(self, package, provider, links):
        from meta_standards_converter.metadata.archive_enrichment import LinkedArchiveEnricher
        service = self.context.service("linked_enricher", lambda: LinkedArchiveEnricher(
            resource_profile=self.context.profile, geo_source=self.context.geo_source(),
            ae_source=self.context.ae_source()))
        return service.enrich_selected(package, accessions=links, session=self.context.session,
            mode=self.context.preparation_policy.enrichment)

    def _standard(self, package, provider):
        from meta_standards_converter.metadata.enrichment import MINiMLEnricher
        service = self.context.service("enricher", lambda: MINiMLEnricher(resource_profile=self.context.profile))
        runs = self.context.preparation_policy.enrichment == "standard" and provider not in {"sra", "ena"}
        result = service.enrich_selected(package, publications=True, run_metadata=runs, session=self.context.session)
        issues = list(getattr(service, "publication_issues", ()))
        if getattr(service, "_sra_failures", 0):
            issues.append("Run evidence retrieval was incomplete")
        return result, issues
