# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Prepare loaded groups once using existing enrichment and merge services."""
from copy import deepcopy
from dataclasses import replace
import hashlib
import json
import re

from meta_standards_converter.miniml import MINiMLCodec
from meta_standards_converter.miniml.codec import MINiMLCompatibilityError
from meta_standards_converter.metadata.preparation_scope import loading_source, convert_source
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

    def prepare(self, loaded):
        if loaded.metadata is None:
            return loaded
        providers = {provider_of(p, loaded.provider) for g in loaded.metadata.groups for p in g.packages}
        if len(providers) == 1:
            loaded.provider = providers.pop()
        loaded = self.families.expand(loaded)
        groups = []
        for group in loaded.metadata.groups:
            packages = []
            for package in group.packages:
                provider = provider_of(package, loaded.provider)
                key = (hashlib.sha256(json.dumps(package.to_mapping(), sort_keys=True).encode()).hexdigest(),
                       provider, self.context.preparation_policy)
                if key not in self.cache:
                    self.cache[key] = self._package(package, provider)
                prepared, records, diagnostics = self.cache[key]
                packages.append(prepared)
                loaded.preparation.extend(deepcopy(records))
                loaded.diagnostics.extend(diagnostics)
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
                # A collaborator gets an isolated, valid candidate. Commit only
                # after it passes the same policy used by strict consumers.
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
                diagnostics.append(Diagnostic("source_partial", f"{identity}: {name} failed ({type(exc).__name__})", "preparation", "warning"))

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
        run("standard", True, lambda p: self._standard(p, provider))
        return package, records, diagnostics

    def _peer(self, package, provider, accession):
        from meta_standards_converter.converters.sra2json import SRA2JSONConverter
        from meta_standards_converter.converters.ena2json import ENA2JSONConverter
        from meta_standards_converter.metadata.archive_enrichment import merge_archive_metadata
        peer = "ena" if provider == "sra" else "sra"
        cls = ENA2JSONConverter if peer == "ena" else SRA2JSONConverter
        primary = self.context.converter.services.get(provider + "2json")
        converter = getattr(primary, "peer_converter", None) or self.context.service(
            peer + "2json", lambda: cls(resource_profile=self.context.profile))
        with loading_source(self.context.preparation_policy):
            result = convert_source(converter, accession, include_peer=False, enrich_from_geo_ae=False,
                                       evidence_dir=self.context.input_options.get("evidence_dir"))
        issues = [issue for outcome in result.studies for issue in outcome.issues]
        if not result.packages:
            issues.append("No peer metadata returned")
        for extra in result.packages:
            if extra == package:
                continue  # An identical peer package adds no evidence or declarations.
            package, found = merge_archive_metadata(package, extra, prefer=False)
            issues.extend(found)
        return package, issues

    def _linked(self, package, provider, links):
        from meta_standards_converter.metadata.archive_enrichment import LinkedArchiveEnricher
        primary = self.context.converter.services.get((provider or "") + "2json")
        service = getattr(primary, "linked_enricher", None) or self.context.service(
            "linked_enricher", lambda: LinkedArchiveEnricher(resource_profile=self.context.profile,
                geo_converter=self.context.converter.services.get("geo2json"),
                ae_converter=self.context.converter.services.get("ae2json")))
        with loading_source(self.context.preparation_policy):
            return service.enrich_selected(package, accessions=links) if hasattr(service, "enrich_selected") else service.enrich(package)

    def _standard(self, package, provider):
        from meta_standards_converter.metadata.enrichment import MINiMLEnricher
        service = self.context.converter.services.get("enricher")
        if service is None:
            primary = self.context.converter.services.get((provider or "") + "2json") or self.context._services.get((provider or "") + "2json")
            if primary is None and provider == "geo":
                primary = self.context.converter.services.get("geo2ae")
            if primary is None and self.context.target == "magetab":
                primary = self.context.converter.services.get("json2ae")
            service = getattr(primary, "enricher", None) or getattr(primary, "publication_enricher", None)
        service = service or self.context.service("enricher", lambda: MINiMLEnricher(resource_profile=self.context.profile))
        runs = self.context.preparation_policy.enrichment == "standard" and provider not in {"sra", "ena"}
        if hasattr(service, "enrich_selected"):
            result = service.enrich_selected(package, publications=True, run_metadata=runs)
        elif self.context.preparation_policy.enrichment == "standard":
            result = service.enrich(package)
        else:
            # A legacy injected enricher cannot promise selective retrieval.
            raise ValueError("Injected enricher must support enrich_selected for curators")
        issues = list(getattr(service, "publication_issues", ()))
        if getattr(service, "_sra_failures", 0):
            issues.append("Run evidence retrieval was incomplete")
        return result, issues
