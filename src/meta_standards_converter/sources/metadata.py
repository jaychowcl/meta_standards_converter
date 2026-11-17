# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Source loading services used by conversion and linked-evidence preparation.

These services own retrieval plus parsing. They never enrich or publish a package.
"""
from pathlib import Path
from typing import Protocol

from meta_standards_converter.miniml import MINiMLPackage
from meta_standards_converter.runtime_contracts import get_resource_profile
from .contracts import request_metrics


class MetadataSource(Protocol):
    def load(self, source: str, **options) -> list[MINiMLPackage]: ...


class GEOMetadataSource:
    def __init__(self, *, fetcher=None, parser=None, resource_profile="standard"):
        from .geo import GEOWebFetcher
        from ..miniml.geo_parser import GEOParser
        self.profile = get_resource_profile(resource_profile)
        self.fetcher = fetcher if fetcher is not None else GEOWebFetcher(resource_profile=self.profile)
        self.parser = parser if parser is not None else GEOParser(resource_profile=self.profile)

    def load(self, accession, *, remove_empty=True):
        text = self.fetcher.fetch_gse_miniml(gse=accession)
        return self.parser.parse(text, remove_empty=remove_empty)

    def metrics(self):
        return request_metrics(self.fetcher)


class MAGETabMetadataSource:
    def __init__(self, *, fetcher=None, parser=None, resource_profile="standard", allowed_hosts=()):
        from .magetab import AEWebFetcher
        from ..magetab.parser import AEParser
        from ..retrieval import RetrievalPolicy
        self.profile = get_resource_profile(resource_profile)
        self.fetcher = fetcher if fetcher is not None else AEWebFetcher(
            resource_profile=self.profile,
            retrieval_policy=RetrievalPolicy(resource_profile=self.profile,
                allowed_hosts=frozenset(allowed_hosts), allowed_schemes=frozenset({"https"})))
        self.parser = parser if parser is not None else AEParser()

    def load(self, source, *, sdrf_sources=None):
        return [self.parser.parse(self.fetcher.resolve(source, sdrf_sources=sdrf_sources))]

    def metrics(self):
        return request_metrics(self.fetcher)


class NativeMetadataSource:
    """Load native records using the existing core/checkpoint recovery boundary."""
    def __init__(self, provider, *, source=None, parser=None, resource_profile="standard"):
        if provider not in {"sra", "ena"}:
            raise ValueError("native provider must be sra or ena")
        if provider == "sra":
            from .sra import SRASource as Source
            from ..miniml.sra_parser import SRAParser as Parser
        else:
            from .ena import ENASource as Source
            from ..miniml.ena_parser import ENAParser as Parser
        self.provider = provider
        self.source = source if source is not None else Source(resource_profile=resource_profile)
        self.parser = parser if parser is not None else Parser()

    def resolve_arrayexpress(self, accession):
        """Discover and verify alternative loading evidence, without loading it."""
        if self.provider != 'ena':
            raise ValueError('ArrayExpress fallback requires ENA')
        return self.source.resolve_arrayexpress(accession)

    def load(self, accession, *, publications=True, evidence_dir=None):
        from ..converters.archive_results import ArchiveImportResult, StudyImportOutcome, _load_native, native_failure
        from ..metadata.archive_diagnostics import consistency_issues
        result = ArchiveImportResult(str(accession))
        if hasattr(self.source, "http"):
            self.source.http.evidence_dir = Path(evidence_dir) if evidence_dir else None
        try:
            resolution = self.source.resolve(accession)
        except Exception as error:
            result.studies.append(StudyImportOutcome(str(accession), str(accession), "failed",
                issues=[type(error).__name__], error=native_failure([error])))
            return result
        if not resolution.studies:
            result.studies.append(StudyImportOutcome(str(accession), str(accession), "failed",
                issues=list(resolution.issues), error=native_failure(getattr(resolution, "failures", ()))))
        seen = set()
        for seed in resolution.studies:
            if seed.study in seen:
                continue
            seen.add(seed.study)
            outcome = StudyImportOutcome(seed.study, seed.primary, "failed", issues=list(resolution.issues))
            result.studies.append(outcome)
            try:
                package, issues = _load_native(self.source, self.parser, seed, _outcome=outcome, publications=publications)
                outcome.issues.extend(issues)
                outcome.issues.extend(consistency_issues(package.to_mapping()))
                outcome.package = package
                outcome.status = "partial" if outcome.issues else "complete"
            except Exception as error:
                outcome.issues.append(f"native load failed: {type(error).__name__}")
                outcome.error = native_failure([outcome.error, error])
        return result

    def metrics(self):
        return request_metrics(self.source)
