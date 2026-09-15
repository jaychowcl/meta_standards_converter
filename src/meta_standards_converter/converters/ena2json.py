# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Independent ena native metadata conversion workflow."""
import logging
from pathlib import Path
from meta_standards_converter.sources.ena import ENASource
from meta_standards_converter.miniml.ena_parser import ENAParser
from meta_standards_converter.runtime_contracts import get_resource_profile
from .archive_results import ArchiveImportResult, StudyImportOutcome, output_name, publish_json

logger = logging.getLogger(__name__)


class ENA2JSONConverter:
    def __init__(self, source=None, parser=None, linked_enricher=None, peer_converter=None,
                 resource_profile='standard', resource_overrides=None, *, publication_enricher=None):
        self.profile = get_resource_profile(resource_profile, overrides=resource_overrides)
        self.source = source or ENASource(resource_profile=self.profile)
        self.parser = parser or ENAParser()
        self.linked_enricher = linked_enricher
        self.peer_converter = peer_converter
        self.publication_enricher = publication_enricher

    def convert(self, accession, *, out=None, enrich_from_geo_ae=False, include_peer=False,
                report_path=None, evidence_dir=None, overwrite=False, seen_studies=None,
                _resolution=None, _filename_seeds=None) -> ArchiveImportResult:
        result = ArchiveImportResult(str(accession))
        seen = seen_studies if seen_studies is not None else set()
        if hasattr(self.source, 'http'):
            self.source.http.evidence_dir = Path(evidence_dir) if evidence_dir else None
        try:
            resolution = _resolution if _resolution is not None else self.source.resolve(accession)
        except Exception as error:
            result.studies.append(StudyImportOutcome(str(accession), str(accession), 'failed', issues=[type(error).__name__]))
            logger.warning('%s: resolution failed (%s)', accession, type(error).__name__)
            resolution = None
        if resolution is not None:
            if not resolution.studies:
                result.studies.append(StudyImportOutcome(str(accession), str(accession), 'failed', issues=resolution.issues))
                logger.warning('%s: no study resolved (%s)', accession, '; '.join(resolution.issues))
            for seed in resolution.studies:
                key = ('ena', seed.study)
                if key in seen:
                    if resolution.issues:
                        result.studies.append(StudyImportOutcome(seed.study, seed.primary, 'partial', issues=list(resolution.issues)))
                        for issue in resolution.issues:
                            logger.warning('%s: %s', accession, issue)
                    continue
                seen.add(key)
                outcome = StudyImportOutcome(seed.study, seed.primary, 'failed', issues=list(resolution.issues))
                result.studies.append(outcome)
                try:
                    records = self.source.fetch(seed)
                    if not records.xml and not any(records.indexed.values()):
                        raise ValueError('No native metadata records retrieved')
                    package = self.parser.parse(records)
                    outcome.issues.extend(records.issues)
                    from meta_standards_converter.metadata.archive_enrichment import LinkedArchiveEnricher, merge_archive_metadata, linked_accessions
                    if include_peer:
                        try:
                            if self.peer_converter is None:
                                from .sra2json import SRA2JSONConverter
                                peer = SRA2JSONConverter(resource_profile=self.profile)
                            else:
                                peer = self.peer_converter
                            peer_result = peer.convert(seed.study, evidence_dir=evidence_dir)
                            for peer_outcome in peer_result.studies:
                                outcome.issues.extend(peer_outcome.issues)
                            for extra in peer_result.packages:
                                package, issues = merge_archive_metadata(package, extra, prefer=False)
                                outcome.issues.extend(issues)
                        except Exception as error:
                            outcome.issues.append(f'peer enrichment unavailable: {type(error).__name__}')
                    if enrich_from_geo_ae:
                        try:
                            enricher = self.linked_enricher or LinkedArchiveEnricher(resource_profile=self.profile)
                            package, issues = enricher.enrich(package)
                            outcome.issues.extend(issues)
                        except Exception as error:
                            outcome.issues.append(f'GEO/ArrayExpress enrichment unavailable: {type(error).__name__}')
                    else:
                        links = linked_accessions(package)
                        if links:
                            logger.warning('%s: linked GEO/ArrayExpress metadata may be richer; use --enrich-from-geo-ae (%s)', seed.primary, ', '.join(links))
                    try:
                        from meta_standards_converter.metadata.enrichment import MINiMLEnricher
                        publication_enricher = self.publication_enricher or MINiMLEnricher(resource_profile=self.profile)
                        package = publication_enricher.enrich(package)
                        outcome.issues.extend(getattr(publication_enricher, 'publication_issues', []))
                    except Exception as error:
                        outcome.issues.append(f'publication enrichment unavailable: {type(error).__name__}')
                    from ..metadata.archive_diagnostics import consistency_issues
                    outcome.issues.extend(consistency_issues(package.to_mapping()))
                    if out is not None:
                        path = Path(out) / output_name(seed, _filename_seeds or resolution.studies)
                        publish_json(path, package.to_mapping(), overwrite)
                        outcome.output = str(path)
                    outcome.package = package
                    outcome.status = 'partial' if outcome.issues else 'complete'
                except Exception as error:
                    outcome.issues.append(type(error).__name__)
                for issue in outcome.issues:
                    logger.warning('%s: %s', seed.study, issue)
                logger.info('%s: %s', seed.primary, outcome.status)
        if report_path is not None:
            publish_json(report_path, result.to_mapping(), overwrite)
        return result
