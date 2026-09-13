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
                 resource_profile='standard', resource_overrides=None):
        self.profile = get_resource_profile(resource_profile, overrides=resource_overrides)
        self.source = source or ENASource(resource_profile=self.profile)
        self.parser = parser or ENAParser()
        self.linked_enricher = linked_enricher
        self.peer_converter = peer_converter

    def convert(self, accession, *, out=None, enrich_from_geo_ae=False, include_peer=False,
                report_path=None, evidence_dir=None, overwrite=False, seen_studies=None):
        result = ArchiveImportResult(str(accession))
        seen = seen_studies if seen_studies is not None else set()
        if hasattr(self.source, 'http'):
            self.source.http.evidence_dir = Path(evidence_dir) if evidence_dir else None
        try:
            resolution = self.source.resolve(accession)
        except Exception as error:
            result.studies.append(StudyImportOutcome(str(accession), str(accession), 'failed', issues=[type(error).__name__]))
            resolution = None
        if resolution is not None:
            if not resolution.studies:
                result.studies.append(StudyImportOutcome(str(accession), str(accession), 'failed', issues=resolution.issues))
            for seed in resolution.studies:
                key = ('ena', seed.study)
                if key in seen:
                    continue
                seen.add(key)
                outcome = StudyImportOutcome(seed.study, seed.primary, 'failed', issues=list(resolution.issues))
                result.studies.append(outcome)
                try:
                    records = self.source.fetch(seed)
                    if not records.xml:
                        raise ValueError('No native metadata records retrieved')
                    package = self.parser.parse(records)
                    outcome.issues.extend(records.issues)
                    from meta_standards_converter.metadata.archive_enrichment import LinkedArchiveEnricher, merge_archive_metadata, linked_accessions
                    if include_peer:
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
                    if enrich_from_geo_ae:
                        enricher = self.linked_enricher or LinkedArchiveEnricher(resource_profile=self.profile)
                        package, issues = enricher.enrich(package)
                        outcome.issues.extend(issues)
                    else:
                        links = linked_accessions(package)
                        if links:
                            logger.warning('%s: linked GEO/ArrayExpress metadata may be richer; use --enrich-from-geo-ae (%s)', seed.primary, ', '.join(links))
                    if out is not None:
                        path = Path(out) / output_name(seed, resolution.studies)
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
