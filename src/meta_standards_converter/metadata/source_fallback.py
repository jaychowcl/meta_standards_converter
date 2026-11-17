# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Portable provenance for incomplete alternative source loading."""
import json
from ..miniml import MINiMLCodec


def annotate(package, requested, mapping, reason):
    data = package.to_mapping()
    evidence = {'requested_accession': requested, 'source': 'ena', 'study': mapping.study,
                'project': mapping.project, 'evidence': list(mapping.evidence),
                'reason': reason, 'complete': False}
    records = data.setdefault('extensions', {}).setdefault('insdc', {'version': '1.0', 'records': []})['records']
    records.append({'provider': 'ena', 'kind': 'source_fallback', 'accession': requested, 'metadata': evidence})
    data['series'].setdefault('comments', []).append({
        'name': 'MSCSourceFallback', 'value': json.dumps(evidence, sort_keys=True, separators=(',', ':'))})
    return MINiMLCodec().decode(data, strict=True).package


def recorded(package):
    """Report persisted incompleteness independently of enrichment settings."""
    for record in package.to_mapping().get('extensions', {}).get('insdc', {}).get('records', []):
        metadata = record.get('metadata')
        if record.get('kind') == 'source_fallback' and isinstance(metadata, dict) and metadata.get('complete') is False:
            yield metadata
