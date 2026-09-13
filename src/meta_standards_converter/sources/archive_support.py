# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Transport and result containers; no shared archive discovery workflow."""
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any
import hashlib
import json
import re
import xml.etree.ElementTree as ET

from meta_standards_converter.helpers.request_helper import RateLimitedRequester, RequestSettings, NCBIApplicationIdentity
from meta_standards_converter.runtime_contracts import get_resource_profile, require_disk_headroom
from meta_standards_converter.xml_safety import parse_xml, read_limited_response


@dataclass(frozen=True)
class StudySeed:
    study: str
    primary: str


@dataclass
class Resolution:
    studies: list[StudySeed] = field(default_factory=list)
    issues: list[str] = field(default_factory=list)


@dataclass
class StudyRecords:
    seed: StudySeed
    xml: list[ET.Element] = field(default_factory=list)
    indexed: dict[str, list[dict]] = field(default_factory=dict)
    linked: list[dict] = field(default_factory=list)
    issues: list[str] = field(default_factory=list)


def accession_kind(value):
    value = str(value).strip().upper()
    for kind, pattern in [('study', r'[SED]RP\d+'), ('experiment', r'[SED]RX\d+'),
                          ('run', r'[SED]RR\d+'), ('sample', r'[SED]RS\d+'),
                          ('project', r'PRJ(?:NA|EB|DB)\d+'),
                          ('biosample', r'SAM(?:N|EA|D)\d+')]:
        if re.fullmatch(pattern, value):
            return value, kind
    raise ValueError('Expected an INSDC study, project, sample, experiment or run accession')


def identifier(node):
    if node is None:
        return None
    return node.get('accession') or node.findtext('IDENTIFIERS/PRIMARY_ID')


def chunks(values, size=100):
    for start in range(0, len(values), size):
        yield values[start:start + size]


class ArchiveHTTP:
    def __init__(self, service, requester=None, resource_profile='standard', evidence_dir=None):
        self.profile = get_resource_profile(resource_profile)
        self.requester = requester or RateLimitedRequester(service=service,
            settings=RequestSettings.from_resource_profile(self.profile,
                request_delay=0.5 if service == 'ncbi_eutils' else 1.0))
        self.identity = NCBIApplicationIdentity()
        self.downloaded_bytes = 0
        self.evidence_dir = Path(evidence_dir) if evidence_dir else None

    def get(self, url, params=None, fmt='xml'):
        params = dict(params or {})
        if url.startswith('https://eutils.ncbi.nlm.nih.gov/'):
            params.update(self.identity.params())
        remaining = self.profile.max_aggregate_download_bytes - self.downloaded_bytes
        if remaining <= 0:
            raise ValueError("Archive aggregate download budget exhausted")
        response = self.requester.get(url, params=params, stream=True)
        try:
            response.raise_for_status()
            raw = read_limited_response(response, max_bytes=min(self.profile.max_xml_bytes, remaining))
            self.downloaded_bytes += len(raw)
        finally:
            response.close()
        if self.evidence_dir:
            require_disk_headroom(self.evidence_dir, required_bytes=max(1, len(raw)),
                                  headroom_fraction=self.profile.disk_headroom_fraction)
            self.evidence_dir.mkdir(parents=True, exist_ok=True)
            # No credentials, request URLs or field-level provenance in the evidence names.
            digest = hashlib.sha256(raw).hexdigest()
            path = self.evidence_dir / (digest + ('.xml' if fmt == 'xml' else '.json' if fmt == 'json' else '.txt'))
            if not path.exists():
                path.write_bytes(raw)
        if fmt == 'xml':
            root = parse_xml(raw, max_bytes=self.profile.max_xml_bytes)
            if root.tag == 'ERROR' or root.find('.//ERROR') is not None:
                raise ValueError('Provider returned an XML error')
            return root
        if fmt == 'json':
            result = json.loads(raw)
            if isinstance(result, dict) and result.get('error'):
                raise ValueError('Provider returned a JSON error')
            return result
        return raw.decode('utf-8-sig')


def attempt(records, label, call):
    """Independent retrieval failure, with no request secrets in diagnostics."""
    try:
        return call()
    except Exception as error:
        records.issues.append(f'{label}: {type(error).__name__}')
        return None
