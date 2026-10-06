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
import threading
from types import SimpleNamespace
import xml.etree.ElementTree as ET

from meta_standards_converter.helpers.request_helper import RateLimitedRequester, RequestSettings, NCBIApplicationIdentity
from meta_standards_converter.runtime_contracts import get_resource_profile, require_disk_headroom
from meta_standards_converter.xml_safety import parse_xml, read_limited_response, XMLSizeLimitError


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
                          ('project', r'PRJ(?:NA|EB|DB|DA)\d+'),
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


def publication_ids(records):
    """All explicitly linked PMIDs, including entity-scoped retrieval targets."""
    from .archive_publications import references
    return {r['pubmed_id'] for r in references(records) if r.get('pubmed_id')}


class ArchiveHTTP:
    def __init__(self, service, requester=None, resource_profile='standard', evidence_dir=None):
        self.profile = get_resource_profile(resource_profile)
        self.requester = requester or RateLimitedRequester(service=service,
            settings=RequestSettings.from_resource_profile(self.profile))
        self.identity = NCBIApplicationIdentity()
        self.downloaded_bytes = 0
        self._download_lock = threading.RLock()
        self.evidence_dir = Path(evidence_dir) if evidence_dir else None

    def get(self, url, params=None, fmt='xml'):
        from ..metadata.retrieval_budget import _check_cancelled
        _check_cancelled()
        params = dict(params or {})
        if url.startswith(('https://eutils.ncbi.nlm.nih.gov/', 'https://pmc.ncbi.nlm.nih.gov/tools/idconv/')):
            params.update(self.identity.params())
        with self._download_lock:
            remaining = self.profile.max_aggregate_download_bytes - self.downloaded_bytes
        if remaining <= 0:
            raise ValueError("Archive aggregate download budget exhausted")
        response = self.requester.get(url, params=params, stream=True)
        try:
            response.raise_for_status()
            def counted_chunks(**kwargs):
                for chunk in response.iter_content(**kwargs):
                    _check_cancelled()
                    with self._download_lock:
                        if self.downloaded_bytes + len(chunk) > self.profile.max_aggregate_download_bytes:
                            raise XMLSizeLimitError('Archive aggregate download budget exhausted')
                        self.downloaded_bytes += len(chunk)
                    yield chunk
            bounded = SimpleNamespace(headers=response.headers, iter_content=counted_chunks)
            raw = read_limited_response(bounded, max_bytes=min(self.profile.max_xml_bytes, remaining), chunk_size=65536)
        finally:
            response.close()
        if self.evidence_dir:
            require_disk_headroom(self.evidence_dir, required_bytes=max(1, len(raw)),
                                  headroom_fraction=self.profile.disk_headroom_fraction)
            self.evidence_dir.mkdir(parents=True, exist_ok=True)
            # No credentials, request URLs or field-level provenance in the evidence names.
            digest = hashlib.sha256(raw).hexdigest()
            path = self.evidence_dir / (digest + ('.xml' if fmt == 'xml' else '.json' if fmt == 'json' else '.txt'))
            with self._download_lock:
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
        from ..metadata.retrieval_budget import _BudgetExhausted
        if isinstance(error, _BudgetExhausted): raise
        records.issues.append(f'{label}: {type(error).__name__}')
        return None
