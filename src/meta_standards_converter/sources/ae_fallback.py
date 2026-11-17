# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Verify public ENA evidence before substituting an unavailable AE source."""
from dataclasses import dataclass
import re
from urllib.parse import urlsplit
from xml.etree import ElementTree as ET

from ..provider_issues import provider_boundary
from .errors import MetadataSourceError


@dataclass(frozen=True)
class VerifiedENAStudy:
    study: str
    project: str | None
    evidence: tuple[str, ...]


def _identifiers(node):
    return {value for value in [node.get('accession'),
        *(n.text for n in node.findall('IDENTIFIERS/PRIMARY_ID')),
        *(n.text for n in node.findall('IDENTIFIERS/SECONDARY_ID')),
        *(n.text for n in node.findall('IDENTIFIERS/EXTERNAL_ID') if n.get('namespace') == 'BioProject')]
        if value and re.fullmatch(r'(?:[SED]RP\d+|PRJ[A-Z]+\d+)', value)}


def _links(node, accession):
    path = f'{node.tag}_LINKS/{node.tag}_LINK/'
    for link in node.findall(path + 'URL_LINK'):
        literal = link.findtext('URL', '')
        url = urlsplit(literal)
        if (url.scheme in {'http', 'https'} and url.hostname in {'www.ebi.ac.uk', 'ebi.ac.uk'}
                and not url.username and not url.password and not url.query and not url.fragment
                and url.path.rstrip('/') in {
                    f'/biostudies/arrayexpress/studies/{accession}',
                    f'/arrayexpress/experiments/{accession}'}):
            yield literal
    for link in node.findall(path + 'XREF_LINK'):
        if link.findtext('DB', '').casefold() in {'arrayexpress', 'biostudies'} and link.findtext('ID') == accession:
            yield f'https://www.ebi.ac.uk/biostudies/arrayexpress/studies/{accession}'


def _response(call, accession):
    from ..xml_safety import XMLSafetyError
    try:
        return call()
    except (MetadataSourceError, XMLSafetyError):
        raise
    except (ValueError, ET.ParseError) as error:
        raise MetadataSourceError('provider_response_invalid', accession, provider='ena') from error


@provider_boundary(provider='ena', operation='source_fallback_verification')
def resolve_arrayexpress(source, accession):
    """An alias query discovers candidates; only exact public record links verify them.

    All candidates must be inspected before selecting one. Transport/budget errors
    propagate, so a partially inspected candidate set cannot become a unique match.
    """
    if not re.fullmatch(r'E-[A-Z]+-\d+', accession):
        raise MetadataSourceError('fallback_unverified', accession, provider='ena')
    rows = _response(lambda: source.search('study', f'study_alias="{accession}"',
                         'study_accession,secondary_study_accession,study_alias'), accession)
    if not isinstance(rows, list):
        raise MetadataSourceError('provider_response_invalid', accession, provider='ena')
    candidates = set()
    for row in rows:
        if not isinstance(row, dict) or not isinstance(row.get('study_alias'), str):
            raise MetadataSourceError('provider_response_invalid', accession, provider='ena')
        identity = row.get('study_accession') or row.get('secondary_study_accession')
        if not isinstance(identity, str) or not re.fullmatch(r'(?:[SED]RP\d+|PRJ[A-Z]+\d+)', identity):
            raise MetadataSourceError('provider_response_invalid', accession, provider='ena')
        if row.get('study_alias') == accession:
            candidates.add(identity)
    verified = []
    seen = set()
    for candidate in sorted(candidates):
        if candidate in seen:
            continue
        root = _response(lambda: source.xml([candidate]), accession)
        if not ET.iselement(root) or root.tag not in {'PROJECT', 'PROJECT_SET', 'STUDY', 'STUDY_SET'}:
            raise MetadataSourceError('provider_response_invalid', accession, provider='ena')
        for node in root.iter():
            if node.tag not in {'PROJECT', 'STUDY'}:
                continue
            ids = _identifiers(node)
            if candidate not in ids:
                continue
            status = [a.findtext('VALUE', '').casefold() for a in node.findall(f'{node.tag}_ATTRIBUTES/{node.tag}_ATTRIBUTE')
                      if a.findtext('TAG') == 'ENA-STATUS']
            evidence = tuple(sorted(set(_links(node, accession))))
            studies = sorted(a for a in ids if re.fullmatch(r'[SED]RP\d+', a))
            projects = sorted(a for a in ids if a.startswith('PRJ'))
            if status != ['public'] or not evidence or len(studies) != 1 or len(projects) > 1:
                continue
            seen.update(ids)
            match = VerifiedENAStudy(studies[0], projects[0] if projects else None, evidence)
            if match not in verified:
                verified.append(match)
    if len(verified) != 1:
        raise MetadataSourceError('fallback_ambiguous' if verified else 'fallback_unverified', accession, provider='ena')
    return verified[0]
