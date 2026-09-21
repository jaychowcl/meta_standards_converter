"""Study-family traversal follows typed provider relationships, never citations."""
from collections import deque
from dataclasses import replace
import re
from urllib.parse import parse_qs, urlsplit

from .contracts import Diagnostic, InputSpec
from meta_standards_converter.metadata.preparation_scope import loading_source

PROJECT = re.compile(r'PRJ(?:NA|EB|DB|DA)\d+')


def nodes(tree):
    if isinstance(tree, dict):
        yield tree
        for child in tree.get('children', []):
            yield from nodes(child)


def project_records(package, provider):
    from .preparation import study_ids
    own = set(study_ids(package))
    for record in package.to_mapping().get('extensions', {}).get('insdc', {}).get('records', []):
        if record.get('provider', '').lower() != provider:
            continue
        root = record.get('metadata', {})
        if root.get('tag') == 'PROJECT' and root.get('attributes', {}).get('accession') in own:
            yield root
        elif root.get('tag') == 'Project':
            ids = [n.get('attributes', {}) for n in nodes(root) if n.get('tag') == 'ArchiveID']
            if any(i.get('accession') in own for i in ids):
                yield root


def recorded_neighbors(package, provider):
    found = set()
    if provider == 'geo':
        for relation in package.to_mapping().get('series', {}).get('relation', []):
            if str(relation.get('type', '')).casefold() not in {'superseries', 'subseries', 'superseries of', 'subseries of'}:
                continue
            target = relation.get('target', '')
            candidates = [target]
            if isinstance(target, str) and target.startswith(('http://', 'https://')):
                url = urlsplit(target)
                if url.hostname in {'www.ncbi.nlm.nih.gov', 'ncbi.nlm.nih.gov'}:
                    candidates = parse_qs(url.query).get('acc', [])
            found.update(a.upper() for a in candidates if isinstance(a, str) and re.fullmatch(r'GSE\d+', a, re.I))
    if provider == 'ena':
        for root in project_records(package, provider):
            for node in nodes(root):
                if node.get('tag') in {'CHILD_PROJECT', 'PARENT_PROJECT'}:
                    accession = node.get('attributes', {}).get('accession', '')
                    if PROJECT.fullmatch(accession):
                        found.add(accession)
    return sorted(found)


class StudyFamilies:
    def __init__(self, context):
        self.context = context
        self.cache = {}
        self.projects = {}

    def _neighbors(self, package, provider):
        """ENA retains explicit XML edges; NCBI exposes typed, verified ELink edges."""
        found = set(recorded_neighbors(package, provider))
        if provider != 'sra':
            return sorted(found), []
        from meta_standards_converter.converters.sra2json import SRA2JSONConverter
        from meta_standards_converter.sources.archive_support import Resolution
        issues = []
        for root in project_records(package, provider):
            for node in nodes(root):
                if node.get('tag') != 'ArchiveID':
                    continue
                attrs = node.get('attributes', {})
                uid, accession = attrs.get('id'), attrs.get('accession', '')
                if not uid or not uid.isdigit() or not PROJECT.fullmatch(accession):
                    continue
                if accession not in self.projects:
                    result, linked = Resolution(), set()
                    try:
                        source = self.context.service('sra2json', lambda: SRA2JSONConverter(resource_profile=self.context.profile)).source
                        for direction in ('u2d', 'd2u'):
                            name = 'bioproject_bioproject_' + direction
                            ids = source.publication_links('bioproject', 'bioproject', [uid], name)[uid]
                            for target in ids:
                                doc = source.project_xml(target, result)
                                if doc is not None:
                                    linked.update(n.get('accession') for n in doc.findall('.//ProjectID/ArchiveID') if PROJECT.fullmatch(n.get('accession', '')))
                    except Exception as exc:
                        result.issues.append(f'{accession}: project hierarchy unavailable ({type(exc).__name__})')
                    self.projects[accession] = sorted(linked), result.issues
                children, warnings = self.projects[accession]
                found.update(children)
                issues.extend(warnings)
        return sorted(found), issues

    def expand(self, loaded):
        from .preparation import provider_of, study_ids
        context = self.context
        if not context.preparation_policy.expand_studies:
            loaded.preparation.append({'operation': 'study_expansion', 'status': 'skipped', 'reason': 'disabled'})
            return loaded
        queue = deque(loaded.metadata.groups)
        groups = list(loaded.metadata.groups)
        seen = {(provider_of(p, loaded.provider), a) for g in groups for p in g.packages for a in study_ids(p)}
        known_groups = {g.dataset_id for g in groups}
        while queue:
            group = queue.popleft()
            for package in group.packages:
                provider = provider_of(package, loaded.provider)
                neighbors, issues = self._neighbors(package, provider)
                loaded.diagnostics.extend(Diagnostic('source_partial', issue, 'expansion', 'warning') for issue in issues)
                for accession in neighbors:
                    identity = (provider, accession)
                    if identity in seen:
                        continue
                    seen.add(identity)
                    record = {'operation': 'study_expansion', 'dataset': group.dataset_id,
                              'provider': provider, 'target': accession, 'status': 'attempted'}
                    loaded.preparation.append(record)
                    key = (*identity, context.preparation_policy, repr(context.input_options))
                    try:
                        if key not in self.cache:
                            # Same provider throughout a family, including curators/off.
                            handler = context.converter.handlers[provider + '_accession']
                            with loading_source(context.preparation_policy):
                                candidate = handler.load(InputSpec(accession), context)
                            if provider == "geo" and not any(accession in study_ids(p) for g in candidate.metadata.groups for p in g.packages):
                                raise ValueError('Retrieved family member does not match the requested identity')
                            if any(d.severity == 'error' for d in candidate.metadata.diagnostics):
                                raise ValueError('Retrieved family member failed metadata validation')
                            self.cache[key] = candidate
                        candidate = self.cache[key]
                        if isinstance(candidate, Exception):
                            raise candidate
                        for new in candidate.metadata.groups:
                            if new.dataset_id not in known_groups:
                                groups.append(new)
                                queue.append(new)
                                known_groups.add(new.dataset_id)
                            seen.update((provider, a) for p in new.packages for a in study_ids(p))
                        loaded.diagnostics.extend(candidate.diagnostics)
                        record['status'] = 'partial' if any(d.code in {'source_partial', 'source_failed'} for d in candidate.diagnostics) else 'completed'
                    except Exception as exc:
                        self.cache.setdefault(key, exc)
                        record['status'] = 'failed'
                        loaded.diagnostics.append(Diagnostic('source_partial', f'{accession}: family expansion failed ({type(exc).__name__})', 'expansion', 'warning'))
        loaded.metadata = replace(loaded.metadata, groups=tuple(groups))
        return loaded
