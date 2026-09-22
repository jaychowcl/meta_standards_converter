# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Conservative biological modality from scoped, positive library evidence.

Assay strategy, manufacturer identity and renderer defaults are not modality.
Evidence and conflicts remain inspectable even when a decision is unknown.
"""
from dataclasses import dataclass
import re
from .preparation import (library_evidence_levels, clauses, NON_PREP, BULK,
                          single_cell_signal, CHEMISTRY_IDENTIFIERS, preparation_operation, package_mapping)


@dataclass(frozen=True)
class ModalityEvidence:
    path: str
    text: str
    candidates: tuple[str, ...]


@dataclass(frozen=True)
class ModalityDiagnostic:
    code: str
    paths: tuple[str, ...]


@dataclass(frozen=True)
class ModalityDecision:
    value: str
    evidence: tuple[ModalityEvidence, ...] = ()
    diagnostics: tuple[ModalityDiagnostic, ...] = ()


_SPATIAL = re.compile(r'\b(?:visium|geomx|slide[- ]seq|stereo[- ]seq|seq[- ]scope)\b|\bspatial\s+(?:transcriptom|rna|gene expression|sequenc)', re.I)
_NUCLEUS = re.compile(r'\bsingle[- ](?:nucleus|nuclei)\b|\bsn[- ]?(?:rna|atac)(?:[- ]?seq)?\b', re.I)
_BULK_NUCLEAR = re.compile(r'\bbulk\s+(?:nuclear\s+)?(?:rna|dna|librar|sequenc)', re.I)
_NUCLEAR_MATERIAL = re.compile(r'\b(?:nuclei|nucleus)\s+(?:(?:were|was|are|is)\s+)?(?:extract|isolat|recover)|\b(?:extraction|isolation|recovery)\s+of\s+(?:nuclei|nucleus)|\b(?:isolated|extracted)\s+nuclei', re.I)
_WHOLE_CELL = re.compile(r'\b(?:whole|intact)\s+cells?\b', re.I)
_NEGATIVE = re.compile(r'\b(?:no|not|without|never)\b', re.I)
_CELL = re.compile(r'\bsingle[- ]cell\b|\b(?:split[- ]seq|evercode)\b|\bcombinatorial indexing\b', re.I)


def _signals(text):
    found = set()
    for clause in clauses(text.replace("_", " ")):
        # A bulk-control description can mention the study it controls.
        if BULK.search(clause) or _BULK_NUCLEAR.search(clause) or re.match(r'^\s*bulk\b', clause, re.I):
            found.add('bulk')
            clause = re.split(r'\bcontrol\s+(?:for|in|of)\b', clause, flags=re.I)[0]
        if NON_PREP.search(clause):
            continue
        if _SPATIAL.search(clause): found.add('spatial')
        if _NUCLEUS.search(clause): found.add('single_nucleus')
        elif _CELL.search(clause) or single_cell_signal(clause): found.add('single_cell')
    if text.strip().casefold() in CHEMISTRY_IDENTIFIERS:
        found.add('single_cell')
    if text.strip().casefold() in {'single cell', 'transcriptomic single cell'}:
        found.add('single_cell')
    return tuple(sorted(found))


def _material(text, pattern):
    return any(pattern.search(clause) and not _NEGATIVE.search(clause)
               and not NON_PREP.search(clause) for clause in clauses(text))


def _library(sample, run, data):
    levels = library_evidence_levels(sample, run=run, series=(data or {}).get('series'))
    # Library identity can exclude a shared sample description, but not a
    # positively bound library preparation. Archive categories remain broad.
    identity = levels[0] if any(_signals(text) for _, text in levels[0]) else levels[1]
    bound = [(path, text) for path, text in levels[2]
             if not identity or identity is levels[1]
             or '.run[' in path or path.startswith('series.protocols[')]
    facts = identity + bound
    evidence = tuple(ModalityEvidence(path, text, _signals(text))
                     for path, text in facts if _signals(text))
    candidates = {value for fact in evidence for value in fact.candidates}
    broad = tuple(ModalityEvidence(path, text, _signals(text))
                  for path, text in levels[3] if _signals(text))
    broad_values = {value for fact in broad for value in fact.candidates}
    # The archive's single-cell umbrella is compatible with spatial and
    # single-nucleus assays, but contradicts an explicit bulk preparation.
    if candidates and candidates <= {'spatial', 'single_nucleus'}:
        broad_values.discard('single_cell')
    candidates.update(broad_values)
    evidence += broad
    nuclear = tuple(ModalityEvidence(path, text, ('single_nucleus',))
                    for path, text in bound if _material(text, _NUCLEAR_MATERIAL))
    whole = any(_material(text, _WHOLE_CELL) for _, text in facts)
    if candidates == {'single_cell'} and nuclear:
        candidates = {'single_nucleus', 'single_cell'} if whole else {'single_nucleus'}
        evidence += nuclear
    if candidates == {'single_cell', 'single_nucleus'} and not whole:
        candidates.remove('single_cell')
    if len(candidates) > 1:
        return ModalityDecision('unknown', evidence,
            (ModalityDiagnostic('conflicting_modality', tuple(f.path for f in evidence)),))
    return ModalityDecision(next(iter(candidates)) if candidates else 'unknown', evidence)


@preparation_operation
def resolve_modality(sample, *, data=None):
    """Resolve libraries independently, then aggregate without treating runs as libraries."""
    data = package_mapping(data)
    runs = sample.get('sra_run') or [None]
    libraries = {}
    for index, run in enumerate(runs):
        run = run or {}
        key = run.get('experiment') or run.get('library_name') or run.get('run') or f'sample:{index}'
        libraries.setdefault(key, []).append(_library(sample, run, data))
    resolved = []
    for decisions in libraries.values():
        evidence = tuple(dict.fromkeys(f for d in decisions for f in d.evidence))
        diagnostics = tuple(dict.fromkeys(f for d in decisions for f in d.diagnostics))
        values = {d.value for d in decisions}
        if len(values) != 1:
            diagnostics += (ModalityDiagnostic('conflicting_library_modality', tuple(f.path for f in evidence)),)
        resolved.append(ModalityDecision(next(iter(values)) if len(values) == 1 else 'unknown', evidence, diagnostics))
    evidence = tuple(dict.fromkeys(f for d in resolved for f in d.evidence))
    diagnostics = tuple(dict.fromkeys(f for d in resolved for f in d.diagnostics))
    values = {d.value for d in resolved}
    value = 'unknown' if 'unknown' in values else next(iter(values)) if len(values) == 1 else 'mixed'
    return ModalityDecision(value, evidence, diagnostics)
