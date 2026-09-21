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
                          single_cell_signal, CHEMISTRY_IDENTIFIERS, preparation_operation)


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
_CELL = re.compile(r'\bsingle[- ]cell\b|\b(?:split[- ]seq|evercode)\b|\bcombinatorial indexing\b', re.I)


def _signals(text):
    found = set()
    for clause in clauses(text.replace("_", " ")):
        # A bulk-control description can mention the study it controls.
        if BULK.search(clause) or re.match(r'^\s*bulk\b', clause, re.I):
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


def _library(sample, run, data):
    levels = library_evidence_levels(sample, run=run, series=(data or {}).get('series'))
    nucleus_material = tuple(ModalityEvidence(path, text, ('single_nucleus',))
        for path, text in levels[2] if re.search(r'\b(?:nuclei|nucleus)\s+isolation\b', text, re.I))
    for facts in levels:
        evidence = tuple(ModalityEvidence(path, text, _signals(text)) for path, text in facts if _signals(text))
        candidates = {value for fact in evidence for value in fact.candidates}
        if candidates:
            if len(candidates) == 1:
                value = next(iter(candidates))
                if value == 'single_cell' and nucleus_material:
                    return ModalityDecision('single_nucleus', evidence + nucleus_material)
                return ModalityDecision(value, evidence)
            return ModalityDecision('unknown', evidence,
                (ModalityDiagnostic('conflicting_modality', tuple(f.path for f in evidence)),))
    return ModalityDecision('unknown')


@preparation_operation
def resolve_modality(sample, *, data=None):
    """Resolve libraries independently, then aggregate without treating runs as libraries."""
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
