# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from __future__ import annotations
from dataclasses import dataclass, replace
import re
from .chemistry import resolve_chemistry
from .preparation import single_cell_signal, scoped_method, methods, clauses, NON_PREP, control_role, incompatible_preparation, bound_protocols, control_preparation
from meta_standards_converter.magetab.protocols import ProtocolRegistry
import os
from urllib.parse import urlparse
from meta_standards_converter.helpers.json_helper import JSONHandler


def normalized_extension(path: str) -> str:
    parsed = urlparse(str(path))
    basename = os.path.basename(parsed.path or str(path)).lower()
    for suffix in (".gz", ".zip", ".bz2", ".xz"):
        if basename.endswith(suffix):
            basename = basename[: -len(suffix)]
            break
    return os.path.splitext(basename)[1]

def has_array_files(data: dict) -> bool:
    handler = JSONHandler()
    values = []
    for path in (
        "platform.*.supplementary_data.*.value",
        "sample.*.supplementary_data.*.value",
        "sample.*.raw_data.*.value",
        "series.supplementary_data.*.value",
    ):
        values.extend(x for x in handler._from_path(data, path) if x)
    extensions = (".cel", ".gpr", ".idat", ".chp", ".txt", ".tif", ".tiff", ".exp", ".rpt", ".cab")
    return any(normalized_extension(value) in extensions for value in values)

def series_identity(data: dict) -> str | None:
    """Return the model-authoritative series iid, with accession fallback."""
    series = data.get("series") if isinstance(data, dict) else None
    for item in series if isinstance(series, list) else [series]:
        if not isinstance(item, dict):
            continue
        iid = ProtocolRegistry.clean(item.get("iid"))
        if iid:
            validated = _validated_study_identity(iid)
            if validated:
                return validated
        accessions = item.get("accession")
        for accession in accessions if isinstance(accessions, list) else [accessions]:
            value = accession.get("value") if isinstance(accession, dict) else accession
            cleaned = ProtocolRegistry.clean(value)
            if not cleaned:
                continue
            validated = _validated_study_identity(cleaned)
            if validated:
                return validated
    return None

def _validated_study_identity(value: str) -> str | None:
    upper = value.upper()
    if upper.startswith("GSE"):
        return upper if upper[3:].isdigit() else None
    return value

def _detect_base_technology(data: dict) -> str:
    """Identify sequencing/array evidence without interpreting free-text methods."""
    handler = JSONHandler()
    values = lambda path: (str(x).casefold() for x in handler._from_path(data, path) if x)
    platform_tech = " ".join(values("platform.*.technology"))
    relations = handler._from_path(data, "sample.*.relation.*")
    has_sra = any(isinstance(r, dict) and str(r.get('type', '')).casefold() == 'sra' for r in relations)
    if ("high-throughput sequencing" in platform_tech or has_sra
            or any(values("sample.*.library_strategy")) or 'sra' in values("sample.*.type")):
        return "bulk_sequencing"
    geo_array_categories = {
        "in situ oligonucleotide", "spotted oligonucleotide",
        "mixed spotted oligonucleotide", "spotted dna/cdna",
        "spotted peptide or protein", "antibody", "tissue", "oligonucleotide beads",
    }
    platform_categories = {" ".join(value.split()) for value in values("platform.*.technology")}
    if platform_categories & geo_array_categories or "array" in platform_tech or has_array_files(data):
        return "array"
    return "generic"


# Routing is separate from chemistry: a technology decision never supplies a
# library version or a sequencing recipe.
@dataclass(frozen=True)
class TechnologyEvidence:
    path: str
    text: str
    candidates: tuple[str, ...]


@dataclass(frozen=True)
class TechnologyDiagnostic:
    code: str
    paths: tuple[str, ...]


@dataclass(frozen=True)
class TechnologyDecision:
    handler: str
    evidence: tuple[TechnologyEvidence, ...] = ()
    diagnostics: tuple[TechnologyDiagnostic, ...] = ()
    control_role: str | None = None


def _items(value):
    return value if isinstance(value, list) else [value] if value else []


_SPATIAL = re.compile(r'\bvisium\b|\bspatial\s+(?:transcriptom|rna|gene expression|sequenc)', re.I)


def _signals(text):
    signals = set()
    for clause in clauses(text):
        if NON_PREP.search(clause):
            continue
        spatial = bool(_SPATIAL.search(clause))
        if spatial:
            signals.add('spatial')
        if single_cell_signal(clause):
            signals.add('single_cell')
    return tuple(sorted(signals))


def resolve_technology(sample: dict, channel: dict | None = None,
                       run: dict | None = None, *, data: dict | None = None) -> TechnologyDecision:
    role = control_role(sample, channel, run)
    incompatible = incompatible_preparation(sample, channel, run, (data or {}).get('series'))
    if incompatible:
        evidence = tuple(TechnologyEvidence(p, t, ('incompatible_rna_preparation',)) for p, t in incompatible)
        identity = str((run or {}).get('library_name') or sample.get('title') or '')
        explicit_chip = bool(re.search(r'(?:single[- ]cell\s+|\bsc[- ]?)chip', identity, re.I))
        handler = 'single_cell_sequencing' if explicit_chip else 'sequencing'
        return TechnologyDecision(handler, evidence, (TechnologyDiagnostic('incompatible_preparation_evidence',
                                  tuple(p for p, _ in incompatible)),), role)
    decision = _resolve_technology(sample, channel, run, data=data)
    return replace(decision, control_role=role)


def _resolve_technology(sample: dict, channel: dict | None = None,
                       run: dict | None = None, *, data: dict | None = None) -> TechnologyDecision:
    """Route one sample/library using identity before shared method descriptions.

    Input dictionaries and source objects are never modified. Equally applicable
    conflicting identities produce generic sequencing and an auditable warning.
    """
    data = data or {}
    # Preserve array/platform routing before interpreting sequencing methods.
    scoped = dict(data, sample=[sample], series={})
    if sample.get('platform_ref'):
        refs = {x.get('ref') for x in _items(sample['platform_ref']) if isinstance(x, dict)}
        scoped['platform'] = [p for p in _items(data.get('platform')) if isinstance(p, dict) and p.get('iid') in refs]
    base_technology = _detect_base_technology(scoped)
    if base_technology == 'array':
        return TechnologyDecision('array')
    prefix = f"sample[{sample.get('iid') or 'unknown'}]"
    channels = [channel] if channel is not None else _items(sample.get('channel'))
    levels = [[], [], [], [], []]  # library/channel identity, sample identity, preparation, source, shared fallback

    def add(level, path, text):
        if text:
            if level == 2:
                text = control_preparation(str(text), control_role(sample, channel, run))
            signals = ('single_cell',) if level == 3 and str(text).casefold() in {'single cell', 'transcriptomic single cell'} else _signals(str(text))
            levels[level].append(TechnologyEvidence(path, str(text), signals))

    from meta_standards_converter.metadata.preparation import library_evidence_levels
    for level, facts in enumerate(library_evidence_levels(sample, channel, run, series=data.get('series'))):
        for path, text in facts:
            add(level, path, text)
    method, method_evidence = scoped_method(sample, channel, run)
    if method == 'bulk':
        return TechnologyDecision('bulk_sequencing', tuple(TechnologyEvidence(p, t, ('bulk',)) for p, t, _ in method_evidence))
    if method == 'ambiguous' and any('bulk' in values for _, _, values in method_evidence):
        return TechnologyDecision('sequencing', tuple(TechnologyEvidence(p, t, values) for p, t, values in method_evidence),
            (TechnologyDiagnostic('ambiguous_technology', tuple(p for p, _, _ in method_evidence)),))
    chemistry = resolve_chemistry(sample, channel=channel, run=run, series=data.get('series'))
    for fact in chemistry.evidence:
        if fact.field == 'identifier' and any(f.path == fact.path and f.field == 'manufacturer' for f in chemistry.evidence):
            levels[0].append(TechnologyEvidence(fact.path, fact.text, ('single_cell',)))
    for i, series in enumerate(_items(data.get('series'))):
        if control_role(sample, channel, run) == 'empty control':
            break
        if isinstance(series, dict):
            for key in ('title','summary','overall_design'):
                text = str(series.get(key) or '')
                # Explicitly assigned or subset protocols cannot become a
                # universal fallback simply because they mention one method.
                for sentence in re.split(r'(?<!\d)\.(?!\d)|\n', text):
                    accessions = set(re.findall(r'\bGSM\d+\b', sentence, re.I))
                    if accessions and str(sample.get('iid', '')).upper() not in {a.upper() for a in accessions}:
                        continue
                    if re.search(r'\b(?:other|some|subset of) (?:samples|libraries)\b', sentence, re.I):
                        continue
                    add(4, f'series[{i}].{key}', sentence)

    for level, evidence in enumerate(levels):
        candidates = {c for e in evidence for c in e.candidates}
        if not candidates:
            continue
        if len(candidates) > 1:
            relevant = tuple(e for e in evidence if e.candidates)
            return TechnologyDecision('sequencing', relevant,
                (TechnologyDiagnostic('ambiguous_technology', tuple(e.path for e in relevant)),))
        selected = next(iter(candidates))
        # Lower-level text supports identity only if it does not describe a
        # competing method. Keep full original evidence for auditability.
        supporting = tuple(e for group in levels[level:] for e in group if e.candidates == (selected,))
        if selected == 'spatial':
            return TechnologyDecision('spatial_sequencing', supporting)
        method_support = tuple(TechnologyEvidence(p, t, (selected,)) for p, t, _ in method_evidence)
        if method == 'ambiguous':
            return TechnologyDecision('single_cell_sequencing', supporting + method_support,
                (TechnologyDiagnostic('ambiguous_preparation', tuple(p for p, _, _ in method_evidence)),))
        if method in {'dropseq', '10x', 'droplet'} or chemistry.manufacturer == '10x Genomics':
            return TechnologyDecision('droplet_single_cell_sequencing', supporting + method_support)
        # A library-source field can establish single-cell identity after its
        # bound preparation supplied the format without an explicit SC keyword.
        # Applicable preparation precedes unrelated study-wide method fallback.
        preparations = tuple(e for e in levels[2] if methods(e.text) - {'bulk'})
        format_evidence = preparations or tuple(e for group in levels[level:] for e in group)
        formats = {m for e in format_evidence for m in methods(e.text) if m != 'bulk'}
        supporting += tuple(e for e in preparations if e not in supporting)
        if method == 'not_10x':
            formats.discard('10x')
        if method == 'plate' or formats == {'plate'}:
            return TechnologyDecision('plate_single_cell_sequencing', supporting + method_support)
        if formats and formats <= {'10x', 'dropseq', 'droplet'}:
            return TechnologyDecision('droplet_single_cell_sequencing', supporting)
        diagnostics = ((TechnologyDiagnostic('ambiguous_preparation', tuple(e.path for e in format_evidence)),)
                       if len(formats) > 1 else ())
        return TechnologyDecision('single_cell_sequencing', supporting, diagnostics)

    return TechnologyDecision('sequencing' if control_role(sample, channel, run) == 'empty control' else base_technology)


def detect_ae_technology(data: dict) -> str:
    """Return an IDF summary key; SDRF dispatch resolves each sample/library."""
    decisions = {resolve_technology(sample, data=data).handler for sample in _items(data.get('sample')) if isinstance(sample, dict)}
    if len(decisions) == 1:
        return next(iter(decisions))
    if decisions and decisions.isdisjoint({'array','generic'}):
        return 'sequencing'
    return _detect_base_technology(data) if not decisions else 'generic'
