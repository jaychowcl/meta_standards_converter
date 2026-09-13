# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Conservative, network-free projection of explicit archive library methods."""
import re
import logging

_PREPARATION = re.compile(r'^(?:The\s+)?(?:library\s+was|libraries\s+were)\s+'
                          r'(?:made|prepared|constructed)\s+(?:using|with)\s+\S', re.I)
_SCOPED = re.compile(
    r'\b(?:[SED]R[APRSX]\d+|SAM[END][A-Z]?\d+|PRJ[A-Z]+\d+)\b|'
    r'\b(?:sample|study|run)\s+(?:accession|identifier|id|name)\b|'
    r'\b(?:sample|study|run)\s+[\w.-]*\d[\w.#-]*\b|'
    r'\b(?:for|from|of)\s+(?:the\s+)?(?:sample|study|run)\b|'
    r'\b(?:multiplex\w*|barcod\w*|index(?:es|ed)?\s+(?:read|sequence|assignment)\w*)\b|'
    r'\bindex\s+[ACGTN]{4,}\b|'
    r'\b(?:reads?\s+tagged|tagged\s+with)\b', re.I)
_ABBREVIATION = re.compile(r'\b(?:e\.g|i\.e|etc|vs|Inc|Co|Dr|Mr|Mrs)\.$')


def method_sentences(description):
    """Select complete explicitly labelled preparation sentences; never paraphrase."""
    if not isinstance(description, str):
        return []
    sentences, start = [], 0
    for end in re.finditer(r'[.!?](?=\s|$)', description):
        candidate = description[start:end.end()].strip()
        if end[0] == '.' and _ABBREVIATION.search(candidate):
            continue
        sentences.append(candidate)
        start = end.end()
    if description[start:].strip():
        sentences.append(description[start:].strip())
    return [s for s in sentences if _PREPARATION.match(s) and not _SCOPED.search(s)]


def library_description(original, experiment_description):
    """Append supplied method text while preserving the original protocol verbatim."""
    result = original or ''
    for sentence in method_sentences(experiment_description):
        if ' '.join(sentence.split()) not in ' '.join(result.split()):
            result = result + '\n\n' + sentence if result else sentence
    return result or None


def prepare_native_protocols(data):
    """Upgrade generated native definitions on a mutable export copy only.

    A native generated name and an exact assay/reference join are required for
    method recovery. Incoming GEO/AE definitions and shared ambiguous workflows
    are never interpreted as native experiment evidence.
    """
    if data.get('source', {}).get('format') not in {'ENA', 'SRA'}:
        return
    series = data.get('series', {})
    instruments = {}
    for sample in data.get('sample', []):
        for run in sample.get('sra_run', []):
            if run.get('experiment') and run.get('instrument_model'):
                instruments.setdefault(run['experiment'], set()).add(run['instrument_model'])
    applications = {}
    for path in series.get('assay_paths', []):
        assays = [s for s in path.get('steps', []) if s.get('kind') == 'assay']
        identity = (assays[0].get('name'), assays[0].get('description')) if len(assays) == 1 else (None, None)
        for step in path.get('steps', []):
            if step.get('kind') == 'protocol_application' and step.get('protocol_ref'):
                applications.setdefault(step['protocol_ref'], set()).add(identity)
    for protocol in series.get('protocols', []):
        match = re.fullmatch(r'([SED]RX\d+):library', protocol.get('name', ''))
        if not match or protocol.get('type', {}).get('value') != 'library construction protocol':
            continue
        experiment = match[1]
        hardware = protocol.get('hardware')
        if isinstance(hardware, list) and len(hardware) == 1 and hardware[0] in instruments.get(experiment, set()):
            protocol.pop('hardware')
        candidates = applications.get(protocol['name'], set())
        if not candidates:
            continue
        if len(candidates) != 1 or next(iter(candidates))[0] != experiment:
            logging.getLogger(__name__).warning('%s: ambiguous assay association; protocol method recovery skipped', protocol['name'])
            continue
        protocol['description'] = library_description(protocol.get('description'), next(iter(candidates))[1])
