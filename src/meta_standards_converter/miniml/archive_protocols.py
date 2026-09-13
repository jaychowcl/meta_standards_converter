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

_PREPARATION = re.compile(r'^(?:The\s+)?(?:library\s+was|libraries\s+were)\s+'
                          r'(?:made|prepared|constructed)\s+(?:using|with)\s+\S', re.I)
_SCOPED = re.compile(
    r'\b(?:[SED]R[APRSX]\d+|SAM[END][A-Z]?\d+|PRJ[A-Z]+\d+)\b|'
    r'\b(?:sample|study|run)\s+(?:accession|identifier|id|name)\b|'
    r'\b(?:sample|study|run)\s+[\w.-]*\d[\w.#-]*\b|'
    r'\b(?:multiplex\w*|barcod\w*|index(?:es|ed)?\s+(?:read|sequence|assignment)\w*)\b|'
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
