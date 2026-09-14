# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Shared affirmative preparation evidence for technology and chemistry.

The bounded grammar never rewrites source text. Cell suspensions, sequencers
and generic vendor names do not establish a single-cell assay.
"""
import re

DROPSEQ = re.compile(r'(?<![a-z0-9])drop[- ]?seq(?![a-z0-9])', re.I)
TENX = re.compile(r'(?<![a-z0-9])(?:10[x×](?:\s+genomics)?|chromium)(?![a-z0-9])', re.I)
SINGLE = re.compile(r"(?<![a-z])(?:sc|sn)[- ]?(?:rna|atac)(?:[-_ ]?seq)?(?![a-z])|(?<![a-z])cite[- ]seq(?![a-z0-9])|single[- ](?:cell|nucleus|nuclei)\s+(?:(?:plate|droplet)[- ](?:based\s+)?)?(?:rna|transcriptom|sequenc|atac|librar|assay|[35]['′’])", re.I)
PLATE = re.compile(r'\bplate[- ]based\b|\bsingle[- ]cell\s+plate\s+(?:assay|sequenc|librar)|\b(?:single[- ]cells?|individual cells?)\b.{0,70}\b(?:sorted|deposited|dispensed)\b.{0,55}\b(?:wells?|plates?)\b', re.I)
NON_PREP = re.compile(r'\b(?:not|without|compatible|compatibility|software|sequencer|recommended)\b', re.I)
BULK = re.compile(r'\bbulk\s+(?:rna|dna|sequenc|librar)|\bnot\s+(?:a\s+)?single[- ]cell|\bwithout\s+single[- ]cell', re.I)
GENOMIC = re.compile(r'\b(?:genome|genomic|wgs|whole[- ]genome)\b', re.I)


def clauses(text):
    if '://' in text and not re.search(r'\s', text):
        return []
    return re.split(r'(?<!\d)\.(?!\d)|\n|;', text)


def single_cell_signal(text):
    return any(not NON_PREP.search(c) and (SINGLE.search(c) or DROPSEQ.search(c)
               or (TENX.search(c) and re.search(r'single[- ]cell|\bgem[- ]x\s+flex\b|\bflex\s+(?:gene expression|kit|reagent|protocol)|\bmultiome\b', c, re.I))) for c in clauses(text))


def methods(text):
    found = set()
    for c in clauses(text):
        bulk = bool(BULK.search(c) or (TENX.search(c) and GENOMIC.search(c) and not single_cell_signal(c)))
        if bulk:
            found.add('bulk')
        if NON_PREP.search(c):
            if TENX.search(c) and re.search(r'\b(?:not|without)\b', c, re.I):
                found.add('not_10x')
            continue
        if DROPSEQ.search(c):
            found.add('dropseq')
        if TENX.search(c) and not bulk:
            found.add('10x')
        if PLATE.search(c):
            found.add('plate')
        if re.search(r'\bdroplet[- ]based\b|\bdroplet\s+single[- ]cell\b|\b(?:encapsulated|partitioned)\b.{0,50}\bdroplets?\b', c, re.I):
            found.add('droplet')
    return found


def scoped_method(sample, channel=None, run=None):
    """Explicit run/channel identity precedes sample identity and description."""
    prefix = f"sample[{sample.get('iid') or 'unknown'}]"
    levels = [[], [], []]
    for key in ('library_name', 'description'):
        if (run or {}).get(key):
            levels[0].append((prefix + '.run.' + key, str(run[key])))
    channels = [channel] if channel is not None else sample.get('channel', [])
    for i, c in enumerate(channels):
        for j, characteristic in enumerate(c.get('characteristics', [])):
            tag = re.sub(r'[\s_-]+', '_', str(characteristic.get('name') or characteristic.get('tag', '')).strip().casefold())
            if tag in {'assay', 'assay_type', 'library_type', 'library_name', 'technology'}:
                levels[0].append((f'{prefix}.channel[{i}].characteristics[{j}].value', str(characteristic.get('value') or '')))
    for key in ('title', 'library_name'):
        if sample.get(key):
            levels[1].append((prefix + '.' + key, str(sample[key])))
    if sample.get('description'):
        levels[2].append((prefix + '.description', str(sample['description'])))
    for group in levels:
        evidence = tuple((path, text, tuple(sorted(methods(text)))) for path, text in group if methods(text))
        values = {m for _, _, values in evidence for m in values}
        if values:
            if 'droplet' in values and values & {'dropseq', '10x'}:
                values.remove('droplet')
            if 'not_10x' in values and values & {'dropseq', 'plate', 'droplet'} and '10x' not in values:
                values.remove('not_10x')
            return (next(iter(values)) if len(values) == 1 else 'ambiguous'), evidence
        if any(single_cell_signal(text) for _, text in group):
            return 'single_cell', tuple((path, text, ('single_cell',)) for path, text in group if single_cell_signal(text))
    return None, ()
