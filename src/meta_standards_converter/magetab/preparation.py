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
SINGLE = re.compile(r"(?<![a-z])(?:sc|sn|sci)[- ]?(?:rna|atac|chip)(?:[-_ ]?seq)?(?![a-z])|(?<![a-z])cite[- ]seq(?![a-z0-9])|single[- ](?:cell|nucleus|nuclei)\s+(?:(?:plate|droplet)[- ](?:based\s+)?)?(?:rna|transcriptom|transcriptional\s+profiling|sequenc|atac|chip|librar|assay|auto\s+prep|[35]['′’])|single[- ]cell\s+combinatorial\s+indexing", re.I)
PLATE = re.compile(r'\bplate[- ]based\b|\bsingle[- ]cell\s+plate\s+(?:assay|sequenc|librar)|\b(?:single[- ](?:[a-z0-9+-]+\s+){0,3}cells?|individual cells?)\b.{0,70}\b(?:sorted|deposited|dispensed)\b.{0,75}\b(?:wells?|plates?)\b', re.I)
NON_PREP = re.compile(r'\b(?:not|without|compatible|compatibility|software|sequencer|recommended)\b', re.I)
BULK = re.compile(r'\bbulk\s+(?:rna|dna|sequenc|librar|control)|\bnot\s+(?:a\s+)?single[- ]cell|\bwithout\s+single[- ]cell', re.I)
GENOMIC = re.compile(r'\b(?:genome|genomic|wgs|whole[- ]genome)\b', re.I)


# Exact identifiers, reviewed against Cell Ranger's documented chemistry options.
# These identify chemistry only, never a sequencing recipe.
CHEMISTRY_IDENTIFIERS = {
    **{f'sc3pv{v}': ('3 prime', (str(v),)) for v in range(1, 5)},
    'sc3pv3ht': ('3 prime', ('3.1',)),
    'sc5p-pe': ('5 prime', ()), 'sc5p-r2': ('5 prime', ()),
    'sc5p-pe-v3': ('5 prime', ('3',)), 'sc5p-r2-v3': ('5 prime', ('3',)),
    'sc5pht': ('5 prime', ('2',)),
}


def clauses(text):
    if '://' in text and not re.search(r'\s', text):
        return []
    # A supplier/device parenthesis can contain a semicolon inside one method
    # sentence (e.g. FACS sorted (BD Influx; BD Biosciences) into wells).
    result = []
    for sentence in re.split(r'(?<!\d)\.(?!\d)|\n', text):
        depth, start = 0, 0
        for i, char in enumerate(sentence):
            depth += char == '('
            depth = max(0, depth - (char == ')'))
            if char == ';' and depth == 0:
                result.append(sentence[start:i]); start = i + 1
        result.append(sentence[start:])
    return result


def control_role(sample, channel=None, run=None):
    """Control identity is distinct from the preparation used for that control."""
    groups = [[str((run or {}).get(k) or '') for k in ('library_name', 'description')],
              [str(sample.get(k) or '') for k in ('title', 'description')]]
    channels = [channel] if channel is not None else sample.get('channel', [])
    groups.insert(1, [str(a.get('value') or '') for c in channels for a in c.get('characteristics', [])
                      if re.sub(r'[ _-]+', ' ', str(a.get('name') or a.get('tag') or '')).casefold()
                      in {'control', 'control type', 'sample type'}])
    for group in groups:
        for value in group:
            if re.search(r'\bbulk control\b', value, re.I):
                return 'bulk control'
            if re.search(r'\bno[- ]cell\s+control\b|\bempty\s+control\b|\bnegative\s+control\b.*(?<![a-z])empty(?![a-z])', value, re.I):
                return 'empty control'
    return None


def incompatible_preparation(sample, channel=None, run=None):
    """Recognize explicit RNA preparation attached to a genomic ChIP assay.

    This is a narrow compatibility check, not a biological target guesser.
    """
    strategy = str((run or {}).get('library_strategy') or sample.get('library_strategy') or '').casefold()
    source = str((run or {}).get('library_source') or sample.get('library_source') or '').casefold()
    if strategy not in {'chip-seq', 'chip seq'} or source != 'genomic':
        return ()
    prefix = f"sample[{sample.get('iid') or 'unknown'}]"
    values = [(prefix+'.run.'+k, str((run or {}).get(k) or ''))
              for k in ('library_protocol', 'library_construction_protocol', 'description')]
    channels = [channel] if channel is not None else sample.get('channel', [])
    values += [(prefix+f'.channel[{i}].extract_protocol', str(c.get('extract_protocol') or ''))
               for i, c in enumerate(channels)]
    rna = re.compile(r'single[- ]cell\s+rna|whole\s+transcriptome\s+amplif|reverse[- ]transcription\s+of\s+mrna', re.I)
    return tuple((path, value) for path, value in values
                 if any(not NON_PREP.search(c) and rna.search(c) for c in clauses(value)))


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
    structured_paths = set()
    for key in ('library_name', 'description'):
        if (run or {}).get(key):
            levels[0].append((prefix + '.run.' + key, str(run[key])))
    channels = [channel] if channel is not None else sample.get('channel', [])
    for i, c in enumerate(channels):
        i = next((j for j, item in enumerate(sample.get('channel', [])) if item is c), i)
        for j, characteristic in enumerate(c.get('characteristics', [])):
            tag = re.sub(r'[\s_-]+', '_', str(characteristic.get('name') or characteristic.get('tag', '')).strip().casefold())
            path = f'{prefix}.channel[{i}].characteristics[{j}].value'
            value = str(characteristic.get('value') or '')
            if tag in {'assay', 'assay_type', 'library_type', 'library_name', 'technology'}:
                levels[0].append((path, value))
            if tag in {'singlecell_type', 'chemistry', 'library_chemistry'} and value.strip().casefold() in CHEMISTRY_IDENTIFIERS:
                levels[0].append((path, value))
                structured_paths.add(path)
    for key in ('title', 'library_name'):
        if sample.get(key):
            levels[1].append((prefix + '.' + key, str(sample[key])))
    if sample.get('description'):
        levels[2].append((prefix + '.description', str(sample['description'])))
    for group in levels:
        evidence = tuple((path, text, ('10x',) if path in structured_paths else tuple(sorted(methods(text))))
                         for path, text in group if path in structured_paths or methods(text))
        values = {m for _, _, values in evidence for m in values}
        if values:
            if 'bulk' in values and any(single_cell_signal(text) for _, text in group):
                values.add('single_cell')
            if 'droplet' in values and values & {'dropseq', '10x'}:
                values.remove('droplet')
            if 'not_10x' in values and values & {'dropseq', 'plate', 'droplet'} and '10x' not in values:
                values.remove('not_10x')
            return (next(iter(values)) if len(values) == 1 else 'ambiguous'), evidence
        if any(single_cell_signal(text) for _, text in group):
            return 'single_cell', tuple((path, text, ('single_cell',)) for path, text in group if single_cell_signal(text))
    return None, ()
