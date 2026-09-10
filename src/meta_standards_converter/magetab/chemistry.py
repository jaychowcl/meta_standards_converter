# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Conservative, network-free interpretation of library preparation evidence.

Matching text is normalized, but every fact retains its original source. Kit
versions never supply a sequencing recipe. Unresolved fields are not rendered.
"""
from __future__ import annotations

from dataclasses import dataclass
import re


@dataclass(frozen=True)
class ChemistryEvidence:
    field: str
    value: str
    path: str
    text: str


@dataclass(frozen=True)
class ChemistryDiagnostic:
    code: str
    field: str
    paths: tuple[str, ...]


@dataclass(frozen=True)
class ChemistryResult:
    manufacturer: str | None
    family: str | None
    versions: tuple[str, ...]
    library_role: str | None
    index_configuration: str | None
    attributes: tuple[tuple[str, str], ...]
    evidence: tuple[ChemistryEvidence, ...]
    diagnostics: tuple[ChemistryDiagnostic, ...]


def _normalize(text: str) -> str:
    text = text.casefold().translate(str.maketrans({c: "'" for c in '′’ʹ‘`´'}))
    text = re.sub(r'\b(three|five)[ -]+prime\b', lambda m: "3'" if m[1] == 'three' else "5'", text)
    return re.sub(r"\b([35])[ -]+prime\b", r"\1'", text)


def _list(value):
    return value if isinstance(value, list) else [value] if value else []


_VENDOR = re.compile(r'\b(?:10[x×](?:\s+genomics)?|chromium)\b')
# A family marker must be part of a library/kit phrase, not a molecular end.
_KIT = re.compile(r"(?:single[- ]cell\s+|10[x×](?:\s+genomics)?\s+)([35])'|\b(?:gem[- ]x\s+)?flex\b|\bmultiome\b")
_VERSION = re.compile(r'\bv\.?\s*(\d+(?:\.\d+)*)\b')
_NON_PREP = re.compile(r'\b(?:not|without|compatible|compatibility|fixation|software|sequencer|recommended)\b')
_ROLES = {'gene expression': r'gene expression|\bgex\b|\bcdna\b', 'vdj': r'v\(d\)j|\bvdj\b', 'feature barcode': r'feature barcod|antibody capture|cell surface protein|crispr'}


# Exact identifiers, reviewed against Cell Ranger's documented chemistry options.
# These identify chemistry only, never a sequencing recipe.
_IDENTIFIERS = {
    **{f'sc3pv{v}': ('3 prime', (str(v),)) for v in range(1, 5)},
    'sc3pv3ht': ('3 prime', ('3.1',)),
    'sc5p-pe': ('5 prime', ()), 'sc5p-r2': ('5 prime', ()),
    'sc5p-pe-v3': ('5 prime', ('3',)), 'sc5p-r2-v3': ('5 prime', ('3',)),
    'sc5pht': ('5 prime', ('2',)),
}


@dataclass(frozen=True)
class _ChemistryCandidate:
    family: str
    versions: tuple[str, ...]
    path: str
    structured: bool = False


def _structured_sources(sample, channel):
    channels = [channel] if channel is not None else _list(sample.get('channel'))
    for i, item in enumerate(channels):
        if not isinstance(item, dict):
            continue
        index = next((j for j, c in enumerate(_list(sample.get('channel'))) if c is item), i)
        for j, characteristic in enumerate(_list(item.get('characteristics'))):
            if not isinstance(characteristic, dict):
                continue
            tag = re.sub(r'[\s_-]+', '_', str(characteristic.get('name') or characteristic.get('tag', '')).strip().casefold())
            if tag in {'singlecell_type', 'chemistry', 'library_chemistry'}:
                path = f"sample[{sample.get('iid') or 'unknown'}].channel[{index}].characteristics[{j}].value"
                yield path, str(characteristic.get('value') or '')


def _roles(text):
    return {role for role, pattern in _ROLES.items() if re.search(pattern, text)}


def _sources(sample, channel, run, series):
    sources = []
    iid = str(sample.get('iid') or 'unknown')
    channels = [channel] if channel is not None else _list(sample.get('channel'))
    for i, item in enumerate(channels):
        if not isinstance(item, dict):
            continue
        original_index = next((j for j,c in enumerate(_list(sample.get('channel'))) if c is item), i)
        for key in ('extract_protocol',):
            if item.get(key):
                sources.append((f'sample[{iid}].channel[{original_index}].{key}', str(item[key])))
    for key in ('title', 'description'):
        if sample.get(key):
            sources.append((f'sample[{iid}].{key}', str(sample[key])))
    if run:
        for key in ('library_construction_protocol', 'library_protocol', 'description'):
            if run.get(key):
                sources.append((f'sample[{iid}].run.{key}', str(run[key])))
    for i, item in enumerate(_list(series)):
        if not isinstance(item, dict):
            continue
        for key in ('overall_design',):
            text = str(item.get(key) or '')
            for sentence in re.split(r'(?<!\d)\.(?!\d)|\n', text):
                if re.search(r'\ball (?:samples|libraries)\b.*\b(?:prepared|generated|constructed)\b', sentence, re.I) or (iid != 'unknown' and re.search(r'(?<!\w)'+re.escape(iid)+r'(?!\w)', sentence)):
                    sources.append((f'series[{i}].{key}', sentence))
    return sources


def resolve_chemistry(sample: dict, channel: dict | None = None, run: dict | None = None,
                      series: dict | list | None = None) -> ChemistryResult:
    """Resolve applicable facts for one sample/channel and optional library/run.

    Shared study text is considered only with an explicit sample accession or
    universal library-preparation statement. An unlabelled read_lengths list is
    deliberately ignored. Supplied imported SDRF attributes remain source data.
    """
    sources = _sources(sample, channel, run, series)
    facts: list[ChemistryEvidence] = []
    diagnostics: list[ChemistryDiagnostic] = []

    def add(field, value, path, text):
        fact = ChemistryEvidence(field, value, path, text)
        if fact not in facts:
            facts.append(fact)

    candidates: list[_ChemistryCandidate] = []
    for path, original in _structured_sources(sample, channel):
        identifier = original.strip().casefold()
        add('identifier', identifier, path, original)
        identity = _IDENTIFIERS.get(identifier)
        if identity is None:
            if identifier and identifier != 'auto':
                diagnostics.append(ChemistryDiagnostic('unknown_identifier', 'identifier', (path,)))
            continue
        family, versions = identity
        candidates.append(_ChemistryCandidate(family, versions, path, structured=True))
        add('manufacturer', '10x Genomics', path, original)
        add('family', family, path, original)
        for version in versions:
            add('version', version, path, original)

    # A library label supplies scope, never a kit version or read recipe.
    selected_roles = _roles(_normalize(str((run or {}).get('library_name') or sample.get('library_name') or sample.get('title') or '')))
    selected_role = next(iter(selected_roles)) if len(selected_roles) == 1 else None
    recipe_segments = []
    for path, original in sources:
        normalized = _normalize(original)
        if '://' in normalized and not re.search(r'\s', normalized):
            continue
        # Explicit preparation scopes end before unrelated fixation clauses.
        clauses = re.split(r'(?<!\d)\.(?!\d)|\n|;|\bor fixed\b|\bbut\b', normalized)
        active_roles = set()
        for clause in clauses:
            clause = clause.strip()
            if not clause or _NON_PREP.search(clause):
                continue
            roles = _roles(clause)
            if re.search(r'(?:libraries|library)\s*:', clause):
                active_roles = roles
            if selected_role and active_roles and selected_role not in active_roles:
                continue
            vendor = bool(_VENDOR.search(clause))
            if re.search(r'\b(?:chromium|10[x×]\s+genomics)\b', clause):
                add('manufacturer', '10x Genomics', path, original)
            matches = list(_KIT.finditer(clause)) if vendor else []
            for i, match in enumerate(matches):
                # Stop before the next named chemistry; don't borrow its version.
                end = matches[i+1].start() if i+1 < len(matches) else len(clause)
                phrase = clause[match.start():end]
                if not re.search(r'\b(?:kit|reagent|protocol|gene expression|library|libraries)\b', phrase):
                    continue
                family = f'{match[1]} prime' if match[1] else ('flex' if 'flex' in match[0] else 'multiome')
                add('manufacturer', '10x Genomics', path, original)
                add('family', family, path, original)
                phrase_versions = tuple(sorted(set(_VERSION.findall(phrase))))
                candidates.append(_ChemistryCandidate(family, phrase_versions, path))
                for version in phrase_versions:
                    add('version', version, path, original)
                for role in sorted(roles):
                    add('library_role', role, path, original)
                if 'dual index' in phrase or 'dual-index' in phrase:
                    add('index_configuration', 'dual', path, original)
                elif 'single index' in phrase or 'single-index' in phrase:
                    add('index_configuration', 'single', path, original)
            recipe_segments.append((path, original, clause, active_roles.copy()))

    # Narrow a compatible alternative within one preparation phrase, not the
    # union of independent, potentially contradictory preparation statements.
    identifiers = [c for c in candidates if c.structured]
    applicable = []
    for candidate in candidates:
        matches = [c for c in identifiers if c.family == candidate.family
                   and c.versions and set(c.versions) <= set(candidate.versions)
                   and ('.channel[' not in candidate.path
                        or c.path.split('.characteristics[')[0] == candidate.path.rsplit('.', 1)[0])]
        if not candidate.structured and len(candidate.versions) > 1 and matches:
            applicable.extend(matches)
        else:
            applicable.append(candidate)

    def values(field):
        if field == 'family':
            return sorted({c.family for c in applicable})
        if field == 'version':
            return sorted({v for c in applicable for v in c.versions})
        return sorted({f.value for f in facts if f.field == field})

    def unique(field):
        options = values(field)
        if len(options) > 1:
            diagnostics.append(ChemistryDiagnostic('ambiguous_'+field.replace(' ', '_'), field, tuple(sorted({f.path for f in facts if f.field == field}))))
        return options[0] if len(options) == 1 else None

    manufacturer = unique('manufacturer')
    family = unique('family')
    versions = tuple(values('version'))
    unique('version')
    library_role = selected_role or unique('library_role')
    # A recipe must describe a specific read. Only transcript-labelled reads
    # become cDNA attributes; cycle counts alone cannot establish read roles.
    if manufacturer:
        for path, original, clause, roles in recipe_segments:
            if roles and not selected_role and len({r for _,_,_,rs in recipe_segments for r in rs}) > 1:
                continue
            for field, size in re.findall(r'\b(cell barcode|umi barcode|cdna read) offset\s*:\s*(\d+)\b', clause):
                add(field+' offset', size, path, original)
            matches = list(re.finditer(r'\b(read\s*[12]|r[12]|i[57])(?:\s+index)?\s*[:=]\s*(\d+)\s*(?:cycles|bp)|\b(\d+)\s*(?:cycles|bp)\s*(read\s*[12]|r[12]|i[57])\b', clause))
            for i, m in enumerate(matches):
                read = re.sub(r'\s+', '', m[1] or m[4]).replace('read','r')
                size = m[2] or m[3]
                tail = clause[m.end():matches[i+1].start() if i+1<len(matches) else len(clause)]
                if read in ('i5','i7'):
                    add(read+'_size', size, path, original)
                if re.search(r'\b(?:transcript|cdna)\b', tail):
                    add('cdna read', read.replace('r','read'), path, original)
                    add('cdna read size', size, path, original)
                barcode = re.search(r'(\d+)\s*bp\s*(?:cell\s+)?barcode', tail)
                umi = re.search(r'(\d+)\s*bp\s*umi', tail)
                if barcode and umi and read in ('r1','r2'):
                    add('cell barcode read', read.replace('r','read'), path, original)
                    add('cell barcode size', barcode[1], path, original)
                    add('umi barcode read', read.replace('r','read'), path, original)
                    add('umi barcode size', umi[1], path, original)
    if values('i5_size') and values('i7_size'):
        for fact in tuple(facts):
            if fact.field in ('i5_size', 'i7_size'):
                add('index_configuration', 'dual', fact.path, fact.text)
    index_configuration = unique('index_configuration')
    attributes = []
    if family in ('3 prime','5 prime'):
        attributes.append(('end bias', family+' tag'))
    if manufacturer:
        attributes.append(('single cell isolation','10x technology'))
        if family:
            construction = f'10x {family}' + (f' v{versions[0]}' if len(versions)==1 else '')
            attributes.append(('library construction', construction))
        else:
            diagnostics.append(ChemistryDiagnostic('unresolved_family','family',tuple(p for p,_ in sources)))
        if not versions:
            diagnostics.append(ChemistryDiagnostic('unresolved_version','version',tuple(p for p,_ in sources)))
    for field in ('cdna read','cdna read size','cdna read offset','cell barcode read','cell barcode size','cell barcode offset','umi barcode read','umi barcode size','umi barcode offset'):
        value=unique(field)
        if value is not None:
            attributes.append((field,value))
    for read,index in (('i7','index1'),('i5','index2')):
        size=unique(read+'_size')
        if size:
            attributes.extend((('sample barcode read',index),('sample barcode size',size)))
    return ChemistryResult(manufacturer,family,versions,library_role,index_configuration,tuple(attributes),tuple(facts),tuple(diagnostics))
