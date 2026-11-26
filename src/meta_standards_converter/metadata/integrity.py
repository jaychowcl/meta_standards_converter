# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Pure semantic checks, separate from advisory scientific interpretation."""
from collections import Counter, defaultdict
import re

from ..converters.unified.contracts import Diagnostic, CurationContext
from .file_resolution import classify, structural_uri, download_assets


def _data(package):
    return package.to_mapping() if hasattr(package, 'to_mapping') else package


def _comment_groups(value):
    if isinstance(value, dict):
        if value.get('comments'):
            group = defaultdict(list)
            for comment in value['comments']:
                group[comment['name'].casefold()].append(str(comment.get('value', '')).casefold())
            yield group
        for key, child in value.items():
            if key != 'comments':
                yield from _comment_groups(child)
    elif isinstance(value, list):
        for child in value:
            yield from _comment_groups(child)


def qualifier_ownership(before, after):
    """Compare supplied annotation owners, retaining repeated occurrences.

    Call only for the same graph before/after a deterministic preparation, not
    across independently curated studies or merges of different source scopes.
    Inserted graph steps do not shift an existing owner's identity.
    """
    def owners(value):
        found = {}
        for pi, path in enumerate(_data(value).get('series', {}).get('assay_paths', [])):
            occurrences = Counter()
            for step in path.get('steps', []):
                identity = (step['kind'], step.get('name', step.get('protocol_ref')), step.get('sample_ref'))
                occurrences[identity] += 1
                for field in ('characteristics', 'factor_values', 'parameter_values'):
                    fields = Counter()
                    for item in step.get(field, []):
                        name = (item.get('name'), item.get('value'))
                        fields[name] += 1
                        found[pi, identity, occurrences[identity], field, name, fields[name]] = item
        return found
    data = _data(before)
    taxa = {}
    for sample in data.get('sample', []):
        channels = sample.get('channel', [])
        if len(channels) == 1:
            for organism in channels[0].get('organism', []):
                if isinstance(organism, dict) and organism.get('taxid'):
                    taxa.setdefault((sample.get('iid'), organism.get('value')), set()).add(str(organism['taxid']))

    def taxon_number(identifier):
        match = re.fullmatch(r'(?:(?:https?://purl\.obolibrary\.org/obo/NCBITaxon_)|NCBITaxon:)?([0-9]+)',
                             str(identifier))
        return match.group(1) if match else None

    def retained(key, original, candidate):
        if candidate == original:
            return True
        if candidate is None or any(candidate.get(k) != v or k not in candidate for k, v in original.items()):
            return False
        _, (kind, _, owner), _, field, (name, label), _ = key
        return (kind in {'sample', 'source'} and field == 'characteristics'
                and str(name).casefold() == 'organism'
                and candidate.keys() - original.keys() <= {'term_source_ref', 'term_accession_number'}
                and candidate.get('term_source_ref') == 'NCBITaxon'
                and taxa.get((owner, label)) == {taxon_number(candidate.get('term_accession_number'))})

    original, generated = owners(before), owners(after)
    # Dict iteration preserves input order. Compare each original occurrence's
    # ordinal among its surviving siblings while allowing genuinely new rows.
    last_ordinal = {}
    reordered = set()
    ordinal = {key: i for i, key in enumerate(original)}
    for key in generated:
        if key not in original:
            continue
        owner = key[:4]
        if ordinal[key] < last_ordinal.get(owner, -1):
            reordered.add(key)
        last_ordinal[owner] = max(last_ordinal.get(owner, -1), ordinal[key])

    return [Diagnostic('integrity_qualifier_ownership',
        'A supplied named-value occurrence or its qualifiers changed owner or value.',
        'preparation', 'error', location='/series/assay_paths/' + str(key[0]))
        for key, value in original.items() if key in reordered or not retained(key, value, generated.get(key))]


def assess_integrity(package):
    """Detect provable collisions; report source ambiguities without invalidating."""
    data = _data(package)
    issues, identities, files = [], defaultdict(set), defaultdict(list)
    samples = {s['iid']: s for s in data.get('sample', [])}
    from ..magetab.technology import series_identity
    dataset = series_identity(data) or 'input'
    def issue(code, message, location, error=False):
        issues.append(Diagnostic(('integrity_' if error else 'curation_') + code,
                                 message, 'preparation' if error else 'curation',
                                 'error' if error else 'warning', location=location,
                                 action=('Correct the generated workflow binding without changing source evidence.'
                                         if error else 'Review the supplied evidence and entity scope; retain unresolved assertions separately.'),
                                 curation=None if error else CurationContext(
                                     dataset, 'workflow', dataset, code, 'supplied_evidence_requires_review')))
    for pi, path in enumerate(data.get('series', {}).get('assay_paths', [])):
        generated = any(c.get('name') == 'msc_workflow_origin' for c in path.get('comments', []))
        groups = list(_comment_groups(path))
        assignments = {role: {value for group in groups for value in group.get(role.replace('_', ' ') + ' read', []) if value}
                       for role in ('cell_barcode', 'sample_barcode', 'umi_barcode')}
        for group in groups:
            for prefix in ('', 'fastq_'):
                for role, index in zip(group.get(prefix + 'read_type', []), group.get(prefix + 'read_index', [])):
                    expected = assignments.get(role, set())
                    if len(expected) == 1 and index and index not in expected:
                        issue('barcode_assignment',
                              'File read role/index contradicts the supplied library barcode assignment; review the applicable chemistry without inferring a swap.',
                              f'/series/assay_paths/{pi}')
        for si, step in enumerate(path.get('steps', [])):
            location = f'/series/assay_paths/{pi}/steps/{si}'
            comments = defaultdict(list)
            for c in step.get('comments', []):
                comments[c.get('name', '').casefold()].append(str(c.get('value', '')))
            if step['kind'] in {'source', 'sample', 'extract', 'labeled_extract'}:
                marker = comments.get('msc_channel')
                if marker:
                    identities[step['kind'], step.get('name'), step.get('sample_ref')].add(tuple(marker))
            for name, values in comments.items():
                if 'barcode' in name and ('read' in name or 'index' in name) and len(set(values) - {''}) > 1:
                    issue('barcode_assignment', 'Conflicting supplied barcode/read assignments require chemistry-specific review.', location)
            link = step.get('link', {})
            uri = structural_uri(link.get('value'))
            if uri:
                digest = comments.get('md5', comments.get('fastq_md5', []))
                files[uri].extend(digest)
            sample = samples.get(step.get('sample_ref'))
            if sample and step['kind'] in {'source', 'sample'} and len(sample.get('channel', [])) == 1:
                expected = {str(o.get('taxid')) for o in sample['channel'][0].get('organism', []) if o.get('taxid')}
                supplied = {re.sub(r'^(?:NCBITaxon[:_]|taxon:)', '', str(v.get('term_accession_number')).rsplit('/', 1)[-1], flags=re.I)
                            for v in step.get('characteristics', [])
                            if v.get('name', '').casefold() == 'organism' and v.get('term_accession_number')}
                supplied = {value for value in supplied if value.isdigit()}
                if expected and supplied and expected != supplied:
                    issue('taxonomy_assertions', 'Workflow and selected sample taxonomy differ; retain alternative assertions and review their scope.',
                          location, error=generated)
        if generated and sum(s['kind'] == 'protocol_application' for s in path['steps']) > 1:
            issue('protocol_chronology', 'Generated workflow uses provider field order; explicit chronology remains in the supplied protocol text.',
                  f'/series/assay_paths/{pi}')
    for identity, owners in identities.items():
        if len(owners) > 1:
            issue('material_collision', f'Material identity {identity[1]} is shared by distinct channel occurrences.',
                  '/series/assay_paths', error=True)
    for uri, digests in files.items():
        if len(set(digests) - {''}) > 1:
            issue('file_conflict', 'Different checksums are supplied for the same location; assets must remain separate.',
                  '/series/assay_paths')
    for i, sample in enumerate(data.get('sample', [])):
        for field in ('raw_data', 'supplementary_data'):
            for j, record in enumerate(sample.get(field, [])):
                if classify(record) == 'supplementary' and (record.get('value') or record.get('uri')):
                    issue('file_role', 'File link retained with unknown role; raw/derived status requires evidence.',
                          f'/sample/{i}/{field}/{j}')
        for j, run in enumerate(sample.get('sra_run', [])):
            for record in run.get('fastq_files', []):
                if str(run.get('library_layout', '')).upper() == 'SINGLE' and record.get('read_role') == 'paired':
                    issue('read_layout', 'Paired-file role conflicts with supplied single-read library layout.',
                          f'/sample/{i}/sra_run/{j}')
    # Repeated file branches should not repeat the identical advisory.
    unique = {}
    for diagnostic in issues:
        unique.setdefault((diagnostic.code, diagnostic.message, diagnostic.location), diagnostic)
    return list(unique.values())


def entity_counts(package):
    """Structural counts only; sample records cannot establish biological units."""
    data = _data(package)
    samples = data.get('sample', [])
    assays, runs, records = set(), set(), []
    for path in data.get('series', {}).get('assay_paths', []):
        for step in path.get('steps', []):
            if step['kind'] in {'assay', 'hybridization'}:
                assays.add(step.get('name'))
            if step['kind'] == 'scan' and re.fullmatch(r'[SED]RR\d+', step.get('name', '')):
                runs.add(step['name'])
    for sample in samples:
        for run in sample.get('sra_run', []):
            if run.get('run'):
                runs.add(run['run'])
        for field in ('raw_data', 'supplementary_data'):
            records.extend(sample.get(field, []))
    return {'biological_units': None, 'biological_units_reason': 'Requires independent biological identity evidence',
            'sample_records': len(samples),
            'materials': sum(max(1, len(s.get('channel', []))) for s in samples),
            'assays': len(assays) if assays else None, 'runs': len(runs),
            'file_source_occurrences': len(records), 'resolved_file_assets': len(download_assets(records))}
