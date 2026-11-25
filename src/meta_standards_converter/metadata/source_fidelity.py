# =============================================================================
# Authors
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Project supplied identifiers only through explicit graph ownership, offline."""
from copy import deepcopy
import re


def organism_annotations(organism, warnings):
    """Preserve occurrence-local qualifiers and an explicit conflicting taxid."""
    record = {'name': 'organism', **{k: deepcopy(v) for k, v in organism.items()
        if k in {'value', 'term_source_ref', 'term_accession_number', 'unit'}}}
    taxid = str(organism.get('taxid') or '').strip()
    if not taxid:
        return [record]
    explicit = str(record.get('term_accession_number') or '')
    if not explicit and not record.get('term_source_ref'):
        return [{**record, 'term_source_ref': 'NCBITaxon', 'term_accession_number': taxid}]
    if (record.get('term_source_ref') == 'NCBITaxon' and explicit in {taxid, f'NCBITaxon:{taxid}', f'NCBITaxon_{taxid}',
        f'http://purl.obolibrary.org/obo/NCBITaxon_{taxid}',
        f'https://purl.obolibrary.org/obo/NCBITaxon_{taxid}'}):
        return [record]
    warnings.append(f"Organism {organism.get('value')!r} has conflicting taxonomy qualifiers and taxid {taxid}; both retained.")
    return [record, {'name': 'organism', 'value': organism.get('value', ''),
        'term_source_ref': 'NCBITaxon', 'term_accession_number': taxid}]


def project_supplied_identifiers(data):
    """Add sample experiment evidence and taxids without fabricating run joins."""
    samples = {s['iid']: s for s in data.get('sample', []) if s.get('iid')}
    warnings = []
    for path in data.get('series', {}).get('assay_paths', []):
        steps = path['steps']
        owners = {s['sample_ref'] for s in steps if s.get('sample_ref')}
        for node in steps:
            sample = samples.get(node.get('sample_ref'))
            if not sample or node['kind'] not in {'source', 'sample'}:
                continue
            channels = sample.get('channel', [])
            if len(channels) == 1:
                characteristics = node.setdefault('characteristics', [])
                for organism in channels[0].get('organism', []):
                    if not isinstance(organism, dict) or not organism.get('value') or not organism.get('taxid'):
                        continue
                    for annotation in organism_annotations(organism, warnings):
                        same = [a for a in characteristics if str(a.get('name', '')).casefold() == 'organism' and a.get('value') == annotation.get('value')]
                        if any({**a, 'name': 'organism'} == annotation for a in same):
                            continue
                        bare = next((a for a in same if not a.get('term_source_ref') and not a.get('term_accession_number')), None)
                        if bare is not None and annotation.get('term_source_ref'):
                            bare.update({k: v for k, v in annotation.items() if k != 'name'})
                        elif not same or annotation.get('term_source_ref'):
                            characteristics.append(annotation)
            relations = []
            for relation in sample.get('relation', []):
                target = str(relation.get('target') or '')
                if str(relation.get('type', '')).upper() != 'SRA':
                    continue
                matches = re.findall(r'(?<![A-Z0-9])(?:SRX|ERX|DRX)\d+(?![A-Z0-9])', target, re.I)
                for identifier in matches:
                    if identifier.upper() not in relations:
                        relations.append(identifier.upper())
            if not relations:
                continue
            assays = [s for s in steps if s['kind'] == 'assay']
            known = {c['value'] for s in assays for c in s.get('comments', []) if c['name'] == 'ENA_EXPERIMENT'}
            known.update(s['name'] for s in assays if re.fullmatch(r'(SRX|ERX|DRX)\d+', s['name']))
            if len(relations) > 1 or (known and not set(relations) <= known):
                warnings.append(f"Sample {sample['iid']} has multiple or conflicting experiment relations; sample evidence retained without new run/file associations.")
            # Only one explicitly owned assay without a competing experiment can
            # receive an identifier. Ambiguous relations remain on their source.
            run_or_file = any(s['kind'].endswith('file') or (s['kind'] == 'scan' and re.fullmatch(r'(SRR|ERR|DRR)\d+', s['name'])) or s.get('link') for s in steps)
            unambiguous = (not run_or_file or bool(known)) and len(owners) == 1 and len(assays) == 1 and len(relations) == 1 and (not known or known == set(relations))
            target = assays[0] if unambiguous else node
            label = 'ENA_EXPERIMENT' if unambiguous else 'SRA_EXPERIMENT_RELATION'
            comments = target.setdefault('comments', [])
            for identifier in relations:
                comment = {'name': label, 'value': identifier}
                if comment not in comments:
                    comments.append(comment)
    # Only the taxonomy namespace introduced by this projection is declared.
    # Authored unresolved/conflicting namespace declarations remain untouched.
    has_taxid = any(o.get('taxid') for sample in samples.values()
        for channel in sample.get('channel', []) for o in channel.get('organism', []) if isinstance(o, dict))
    if has_taxid and 'NCBITaxon' not in {d.get('iid') for d in data.get('database', [])}:
        data.setdefault('database', []).append({'iid': 'NCBITaxon', 'name': 'NCBITaxon'})
    if warnings:
        stored = data.setdefault('extensions', {}).setdefault('msc_workflow_preparation', {}).setdefault('warnings', [])
        stored.extend(w for w in dict.fromkeys(warnings) if w not in stored)
