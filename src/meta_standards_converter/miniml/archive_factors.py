# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Explicit experiment factors, scoped to verified native assay identities."""
from copy import deepcopy
import json
import logging
import re
from .factor_terms import normalized_factor_category


def declared_factor(name, value, unit=None):
    prefix = 'Experimental Factor:'
    if not isinstance(name, str) or not name.startswith(prefix):
        return None
    label = name[len(prefix):].strip()
    if not label or value is None or not str(value).strip():
        return None
    result = {'name': label, 'value': str(value)}
    if unit is not None and str(unit).strip():
        result['unit'] = {'value': str(unit)}
    return result


class ExperimentBindings:
    """Operation-local index shared by factor restoration and residual matching."""

    def __init__(self, data):
        self.members = {}
        self.runs = {}
        for sample in data.get('sample', []):
            for run in sample.get('sra_run', []):
                experiment = run.get('experiment')
                self.members.setdefault(experiment, set()).add(sample['iid'])
                self.runs.setdefault((experiment, run.get('run')), set()).add(sample['iid'])

    def experiment(self, path):
        steps = path.get('steps', [])
        assays = [s for s in steps if s.get('kind') == 'assay']
        if len(assays) != 1:
            return None
        accession = assays[0].get('name')
        refs = {s['sample_ref'] for s in steps if s.get('sample_ref')}
        if not refs or not refs <= self.members.get(accession, set()):
            return None
        for scan in (s for s in steps if s.get('kind') == 'scan'):
            # A descriptive alias requires an explicit run comment. Every supplied
            # run identity must agree, including a conflicting accession as name.
            identities = {str(c.get('value') or '') for c in scan.get('comments', [])
                          if c.get('name') in {'ENA_RUN', 'SRA_RUN', 'DRA_RUN'}}
            name = str(scan.get('name') or '')
            if re.fullmatch(r'[SED]RR\d+', name):
                identities.add(name)
            if not identities or any(not refs <= self.runs.get((accession, run), set()) for run in identities):
                return None
        return accession


def restore_native_factors(data, records=None):
    """Complete missing assay factor groups; never infer factors from characteristics."""
    provider = data.get('source', {}).get('format', '').lower()
    if provider not in {'ena', 'sra'}:
        return
    from ..metadata.archive_enrichment import informative
    from .archive_residuals import children, child_text
    records = records if records is not None else data.get('extensions', {}).get('insdc', {}).get('records', [])
    candidates = {}
    for record in records:
        # Standalone entities avoid counting the same experiment's package wrapper.
        if record.get('kind') != 'EXPERIMENT' or record.get('provider') not in {'ena', 'sra'}:
            continue
        metadata = record.get('metadata', {})
        groups = [declared_factor(child_text(a, 'TAG'), child_text(a, 'VALUE'), child_text(a, 'UNITS'))
                  for container in children(metadata, 'EXPERIMENT_ATTRIBUTES')
                  for a in children(container, 'EXPERIMENT_ATTRIBUTE')]
        groups = [g for g in groups if g is not None]
        if groups:
            variants = candidates.setdefault(record['accession'], {}).setdefault(record['provider'], [])
            if groups not in variants:
                variants.append(groups)
    accepted = {}
    for accession, sources in candidates.items():
        variants = sources.get(provider) or [v for source in sources.values() for v in source]
        unique = {json.dumps(v, sort_keys=True): v for v in variants}
        if len(unique) == 1:
            accepted[accession] = next(iter(unique.values()))
        else:
            logging.getLogger(__name__).warning('%s: conflicting explicit experiment-factor records remain residual', accession)
    bindings = ExperimentBindings(data)
    for path in data.get('series', {}).get('assay_paths', []):
        accession = bindings.experiment(path)
        if accession not in accepted:
            continue
        steps = path['steps']
        assay = next(s for s in steps if s['kind'] == 'assay')
        groups = {}
        for value in accepted.get(accession, []):
            groups.setdefault(value['name'], []).append(value)
        for name, values in groups.items():
            existing = [v for step in steps for v in step.get('factor_values', []) if v['name'] == name]
            if existing and (informative(existing) or not informative(values)):
                continue
            for step in steps:
                if 'factor_values' in step:
                    step['factor_values'] = [v for v in step['factor_values'] if v['name'] != name]
            assay.setdefault('factor_values', []).extend(deepcopy(values))
    series = data.get('series', {})
    declared = {v.get('name') or v.get('factor') for v in series.get('variable', [])}
    for path in series.get('assay_paths', []):
        for step in path['steps']:
            for value in step.get('factor_values', []):
                name = value['name']
                if name not in declared:
                    series.setdefault('variable', []).append({'name': name, 'factor': normalized_factor_category(name)})
                    declared.add(name)
