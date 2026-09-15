# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Narrow source-consistency diagnostics; never repair experimental evidence."""
from ..magetab.preparation import incompatible_preparation, preparation_operation


@preparation_operation
def consistency_issues(data):
    issues = []
    from ..miniml.archive_results import normalize_index_companions
    normalize_index_companions(data, issues=issues, check_only=True)
    for sample in data.get('sample', []):
        for run in sample.get('sra_run', []):
            conflicts = incompatible_preparation(sample, run=run, series=data.get('series'))
            if conflicts:
                issues.append(f"{sample['iid']}/{run.get('experiment')}/{run.get('run')}: "
                    f"library_strategy={run.get('library_strategy')!r}, library_source={run.get('library_source')!r} "
                    f"conflicts with RNA-specific preparation at {', '.join(p for p,_ in conflicts)}; source statements retained.")
        for index, channel in enumerate(sample.get('channel', [])):
            attributes = channel.get('characteristics', [])
            antibodies = {str(a['value']).strip().casefold() for a in attributes
                          if a.get('name', '').strip().casefold() == 'antibody target' and a.get('value')}
            targets = {str(a['value']).strip().casefold() for a in attributes
                       if a.get('name', '').strip().casefold() == 'chip target' and a.get('value')}
            if len(antibodies) == len(targets) == 1 and antibodies != targets:
                issues.append(f"{sample['iid']}/channel[{index}]: explicit antibody target {sorted(antibodies)} "
                              f"differs from explicit ChIP target {sorted(targets)}; source statements retained.")
    for occurrence in data.get('extensions', {}).get('magetab', {}).get('unbound_annotations', []):
        issues.append(f"Unbound {occurrence.get('header')} at {occurrence.get('sdrf')} row "
                      f"{occurrence.get('row_index')}, column {occurrence.get('column_index')}; source occurrence retained.")
    return issues
