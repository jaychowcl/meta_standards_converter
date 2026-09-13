# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Pure parser for Entrez experiment packages and their linked metadata."""
from . import insdc_support as m
from ..sources.archive_support import identifier


class SRAParser:
    def parse(self, records):
        packages = [p for root in records.xml for p in root.findall('.//EXPERIMENT_PACKAGE')]
        packages = [p for p in packages if identifier(p.find('EXPERIMENT/STUDY_REF')) == records.seed.study]
        study = next((p.find('STUDY') for p in packages if identifier(p.find('STUDY')) == records.seed.study), None)
        series = m.study_record(study, records.seed)
        samples, experiments = {}, {}
        for package in packages:
            node = package.find('SAMPLE')
            if node is not None and identifier(node) not in samples:
                samples[identifier(node)] = m.sample_record(node, identifier(node))
            experiment = package.find('EXPERIMENT')
            if experiment is not None:
                experiments[identifier(experiment)] = experiment
        for root in records.xml:
            if root.tag == 'SAMPLE_SET':
                for node in root.findall('SAMPLE'):
                    samples.setdefault(identifier(node), m.sample_record(node, identifier(node)))
        protocols, paths, seen = [], [], set()
        for package in packages:
            experiment = package.find('EXPERIMENT')
            if experiment is None:
                continue
            expt_id = identifier(experiment)
            if identifier(experiment.find('STUDY_REF')) != records.seed.study:
                continue
            protocol = m.protocol_for(experiment)
            if protocol and not any(p['name'] == protocol['name'] for p in protocols):
                protocols.append(protocol)
            for run in package.findall('RUN_SET/RUN'):
                if identifier(run) in seen:
                    continue
                seen.add(identifier(run))
                # Explicit pools may name multiple biological members.
                members = [identifier(n) for n in run.findall('Pool/Member')]
                members = list(dict.fromkeys(v for v in members if v))
                if not members:
                    members = [identifier(n) for n in experiment.findall('DESIGN/SAMPLE_DESCRIPTOR/POOL/*') if identifier(n)]
                if not members:
                    members = [identifier(experiment.find('DESIGN/SAMPLE_DESCRIPTOR'))]
                files = m.files_from_sra(run)
                for member in members:
                    if member not in samples:
                        records.issues.append(f'{identifier(run)}: unresolved sample {member}')
                        continue
                    sample = samples[member]
                    value = m.attach_run(sample, experiment, run, records.seed.study, files)
                    paths.extend(m.assay_paths(sample, value, experiment, files, protocol))
        return m.finish('sra', records, series, list(samples.values()), protocols, paths)
