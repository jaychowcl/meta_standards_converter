"""Pure parser for ENA Browser objects and indexed Portal records."""
from . import insdc_support as m
from ..sources.archive_support import identifier


class ENAParser:
    def parse(self, records):
        entities = {}
        for root in records.xml:
            for node in root:
                if node.tag in ('STUDY', 'PROJECT', 'SAMPLE', 'EXPERIMENT', 'RUN'):
                    entities.setdefault(node.tag, {})[identifier(node)] = node
        study = entities.get('STUDY', {}).get(records.seed.study)
        if study is None:
            study = entities.get('PROJECT', {}).get(records.seed.primary)
        series = m.study_record(study, records.seed)
        samples, aliases = {}, {}
        for accession, node in entities.get('SAMPLE', {}).items():
            ids = m.accessions(node, accession)
            primary = next((a['value'] for a in ids if a['database'] == 'BioSample'), accession)
            for a in ids:
                aliases[a['value']] = primary
            samples[primary] = m.sample_record(node, primary)
        rows = {}
        for row in records.indexed.get('read_run', []):
            rows.setdefault(row.get('run_accession'), []).append(row)
        protocols, paths = [], []
        for run_id, run in entities.get('RUN', {}).items():
            experiment = entities.get('EXPERIMENT', {}).get(identifier(run.find('EXPERIMENT_REF')))
            if experiment is None:
                records.issues.append(f'{run_id}: unresolved experiment')
                continue
            ref = identifier(experiment.find('DESIGN/SAMPLE_DESCRIPTOR'))
            members = [identifier(n) for n in experiment.findall('DESIGN/SAMPLE_DESCRIPTOR/POOL/MEMBER')]
            members = list(dict.fromkeys(v for v in members if v)) or [ref]
            files = [f for row in rows.get(run_id, []) for f in m.files_from_ena(row)]
            if not files:
                files = m.files_from_sra(run)
            protocol = m.protocol_for(experiment)
            if protocol and not any(p['name'] == protocol['name'] for p in protocols):
                protocols.append(protocol)
            for member in members:
                sample = samples.get(aliases.get(member, member))
                if sample is None:
                    records.issues.append(f'{run_id}: unresolved sample {member}')
                    continue
                value = m.attach_run(sample, experiment, run, records.seed.study, files)
                paths.extend(m.assay_paths(sample, value, experiment, files, protocol))
        return m.finish('ena', records, series, list(samples.values()), protocols, paths)
