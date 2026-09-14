# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
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
        _indexed_entities(records, entities)
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
            if identifier(experiment.find('STUDY_REF')) not in (records.seed.study, records.seed.primary):
                records.issues.append(f'{run_id}: experiment outside resolved study')
                continue
            ref = identifier(experiment.find('DESIGN/SAMPLE_DESCRIPTOR'))
            members = [identifier(n) for n in experiment.findall('DESIGN/SAMPLE_DESCRIPTOR/POOL/*')]
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
                indexed = [r for r in records.indexed.get('read_experiment', [])
                           if r.get('experiment_accession') == value['experiment']] + rows.get(run_id, [])
                for key in ('instrument_model', 'instrument_platform'):
                    candidates = {r[key] for r in indexed if r.get(key)}
                    if not value.get(key) and len(candidates) == 1:
                        value[key] = candidates.pop()
                statistics = {k: row[k] for row in rows.get(run_id, [])
                              for k in ('read_count', 'base_count') if row.get(k) not in (None, '')}
                if statistics:
                    value['indexed_statistics'] = statistics
                paths.extend(m.assay_paths(sample, value, experiment, files, protocol))
        return m.finish('ena', records, series, list(samples.values()), protocols, paths)


def _indexed_entities(records, entities):
    """Use indexed fields when Browser records failed; retain the original rows separately."""
    import xml.etree.ElementTree as ET
    def put(node, path, value):
        if value in (None, ''):
            return
        for part in path.split('/'):
            child = node.find(part)
            if child is None:
                child = ET.SubElement(node, part)
            node = child
        node.text = str(value)
    for kind, tag, key in [('sample', 'SAMPLE', 'sample_accession'),
                           ('read_experiment', 'EXPERIMENT', 'experiment_accession'),
                           ('read_run', 'RUN', 'run_accession')]:
        for row in records.indexed.get(kind, []):
            acc = row.get(key)
            if not acc:
                continue
            present = entities.setdefault(tag, {})
            if acc in present or any(acc == i.text for n in present.values() for i in n.findall('IDENTIFIERS/*')):
                continue
            node = ET.Element(tag, accession=acc)
            put(node, 'IDENTIFIERS/PRIMARY_ID', acc)
            put(node, 'TITLE', row.get(kind.removeprefix('read_') + '_title'))
            if kind == 'sample':
                put(node, 'IDENTIFIERS/SECONDARY_ID', row.get('secondary_sample_accession'))
                put(node, 'SAMPLE_NAME/SCIENTIFIC_NAME', row.get('scientific_name'))
                put(node, 'SAMPLE_NAME/TAXON_ID', row.get('tax_id'))
                put(node, 'DESCRIPTION', row.get('sample_description'))
                for field in ('isolation_source', 'host', 'host_tax_id', 'sex', 'tissue_type', 'cell_type', 'disease', 'developmental_stage'):
                    if row.get(field):
                        attrs = node.find('SAMPLE_ATTRIBUTES')
                        if attrs is None:
                            attrs = ET.SubElement(node, 'SAMPLE_ATTRIBUTES')
                        attr = ET.SubElement(attrs, 'SAMPLE_ATTRIBUTE')
                        put(attr, 'TAG', field); put(attr, 'VALUE', row[field])
            elif kind == 'read_experiment':
                study = row.get('secondary_study_accession') or row.get('study_accession')
                if study:
                    ET.SubElement(node, 'STUDY_REF', accession=study)
                design = ET.SubElement(node, 'DESIGN')
                sample = row.get('sample_accession') or row.get('secondary_sample_accession')
                if sample:
                    ET.SubElement(design, 'SAMPLE_DESCRIPTOR', accession=sample)
                for field in ('library_strategy', 'library_source', 'library_selection'):
                    put(design, 'LIBRARY_DESCRIPTOR/' + field.upper(), row.get(field))
                layout = row.get('library_layout')
                if layout in ('SINGLE', 'PAIRED'):
                    put(design, 'LIBRARY_DESCRIPTOR/LIBRARY_LAYOUT/' + layout, 'indexed')
                put(node, 'PLATFORM/INDEXED/INSTRUMENT_MODEL', row.get('instrument_model'))
                put(design, 'LIBRARY_DESCRIPTOR/LIBRARY_CONSTRUCTION_PROTOCOL', row.get('library_construction_protocol'))
            elif row.get('experiment_accession'):
                ET.SubElement(node, 'EXPERIMENT_REF', accession=row['experiment_accession'])
            present[acc] = node
            records.issues.append(f'{acc}: using indexed metadata without full Browser record')
