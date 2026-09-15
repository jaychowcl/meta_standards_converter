# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Small MINiML mapping primitives shared by the two native XML parsers."""
from copy import deepcopy
from pathlib import PurePosixPath
import re
import xml.etree.ElementTree as ET
from urllib.parse import urlsplit, quote, unquote

from .codec import MINiMLCodec
from ..sources.archive_support import identifier


def text(node, path, default=None):
    if node is None:
        return default
    item = node.find(path)
    if item is None:
        return default
    value = ''.join(item.itertext()).strip()
    return value or default


def tree(node):
    result = {'tag': node.tag, 'attributes': dict(node.attrib), 'children': [tree(c) for c in node]}
    if node.text is not None:
        result['text'] = node.text
    if node.tail and node.tail.strip():
        result['tail'] = node.tail
    return result


def retained(provider, records):
    values = []
    for root in records.xml:
        for node in list(root) if root.tag.endswith('_SET') or root.tag in ('BioSampleSet', 'RecordSet', 'PubmedArticleSet') else [root]:
            accession = node.findtext('MedlineCitation/PMID') if node.tag == 'PubmedArticle' else identifier(node) or identifier(node.find('EXPERIMENT'))
            values.append({'provider': provider, 'kind': node.tag, 'accession': accession, 'metadata': tree(node)})
            if node.tag == 'EXPERIMENT_PACKAGE':
                for entity in node.iter():
                    if entity.tag in ('STUDY', 'SAMPLE', 'EXPERIMENT', 'RUN', 'SUBMISSION'):
                        values.append({'provider': provider, 'kind': entity.tag, 'accession': identifier(entity), 'metadata': tree(entity)})
    for kind, rows in records.indexed.items():
        for row in rows:
            values.append({'provider': provider, 'kind': kind, 'accession': row.get(kind.removeprefix('read_') + '_accession') or row.get('accession'), 'metadata': deepcopy(row)})
    values.extend(deepcopy(records.linked))
    return {'version': '1.0', 'records': values}


def database_for(value):
    for prefix, db in [('GSE', 'GEO'), ('GSM', 'GEO'), ('E-', 'ArrayExpress'),
                       ('PRJ', 'BioProject'), ('SAM', 'BioSample'), ('SR', 'SRA'), ('ER', 'ENA'), ('DR', 'DRA')]:
        if value.startswith(prefix):
            return db
    return 'INSDC'


def accessions(node, primary):
    values = [primary] if primary else []
    if node is not None:
        for item in node.findall('IDENTIFIERS/*'):
            value = (item.text or '').strip()
            if re.fullmatch(r'(?:[SED]R[PSXR]\d+|PRJ(?:NA|EB|DB|DA)\d+|SAM(?:N|EA|D)\d+|GS[EM]\d+|E-[A-Z]+-\d+)', value):
                if value not in values:
                    values.append(value)
    return [{'value': v, 'database': database_for(v)} for v in values]


def relations(node):
    from .reference_targets import parse_reference_targets
    values = []
    if node is not None:
        for link in node.findall('.//XREF_LINK'):
            db, target = text(link, 'DB'), text(link, 'ID')
            if db and target:
                values.extend({'type': db, 'target': part} for part in parse_reference_targets(db, target))
        for link in node.findall('.//URL_LINK'):
            target = text(link, 'URL')
            if target:
                values.append({'type': text(link, 'LABEL', 'external'), 'target': target})
    return values


def attributes(node, kind='SAMPLE'):
    result = []
    if node is not None:
        for attr in node.findall(f'{kind}_ATTRIBUTES/{kind}_ATTRIBUTE'):
            name = text(attr, 'TAG')
            value = text(attr, 'VALUE', '')
            if name:
                item = {'name': name, 'value': value}
                if text(attr, 'UNITS'):
                    item['unit'] = {'value': text(attr, 'UNITS')}
                result.append(item)
    return result


def sample_record(node, primary):
    organism = text(node, 'SAMPLE_NAME/SCIENTIFIC_NAME')
    taxid = text(node, 'SAMPLE_NAME/TAXON_ID')
    attrs = attributes(node)
    channel = {'characteristics': attrs}
    explicit = {}
    for attr in attrs:
        explicit.setdefault(attr['name'].casefold(), []).append(attr['value'])
    molecules = set(explicit.get('molecule', []))
    if len(molecules) == 1:
        channel['molecule'] = {'value': next(iter(molecules))}
    host_ids = set(explicit.get('host_taxid', []) + explicit.get('host_tax_id', []))
    if len(host_ids) == 1:
        for attr in attrs:
            if attr['name'].casefold() == 'host':
                attr.update(term_source_ref='NCBITaxon', term_accession_number=next(iter(host_ids)))
    if organism or taxid:
        channel['organism'] = [{'value': organism or '', 'taxid': taxid}]
    statuses = []
    for a in attrs:
        field = {'ENA-FIRST-PUBLIC': 'release_date', 'ENA-LAST-UPDATE': 'last_update_date'}.get(a['name'])
        if field and a['value']:
            statuses.append({field: a['value'], 'database': 'ENA'})
    source = next((a['value'] for a in attrs if a['name'] == 'source_name'), None)
    if source is None:
        source = next((a['value'] for a in attrs if a['name'] == 'isolation_source'), None)
    if source:
        channel['source'] = {'value': source}
    return {'iid': primary, 'accession': accessions(node, primary), 'title': text(node, 'TITLE'),
            'description': text(node, 'DESCRIPTION'), 'type': 'SRA', 'channel_count': '1',
            'channel': [channel], 'status': statuses, 'relation': relations(node), 'sra_run': []}


def library(experiment):
    descriptor = experiment.find('DESIGN/LIBRARY_DESCRIPTOR') if experiment is not None else None
    result = {}
    for name in ('strategy', 'source', 'selection'):
        value = text(descriptor, 'LIBRARY_' + name.upper())
        if value:
            result['library_' + name] = value
    layout = descriptor.find('LIBRARY_LAYOUT') if descriptor is not None else None
    if layout is not None and len(layout):
        result['library_layout'] = layout[0].tag
    instrument = text(experiment, 'PLATFORM/*/INSTRUMENT_MODEL')
    if instrument:
        result['instrument_model'] = instrument
    platform = experiment.find('PLATFORM') if experiment is not None else None
    if platform is not None and len(platform) == 1 and platform[0].tag != 'INDEXED':
        result['instrument_platform'] = platform[0].tag
    return result


def files_from_ena(row):
    files = []
    for family in ('fastq', 'submitted', 'sra', 'bam'):
        columns = {key: str(row.get(family + '_' + key) or '').split(';')
                   for key in ('ftp', 'md5', 'bytes', 'file_role', 'format', 'aspera', 'galaxy')}
        if not any(any(v) for v in columns.values()):
            continue
        for i in range(max(map(len, columns.values()))):
            value = lambda key: columns[key][i] if i < len(columns[key]) else ''
            path = value('ftp')
            # Portal FTP locations are paths: a literal # is not a URI fragment.
            uri = (path if '://' in path else 'ftp://' + path) if path else None
            if uri:
                scheme, location = uri.split('://', 1)
                host, separator, file_path = location.partition('/')
                uri = scheme + '://' + host + separator + quote(unquote(file_path), safe='/')
            item = {'uri': uri, 'format': value('format') or family, 'role': value('file_role')}
            if uri:
                item['filename'] = PurePosixPath(unquote(urlsplit(uri).path)).name
            for key in ('md5', 'bytes', 'aspera', 'galaxy'):
                if value(key):
                    item[key] = value(key)
            files.append(item)
    return files


def _filename_uri(uri, filename):
    """A supplied literal filename can prove that a trailing # is a path byte."""
    if (isinstance(uri, str) and isinstance(filename, str) and '#' in filename
            and '/' not in filename and '://' in uri and uri.endswith('/' + filename)
            and not any(c in uri[:-len(filename)].split('://', 1)[1] for c in '?#')):
        return uri[:-len(filename)] + filename.replace('#', '%23')
    return uri


def files_from_sra(run):
    files = []
    for node in run.findall('SRAFiles/SRAFile'):
        base = {'filename': node.get('filename'), 'bytes': node.get('size'), 'md5': node.get('md5'),
                'format': node.get('semantic_name', ''), 'role': node.get('supertype', '')}
        urls = ([{'url': node.get('url')}] if node.get('url') else []) + [dict(x.attrib) for x in node.findall('Alternatives')]
        for alternative in urls or [{}]:
            files.append({**base, **alternative, 'uri': _filename_uri(alternative.get('url'), base['filename'])})
    # Some partners retain submitted file descriptors instead of SRAFiles.
    for node in run.findall('.//DATA_BLOCK/FILES/FILE'):
        item = {'filename': node.get('filename'), 'format': node.get('filetype', ''), 'role': 'submitted', 'checksum_method': node.get('checksum_method'), 'checksum': node.get('checksum')}
        if node.get('checksum_method', '').upper() == 'MD5':
            item['md5'] = node.get('checksum')
        if '://' in (node.get('filename') or ''):
            item['uri'] = node.get('filename')
        files.append(item)
    return files


def attach_run(sample, experiment, run, study, files):
    expt_id, run_id = identifier(experiment), identifier(run)
    biosample = next((a['value'] for a in sample['accession'] if a['database'] == 'BioSample'), None)
    archive_sample = next((a['value'] for a in sample['accession'] if a['value'].startswith(('SRS', 'ERS', 'DRS'))), sample['iid'])
    record = {'run': run_id, 'study': study, 'experiment': expt_id, 'sample': archive_sample,
              'biosample': biosample, 'scan_name': run.get('alias') or run_id, **library(experiment),
              'read_lengths': [r.get('average') for r in run.findall('Statistics/Read') if r.get('average')]}
    fastqs = [f for f in files if 'fastq' in str(f.get('format', '')).lower()]
    if fastqs:
        record['fastq_files'] = fastqs
    record['files'] = deepcopy(files)
    if run.find('Statistics') is not None:
        record['statistics'] = tree(run.find('Statistics'))
    sample['sra_run'].append(record)
    for f in files:
        if f.get('uri'):
            link = {'value': f['uri'], 'type': f.get('format') or 'raw'}
            if re.fullmatch('[a-fA-F0-9]{32}', f.get('md5') or ''):
                link['checksum'] = f['md5']
            sample.setdefault('raw_data', []).append(link)
    return record


def protocol_for(experiment):
    from .archive_protocols import library_description
    description = library_description(
        text(experiment, 'DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_CONSTRUCTION_PROTOCOL'),
        text(experiment, 'DESIGN/DESIGN_DESCRIPTION'))
    if not description:
        return None
    result = {'name': identifier(experiment) + ':library', 'description': description,
              'type': {'value': 'library construction protocol'}}
    return result


def assay_paths(sample, run, experiment, files, protocol):
    # Repeating attributes remain occurrences rather than a dict keyed by name.
    channel = sample['channel'][0]
    attrs = deepcopy(channel.get('characteristics', []))
    for organism in channel.get('organism', []):
        value = {'name': 'organism', 'value': organism['value']}
        if organism.get('taxid'):
            value.update(term_source_ref='NCBITaxon', term_accession_number=organism['taxid'])
        attrs.append(value)
    source = {'kind': 'source', 'name': sample['iid'], 'sample_ref': sample['iid'], 'characteristics': attrs}
    if channel.get('source'):
        source['description'] = channel['source']['value']
    assay = {'kind': 'assay', 'name': run['experiment'], 'sample_ref': sample['iid'],
             'technology_type': {'value': 'sequencing assay'}, 'comments': []}
    if text(experiment, 'DESIGN/DESIGN_DESCRIPTION'):
        assay['description'] = text(experiment, 'DESIGN/DESIGN_DESCRIPTION')
    for key in ('library_strategy', 'library_source', 'library_selection', 'library_layout', 'instrument_model', 'instrument_platform'):
        if run.get(key):
            assay['comments'].append({'name': key.upper(), 'value': run[key]})
    scan = {'kind': 'scan', 'name': run['run'], 'comments': [
        {'name': 'ENA_RUN', 'value': run['run']}, {'name': 'ENA_EXPERIMENT', 'value': run['experiment']},
        {'name': 'ENA_SAMPLE', 'value': run['sample']}]}
    result = []
    for file in files or [None]:
        steps = [deepcopy(source)]
        if protocol:
            steps.append({'kind': 'protocol_application', 'protocol_ref': protocol['name']})
        steps.extend([deepcopy(assay), deepcopy(scan)])
        if file and (file.get('uri') or file.get('filename')):
            node = {'kind': 'array_data_file', 'name': file.get('filename') or file['uri'], 'comments': []}
            if file.get('uri'):
                node['link'] = {'value': file['uri'], 'type': file.get('format') or 'raw'}
            for key in ('md5', 'bytes', 'format', 'role', 'checksum_method', 'checksum'):
                if file.get(key):
                    node['comments'].append({'name': key.upper(), 'value': str(file[key])})
            steps.append(node)
        result.append({'steps': steps})
    return result


def finish(provider, records, series, samples, protocols, paths):
    fill_linked_metadata(records, series, samples, protocols, paths)
    for sample in samples:
        runs = sample['sra_run']
        sample['sra_accession'] = list(dict.fromkeys(r['experiment'] for r in runs))
        sample['ena_accession'] = list(dict.fromkeys(r['study'] for r in runs))
        for key in ('library_strategy', 'library_source', 'library_selection'):
            values = {r[key] for r in runs if r.get(key)}
            if len(values) == 1 and all(r.get(key) for r in runs):
                sample[key] = values.pop()
    for record in records.linked:
        if record['kind'] == 'accession_range':
            target = next((s for s in samples if record['accession'] in {a['value'] for a in s['accession']}), None)
            if target is not None:
                value = record['metadata']
                target['relation'] = [r for r in target.get('relation', []) if r != {'type':value['database'], 'target':value['literal']}]
                target['relation'].extend({'type':value['database'], 'target':a} for a in value['accessions'])
    by_sample = {s['iid']: s for s in samples}
    for path in paths:
        for step in path['steps']:
            if step.get('kind') == 'source' and step.get('sample_ref') in by_sample:
                channel = by_sample[step['sample_ref']]['channel'][0]
                organism = []
                for value in channel.get('organism', []):
                    organism.append({'name': 'organism', 'value': value['value'], **({'term_source_ref': 'NCBITaxon', 'term_accession_number': value['taxid']} if value.get('taxid') else {})})
                step['characteristics'] = deepcopy(channel['characteristics']) + organism
    for row in records.indexed.get('study', []):
        if row.get('study_accession') not in (records.seed.primary, records.seed.study):
            continue
        status = {key: row[source] for key, source in (('release_date', 'first_public'), ('last_update_date', 'last_updated')) if row.get(source)}
        if status:
            series.setdefault('status', []).append({**status, 'database': provider.upper()})
        for key, source in (('title', 'study_title'), ('summary', 'study_description')):
            if not series.get(key) and row.get(source):
                series[key] = row[source]
    series.update(sample_ref=[{'ref': s['iid']} for s in samples], protocols=protocols, assay_paths=paths)
    data = {'miniml_schema_version': '3.0', 'source': {'format': provider.upper()}, 'series': series,
            'sample': samples, 'extensions': {'insdc': retained(provider, records)}}
    dbs = {a['database'] for entity in [series, *samples] for a in entity.get('accession', [])}
    data['database'] = [{'iid': db, 'name': db} for db in sorted(dbs)]
    from .archive_entities import actors, declare_ontologies
    data['organization'], data['contributor'] = actors(records, provider)
    sample_actors = {o['iid']: o.get('sample_accession') for o in data['organization']}
    if data['contributor']:
        series['contributor_ref'] = [{'ref': c['iid']} for c in data['contributor'] if not sample_actors.get(c.get('organization_ref', {}).get('ref'))]
        for sample in samples:
            aliases = {sample['iid'], *[a['value'] for a in sample.get('accession', [])]}
            refs = [{'ref': c['iid']} for c in data['contributor'] if sample_actors.get(c.get('organization_ref', {}).get('ref')) in aliases]
            if refs: sample['contact_ref'] = refs
    declare_ontologies(data)
    from .archive_paths import complete_native_paths
    complete_native_paths(data)
    from .archive_residuals import finalize
    return finalize(data)


def study_record(node, seed):
    result = {'iid': seed.primary, 'accession': accessions(node, seed.primary),
              'title': text(node, 'DESCRIPTOR/STUDY_TITLE') or text(node, 'TITLE'),
              'summary': text(node, 'DESCRIPTOR/STUDY_ABSTRACT') or text(node, 'DESCRIPTION'),
              'relation': relations(node)}
    if seed.study not in [a['value'] for a in result['accession']]:
        result['accession'].append({'value': seed.study, 'database': database_for(seed.study)})
    kind = node.find('DESCRIPTOR/STUDY_TYPE') if node is not None else None
    if kind is not None and kind.get('existing_study_type'):
        result['type'] = [{'value': kind.get('existing_study_type')}]
    result['pubmed_id'] = [r['target'] for r in result['relation'] if r['type'].lower() == 'pubmed']
    return result


def biosample_characteristic(name, item):
    """Project the declared BioSamples value/unit and singleton ontology group."""
    value = {'name': name, 'value': item.get('text', '')}
    if item.get('unit'):
        value['unit'] = {'value': item['unit']}
    terms = item.get('ontologyTerms', [])
    if len(terms) == 1 and terms[0]:
        value['term_accession_number'] = terms[0]
        if re.fullmatch(r'https?://www\.ebi\.ac\.uk/efo/EFO_\d+', terms[0]):
            value['term_source_ref'] = 'EFO'
    return value


def fill_linked_metadata(records, series, samples, protocols, paths):
    """Project linked fields only when their entity binding is explicit."""
    from .archive_publications import project_publications
    project_publications(records, series, samples)
    by_id = {a['value']: s for s in samples for a in s['accession']}
    seen_contacts = set()
    for root in records.xml:
        for bio in root.findall('.//BioSample'):
            sample = by_id.get(bio.get('accession'))
            if sample is None:
                continue
            channel = sample['channel'][0]
            names = {a['name'] for a in channel['characteristics']}
            for a in bio.findall('Attributes/Attribute'):
                name = a.get('attribute_name') or a.get('display_name')
                if name and name not in names:
                    value = {'name': name, 'value': ''.join(a.itertext())}
                    if a.get('unit'):
                        value['unit'] = {'value': a.get('unit')}
                    channel['characteristics'].append(value)
            sample.setdefault('status', []).append({'database': 'BioSample', **{k: bio.get(v) for k, v in [
                ('submission_date', 'submission_date'), ('release_date', 'publication_date'),
                ('last_update_date', 'last_update')] if bio.get(v)}})
            for name, literal in [('access', bio.get('access')), ('record status', bio.find('Status').get('status') if bio.find('Status') is not None else None)]:
                if literal:
                    sample['status'].append({'database': 'BioSample', 'comment': [{'name':name, 'value':literal}]})
        for project in root.findall('.//Project'):
            desc = project.find('ProjectDescr')
            if not series.get('title'):
                series['title'] = text(desc, 'Title')
            if not series.get('summary'):
                series['summary'] = text(desc, 'Description')
    # EBI BioSamples supplies repeated values and optional ontology URLs.
    for record in records.linked:
        if record['provider'] != 'biosamples':
            continue
        sample = by_id.get(record['accession'])
        if sample is None:
            continue
        if record['metadata'].get('status'):
            sample.setdefault('status', []).append({'database':'BioSamples', 'comment':[{'name':'status','value':record['metadata']['status']}]})
        attrs = sample['channel'][0]['characteristics']
        for name, values in record['metadata'].get('characteristics', {}).items():
            if name.casefold() in ('organism', 'title', 'description'):
                continue
            supplied = [biosample_characteristic(name, item) for item in values]
            existing = [a for a in attrs if a['name'].replace('_', ' ').casefold() == name.replace('_', ' ').casefold()]
            if not existing:
                attrs.extend(supplied)
                continue
            for value in supplied:
                # An occurrence is completed only when both sources identify a
                # unique literal value. Repetitions/conflicts remain residual.
                candidates = [a for a in existing if a.get('value') == value['value']]
                if len(candidates) != 1 or sum(v['value'] == value['value'] for v in supplied) != 1:
                    continue
                current = candidates[0]
                group = {k:v for k,v in value.items() if k not in ('name', 'value')}
                if all(not current.get(k) or current[k] == v for k,v in group.items()):
                    for key, literal in group.items():
                        current.setdefault(key, deepcopy(literal))


    project_results(records, series, samples, protocols, paths)


def project_results(records, series, samples, protocols, paths):
    """Analysis/assembly links use explicit sample/run bindings; directories stay relations."""
    by_sample = {a['value']: s for s in samples for a in s['accession']}
    by_run = {}
    for sample in samples:
        for run in sample['sra_run']:
            by_run.setdefault(run['run'], []).append(sample)
    entries = {}
    for root in records.xml:
        for node in root.iter():
            if node.tag in ('ANALYSIS', 'ASSEMBLY'):
                entries.setdefault(identifier(node), {'rows': [], 'node': None})['node'] = node
    for kind in ('analysis', 'assembly'):
        for row in records.indexed.get(kind, []):
            acc = (row.get('assembly_set_accession') if kind == 'assembly' else None) or row.get(kind + '_accession') or row.get('accession')
            if acc:
                entries.setdefault(acc, {'rows': [], 'node': None})['rows'].append(row)
    for record in records.linked:
        if record['kind'] == 'assembly':
            metadata = record['metadata']
            ids = list(dict.fromkeys([record['accession'], *[metadata.get('synonym', {}).get(k) for k in ('genbank', 'refseq')]]))
            sample = by_sample.get(metadata.get('biosampleaccn'))
            for acc in filter(None, ids):
                for target in ([series, sample] if sample else [series]):
                    target.setdefault('relation', []).append({'type': 'assembly', 'target': acc})
            for field in ('ftppath_genbank', 'ftppath_refseq'):
                if metadata.get(field): series['relation'].append({'type': 'assembly directory', 'target': metadata[field]})
            for field in ('ftppath_stats_rpt', 'ftppath_regions_rpt', 'ftppath_assembly_rpt'):
                if metadata.get(field):
                    result_file(sample or series, {'uri': metadata[field], 'format': 'assembly report',
                                                   'assembly_accession': record['accession']}, paths)
    for acc, entry in entries.items():
        node, rows = entry['node'], entry['rows']
        sample_refs, run_refs, files = [], [], []
        study_ref = identifier(node.find('STUDY_REF')) if node is not None else None
        if study_ref and study_ref not in (records.seed.study, records.seed.primary) and not (node.findall('SAMPLE_REF') or node.findall('RUN_REF')):
            continue
        if node is not None:
            sample_refs.extend(identifier(n) for n in node.findall('SAMPLE_REF') if identifier(n))
            run_refs.extend(identifier(n) for n in node.findall('RUN_REF') if identifier(n))
            for f in node.findall('.//FILES/FILE'):
                uri = f.get('filename', '')
                if urlsplit(uri).scheme:
                    files.append({'uri': uri, 'format': f.get('filetype', 'analysis'),
                        'checksum_method': f.get('checksum_method'), 'checksum': f.get('checksum'),
                        'role': f.get('role'), 'bytes': f.get('bytes')})
        for row in rows:
            sample_refs.extend(v for v in (row.get('sample_accession') or '').split(';') if v)
            run_refs.extend(v for v in (row.get('run_accession') or '').split(';') if v)
            files.extend(files_from_ena(row))
        selected = {by_sample[r]['iid']: by_sample[r] for r in sample_refs if r in by_sample}
        run_selected = {s['iid']: s for r in run_refs for s in by_run.get(r, [])}
        if sample_refs and run_refs:
            selected = {k: s for k, s in selected.items() if k in run_selected}
        elif run_refs:
            selected = run_selected
        targets = list(selected.values()) if sample_refs or run_refs else [series]
        series.setdefault('relation', []).append({'type': 'analysis' if acc.startswith(('ERZ', 'SRZ', 'DRZ')) else 'assembly', 'target': acc})
        for target in targets:
            if target is not series:
                target.setdefault('relation', []).append({'type': 'analysis' if acc.startswith(('ERZ','SRZ','DRZ')) else 'assembly', 'target': acc})
        description = text(node, 'PROTOCOL') or text(node, 'METHOD') or text(node, 'ANALYSIS_TYPE/SEQUENCE_ASSEMBLY/ASSEMBLY_METHOD')
        software = [n.text for n in node.findall('.//PROGRAM') if n.text] if node is not None else []
        protocol = None
        if description:
            protocol = {'name': acc + ':analysis', 'type': {'value': 'data analysis protocol'},
                        'description': description, 'software': software}
            protocols.append(protocol)
        for f in files:
            if not f.get('uri'):
                continue
            f['analysis_accession' if acc.startswith(('ERZ','SRZ','DRZ')) else 'assembly_accession'] = acc
            link = {'value': f['uri'], 'type': f.get('format') or 'analysis', **{k:f[k] for k in ('role','bytes','aspera','galaxy','checksum_method') if f.get(k)}}
            link.update({k:f[k] for k in ('analysis_accession','assembly_accession') if f.get(k)})
            if f.get('checksum') and (f.get('checksum_method') or '').upper() != 'MD5': link['file_checksum'] = f['checksum']
            checksum = f.get('md5') or (f.get('checksum') if (f.get('checksum_method') or '').upper() == 'MD5' else None)
            if checksum:
                link['checksum'] = checksum
            for target in targets:
                target.setdefault('supplementary_data', []).append(deepcopy(link))
                if target is series:
                    continue
                result_file(target, f, paths, run_refs, protocol, add_link=False)


def result_file(target, file, paths, run_refs=(), protocol=None, *, add_link=True):
    """A sample result without run evidence is a source-to-file branch."""
    link = {'value': file['uri'], 'type': file.get('format') or 'analysis', **{k:file[k] for k in ('role','bytes','aspera','galaxy','checksum_method','assembly_accession','analysis_accession') if file.get(k)}}
    if file.get('checksum'):
        link['checksum' if (file.get('checksum_method') or '').upper() == 'MD5' else 'file_checksum'] = file['checksum']
    if file.get('md5'): link['checksum'] = file['md5']
    if add_link: target.setdefault('supplementary_data', []).append(deepcopy(link))
    if file.get('format') == 'assembly report':
        return  # Supporting documentation has no asserted experimental processing edge.
    if 'channel' not in target:
        return
    bases = [p for p in paths if any(s.get('sample_ref') == target['iid'] for s in p['steps'])
             and not any(s.get('kind', '').startswith('derived_') for s in p['steps'])
             and (not run_refs or any(s.get('kind') == 'scan' and s.get('name') in run_refs for s in p['steps']))]
    if not bases:
        return
    if run_refs:
        branches = []
        for base in bases:
            end = next((i for i,s in enumerate(base['steps']) if s.get('kind') == 'array_data_file'), len(base['steps']))
            prefix = base['steps'][:end]
            if prefix not in branches:
                branches.append(prefix)
    else:
        branches = [[s for s in bases[0]['steps'] if s['kind'] == 'source']]
    for base in branches:
        steps = deepcopy(base)
        if protocol: steps.append({'kind':'protocol_application', 'protocol_ref':protocol['name']})
        steps.append({'kind':'derived_array_data_file', 'name':unquote(PurePosixPath(urlsplit(file['uri']).path).name), 'link':deepcopy(link),
                      'comments':[{'name':k.upper(),'value':str(file[k])} for k in ('format','bytes','role','checksum_method','checksum','aspera','galaxy') if file.get(k)]})
        paths.append({'steps':steps})
