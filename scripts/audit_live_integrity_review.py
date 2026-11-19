# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Source-backed follow-up checks and semantic review of sparse deposited rows."""
import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path
import re
import statistics
import sys

REPO = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPO), str(REPO / 'src')]
from tests.support.authorship_capture import write_json
from tests.support.authorship_sources import documents, provenance
from tests.support.authorship_fidelity import table
from tests.support.deposited_fidelity import without_derived_files, without_added_protocol_placeholders, NODE_HEADERS
from tests.support.live_integrity_semantics import normalize, compare_rows, row_graph, nonempty, organism_taxids
from scripts.audit_additional_conversions import clean_revision
from scripts.audit_authorship_candidate import verify_seal
from scripts.audit_live_integrity import PRODUCTION


def read(path):
    return json.loads(Path(path).read_text())


def semantic_fields(graph):
    def walk(item, path):
        if isinstance(item, dict):
            if 'value' in item and nonempty(item['value']):
                yield {'field': item['kind'], 'value': item['value'], 'owner_path': path}
            for key in ('graph', 'factors', 'attributes', 'qualifiers'):
                for index, child in enumerate(item.get(key, [])):
                    yield from walk(child, path + '/' + key + '/' + str(index))
    return list(walk(graph, 'row'))


def reference_review(root, output):
    results, fields = [], []
    for acc in dict.fromkeys(s['reference'] for s in read(root/'cohort.json')['inputs']):
        original = normalize(table(root/'references'/f'{acc}.sdrf.txt'))
        source = without_derived_files(original)
        result = {'accession': acc, 'routes': {}, 'excluded_derived_cells':
            sum(nonempty(v) for row in original[1:] for v in row) -
            sum(nonempty(v) for row in source[1:] for v in row)}
        for route in ('direct', 'via-json'):
            document = read(root/'deposited'/acc/route/'payload-0.json')[0]
            actual = normalize(without_added_protocol_placeholders(table(root/'references'/f'{acc}.idf.txt'), document))
            detail = compare_rows(source, actual)
            result['routes'][route] = detail
            if route == 'via-json':
                for pair in detail['matches']:
                    left = semantic_fields(row_graph(source[0], source[pair['source_row']]))
                    right = semantic_fields(row_graph(actual[0], actual[pair['generated_row']]))
                    assert left == right
                    fields.extend(dict(accession=acc, **pair, **cell, source=cell['value'],
                        generated=cell['value'], reference=cell['value'],
                        classification='equivalent-or-documented-transformation',
                        evidence=f'references/{acc}.sdrf.txt',
                        reason='Full ordered row graph and duplicate-aware qualifier ownership match; empty slots/header syntax normalized.') for cell in left)
        results.append(result)
    write_json(output/'deposited-semantic-review.json', results)
    with (output/'deposited-semantic-fields.jsonl').open('w') as stream:
        for field in fields:
            stream.write(json.dumps(field, ensure_ascii=False)+'\n')
    return results


def native_files(docs):
    """Independent file assertions from provider XML and aligned ENA columns."""
    facts = {}
    for doc in docs:
        if doc['reference_assisted']:
            continue
        if doc['kind'] == 'xml':
            for run in doc['xml'].iter('RUN'):
                for file in run.findall('SRAFiles/SRAFile'):
                    for location in [file, *file.findall('Alternatives')]:
                        uri = location.get('url')
                        if not uri:
                            continue
                        item = {'run': run.get('accession'), 'uri': uri, 'filename': file.get('filename'),
                            'md5': file.get('md5'), 'bytes': file.get('size'),
                            'evidence': provenance(doc, f"RUN[@accession={run.get('accession')}]/SRAFiles/SRAFile")}
                        facts[tuple(item[k] for k in ('run','uri','filename','md5','bytes'))] = item
        elif doc['kind'] == 'json' and isinstance(doc['data'], list) and 'filereport' in doc['url']:
            for row in doc['data']:
                if not isinstance(row, dict) or not row.get('run_accession'):
                    continue
                for prefix in ('fastq', 'submitted'):
                    if not row.get(prefix+'_ftp'):
                        continue
                    links = row[prefix+'_ftp'].split(';')
                    md5s = row.get(prefix+'_md5', '').split(';')
                    sizes = row.get(prefix+'_bytes', '').split(';')
                    assert len(md5s) == len(sizes) == len(links), 'Unaligned provider file columns'
                    for uri, md5, size in zip(links, md5s, sizes):
                        if '://' not in uri:
                            uri = 'ftp://' + uri
                        item = {'run': row['run_accession'], 'uri': uri, 'filename': uri.rsplit('/',1)[-1],
                            'md5': md5, 'bytes': size, 'evidence': provenance(doc, f"read_run[{row['run_accession']}]/{prefix}")}
                        facts[tuple(item[k] for k in ('run','uri','filename','md5','bytes'))] = item
    return list(facts.values())


def exported_files(sdrf):
    result = defaultdict(list)
    for row in sdrf[1:]:
        runs = {v for h,v in zip(sdrf[0],row) if h in ('Comment[ENA_RUN]', 'Scan Name') and re.fullmatch(r'[SED]RR\d+',str(v))}
        if len(runs) != 1:
            continue
        groups = defaultdict(lambda: defaultdict(list))
        fastq_owner = False
        for h,v in zip(sdrf[0],row):
            if h in NODE_HEADERS:
                fastq_owner = False
            if h == 'Comment[FASTQ_URI]':
                fastq_owner = True
            elif h == 'Comment[MD5]' and fastq_owner:
                groups['FASTQ']['MD5'].append(v)
            match = re.fullmatch(r'Comment\[(ARCHIVE_FILE|SUBMITTED_FILE|FASTQ)_(.+)\]',h)
            if match:
                key = {'FILE_NAME':'NAME','FILE_FORMAT':'FORMAT'}.get(match[2],match[2])
                groups[match[1]][key].append(v)
        for group in groups.values():
            for i in range(max(map(len,group.values()),default=0)):
                item = {k: v[i] for k,v in group.items() if i<len(v)}
                if item.get('URI'):
                    # The complete URI itself retains its basename. Explicit
                    # supplied filenames take precedence and cannot be replaced.
                    item.setdefault('NAME', item['URI'].rsplit('/',1)[-1])
                    result[next(iter(runs)),item['URI']].append(item)
    return result


def source_checks(root, output):
    checks, counts, missing_files = [], [], []
    def record(case, mode, scope, entity, field, expected, actual, evidence, issue=None):
        checks.append({'case':case,'mode':mode,'scope':scope,'entity':entity,'field':field,
            'source':expected,'generated':actual,'reference':None,'evidence':evidence,
            'classification':'equivalent-or-documented-transformation' if expected==actual else
                'source-supported-loss-alteration-or-misassociation','issue':issue if expected!=actual else None})
    for spec in read(root/'cohort.json')['inputs']:
        case_name=spec['provider']+'-'+spec['accession']; case=root/'cases'/case_name
        for mode in ('off','curators','standard'):
            package=read(case/mode/'source-json/payload-0.json')[0]
            document=read(case/mode/'replayed-magetab/payload-0.json')[0]
            sdrf=next(r[1] for r in document if r[0]=='SDRF File')
            docs=documents(case,mode); paths=package['series'].get('assay_paths',[])
            sample_ids={s['iid'] for s in package['sample']}
            counts.append({'case':case_name,'mode':mode,'sample_records':len(sample_ids),
                'sdrf_rows':len(sdrf)-1,'canonical_paths':len(paths),
                'source_materials':len({n['name'] for p in paths for n in p['steps'] if n['kind']=='source'}),
                'assays':len({n['name'] for p in paths for n in p['steps'] if n['kind']=='assay'}),
                'runs':len({r['run'] for s in package['sample'] for r in s.get('sra_run',[])}),
                'file_occurrences':sum(len(r.get('files',[])) for s in package['sample'] for r in s.get('sra_run',[])),
                'biological_units':None})
            if spec['provider']=='geo':
                source=next(d for d in docs if d['kind']=='xml' and d['xml'].tag=='MINiML')
                for sample in source['xml'].findall('Sample'):
                    sid=sample.get('iid'); channels=sample.findall('Channel')
                    owned=[p for p in paths if any(n.get('sample_ref')==sid for n in p['steps'])]
                    for kind in ('source','extract'):
                        values={n['name'] for p in owned for n in p['steps'] if n['kind']==kind}
                        record(case_name,mode,'JSON',sid,kind+' material count',len(channels),len(values),provenance(source,f'Sample[{sid}]/Channel'))
                    source_index=sdrf[0].index('Source Name')
                    owned_rows=[row for row in sdrf[1:] if row[source_index] in {sid,*[f'{sid}:channel{i+1}' for i in range(len(channels))]}]
                    for i, channel in enumerate(channels):
                        material = sid if len(channels)==1 else f'{sid}:channel{i+1}'
                        channel_rows = [r for r in owned_rows if r[source_index] == material]
                        for organism in channel.findall('Organism'):
                            taxid = organism.get('taxid')
                            if not taxid: continue
                            preserved = bool(channel_rows) and all(taxid in organism_taxids(sdrf[0],r) for r in channel_rows)
                            record(case_name,mode,'SDRF',[sid,i],'owned organism taxid',taxid,
                                taxid if preserved else None,
                                provenance(source,f'Sample[{sid}]/Channel[{i+1}]/Organism/@taxid'),'MSC-LIVE-003')
                    for rel in sample.findall('Relation'):
                        if rel.get('type')!='SRA':continue
                        target=rel.get('target'); accession=target.rsplit('=',1)[-1]
                        record(case_name,mode,'SDRF',sid,'source SRA experiment identity',accession,
                            accession if any(target in row or accession in row for row in owned_rows) else None,
                            provenance(source,f'Sample[{sid}]/Relation[type=SRA]'),'MSC-LIVE-001')
                    for file in sample.findall('Supplementary-Data'):
                        uri=(file.text or '').strip()
                        if not uri:continue
                        if uri == 'NONE':
                            missing_files.append({'case':case_name,'mode':mode,'sample':sid,
                                'source':file.text,'generated':None,
                                'reason':'Exact GEO missing-file sentinel; not a populated file relationship.',
                                'evidence':provenance(source,f'Sample[{sid}]/Supplementary-Data')})
                            continue
                        record(case_name,mode,'SDRF',sid,'owned supplied file URI',uri,
                            uri if any(uri in row for row in owned_rows) else None,
                            provenance(source,f'Sample[{sid}]/Supplementary-Data'))
            inventory=defaultdict(list)
            for sample in package['sample']:
                for run in sample.get('sra_run',[]):
                    for file in [*run.get('files',[]),*run.get('fastq_files',[])]:
                        inventory[run['run'],str(file.get('uri') or file.get('url') or '').strip()].append(file)
            exported=exported_files(sdrf)
            for source in native_files(docs):
                key=(source['run'],source['uri'])
                for field,label in [('filename','NAME'),('md5','MD5'),('bytes','BYTES')]:
                    if not source[field]:continue
                    expected=[str(source[field])]
                    observed=sorted({str(f[field]) for f in inventory[key] if f.get(field)})
                    record(case_name,mode,'JSON',list(key),field,expected,observed,source['evidence'])
                    # Explicit owned export blocks only; unsupported layouts are limitations.
                    if key in exported:
                        actual=sorted({str(f[label]) for f in exported[key] if f.get(label)})
                        record(case_name,mode,'SDRF',list(key),field,expected,actual,source['evidence'])
    write_json(output/'repair-target-checks.json',checks)
    write_json(output/'entity-counts.json',counts)
    write_json(output/'missing-file-sentinels.json',missing_files)
    return checks


def protocol_projection_review(root, output):
    """Contrast supplied processing semantics with generated JSON scalar copies."""
    idf = {r[0]: r[1:] for r in table(root/'references/E-GEOD-37703.idf.txt')}
    index = idf['Protocol Type'].index('feature_extraction')
    name, text = idf['Protocol Name'][index], idf['Protocol Description'][index]
    assert 'TopHat' in text and 'Cuffdiff' in text
    values = [('deposited/E-GEOD-37703', 'off', root/'deposited/E-GEOD-37703/saved.json')]
    values += [(case, 'standard', root/'cases'/case/'standard/saved-0.json')
               for case in ('sra-SRP012607', 'ena-SRP012607')]
    records = []
    for case, mode, path in values:
        package = read(path)
        for sample in package['sample']:
            for i, channel in enumerate(sample.get('channel', [])):
                if channel.get('extract_protocol') != text:
                    continue
                records.append({'case':case, 'mode':mode, 'entity':sample['iid'],
                    'field':f'channel/{i}/extract_protocol', 'source':None, 'generated':text,
                    'reference':{'protocol_name':name,'type':'feature_extraction','description':text},
                    'classification':'unsupported-information-introduced', 'issue':'MSC-LIVE-002',
                    'evidence':'references/E-GEOD-37703.idf.txt:Protocol Type/Protocol Description',
                    'generated_extraction_definitions':[p['name'] for p in package['series'].get('protocols', [])
                        if p.get('description') == text and p.get('type', {}).get('value') == 'nucleic acid extraction protocol'
                        and ':'+sample['iid']+':extract_protocol' in p['name']]})
    write_json(output/'protocol-projection-review.json',records)
    return records


def json_inventory(root, output):
    """Retain every populated scalar and its occurrence path, without inventing verdicts."""
    def leaves(value, path=''):
        if isinstance(value, dict):
            for key, child in value.items():
                yield from leaves(child, path+'/'+str(key).replace('~','~0').replace('/','~1'))
        elif isinstance(value, list):
            for i, child in enumerate(value): yield from leaves(child, path+'/'+str(i))
        elif nonempty(value):
            yield path, value
    counts = Counter()
    with (output/'json-populated-fields.jsonl').open('w') as stream:
        for path in sorted((root/'cases').glob('*/*/saved-0.json')):
            case, mode = path.parts[-3:-1]
            for pointer, value in leaves(read(path)):
                stream.write(json.dumps({'case':case,'mode':mode,'path':pointer,'generated':value,
                    'assessment':'inventory; consult source checks and review records for verdicts'},ensure_ascii=False)+'\n')
                counts[case+'/'+mode] += 1
    write_json(output/'json-field-counts.json',dict(counts))


def paired_reference_review(root, output):
    """Review every supplied reference cell after source-backed entity grouping.

    File-expanded rows remain separate. This is agreement, not a source-loss
    oracle: unmatched annotations stay review requirements, never inferred loss.
    """
    aliases, evidence = {}, []
    for doc in documents(root/'cases/sra-SRP012607', 'off'):
        if doc['kind'] != 'xml': continue
        for sample in doc['xml'].findall('.//SAMPLE'):
            geo = sample.find("IDENTIFIERS/EXTERNAL_ID[@namespace='GEO']")
            if geo is None or not re.fullmatch(r'GSM\d+',geo.text or ''): continue
            for value in [sample.get('accession'), *(n.text for n in sample.findall('IDENTIFIERS/*'))]:
                if value:
                    assert value not in aliases or aliases[value] == geo.text
                    aliases[value] = geo.text
            evidence.append(provenance(doc,f"SAMPLE[{sample.get('accession')}]/IDENTIFIERS"))
    channel_evidence = {}
    for case in sorted((root/'cases').glob('geo-*')):
        source = next(d for d in documents(case,'off') if d['kind']=='xml' and d['xml'].tag=='MINiML')
        for sample in source['xml'].findall('Sample'):
            aliases[sample.get('iid')] = sample.get('iid')
            channel_evidence[sample.get('iid')] = {'count':len(sample.findall('Channel')),
                'evidence':provenance(source,f"Sample[{sample.get('iid')}]/Channel")}
    def identity(headers, row):
        name = row[headers.index('Source Name')]
        match = re.fullmatch(r'(GSM\d+)(?::channel|\s+)(\d+)',name)
        sample, channel = (match[1],int(match[2])) if match else (aliases.get(name,name),1)
        runs = {v for h,v in zip(headers,row) if h in ('Scan Name','Comment[ENA_RUN]') and re.fullmatch(r'[SED]RR\d+',v)}
        return sample, channel, runs
    def fields(headers,row):
        graph=row_graph(headers,row)
        result=defaultdict(list)
        def add(item, owner):
            kind=item['kind'].casefold()
            if nonempty(item['value']): result[owner+(kind,)].append(item['value'])
            for child in item.get('attributes',[]): add(child,owner+(kind,))
            for child in item.get('qualifiers',[]): add(child,owner+(kind,))
        for node in graph['graph']:add(node,())
        for factor in graph['factors']:add(factor,('row-factor',))
        return result
    records, groups=[],[]
    for spec in read(root/'cohort.json')['inputs']:
        case=spec['provider']+'-'+spec['accession']
        reference=normalize(table(root/'references'/(spec['reference']+'.sdrf.txt')))
        assessment=read(root/'cases'/case/'assessment.json')
        for mode in ('off','curators','standard'):
            document=read(root/'cases'/case/mode/'replayed-magetab/payload-0.json')[0]
            actual=normalize(next(r[1] for r in document if r[0]=='SDRF File'))
            actual_ids=[identity(actual[0],row) for row in actual[1:]]
            actual_fields=[fields(actual[0],row) for row in actual[1:]]
            for ri,row in enumerate(reference[1:],1):
                sid,channel,runs=identity(reference[0],row)
                assert sid in channel_evidence and channel <= channel_evidence[sid]['count']
                matches=[i for i,(gs,gc,gr) in enumerate(actual_ids,1)
                    if gs==sid and gc==channel and (not runs or not gr or bool(runs & gr))]
                group={'case':case,'mode':mode,'reference':spec['reference'],'reference_row':ri,
                    'sample':sid,'channel':channel,'runs':sorted(runs),'generated_rows':matches,
                    'reference_assisted':assessment['modes'][mode]['assistance']['reference_assisted'],
                    'identity_evidence':channel_evidence[sid],
                    'grouping':'Explicit source sample aliases and channel occurrence; run must agree when both supplied. File-expanded rows remain separate.'}
                groups.append(group)
                for field,expected in fields(reference[0],row).items():
                    generated=[{'row':i,'values':actual_fields[i-1].get(field,[])} for i in matches]
                    clean=lambda values:[' '.join(str(v).split()) for v in values]
                    equivalent=bool(generated) and all(clean(g['values'])==clean(expected) for g in generated)
                    classification='equivalent-or-documented-transformation' if equivalent else (
                        'reference-only-annotation-requiring-assessment' if generated and all(not g['values'] for g in generated)
                        else 'source-version-reference-defect-or-unresolved')
                    records.append({**group,'owner_and_field':list(field),'source':None,'reference_values':expected,
                        'generated':generated,'classification':classification,
                        'reason':'Agreement compares exact owned occurrences after whitespace normalization; other differences require source/curation review, not an automatic converter-failure verdict.'})
    write_json(output/'paired-reference-groups.json',{'native_aliases':aliases,'native_alias_evidence':evidence,'groups':groups})
    with (output/'paired-reference-fields.jsonl').open('w') as stream:
        for record in records:stream.write(json.dumps(record,ensure_ascii=False)+'\n')
    return {'groups':len(groups),'matched_groups':sum(bool(g['generated_rows']) for g in groups),
        'field_groups':len(records),'classifications':dict(Counter(r['classification'] for r in records)),
        'source_fidelity_warning':'Agreement/inventory denominator, not converter-loss counts.'}


def timings(root, output):
    result=[]
    for path in sorted((root/'cases').glob('*/*/timing-*.json')):
        value=read(path); stages=defaultdict(list)
        for run in value['records'][1:]:
            totals=Counter()
            for event in run['events']:
                if event['event_type']=='stage_finished':
                    totals[event['stage']]+=event['duration_seconds']
            for name,duration in totals.items():stages[name].append(duration)
        result.append({k:value[k] for k in ('case','mode','route','revision','summary')} |
            {'stage_seconds':{name:{'median':statistics.median(v),'min':min(v),'max':max(v)} for name,v in stages.items()},
             'stage_warning':'Nested inclusive stages overlap; do not sum across stage names.'})
    write_json(output/'runtime-stages.json',result)
    lines = ['# Runtime measurements', '',
        'Live values are one acquisition per input across off, curators and standard, in that order. '
        'Later phases reuse earlier successful responses. Pacing remains enabled, including on reused requests. '
        'Conversion-call time includes capture I/O; output writing is separate. All values are seconds except memory.', '',
        '| Input | Total acquisition | Physical | Reused | Retry identities |',
        '|---|---:|---:|---:|---:|']
    captures = [read(p) for p in sorted((root/'cases').glob('*/acquisition.json'))]
    for value in captures:
        n=value['request_counts']
        lines.append(f"| {value['provider']}-{value['accession']} | {value['elapsed_seconds']:.3f} | {n['physical_requests']} | {n['reused_responses']} | {n['repeated_physical_identities']} |")
    lines += ['', '| Input/preset | Call | CPU | Artifact write | Pacing wait | Physical/reused | Peak RSS MiB |',
        '|---|---:|---:|---:|---:|---:|---:|']
    for value in captures:
        for mode, phase in value['phases'].items():
            n=phase['request_counts']
            lines.append(f"| {value['provider']}-{value['accession']}/{mode} | {phase['conversion_call_seconds']:.3f} | {phase['cpu_seconds']:.3f} | {phase['artifact_write_seconds']:.3f} | {phase['gate_wait_seconds']:.3f} | {n['physical_requests']}/{n['reused_responses']} | {phase['peak_rss_bytes']/2**20:.1f} |")
    lines += ['', 'Offline cells run sequentially in separate processes: one warm-up and three measured repetitions. '
        'No external network fallback is allowed, and offline pacing is suppressed. '
        'These measure reused-response compute, not cold live performance. '
        'Saved-JSON cells use local fresh JSON. Imports/manifest setup and artifact writes are excluded from call time. '
        'Peak RSS is the process high-water mark including imports and warm-up. '
        'Stage durations are inclusive and overlap; they must not be summed.', '',
        '| Input | Preset | Route | Median | Range | Median CPU | Peak RSS MiB | Artifact write |',
        '|---|---|---|---:|---:|---:|---:|---:|']
    for value in result:
        s=value['summary']; detail=read(root/'cases'/value['case']/value['mode']/('timing-'+value['route']+'.json'))
        written=next(r.get('artifact_write_seconds',0) for r in detail['records'] if r['iteration']==1)
        lines.append(f"| {value['case']} | {value['mode']} | {value['route']} | {s['median_seconds']:.3f} | {s['min_seconds']:.3f}–{s['max_seconds']:.3f} | {s['median_cpu_seconds']:.3f} | {s['peak_rss_bytes']/2**20:.1f} | {written:.3f} |")
    consistent = sum(v['summary']['payloads_identical'] and v['summary']['outcomes_identical'] for v in result)
    physical = sum(v['request_counts']['physical_requests'] for v in captures)
    lines += ['', f'{consistent}/{len(result)} cells have identical measured payloads/outcomes. '
        'Errors and replay misses are retained per iteration in timing JSON. '
        'The initial mixed-transport harness attempt is preserved separately in '
        '`harness-failures/initial-mixed-transport/timing-json.json`: its four immediate failures were excluded, '
        'then the harness was corrected at 76e01fd before the complete measurement matrix. '
        f'Source metadata needed {physical} physical requests; reference acquisition is recorded separately.', '']
    (output/'RUNTIMES.md').write_text('\n'.join(lines))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bundle',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--revision',required=True)
    args=parser.parse_args(); root=args.bundle.resolve(); output=args.output.resolve()
    revision = clean_revision(args.revision, PRODUCTION)
    if output == root or output in root.parents:
        raise ValueError('Output cannot replace the input or an ancestor')
    sealed = (root/'SHA256SUMS').exists()
    if sealed and root in output.parents:
        raise ValueError('Sealed evidence is immutable; choose a separate output')
    before = verify_seal(root) if sealed else None
    output.mkdir(parents=True,exist_ok=False)
    references=reference_review(root,output)
    checks=source_checks(root,output); timings(root,output); json_inventory(root,output)
    protocols = protocol_projection_review(root,output)
    paired = paired_reference_review(root,output)
    if sealed: assert verify_seal(root) == before
    write_json(output/'summary.json',{'repair_target_checks':len(checks),
        'review_revision':revision,'sealed_input':before,
        'classifications':dict(Counter(c['classification'] for c in checks)),
        'issues':dict(Counter(c['issue'] or 'unadjudicated' for c in checks if c['source']!=c['generated'])),
        'unsupported_protocol_scalar_projections':len(protocols),
        'paired_reference_review':paired,
        'references_preserved':sum(r['routes']['via-json']['equivalent'] for r in references),
        'meaningful_reference_cells':sum(r['routes']['via-json']['meaningful_source_cells'] for r in references),
        'whitespace_only_reference_cells':sum(r['routes']['via-json']['whitespace_only_source_cells'] for r in references),
        'excluded_derived_cells':sum(r['excluded_derived_cells'] for r in references)})


if __name__=='__main__':
    main()
