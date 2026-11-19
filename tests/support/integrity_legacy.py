# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Review bounded legacy-output changes against the unchanged JSON input.

This oracle imports no production normalizers. Original row occurrences and
qualifier ownership remain subject to the deposited-data comparator after each
explicit, independently checked representation transformation below.
"""
from copy import deepcopy

from tests.support.deposited_fidelity import NODE_HEADERS


def reviewed_legacy_tables(package, previous, current, definitions):
    old, new = deepcopy(previous), deepcopy(current)
    review = []
    samples = {s['iid']: s for s in package['sample']}
    def remove(table, columns):
        return [[v for i, v in enumerate(row) if i not in columns] for row in table]
    # An entirely empty protocol column with no owned annotations is not an
    # application. Never remove a repeated populated application or qualifier.
    empty = {i for i, header in enumerate(old[0]) if header == 'Protocol REF'
             and all(row[i] == '' for row in old[1:])
             and (i + 1 == len(old[0]) or old[0][i + 1] in NODE_HEADERS)}
    if empty:
        review.append({'transformation': 'remove empty unqualified protocol columns', 'columns': sorted(empty)})
        old = remove(old, empty)
    new_headers = new[0]
    source_index = new_headers.index('Source Name')
    marker_columns = [i for i, h in enumerate(new_headers) if h in {'Comment[msc_sample]', 'Comment[msc_channel]'}]
    if marker_columns:
        assert not package['series'].get('assay_paths'), 'Authored paths cannot use the legacy generated-identity review'
        sid_index, channel_index = new_headers.index('Comment[msc_sample]'), new_headers.index('Comment[msc_channel]')
        for row in new[1:]:
            if not row[sid_index] and not row[channel_index]: continue
            sid = row[sid_index]
            sample = samples[sid]
            ci = int(row[channel_index])
            assert len(sample['channel']) > 1 and 0 <= ci < len(sample['channel'])
            identity = f'{sid}:channel{ci + 1}'
            assert row[source_index] == identity
            for i, header in enumerate(new_headers):
                if header in {'Source Name', 'Extract Name'}:
                    assert row[i] == identity
                    row[i] = sid
                elif header == 'Labeled Extract Name':
                    assert row[i].startswith(identity + ':')
                    row[i] = sid + row[i][len(identity):]
        review.append({'transformation': 'distinct channel identities verified against original sample/channel occurrences',
                       'columns': marker_columns})
        new = remove(new, set(marker_columns))
    # Added values must be present on the same source sample/file in old JSON.
    added = set(new[0]) - set(old[0])
    allowed = {'Comment[Unclassified File]', 'Comment[FASTQ_BYTES]'}
    remove_new = set()
    for i, header in enumerate(new[0]):
        if header not in added & allowed: continue
        assert i + 1 == len(new[0]) or new[0][i + 1] in NODE_HEADERS or new[0][i + 1].startswith(('Comment[', 'Characteristics[', 'Factor Value['))
        for row in new[1:]:
            if not row[i]: continue
            sid = row[new[0].index('Source Name')]
            sample = samples[sid]
            if header == 'Comment[Unclassified File]':
                assert not any(h in NODE_HEADERS for h in new[0][1:i]), 'Unknown link must remain on its source owner'
                assert row[i] in {f.get('value', '').strip(' \t\r\n') for f in sample.get('supplementary_data', [])}
            else:
                run_ids = {row[j] for j, h in enumerate(new[0]) if h in {'Comment[ENA_RUN]', 'Comment[SRA_RUN]'}}
                locations = [j for j in range(i) if new[0][j] == 'Comment[FASTQ_URI]']
                assert locations
                uri = row[locations[-1]]
                matches = [f for r in sample.get('sra_run', []) if r.get('run') in run_ids
                           for f in r.get('fastq_files', []) if f.get('uri') == uri]
                assert matches and {f.get('bytes') for f in matches} == {row[i]}
        remove_new.add(i)
        review.append({'transformation': 'retain previously omitted original JSON file metadata', 'column': header, 'index': i})
    new = remove(new, remove_new)
    # Only the specific generated GEO treatment/growth adjacency is changed;
    # explicitly authored order, repeated stages and owned qualifiers are sacred.
    if not package['series'].get('assay_paths') and package.get('source', {}).get('format') == 'GEO MINiML':
        descriptions = {name: ' '.join((text or '').split()) for name, text in definitions.items()}
        for row in old[1:]:
            sid = row[old[0].index('Source Name')]
            sample = samples.get(sid)
            if not sample or len(sample['channel']) != 1: continue
            channel = sample['channel'][0]
            growth, treatment = (' '.join((channel.get(key) or '').split()) for key in ('growth_protocol', 'treatment_protocol'))
            if not growth or not treatment or growth == treatment: continue
            grow = [i for i, h in enumerate(old[0]) if h == 'Protocol REF' and descriptions.get(row[i]) == growth]
            treat = [i for i, h in enumerate(old[0]) if h == 'Protocol REF' and descriptions.get(row[i]) == treatment]
            if len(grow) == len(treat) == 1 and grow[0] == treat[0] + 1:
                after = grow[0] + 1
                assert after == len(old[0]) or old[0][after] in NODE_HEADERS
                row[grow[0]], row[treat[0]] = row[treat[0]], row[grow[0]]
                review.append({'transformation': 'provider-generated growth before treatment', 'sample': sid,
                               'source_growth': growth, 'source_treatment': treatment})
    return old, new, review
