# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Create evidence-labelled scientific derivatives without editing converter policy."""
import argparse
from collections import Counter
from copy import deepcopy
import json
from pathlib import Path
import re
import shutil

from scripts.audit_authorship_candidate import verify_seal, sha256
from scripts.audit_authorship import policy, save_result
from scripts.audit_recorded_provider import primitive
from tests.support.authorship_capture import write_json
from tests.support.authorship_sources import documents, citation_facts, provenance
from tests.support.recorded_replay import ReplayTransport
from meta_standards_converter.miniml import MINiMLCodec
from meta_standards_converter.converters import Converter
from meta_standards_converter.metadata.integrity import entity_counts


def review(source, candidate, output, corrections):
    source, candidate, output, corrections = (Path(p).resolve() for p in (source, candidate, output, corrections))
    if output.exists() or any(output.is_relative_to(p) or p.is_relative_to(output) for p in (source, candidate, corrections)):
        raise ValueError('Use a separate fresh derivative directory')
    seals = {str(p): verify_seal(p) for p in (source, candidate)}
    output.mkdir(parents=True)
    empty = output / 'empty.manifest.json'
    write_json(empty, {'schema_version': '1.0', 'cohort': 'reviewed-integrity', 'source_revision': 'review',
                      'captured_at': 'frozen', 'evidence_kind': 'recorded-http', 'files': {}, 'requests': []})
    records, pending = [], []
    for case in ('geo-GSE13398', 'geo-GSE31356', 'geo-GSE52778', 'sra-ERP148790', 'sra-ERP185509'):
        study = case.split('-', 1)[1]
        original = candidate / 'cases' / case / 'off/saved-0.json'
        data = json.loads(original.read_text())
        before = deepcopy(data)
        changes = []
        docs = documents(source / 'cases' / case, 'off')
        evidence = [provenance(d, 'source study and sample records') for d in docs if not d['reference_assisted']]
        def change(owner, field, value, location, reason, supplied=None, support=None):
            previous = deepcopy(owner.get(field))
            if previous == value:
                return
            changes.append({'entity': location, 'field': field, 'original_value': previous,
                            'proposed_value': deepcopy(value), 'source_value': supplied,
                            'evidence': support or evidence, 'rationale': reason,
                            'status': 'source-reviewed derivative; not an automatic converter rule'})
            owner[field] = deepcopy(value)
        def annotate(owner, field, item, location, reason, supplied=None):
            change(owner, field, [*owner.get(field, []), item], location, reason, supplied)
        def note(owner, location, value, supplied=None):
            annotate(owner, 'comments', {'name': 'reviewed interpretation', 'value': value},
                     location, 'Preserve source text and record this interpretation separately.', supplied)
        samples = {s['iid']: s for s in data['sample']}
        facts = {}
        for sid, sample in samples.items():
            facts[sid] = [{c['name']: c['value'].strip() for c in channel.get('characteristics', [])}
                          for channel in sample.get('channel', [])]
        reviewed_counts = entity_counts(data)
        if study == 'GSE13398':
            assert len(samples) == 16 and all(len(s['channel']) == 2 for s in samples.values())
            assert 'common reference' in data['series']['overall_design']
            classifications = {}
            for sid, sample in samples.items():
                assert 'Stratagene' in sample['channel'][0]['source']['value']
                text = ' '.join(v['value'] for v in sample['channel'][1]['characteristics']).lower()
                description = (sample.get('description') or '').lower()
                contrast = 'tumor' if 'tumor' in text or 'carcinoma' in description else 'normal' if 'normal' in text + description else None
                assert contrast
                classifications[sid] = contrast
                for index, role in enumerate(('common reference', 'biological comparison sample')):
                    annotate(sample['channel'][index], 'characteristics', {'name': 'reviewed channel role', 'value': role},
                             f'/sample/{sid}/channel/{index}', 'Design and channel source jointly establish role; label alone is insufficient.',
                             sample['channel'][index]['source']['value'])
            assert Counter(classifications.values()) == {'tumor': 8, 'normal': 8}
            note(data['series'], '/series', '16 assays, 32 channel materials: eight tumor and eight normal comparisons against a modified universal reference. Donor identities/count remain unspecified.',
                 data['series']['overall_design'])
        elif study == 'GSE31356':
            assert '3 patients' in data['series']['overall_design'] and 'Primary cell cultures' in data['series']['overall_design']
            pairs = {}
            for sid, sample in samples.items():
                match = re.fullmatch(r'(ND|D)(\d+)', sample['title'])
                assert match
                pairs.setdefault(match[2], set()).add(match[1])
                for name, value in [('reviewed patient pair', match[2]), ('reviewed specimen context', 'cultured primary cells')]:
                    annotate(sample['channel'][0], 'characteristics', {'name': name, 'value': value},
                             f'/sample/{sid}/channel/0', 'Patient-matched design supports the supplied D/ND sample code; cultured-cell context is explicit.',
                             {'title': sample['title'], 'design': data['series']['overall_design']})
            assert pairs == {'2': {'D', 'ND'}, '6': {'D', 'ND'}, '18': {'D', 'ND'}}
            reviewed_counts['biological_units'] = 3
            reviewed_counts['biological_units_reason'] = 'Three patient pairs explicitly stated in source design'
            note(data['series'], '/series', 'Three patient-matched D/ND pairs; the tissue-origin contrast is measured after primary cell culture. Six CEL relationships are retained.', data['series']['overall_design'])
        elif study == 'GSE52778':
            crossing = Counter((facts[sid][0]['cell line'], facts[sid][0]['treatment']) for sid in samples)
            assert len(crossing) == 16 and set(crossing.values()) == {1}
            assert len({donor for donor, treatment in crossing}) == 4
            reviewed_counts['biological_units'] = 4
            reviewed_counts['biological_units_reason'] = 'Source growth protocol explicitly states four donors; four cell-line identifiers each cross four treatments'
            note(data['series'], '/series', 'Four donor-derived cell lines crossed with untreated, albuterol, dexamethasone and combination conditions. ERCC mix is uneven; modelling and normalization choices remain downstream decisions.',
                 {'crossing': [list(k) for k in crossing], 'ERCC': {sid: facts[sid][0]['ercc_mix'] for sid in samples}})
        elif study == 'ERP148790':
            assert len(samples) == 1 and reviewed_counts['runs'] == 4
            assert 'nuclei isolation' in next(iter(samples.values()))['description'].lower()
            note(data['series'], '/series', 'Supported interpretation: single-nucleus ATAC-seq; one sample and four sequencing runs. File rows are not biological replicates.',
                 next(iter(samples.values()))['description'])
            pending.append({'study': study, 'entity': 'E-MTAB-13131 barcode annotations', 'original_value': {
                'library_cell_barcode_read': 'index1', 'library_sample_barcode_read': 'index2',
                'file_I1_role': 'sample_barcode/index1', 'file_R2_role': 'cell_barcode/index2'},
                'proposed_value': None, 'status': 'review required; no index swap applied',
                'evidence': [{'path': str(source / 'references/E-MTAB-13131.sdrf.txt'),
                              'sha256': sha256(source / 'references/E-MTAB-13131.sdrf.txt')},
                             {'url': 'https://www.10xgenomics.com/support/epi-atac/documentation/steps/sequencing/sequencing-requirements-for-single-cell-atac'}],
                'rationale': 'Library and file-level assignments contradict each other. Verify applicable ATAC chemistry and original read metadata; filenames alone cannot establish a correction.'})
        elif study == 'ERP185509':
            assert len(samples) == 8 and not data['series'].get('pubmed_id') and not data['series'].get('pubmed_publication')
            sequential = [sid for sid in samples if facts[sid][0].get('duration') == '20 followed by 4']
            assert len(sequential) == 1
            note(data['series'], '/series', 'Bacterial RNA-seq with sequential exposure retained: 20 hours with alpibectir followed by 4 hours with ethionamide in one sample. Continuation/washout of the first compound is not specified. No supplied publication is added.',
                 {'sample': sequential[0], 'duration': facts[sequential[0]][0]['duration'], 'protocol': samples[sequential[0]]['description']})
            pending.append({'study': study, 'entity': 'four DMSO vehicle samples',
                'original_value': {'description': '0.1% DMSO', 'deposited_factor': '0.1 microgram per milliliter'},
                'proposed_value': None, 'status': 'review required; percentage basis and dose not corrected',
                'evidence': [{'path': str(source / 'references/E-MTAB-16253.sdrf.txt'),
                              'sha256': sha256(source / 'references/E-MTAB-16253.sdrf.txt')}],
                'rationale': 'Percentage and mass concentration are incompatible without additional basis/density evidence. Preserve sequential timing and original quantities.'})
        # Bind annotations to original sample/channel and protocol application.
        protocols = {p['name']: p for p in data['series'].get('protocols', [])}
        for pi, path in enumerate(data['series'].get('assay_paths', [])):
            refs = {s['sample_ref'] for s in path['steps'] if s.get('sample_ref')}
            if len(refs) != 1:
                continue
            sid = next(iter(refs))
            if sid not in samples:
                continue
            for si, step in enumerate(path['steps']):
                location = f'/series/assay_paths/{pi}/steps/{si}'
                if step['kind'] in {'source', 'sample'}:
                    markers = [c['value'] for c in step.get('comments', []) if c['name'] == 'msc_channel']
                    ci = int(markers[0]) if markers else 0
                    for field in samples[sid]['channel'][ci].get('characteristics', []):
                        if field['name'].startswith('reviewed '):
                            annotate(step, 'characteristics', field, location, 'Carry reviewed channel annotation onto the exact bound material.', sid)
                if study == 'GSE52778' and step['kind'] == 'protocol_application':
                    protocol = protocols.get(step['protocol_ref'], {})
                    if 'treatment' in str(protocol.get('type', '')).lower():
                        text = samples[sid]['channel'][0]['treatment_protocol']
                        assert '1 μM' in text and '18 h' in text
                        condition = facts[sid][0]['treatment'].lower()
                        for compound in ('albuterol', 'dexamethasone'):
                            if compound in condition:
                                for name, value, unit in [(compound + ' dose', '1', 'micromolar'),
                                                          (compound + ' exposure duration', '18', 'hour')]:
                                    annotate(step, 'parameter_values', {'name': name, 'value': value, 'unit': {'value': unit},
                                             'comments': [{'name': 'compound', 'value': compound}]},
                                             location, 'Dose/duration belong to this compound in this sample treatment application.', text)
        # Source-backed publication proposal, independent PMID/DOI/author checks.
        enriched = candidate / 'cases' / case / 'curators/saved-0.json'
        if enriched.exists() and case.startswith('geo-'):
            supplied = json.loads(enriched.read_text())['series'].get('pubmed_publication', [])
            pubfacts = citation_facts(documents(source / 'cases' / case, 'curators'))
            checked = []
            for publication in supplied:
                matches = [p for p in pubfacts.get(publication.get('pubmed_id'), [])
                           if all(p.get(k) == publication.get(k) for k in ('doi', 'author_list', 'title'))]
                assert matches, (study, publication)
                checked.extend(p['evidence'] for p in matches)
            if supplied:
                change(data['series'], 'pubmed_publication', supplied, '/series',
                       'Publication identity and ordered Unicode authors independently match frozen PubMed evidence.', support=checked)
        dest = output / study
        dest.mkdir()
        MINiMLCodec().decode(data, strict=True)
        annotated = dest / 'annotated.json'
        write_json(annotated, data)
        isolation = ReplayTransport(empty)
        with isolation.offline():
            result = Converter().convert(annotated, out_type='magetab', preparation=policy('off'), outdir=dest / 'magetab')
        isolation.assert_complete()
        assert result.items[0].execution == 'succeeded', result.to_dict()
        save_result(result, dest / 'conversion')
        write_json(dest / 'changes.json', changes)
        assert json.loads(original.read_text()) == before
        records.append({'study': study, 'input': str(original), 'input_sha256': sha256(original),
                        'annotated_sha256': sha256(annotated), 'changes': len(changes),
                        'counts': reviewed_counts, 'result': result.to_dict()})
    pending.append({'study': 'contacts', 'entity': 'name-only Himes and detailed Himes contact; repeated Ratkaj reference records',
                    'original_value': 'Separate supplied contributors', 'proposed_value': None,
                    'status': 'identity evidence required; no name-only consolidation',
                    'evidence': [{'path': str(source / 'contact-review.json'), 'sha256': sha256(source / 'contact-review.json')}],
                    'rationale': 'Names alone do not establish identity; JSON contributors and sample references remain intact.'})
    inherited = output / 'earlier-doi-corrections'
    shutil.copytree(corrections, inherited)
    write_json(output / 'earlier-doi-corrections-provenance.json',
               {p.name: {'original_sha256': sha256(p), 'copied_sha256': sha256(inherited / p.name)} for p in corrections.iterdir() if p.is_file()})
    write_json(output / 'review-records.json', pending)
    write_json(output / 'derivatives.json', records)
    assert all(verify_seal(Path(p)) == seal for p, seal in seals.items())
    write_json(output / 'input-seals.json', seals)
    return records


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('source', 'candidate', 'output', 'corrections'):
        parser.add_argument('--' + name, type=Path, required=True)
    review(**vars(parser.parse_args()))
