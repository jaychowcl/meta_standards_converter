# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Semantic errors and scientific advisories have different outcome semantics."""
from copy import deepcopy
from meta_standards_converter.metadata.integrity import assess_integrity, qualifier_ownership, entity_counts
from meta_standards_converter.metadata.workflows import prepare_workflows
from tests.magetab.test_metadata_integrity import geo


def test_distinct_channels_with_identical_material_name_fail():
    data = prepare_workflows(geo('GSE13398')).to_mapping()
    assert not any(d.severity == 'error' for d in assess_integrity(data))
    data['series']['assay_paths'][1]['steps'][0]['name'] = data['series']['assay_paths'][0]['steps'][0]['name']
    assert any(d.code == 'integrity_material_collision' and d.severity == 'error' for d in assess_integrity(data))


def test_qualifier_swap_fails_but_inserted_application_does_not():
    a = {'series': {'assay_paths': [{'steps': [{'kind': 'source', 'name': 'S',
        'characteristics': [{'name': 'dose', 'value': '', 'comments': [{'name': 'compound', 'value': 'A'}]},
                            {'name': 'dose', 'value': '', 'comments': [{'name': 'compound', 'value': 'B'}]}]}]}]}}
    b = deepcopy(a)
    b['series']['assay_paths'][0]['steps'].insert(0, {'kind': 'protocol_application', 'protocol_ref': 'new'})
    assert not qualifier_ownership(a, b)
    b['series']['assay_paths'][0]['steps'][1]['characteristics'].reverse()
    assert len(qualifier_ownership(a, b)) == 2


def test_source_conflicts_are_advisory_and_chemistry_is_not_guessed():
    data = {'series': {'assay_paths': [{'steps': [
        {'kind': 'array_data_file', 'name': 'one', 'link': {'value': 'https://example/one'},
         'comments': [{'name': 'MD5', 'value': 'a' * 32}, {'name': 'cell barcode read', 'value': 'I1'},
                      {'name': 'cell barcode read', 'value': 'I2'}]},
        {'kind': 'array_data_file', 'name': 'two', 'link': {'value': 'https://example/one'},
         'comments': [{'name': 'MD5', 'value': 'b' * 32}]}]}]}}
    issues = assess_integrity(data)
    assert {d.code for d in issues} == {'curation_file_conflict', 'curation_barcode_assignment'}
    assert all(d.severity == 'warning' for d in issues)


def test_counts_do_not_invent_donors_from_channels():
    counts = entity_counts(prepare_workflows(geo('GSE13398')))
    assert counts['sample_records'] == 1
    assert counts['materials'] == 2
    assert counts['assays'] == 1
    assert counts['biological_units'] is None


def test_advisory_findings_do_not_make_conversion_partial():
    from meta_standards_converter.converters import Converter, PreparationPolicy
    result = Converter().convert(geo('GSE52778'), out_type='json',
        preparation=PreparationPolicy(enrichment='off', expand_studies=False))
    assert result.status == 'complete'
    assert any(d.code == 'curation_protocol_chronology' for d in result.items[0].diagnostics)


def test_deposited_barcode_contradiction_is_reported_without_guessing_a_swap():
    from tests.magetab.test_metadata_integrity import SOURCE
    from tests.support.magetab import parse
    source = SOURCE['reference_fastq']
    assert source['Comment[sample barcode read]'] != source['Comment[READ_INDEX]']
    package = parse(list(source), [list(source.values())])
    before = package.to_mapping()
    issues = assess_integrity(package)
    assert any(d.code == 'curation_barcode_assignment' for d in issues)
    assert all(d.severity != 'error' for d in issues)
    assert package.to_mapping() == before


def test_taxonomy_identifier_spellings_and_authored_alternatives():
    data = {'sample': [{'iid': 'S1', 'channel': [{'organism': [{'value': 'H37Rv', 'taxid': '83332'}]}]}],
            'series': {'assay_paths': [{'comments': [{'name': 'msc_workflow_origin', 'value': 'generated:SRA'}],
                'steps': [{'kind': 'source', 'name': 'S1', 'sample_ref': 'S1', 'characteristics': [
                    {'name': 'organism', 'value': 'H37Rv', 'term_accession_number': 'http://purl.obolibrary.org/obo/NCBITaxon_83332'}]}]}]}}
    assert not assess_integrity(data)
    data['series']['assay_paths'][0]['steps'][0]['characteristics'][0]['term_accession_number'] = 'taxon:1773'
    assert any(d.code == 'integrity_taxonomy_assertions' for d in assess_integrity(data))
    data['series']['assay_paths'][0]['comments'] = []
    assert all(d.severity == 'warning' for d in assess_integrity(data))
