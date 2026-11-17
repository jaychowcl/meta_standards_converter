"""GEO literals remain source evidence rather than inferred curation."""
from copy import deepcopy
import csv

import pytest

from meta_standards_converter.miniml.geo_parser import GEOParser
from meta_standards_converter.miniml import MINiMLCodec
from meta_standards_converter.magetab.constructor import AEConstructor
from meta_standards_converter.magetab.writer import MAGETabWriter
from meta_standards_converter.metadata.archive_enrichment import merge_archive_metadata
from tests.converters.test_audit_source_precedence import package as archive_package


def geo_package():
    samples = ''.join(f'''<Sample iid="GSM{i}"><Title>Sample  {i}</Title>
      <Library-Strategy>RNA-Seq</Library-Strategy><Library-Source>transcriptomic</Library-Source>
      <Channel><Source>  whole  embryo </Source><Organism>Drosophila melanogaster</Organism>
      <Molecule>total RNA</Molecule><Characteristics tag="tissue">whole  embryo</Characteristics>
      <Characteristics tag="developmental stage">stage {i}</Characteristics>
      <Extract-Protocol>10x Chromium Single Cell 3' v3.1.  Author's method.</Extract-Protocol>
      </Channel><Supplementary-Data type="TXT">https://example.test/processed{i}.txt</Supplementary-Data>
      </Sample>''' for i in (1, 2))
    return GEOParser().parse(f'''<MINiML version="0.5.4">{samples}
      <Series iid="GSE1"><Title>Study</Title><Summary>Summary.</Summary><Overall-Design>Design.</Overall-Design>
      <Status><Submission-Date>May 01 2023</Submission-Date><Release-Date>May 02 2023</Release-Date></Status>
      <Sample-Ref ref="GSM1"/><Sample-Ref ref="GSM2"/></Series></MINiML>''')[0].to_mapping()


def rows(data, handler='single_cell_sequencing'):
    return AEConstructor().miniml2magetab(MINiMLCodec().decode(data).package, platform_handler=handler)


def row(table, label):
    return next(r for r in table if r[0] == label)


def test_geo_xml_preserves_decoded_text_and_ignores_indentation():
    text = '  Mixed  case\tand\nnewlines & punctuation.  '
    xml = '<MINiML>\n <Series iid="GSE1"><Summary>' + text.replace('&', '&amp;') + '</Summary></Series>\n</MINiML>'
    data = GEOParser().parse(xml)[0].to_mapping()
    assert data['series']['summary'] == text
    assert 'value' not in data['series']


@pytest.mark.parametrize('field', ['title', 'description'])
def test_geo_primary_values_win_over_linked_arrayexpress(field):
    a, b = archive_package('GEO'), archive_package('MAGE-TAB')
    a['sample'][0][field] = 'GEO  original'; b['sample'][0][field] = 'AE alternative'
    b['sample'][0]['description' if field == 'title' else 'title'] = 'Explicit missing fact'
    before = deepcopy(a)
    merged, issues = merge_archive_metadata(MINiMLCodec().decode(a).package,
        MINiMLCodec().decode(b).package, prefer=True, mode='standard')
    output = merged.to_mapping()
    assert output['sample'][0][field] == 'GEO  original'
    assert output['sample'][0]['description' if field == 'title' else 'title'] == 'Explicit missing fact'
    assert a == before
    conflicts = [r for r in output['extensions']['insdc']['records'] if r['kind'] == 'source_conflict']
    assert any(r['metadata']['selected']['value'] == 'GEO  original' for r in conflicts)


def test_geo_encoding_suggestion_does_not_replace_source():
    clean = 'Müller’s protocol'; broken = clean.encode('utf8').decode('cp1252')
    a, b = archive_package('GEO'), archive_package('MAGE-TAB')
    a['sample'][0]['description'] = broken; b['sample'][0]['description'] = clean
    output, _ = merge_archive_metadata(MINiMLCodec().decode(a).package,
        MINiMLCodec().decode(b).package, prefer=True, mode='standard')
    assert output.to_mapping()['sample'][0]['description'] == broken
    assert clean in str(output.to_mapping()['extensions']['insdc']['records'])


def test_geo_placeholder_word_is_still_an_authored_value():
    a, b = archive_package('GEO'), archive_package('MAGE-TAB')
    a['sample'][0]['description'] = 'not provided'
    b['sample'][0]['description'] = 'An alternative description'
    output, _ = merge_archive_metadata(MINiMLCodec().decode(a).package,
        MINiMLCodec().decode(b).package, prefer=True, mode='standard')
    assert output.to_mapping()['sample'][0]['description'] == 'not provided'


@pytest.mark.parametrize('field', ['organism', 'library_source'])
def test_geo_placeholder_organism_and_library_source_survive_conflicting_runs(field):
    a, b = archive_package('GEO'), archive_package('MAGE-TAB')
    a['sample'][0]['channel'][0]['organism'] = [{'value': 'not provided'}]
    a['sample'][0]['library_source'] = 'Authored library source'
    a['sample'][0]['sra_run'][0]['library_source'] = 'TRANSCRIPTOMIC'
    a['sample'][0]['sra_run'].append({'run': 'SRR2', 'experiment': 'SRX2', 'sample': 'SRS1', 'library_source': 'GENOMIC'})
    output, _ = merge_archive_metadata(MINiMLCodec().decode(a).package,
        MINiMLCodec().decode(b).package, prefer=True, mode='standard')
    sample = output.to_mapping()['sample'][0]
    if field == 'organism':
        assert sample['channel'][0]['organism'] == [{'value': 'not provided'}]
    else:
        assert sample['library_source'] == 'Authored library source'


def test_geo_characteristic_conflicts_keep_original_and_add_missing_tags():
    a, b = archive_package('GEO'), archive_package('MAGE-TAB')
    a['sample'][0]['channel'][0]['characteristics'] = [{'name': 'Tissue', 'value': '  lung  '}]
    b['sample'][0]['channel'][0]['characteristics'] = [
        {'name': 'Tissue', 'value': 'heart'}, {'name': 'genotype', 'value': 'WT'}]
    output, _ = merge_archive_metadata(MINiMLCodec().decode(a).package,
        MINiMLCodec().decode(b).package, prefer=True, mode='standard')
    assert output.to_mapping()['sample'][0]['channel'][0]['characteristics'] == [
        {'name': 'Tissue', 'value': '  lung  '}, {'name': 'genotype', 'value': 'WT'}]


def test_geo_protocols_with_distinct_authored_spacing_are_not_merged():
    from meta_standards_converter.magetab.protocol_export import prepare_protocols
    data = geo_package()
    data['sample'][0]['channel'][0]['extract_protocol'] = 'Authored  method'
    data['sample'][1]['channel'][0]['extract_protocol'] = 'Authored method'
    descriptions = row(rows(data), 'Protocol Description')[1:]
    assert 'Authored  method' in descriptions and 'Authored method' in descriptions
    data['series']['protocols'] = [{'name': 'one', 'description': 'Authored  method'},
                                 {'name': 'two', 'description': 'Authored method'}]
    prepare_protocols(data)
    assert [p['description'] for p in data['series']['protocols']] == ['Authored  method', 'Authored method']


def test_explicit_canonical_chemistry_and_barcodes_are_retained():
    data = geo_package()
    data['series']['assay_paths'] = [{'steps': [
        {'kind': 'source', 'name': 'GSM1', 'sample_ref': 'GSM1'},
        {'kind': 'extract', 'name': 'extract', 'comments': [
            {'name': 'library construction', 'value': 'Curated chemistry'},
            {'name': 'cell barcode size', 'value': '12'}]},
        {'kind': 'assay', 'name': 'assay'}]}]
    table = row(rows(data), 'SDRF File')[1]
    assert table[1][table[0].index('Comment[library construction]')] == 'Curated chemistry'
    assert table[1][table[0].index('Comment[cell barcode size]')] == '12'


def test_export_preserves_characteristics_material_and_library_literals():
    table = row(rows(geo_package()), 'SDRF File')[1]
    header, first = table[:2]
    assert first[header.index('Characteristics[tissue]')] == 'whole  embryo'
    assert not first[header.index('Characteristics[organism part]')]
    assert first[header.index('Comment[LIBRARY_SOURCE]')] == 'transcriptomic'
    material = [first[i] for i, label in enumerate(header) if label == 'Material Type']
    assert material == ['', 'total RNA']
    assert not any(label.startswith('Factor Value[') for label in header)


@pytest.mark.parametrize('handler', ['single_cell_sequencing', 'spatial_sequencing'])
def test_geo_prose_and_platform_selection_do_not_manufacture_annotations(handler):
    table = row(rows(geo_package(), handler), 'SDRF File')[1]
    forbidden = ['library construction', 'end bias', 'cell barcode', 'umi barcode', 'cdna read', 'read_type', 'technical replicate group']
    assert not any(any(term in h for term in forbidden) for h in table[0])


def test_dates_descriptions_and_protocols_keep_their_source_meaning():
    data = geo_package(); output = rows(data)
    assert not any(row(output, 'Date of Experiment')[1:])
    assert not any(r[0] == 'Comment[ArrayExpressSubmissionDate]' for r in output)
    assert row(output, 'Comment[GEOSubmissionDate]')[1:] == ['2023-05-01']
    assert row(output, 'Experiment Description')[1:] == ['Summary. Design.']
    descriptions = row(output, 'Protocol Description')[1:]
    assert data['sample'][0]['channel'][0]['extract_protocol'] in descriptions
    assert not any(' | ' in str(v) for v in descriptions)
    protocol_types = row(output, 'Protocol Type')[1:]
    library = protocol_types.index('nucleic acid library construction protocol')
    assert not descriptions[library]
    reference = row(output, 'Protocol Name')[library + 1]
    assert reference in str(row(output, 'SDRF File')[1])


def test_tsv_safety_is_only_applied_when_writing(tmp_path):
    data = geo_package(); literal = '  Authored  text\twith\r\na newline. '
    data['sample'][0]['channel'][0]['extract_protocol'] = literal
    data['sample'][0]['description'] = literal
    data['series']['summary'] = literal; data['series'].pop('overall_design')
    output = rows(data)
    assert literal in row(output, 'Protocol Description')
    assert literal == row(output, 'Experiment Description')[1]
    table = row(output, 'SDRF File')[1]
    assert table[1][table[0].index('Comment[Sample_description]')] == literal
    before = deepcopy(output)
    path = MAGETabWriter().write(output, str(tmp_path))
    with open(path) as f: published = list(csv.reader(f, delimiter='\t'))
    assert row(published, 'Experiment Description')[1] == '  Authored  text with a newline. '
    assert output == before


def test_explicit_curator_annotations_survive():
    data = geo_package()
    data['series']['variable'] = [{'name': 'developmental stage'}]
    for sample in data['sample']:
        sample['channel'][0]['characteristics'].extend([
            {'name': 'organism part', 'value': 'whole organism'}, {'name': 'replicate', 'value': '1'}])
    output = rows(data); table = row(output, 'SDRF File')[1]
    assert table[1][table[0].index('Characteristics[organism part]')] == 'whole organism'
    assert table[1][table[0].index('Characteristics[replicate]')] == '1'
    assert 'Factor Value[developmental stage]' in table[0]
    assert not any(row(output, 'Experimental Factor Type')[1:])


@pytest.mark.parametrize('provider', ['SRA', 'MAGE-TAB'])
def test_export_does_not_invent_a_submission_date_for_other_sources(provider):
    data = archive_package(provider)
    output = rows(data)
    assert not any(r[0] == 'Comment[ArrayExpressSubmissionDate]' for r in output)


def test_reimport_preserves_authored_protocol_and_submission_provenance(tmp_path):
    from meta_standards_converter.magetab.parser import AEParser
    from meta_standards_converter.sources.magetab import AEWebFetcher
    data = geo_package(); literal = '  Author\'s  extraction.  '
    data['sample'][0]['channel'][0]['extract_protocol'] = literal
    output = rows(data)
    output.append(['Comment[ArrayExpressSubmissionDate]', '2020-01-01'])
    path = MAGETabWriter().write(output, str(tmp_path))
    package = AEParser().parse(AEWebFetcher().resolve(path))
    assert package.to_mapping()['sample'][0]['channel'][0]['extract_protocol'] == literal
    again = AEConstructor().miniml2magetab(package)
    assert row(again, 'Comment[ArrayExpressSubmissionDate]')[1:] == ['2020-01-01']
    assert row(again, 'Comment[GEOSubmissionDate]')[1:] == ['2023-05-01']
    assert not any(row(again, 'Date of Experiment')[1:])
