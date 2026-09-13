import copy
import json
from meta_standards_converter.miniml import MINiMLCodec
from meta_standards_converter.miniml.sra_parser import SRAParser
from meta_standards_converter.metadata.enrichment import MINiMLEnricher
from meta_standards_converter.metadata.interpretation import MINiMLMetadataService
from meta_standards_converter.converters.json2ae import JSON2AEConverter
from meta_standards_converter.converters.json2tsv import JSON2TSVConverter
from meta_standards_converter.magetab.constructor import AEConstructor
from tests.test_native_archive_parsers import fixture_records


def package():
    import xml.etree.ElementTree as ET
    records = fixture_records('sra', 'SRX017289')
    ET.SubElement(records.xml[0].find('.//LIBRARY_DESCRIPTOR'), 'LIBRARY_CONSTRUCTION_PROTOCOL').text = 'Explicit library preparation'
    data = SRAParser().parse(records).to_mapping()
    data['series']['accession'].append({'value': 'GSE18729', 'database': 'GEO'})
    data['series']['relation'] = [{'type': 'ArrayExpress', 'target': 'E-MTAB-1'}]
    data['sample'][0]['accession'].append({'value': 'GSM465245', 'database': 'GEO'})
    return data


def test_identity_and_default_enrichment_do_not_replace_native():
    data = package()
    interpreter = MINiMLMetadataService()
    assert interpreter.study_accession([data]) == 'SRP002056'
    assert interpreter.sample_accession(data['sample'][0]) == 'SRS011830'
    class Forbidden:
        def __getattr__(self, name):
            raise AssertionError('Native imports already control retrieval')
    typed = MINiMLCodec().decode(data).package
    assert MINiMLEnricher(pubmed_fetcher=Forbidden(), insdc_fetcher=Forbidden()).enrich(typed) == typed


def test_native_idf_accessions_and_no_implicit_factors():
    data = package()
    second = copy.deepcopy(data['sample'][0]); second['iid'] = 'SRS2'
    data['sample'].append(second)
    data['sample'][0]['channel'][0]['characteristics'].append({'name': 'dose', 'value': '1'})
    second['channel'][0]['characteristics'].append({'name': 'dose', 'value': '2'})
    data['series'].pop('assay_paths', None)
    rendered = AEConstructor().miniml2magetab(MINiMLCodec().decode(data).package, platform_handler='generic')
    rows = dict((row[0], row[1:]) for row in rendered if row[0] != 'SDRF File')
    assert rows['Comment[ArrayExpressAccession]'] == ['E-MTAB-1']
    i = rows['Comment[SecondaryAccession]'].index('GSE18729')
    assert rows['Comment[SecondaryAccessionTermSourceRef]'][i] == 'GEO'
    assert not any(rows['Experimental Factor Name'])
    assert not any('Factor Value[' in h for h in next(r[1] for r in rendered if r[0] == 'SDRF File')[0])


def test_native_json_downstream_exports_keep_biology_protocols_and_files(tmp_path):
    data = package()
    path = tmp_path / 'native.json'; path.write_text(json.dumps(data))
    tsv = tmp_path / 'native.tsv'
    JSON2TSVConverter().convert_source(path, tsv)
    text = tsv.read_text()
    assert 'SRS011830' in text and 'Caenorhabditis elegans' in text
    output = JSON2AEConverter().convert(str(path), platform_handler='generic')[0]
    sdrf = next(r[1] for r in output if r[0] == 'SDRF File'); header = sdrf[0]
    assert 'Characteristics[organism]' in header
    assert 'Protocol REF' in header
    assert 'Array Data File' in header
    assert all(row[0] == 'SRS011830' for row in sdrf[1:])
    assert 'SRR037073' in str(sdrf)


def test_custom_factor_original_name_is_preserved():
    from meta_standards_converter.magetab.parser import AEParser
    from meta_standards_converter.sources.magetab import MAGETabInput, TextResource
    from tests.converters.test_ae2json import resolved_input
    source = resolved_input(idf='Investigation Title\tTest\nInvestigation Accession\tE-MTAB-1\nExperimental Factor Name\tcustom exposure\nExperimental Factor Type\tcompound\n', sdrfs=['Source Name\tCharacteristics[organism]\tFactor Value[custom exposure]\na\tHomo sapiens\tx\n'])
    typed = AEParser().parse(source)
    variable = typed.to_mapping()['series']['variable'][0]
    assert variable['name'] == 'custom exposure'
    assert variable['factor'] == 'other'
    rendered = AEConstructor().miniml2magetab(typed, platform_handler='generic')
    assert next(r[1] for r in rendered if r[0] == 'Experimental Factor Name') == 'custom exposure'
