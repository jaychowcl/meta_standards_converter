# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Native metadata contracts using the independently vendored provider responses."""
from pathlib import Path
import copy
import json
import xml.etree.ElementTree as ET

import pytest

from meta_standards_converter.sources.archive_support import StudyRecords, StudySeed
from meta_standards_converter.miniml.sra_parser import SRAParser
from meta_standards_converter.miniml.ena_parser import ENAParser
from meta_standards_converter.miniml import MINiMLCodec

DOCS = Path(__file__).resolve().parents[1] / "docs"


def fixture_records(provider, experiment="SRX7812918"):
    study, project = (("SRP250911", "PRJNA609050") if experiment == "SRX7812918"
                      else ("SRP002056", "PRJNA123835"))
    records = StudyRecords(StudySeed(study, project if provider == "ena" else study))
    directory = DOCS / provider / "fixtures" / experiment
    if provider == "sra":
        records.xml.append(ET.parse(directory / "sra-efetch.xml").getroot())
        for kind in ("biosample", "bioproject"):
            records.xml.append(ET.parse(directory / f"{kind}-efetch.xml").getroot())
    else:
        records.xml.extend(ET.parse(directory / f"{kind}.xml").getroot()
                           for kind in ("study", "sample", "experiment", "run"))
        records.indexed["read_run"] = json.loads((directory / "file-report.json").read_text())
    return records


@pytest.mark.parametrize("provider,parser", [("sra", SRAParser), ("ena", ENAParser)])
def test_native_identity_biology_and_relationships(provider, parser):
    package = parser().parse(fixture_records(provider))
    data = package.to_mapping()
    assert data["series"]["iid"] == ("SRP250911" if provider == "sra" else "PRJNA609050")
    assert len(data["sample"]) == 1
    sample = data["sample"][0]
    assert sample["iid"] == ("SRS6225446" if provider == "sra" else "SAMN14218700")
    assert sample["channel"][0]["organism"][0]["taxid"] == "408170"
    attrs = sample["channel"][0]["characteristics"]
    assert any(a["name"] == "host" and a["value"] == "Homo sapiens" for a in attrs)
    assert sample["sra_run"][0]["run"] == "SRR11192680"
    assert sample["sra_run"][0]["library_layout"] == "PAIRED"
    assert data["series"]["assay_paths"]
    assert data["extensions"]["insdc"]["version"] == "1.0"
    assert data["extensions"]["insdc"]["records"]
    assert "completeness" not in json.dumps(data)
    MINiMLCodec().decode(data, strict=True)


def test_ena_parallel_files_preserve_duplicates_and_short_hash_columns():
    records = fixture_records("ena")
    records.indexed["read_run"] = [{"run_accession": "SRR11192680",
        "fastq_ftp": "ftp.sra.ebi.ac.uk/a.fastq.gz;ftp.sra.ebi.ac.uk/a.fastq.gz",
        "fastq_bytes": "10;20", "fastq_md5": "a" * 32}]
    files = ENAParser().parse(records).to_mapping()["sample"][0]["sra_run"][0]["fastq_files"]
    assert len(files) == 2
    assert files[0]["uri"] == files[1]["uri"]
    assert files[1].get("md5") is None
    assert files[1]["bytes"] == "20"


def test_repeating_attributes_units_and_missing_literals_survive():
    records = fixture_records("sra")
    attrs = records.xml[0].find(".//SAMPLE/SAMPLE_ATTRIBUTES")
    for value in ("missing: restricted access", "5"):
        attr = ET.SubElement(attrs, "SAMPLE_ATTRIBUTE")
        for tag, text in (("TAG", "dose"), ("VALUE", value), ("UNITS", "mg")):
            ET.SubElement(attr, tag).text = text
    values = SRAParser().parse(records).to_mapping()["sample"][0]["channel"][0]["characteristics"]
    assert [(a["value"], a["unit"]["value"]) for a in values if a["name"] == "dose"] == [
        ("missing: restricted access", "mg"), ("5", "mg")]


def test_multiple_libraries_do_not_create_samples_or_choose_scalar():
    records = fixture_records("sra")
    original = records.xml[0].find("EXPERIMENT_PACKAGE")
    second = copy.deepcopy(original)
    experiment = second.find("EXPERIMENT")
    experiment.set("accession", "SRX999999")
    experiment.find("IDENTIFIERS/PRIMARY_ID").text = "SRX999999"
    experiment.find("DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_STRATEGY").text = "WGS"
    run = second.find("RUN_SET/RUN")
    run.set("accession", "SRR999999")
    run.find("EXPERIMENT_REF").set("accession", "SRX999999")
    records.xml[0].append(second)
    sample = SRAParser().parse(records).to_mapping()["sample"][0]
    assert len(sample["sra_run"]) == 2
    assert "library_strategy" not in sample
    assert {r["library_strategy"] for r in sample["sra_run"]} == {"AMPLICON", "WGS"}


def test_sra_archives_are_not_fastq_and_statistics_remain_scoped():
    data = SRAParser().parse(fixture_records("sra")).to_mapping()
    run = data["sample"][0]["sra_run"][0]
    assert run["read_lengths"]
    assert all(not f["filename"].endswith("sralite") for f in run.get("fastq_files", []))
    assert any("sralite" in str(r) for r in data["extensions"]["insdc"]["records"])


def test_linked_biosample_attributes_and_project_publications_are_mapped():
    records = fixture_records('sra', 'SRX017289')
    biosample = records.xml[1].find('BioSample')
    attr = ET.SubElement(biosample.find('Attributes'), 'Attribute', attribute_name='tissue')
    attr.text = 'explicit tissue'
    records.xml.append(ET.fromstring('''<PubmedArticleSet><PubmedArticle><MedlineCitation><PMID>20133686</PMID><Article><ArticleTitle>Publication title</ArticleTitle><AuthorList><Author><LastName>Example</LastName><ForeName>Alice</ForeName></Author></AuthorList></Article></MedlineCitation><PubmedData><ArticleIdList><ArticleId IdType="doi">10.1/example</ArticleId></ArticleIdList></PubmedData></PubmedArticle></PubmedArticleSet>'''))
    data = SRAParser().parse(records).to_mapping()
    assert any(c['name'] == 'tissue' for c in data['sample'][0]['channel'][0]['characteristics'])
    assert data['series']['pubmed_publication'][0]['doi'] == '10.1/example'


def test_analysis_files_are_sample_scoped_only_with_explicit_references():
    records = fixture_records('ena')
    records.xml.append(ET.fromstring('''<ANALYSIS_SET><ANALYSIS accession="ERZ1"><STUDY_REF accession="SRP250911"/><SAMPLE_REF accession="SRS6225446"/><TITLE>assembly</TITLE><PROTOCOL>Assembly method</PROTOCOL><ANALYSIS_TYPE><SEQUENCE_ASSEMBLY><PROGRAM>assembler</PROGRAM><COVERAGE>30</COVERAGE></SEQUENCE_ASSEMBLY></ANALYSIS_TYPE><FILES><FILE filename="https://example.org/assembly.fasta" filetype="fasta" checksum_method="SHA256" checksum="abc"/></FILES></ANALYSIS></ANALYSIS_SET>'''))
    data = ENAParser().parse(records).to_mapping()
    link = data['sample'][0]['supplementary_data'][0]
    assert link['value'].endswith('assembly.fasta')
    assert 'checksum' not in link
    assert data['series']['protocols'][0]['software'] == ['assembler']


def test_ena_missing_file_position_is_not_dropped_or_shifted():
    from meta_standards_converter.miniml.insdc_support import files_from_ena
    files = files_from_ena({'fastq_ftp': 'host/a;;host/c', 'fastq_md5': 'a;b;c', 'fastq_bytes': '1;2;3;4'})
    assert len(files) == 4
    assert files[1].get('uri') is None and files[1]['md5'] == 'b'
    assert files[2]['uri'].endswith('/c') and files[2]['bytes'] == '3'
    assert files[3].get('uri') is None and files[3]['bytes'] == '4'


def test_parser_scope_excludes_unrelated_study_samples_and_runs():
    records = fixture_records('sra')
    extra = copy.deepcopy(records.xml[0].find('EXPERIMENT_PACKAGE'))
    extra.find('EXPERIMENT/STUDY_REF').set('accession', 'SRP999')
    extra.find('EXPERIMENT').set('accession', 'SRX999')
    extra.find('SAMPLE').set('accession', 'SRS999')
    extra.find('RUN_SET/RUN').set('accession', 'SRR999')
    records.xml[0].append(extra)
    assert len(SRAParser().parse(records).samples) == 1
    records = fixture_records('ena')
    records.xml[2].find('.//EXPERIMENT/STUDY_REF').set('accession', 'SRP999')
    parsed = ENAParser().parse(records)
    assert not any(s.sra_runs for s in parsed.samples)
    assert records.issues


def test_sra_explicit_experiment_pool_members_receive_same_run():
    records = fixture_records('sra')
    package = records.xml[0].find('EXPERIMENT_PACKAGE')
    second = copy.deepcopy(package.find('SAMPLE')); second.set('accession', 'SRS2')
    second.find('IDENTIFIERS/PRIMARY_ID').text = 'SRS2'
    records.xml.append(ET.fromstring('<SAMPLE_SET/>')); records.xml[-1].append(second)
    desc = package.find('EXPERIMENT/DESIGN/SAMPLE_DESCRIPTOR')
    pool = ET.SubElement(desc, 'POOL')
    for accession in ('SRS6225446', 'SRS2'):
        ET.SubElement(pool, 'MEMBER', accession=accession)
    package.find('RUN_SET/RUN').remove(package.find('RUN_SET/RUN/Pool'))
    data = SRAParser().parse(records).to_mapping()
    assert {s['iid'] for s in data['sample']} == {'SRS6225446', 'SRS2'}
    assert all(s['sra_run'][0]['run'] == 'SRR11192680' for s in data['sample'])


def test_analysis_portal_files_and_run_associations_project_to_paths():
    records = fixture_records('ena')
    records.indexed['analysis'] = [{'analysis_accession': 'ERZ1', 'study_accession': 'PRJNA609050', 'run_accession': 'SRR11192680', 'submitted_ftp': 'host/assembly.fa.gz', 'submitted_format': 'fasta'}]
    records.xml.append(ET.fromstring('<ANALYSIS_SET><ANALYSIS accession="ERZ1"><STUDY_REF accession="SRP250911"/><RUN_REF accession="SRR11192680"/><PROTOCOL>Assembly pipeline</PROTOCOL><FILES><FILE filename="relative/assembly.fa.gz" filetype="fasta"/></FILES></ANALYSIS></ANALYSIS_SET>'))
    data = ENAParser().parse(records).to_mapping()
    assert data['sample'][0]['supplementary_data'][0]['value'] == 'ftp://host/assembly.fa.gz'
    paths = data['series']['assay_paths']
    assert any(s.get('link', {}).get('value') == 'ftp://host/assembly.fa.gz' for p in paths for s in p['steps'])
    assert any(s.get('protocol_ref') == 'ERZ1:analysis' for p in paths for s in p['steps'])


def test_retained_records_have_entity_accessions():
    from meta_standards_converter.miniml.insdc_support import retained
    records = fixture_records('sra')
    records.indexed['read_run'] = [{'run_accession': 'SRR1'}]
    data = retained('sra', records)
    assert any(r['kind'] == 'RUN' and r['accession'] == 'SRR11192680' for r in data['records'])
    assert next(r for r in data['records'] if r['kind'] == 'read_run')['accession'] == 'SRR1'


def test_ena_indexed_evidence_survives_missing_browser_records():
    records = fixture_records('ena')
    records.xml = records.xml[:1]
    records.indexed['read_experiment'] = [{'experiment_accession': 'SRX7812918', 'study_accession': 'PRJNA609050', 'secondary_study_accession': 'SRP250911', 'sample_accession': 'SAMN14218700', 'secondary_sample_accession': 'SRS6225446', 'library_layout': 'PAIRED', 'library_strategy': 'AMPLICON', 'instrument_model': 'MiSeq'}]
    records.indexed['sample'] = [{'sample_accession': 'SAMN14218700', 'secondary_sample_accession': 'SRS6225446', 'scientific_name': 'human gut metagenome', 'tax_id': '408170', 'sample_title': 'sample'}]
    records.indexed['read_run'][0]['experiment_accession'] = 'SRX7812918'
    records.issues.append('Browser unavailable')
    data = ENAParser().parse(records).to_mapping()
    assert len(data['sample']) == 1
    assert data['sample'][0]['sra_run'][0]['run'] == 'SRR11192680'
    assert data['sample'][0]['channel'][0]['organism'][0]['taxid'] == '408170'
    assert data['series']['assay_paths']


def test_explicit_molecule_host_taxon_and_indexed_dates():
    records = fixture_records('ena')
    attrs = records.xml[1].find('.//SAMPLE_ATTRIBUTES')
    for name, value in [('molecule', 'total RNA'), ('host_taxid', '9606')]:
        attr = ET.SubElement(attrs, 'SAMPLE_ATTRIBUTE')
        ET.SubElement(attr, 'TAG').text = name
        ET.SubElement(attr, 'VALUE').text = value
    records.indexed['study'] = [{'study_accession': 'PRJNA609050', 'first_public': '2020-03', 'last_updated': '2024-03-01'}]
    data = ENAParser().parse(records).to_mapping()
    channel = data['sample'][0]['channel'][0]
    assert channel['molecule']['value'] == 'total RNA'
    host = next(c for c in channel['characteristics'] if c['name'] == 'host')
    assert host['term_accession_number'] == '9606'
    assert data['series']['status'][0]['release_date'] == '2020-03'


def test_pool_default_member_and_unresolved_member_are_not_silent():
    records = fixture_records('sra')
    package = records.xml[0].find('EXPERIMENT_PACKAGE')
    run = package.find('RUN_SET/RUN'); run.remove(run.find('Pool'))
    pool = ET.SubElement(package.find('EXPERIMENT/DESIGN/SAMPLE_DESCRIPTOR'), 'POOL')
    ET.SubElement(pool, 'MEMBER', accession='SRS999')
    ET.SubElement(pool, 'DEFAULT_MEMBER', accession='SRS6225446')
    data = SRAParser().parse(records).to_mapping()
    assert data['sample'][0]['sra_run'][0]['run'] == 'SRR11192680'
    assert any('SRS999' in issue for issue in records.issues)
