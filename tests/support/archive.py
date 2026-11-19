# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Shared archive builders and deterministic collaborators."""

from pathlib import Path
import json
import xml.etree.ElementTree as ET
from meta_standards_converter.sources.archive_support import StudyRecords, StudySeed
from copy import deepcopy
from meta_standards_converter.miniml import MINiMLCodec
from meta_standards_converter.miniml.sra_parser import SRAParser
from meta_standards_converter.magetab.parser import AEParser
from tests.support.ae import resolved_input
import requests
from meta_standards_converter.miniml.ena_parser import ENAParser


DOCS = Path(__file__).resolve().parents[2] / "docs"


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


def native():
    return SRAParser().parse(fixture_records('sra'))


def linked(package, accession, title):
    data = package.to_mapping()
    data['source'] = {'format': 'MAGE-TAB' if accession.startswith('E-') else 'GEO'}
    data['series']['iid'] = accession
    data['series']['accession'].append({'value': accession, 'database': 'ArrayExpress' if accession.startswith('E-') else 'GEO'})
    data['series']['title'] = title
    data['sample'][0]['title'] = title
    old_iid = data['sample'][0]['iid']
    data['sample'][0]['iid'] = 'GSM1'
    for path in data['series'].get('assay_paths', []):
        for step in path['steps']:
            if step.get('sample_ref') == old_iid: step['sample_ref'] = 'GSM1'
    data['sample'][0]['channel'][0]['characteristics'] = [{'name': 'host', 'value': title, 'term_accession_number': 'ONT:1', 'term_source_ref': 'ONT'}]
    # This fixture clones generated native paths: update their same host projection too.
    for path in data['series'].get('assay_paths', []):
        for step in path['steps']:
            if step['kind'] in ('source', 'sample'):
                step['characteristics'] = [c for c in step.get('characteristics', []) if c['name'] != 'host'] + deepcopy(data['sample'][0]['channel'][0]['characteristics'])
    return MINiMLCodec().decode(data).package


def workflow():
    header=['Source Name','Characteristics[organism]','Characteristics[Age]','Unit[TimeUnit]','Term Source REF','Term Accession Number','Material Type','Protocol REF','Extract Name','Material Type','Protocol REF','Extract Name','Material Type','Assay Name','Comment[ENA_EXPERIMENT]','Scan Name','Comment[ENA_RUN]','Factor Value[AGE]','Unit[TimeUnit]']
    row=['embryo','Danio rerio','1','days','UO','UO:0000033','organism part','P1','RNA','total RNA','P2','mRNA','3 prime end of mRNA','assay','SRX7812918','scan','SRR11192680','1','d']
    idf='Investigation Accession\tE-MTAB-1\nInvestigation Title\tStudy\nComment[SecondaryAccession]\tSRP250911\nExperimental Factor Name\tAGE\nExperimental Factor Type\tage\nProtocol Name\tP1\tP2\nProtocol Type\textraction\tlibrary construction\nProtocol Description\textract RNA\tselect mRNA\nPerson First Name\tJane\nPerson Last Name\tDoe\nPerson Email\tjane@example.org\n'
    return AEParser().parse(resolved_input(idf=idf,sdrfs=['\t'.join(header)+'\n'+'\t'.join(row)+'\n']))


def enriched_native():
    data=native().to_mapping()
    data['sample'][0]['channel'][0]['organism']=[{'value':'Danio rerio','taxid':'7955'}]
    return MINiMLCodec().decode(data).package


def nodes(value):
    if isinstance(value,dict):
        yield value
        for v in value.values(): yield from nodes(v)
    elif isinstance(value,list):
        for v in value: yield from nodes(v)


class PubMed:
    def __init__(self, fail=False): self.calls=[]; self.fail=fail
    def pubmed_summary(self, pubmed_id):
        self.calls.append(pubmed_id)
        if self.fail: raise requests.ConnectionError('unavailable')
        return ('10.1/fetched','Fetched authors','Fetched title','published','EFO','EFO_0001796')


class NoRuns:
    def extract_sra_accessions(self, *args): raise AssertionError('native sequencing lookup')
    def fetch_sra_runs(self, *args, **kwargs): raise AssertionError('native sequencing lookup')


def records():
    return StudyRecords(StudySeed('ERP1', 'PRJEB1'))


def missing_sample_records():
    value = records()
    value.xml = [ET.fromstring('''<STUDY_SET><STUDY accession="ERP1"><DESCRIPTOR><STUDY_TITLE>Recovery fixture</STUDY_TITLE></DESCRIPTOR></STUDY></STUDY_SET>'''),
        ET.fromstring('''<EXPERIMENT_SET><EXPERIMENT accession="ERX1"><STUDY_REF accession="ERP1"/><DESIGN><SAMPLE_DESCRIPTOR accession="ERS1"/><LIBRARY_DESCRIPTOR><LIBRARY_STRATEGY>RNA-Seq</LIBRARY_STRATEGY><LIBRARY_SOURCE>TRANSCRIPTOMIC</LIBRARY_SOURCE><LIBRARY_SELECTION>cDNA</LIBRARY_SELECTION><LIBRARY_LAYOUT><SINGLE/></LIBRARY_LAYOUT></LIBRARY_DESCRIPTOR></DESIGN><PLATFORM><ILLUMINA><INSTRUMENT_MODEL>Illumina HiSeq 2000</INSTRUMENT_MODEL></ILLUMINA></PLATFORM></EXPERIMENT></EXPERIMENT_SET>'''),
        ET.fromstring('<RUN_SET><RUN accession="ERR1"><EXPERIMENT_REF accession="ERX1"/></RUN></RUN_SET>')]
    row = {'study_accession': 'PRJEB1', 'secondary_study_accession': 'ERP1',
        'run_accession': 'ERR1', 'experiment_accession': 'ERX1', 'sample_accession': 'SAMEA1',
        'secondary_sample_accession': 'ERS1', 'sample_title': 'sample one',
        'scientific_name': 'Homo sapiens', 'tax_id': '9606',
        'fastq_ftp': 'ftp.sra.ebi.ac.uk/one.fastq.gz', 'fastq_md5': 'a'*32, 'fastq_bytes': '42'}
    value.indexed = {'read_run': [row], 'read_experiment': [deepcopy(row)]}
    value.linked = [{'provider': 'biosamples', 'kind': 'sample', 'accession': 'SAMEA1',
        'metadata': {'accession': 'SAMEA1', 'characteristics': {'tissue': [{'text': 'pancreas'}]}}},
        {'provider': 'biosamples', 'kind': 'sample', 'accession': 'SAMEA999',
        'metadata': {'accession': 'SAMEA999', 'characteristics': {'tissue': [{'text': 'unrelated'}]}}}]
    return value


def legacy():
    data=ENAParser().parse(fixture_records('ena')).to_mapping()
    sample=data['sample'][0]
    dates=[{'name':'ENA-FIRST-PUBLIC','value':'2020-01-01'},
           {'name':'ENA-FIRST-PUBLIC','value':'2020-01-01'},
           {'name':'ENA-LAST-UPDATE','value':'2020-02-03T12:01:02Z'},
           {'name':'INSDC last update','value':'2020-02-04'},
           {'name':'collection_date','value':'2019-12'}]
    sample['channel'][0]['characteristics'].extend(deepcopy(dates))
    sample['status']=[{'database':'ENA','release_date':'2020-01-01'},
                      {'database':'BioSample','last_update_date':'2020-04-05'}]
    for path in data['series']['assay_paths']:
        path['steps'][0]['characteristics'].extend(deepcopy(dates))
    return data


def package(provider='SRA', layout='SINGLE'):
    accession = {'SRA': 'SRP1', 'GEO': 'GSE1', 'MAGE-TAB': 'E-MTAB-1'}[provider]
    data = {'miniml_schema_version': '3.0', 'source': {'format': provider},
        'database': [{'iid': p, 'name': p} for p in ['SRA','GEO','ArrayExpress']],
        'series': {'iid': accession, 'title': provider, 'accession': [{'value': accession}],
                   'sample_ref': [{'ref': 'SRS1'}], 'relation': [{'type': 'SRA', 'target': 'SRP1'}],
                   'assay_paths': [{'steps': [{'kind': 'source', 'name': 'SRS1', 'sample_ref': 'SRS1'},
                       {'kind': 'assay', 'name': 'SRX1'}, {'kind': 'scan', 'name': 'SRR1',
                        'comments': [{'name': 'LIBRARY_LAYOUT', 'value': layout}]}]}]},
        'sample': [{'iid': 'SRS1', 'accession': [{'value': 'SRS1', 'database': 'SRA'}],
                    'channel': [{'organism': [{'value': 'Homo sapiens', 'taxid': '9606'}]}],
                    'sra_run': [{'run': 'SRR1', 'experiment': 'SRX1', 'sample': 'SRS1',
                                 'library_strategy': 'RNA-Seq', 'library_layout': layout}]}]}
    return data


class NoPublication:
    publication_issues = []
    def enrich_selected(self, package, **kwargs): return package
