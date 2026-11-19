# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Shared ae builders and deterministic collaborators."""

import json
import requests
from meta_standards_converter.helpers.request_helper import HostRequestGate, RateLimitedRequester, RequestSettings
from meta_standards_converter.retrieval import RetrievalPolicy
from meta_standards_converter.sources.magetab import AEWebFetcher
from types import SimpleNamespace
from unittest.mock import Mock
from meta_standards_converter.converters.archive_results import ArchiveImportResult, StudyImportOutcome
from meta_standards_converter.miniml import MINiMLCodec
from meta_standards_converter.sources.errors import MetadataSourceError
from meta_standards_converter.provider_issues import record
from tests.support.conversion import package
from meta_standards_converter.sources.magetab import (  # noqa: E402
    MAGETabInput,
    TextResource,
)


MAPPINGS = [("E-MTAB-16245", "ERP185477", "PRJEB104172"),
            ("E-MTAB-13162", "ERP149729", "PRJEB64550")]


ACCESSION = "E-MTAB-16245"


BASE = "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/245/" + ACCESSION


INFO = {"httpLink": BASE, "isPublic": True}


def response(payload, status=200):
    result = requests.Response()
    result.status_code = status
    result._content = payload if isinstance(payload, bytes) else json.dumps(payload).encode()
    result._content_consumed = True
    result.headers["Content-Length"] = str(len(result.content))
    return result


def fetcher(tmp_path, reply, *, retries=1):
    calls = []
    def get(url, **kwargs):
        calls.append((url, kwargs))
        return reply(url, kwargs)
    requester = RateLimitedRequester("biostudies", get=get,
        settings=RequestSettings(max_retries=retries, request_delay=0,
                                 backoff_base=0, backoff_max=0),
        host_gate=HostRequestGate(tmp_path))
    policy = RetrievalPolicy(resolver=lambda *a, **k: [(None, None, None, None, ("93.184.216.34", 443))])
    return AEWebFetcher(requester=requester, retrieval_policy=policy), calls


def services(ae, study, primary, *, failure='source_inaccessible', optional_failure=False):
    original = MetadataSourceError(failure, ae, http_status=403 if failure == 'source_inaccessible' else None)
    def deny(*args, **kwargs):
        record(original, 'biostudies')
        raise original
    ae_source = SimpleNamespace(load=Mock(side_effect=deny))
    data = package(study, 'SAMEA1')
    data['source']['format'] = 'ENA'
    data['database'] = [{'iid': db, 'name': db} for db in ['ENA', 'ArrayExpress', 'BioProject']]
    data['series']['iid'] = ae
    data['series']['accession'] = [{'value': ae}, {'value': study}, {'value': primary}]
    data['series']['relation'] = [{'type': 'ArrayExpress', 'target': ae}]
    def load(*args, **kwargs):
        if optional_failure:
            import requests
            response = requests.Response(); response.status_code = 500
            record(requests.HTTPError(response=response), 'biosamples')
        return ArchiveImportResult(study, [StudyImportOutcome(study, primary, 'complete', MINiMLCodec().decode(data).package)])
    ena = SimpleNamespace(load=Mock(side_effect=load), resolve_arrayexpress=Mock(return_value=SimpleNamespace(
        study=study, project=primary, evidence=(f'https://www.ebi.ac.uk/biostudies/arrayexpress/studies/{ae}',))))
    class Enricher:
        def enrich_selected(self, data, **kwargs): return data
    return {'ae_source': ae_source, 'ena_source': ena, 'sra_source': ena, 'enricher': Enricher()}, data


IDF = """MAGE-TAB Version\t1.1
Investigation Title\tExample study
Investigation Accession\tE-MTAB-1
Comment[SecondaryAccession]\tGSE123
Comment[SecondaryAccessionTermSourceRef]\tGEO
Experimental Design\tRNA-seq
Experimental Factor Name\tdisease
Experimental Factor Type\tdisease
Person Last Name\tDoe
Person First Name\tJane
Person Email\tjane@example.org
Person Affiliation\tExample Institute
Public Release Date\t2025-01-02
PubMed ID\t12345
Publication DOI\t10.1/example
Publication Author List\tDoe J
Publication Title\tExample paper
Publication Status\tpublished
Status Term Source Ref\tEFO
Status Term Accession Number\tEFO:0000001
Protocol Name\tP-collect\tP-extract
Protocol Type\tsample collection protocol\tnucleic acid extraction protocol
Protocol Description\tCollect samples\tExtract material
Protocol Type Term Source REF\tEFO\tEFO
Protocol Type Term Accession Number\tEFO:0005518\tEFO:0002944
SDRF File\tstudy.sdrf.txt
Term Source Name\tEFO
Term Source File\thttps://www.ebi.ac.uk/efo/
Mystery Row\tkeep me
"""


SDRF_HEADER = [
    "Source Name",
    "Comment[Sample_title]",
    "Characteristics[organism]",
    "Term Source REF",
    "Characteristics[disease]",
    "Factor Value[disease]",
    "Protocol REF",
    "Extract Name",
    "Material Type",
    "Array Design REF",
    "Technology Type",
    "Comment[LIBRARY_SOURCE]",
    "Comment[LIBRARY_STRATEGY]",
    "Comment[ENA_RUN]",
    "Comment[FASTQ_URI]",
    "Comment[MD5]",
    "Mystery Column",
]


def sdrf(rows=None):
    rows = rows or [
        ["GSM1", "Sample one", "Homo sapiens", "EFO", "case", "case", "P-extract", "GSM1", "RNA", "A-TEST-1", "sequencing assay", "TRANSCRIPTOMIC", "RNA-SEQ", "ERR1", "https://example/1.fastq.gz", "aaa", "x"],
        ["GSM1", "Sample one", "Homo sapiens", "EFO", "case", "case", "P-extract", "GSM1", "RNA", "A-TEST-1", "sequencing assay", "TRANSCRIPTOMIC", "RNA-SEQ", "ERR1", "https://example/2.fastq.gz", "bbb", "y"],
    ]
    return "\n".join("\t".join(row) for row in [SDRF_HEADER, *rows]) + "\n"


def resolved_input(idf=IDF, sdrfs=None):
    return MAGETabInput(
        idf=TextResource("study.idf.txt", idf, "memory:idf"),
        sdrfs=tuple(
            TextResource(f"study{index}.sdrf.txt", text, f"memory:sdrf:{index}")
            for index, text in enumerate(sdrfs or [sdrf()], start=1)
        ),
        source="E-MTAB-1",
        source_kind="accession",
    )
