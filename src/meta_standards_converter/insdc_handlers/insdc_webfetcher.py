# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
'''
Handler class for INSDC/SRA metadata.
'''

import re
import xml.etree.ElementTree as ET
import requests
from urllib.parse import urlparse

from meta_standards_converter.helpers.request_helper import RateLimitedRequester


class INSDCWebfetcher():
    def __init__(
        self,
        ncbi_requester=None,
        ena_requester=None,
        ncbi_request_settings=None,
        ena_request_settings=None,
    ):
        self.ncbi_requester = ncbi_requester or RateLimitedRequester(
            service="ncbi_eutils",
            settings=ncbi_request_settings,
        )
        self.ena_requester = ena_requester or RateLimitedRequester(
            service="ena_portal",
            settings=ena_request_settings,
        )

    def _extract_sra(self, sra: str) -> list:
        '''
        Extracts sra accession within substring
        '''
        pattern = r'\b[SED]R[RXSP]\d+'
        matches = re.findall(pattern, sra, flags=re.IGNORECASE)

        return matches

    def _ncbi_nrx(self, nrx: str) -> list:
        '''
        lookup nrx accession to get nrr accessions
        '''
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id={nrx}&retmode=xml"
        response = self.ncbi_requester.get(url)
        response.raise_for_status()

        root = ET.fromstring(response.content)
        

        return root

    def fetch_sra_runs(self, accession: str) -> list:
        root = self._ncbi_nrx(nrx=accession)
        ena_fastqs_by_run = self.fetch_ena_fastq_files(accession=accession)

        records = []
        for package in root.findall(".//EXPERIMENT_PACKAGE"):
            study_id = self._parse_sra_study_id(package=package)
            experiment = package.find("./EXPERIMENT")
            sample = package.find("./SAMPLE")
            library = self._parse_sra_library(experiment=experiment)
            sample_ids = self._parse_sra_sample_ids(sample=sample)
            instrument_model = self._parse_sra_instrument_model(experiment=experiment)
            experiment_id = self._element_accession(node=experiment)

            for run in package.findall(".//RUN_SET/RUN"):
                run_accession = self._element_accession(node=run)
                fastq_files = ena_fastqs_by_run.get(run_accession) or self._parse_sra_fastqs(run=run)
                records.append({
                    **sample_ids,
                    **library,
                    "study": study_id,
                    "experiment": experiment_id,
                    "run": run_accession,
                    "scan_name": run.get("alias") or run_accession,
                    "instrument_model": instrument_model,
                    "fastq_files": fastq_files,
                    "submitted_file_name": fastq_files[0].get("filename") if fastq_files else None,
                    "md5": fastq_files[0].get("md5") if fastq_files else None,
                    "read_lengths": [
                        read.get("average")
                        for read in run.findall("./Statistics/Read")
                        if read.get("average")
                    ],
                })

        return records

    def _parse_sra_study_id(self, package: ET.Element):
        study = package.find("./STUDY")
        study_id = self._element_accession(node=study)
        if study_id:
            return study_id

        study_ref = package.find(".//STUDY_REF")
        study_id = self._element_accession(node=study_ref)
        if study_id:
            return study_id

        if study is not None:
            for external_id in study.findall(".//EXTERNAL_ID"):
                value = self._clean_sdrf_text(external_id.text)
                if value:
                    return value
        return None

    def fetch_ena_file_report(self, accession: str) -> list:
        url = "https://www.ebi.ac.uk/ena/portal/api/filereport"
        params = {
            "accession": accession,
            "result": "read_run",
            "fields": "run_accession,fastq_ftp,fastq_md5,fastq_bytes",
            "format": "json",
        }
        response = self.ena_requester.get(url, params=params)
        response.raise_for_status()
        rows = response.json()
        return rows if isinstance(rows, list) else []

    def fetch_ena_fastq_files(self, accession: str) -> dict:
        try:
            rows = self.fetch_ena_file_report(accession=accession)
        except (requests.RequestException, ValueError):
            return {}

        fastqs_by_run = {}
        for row in rows:
            if not isinstance(row, dict):
                continue
            run_accession = self._clean_sdrf_text(row.get("run_accession"))
            fastqs = self._parse_ena_fastq_report(row=row)
            if run_accession and fastqs:
                fastqs_by_run[run_accession] = fastqs
        return fastqs_by_run

    def _parse_ena_fastq_report(self, row: dict) -> list:
        uris = self._split_ena_file_field(row.get("fastq_ftp"))
        md5s = self._split_ena_file_field(row.get("fastq_md5"))
        bytes_values = self._split_ena_file_field(row.get("fastq_bytes"))

        fastqs = []
        for index, uri in enumerate(uris):
            normalized_uri = self._normalize_ena_ftp_uri(uri)
            fastqs.append({
                "filename": self._filename_from_uri(normalized_uri),
                "uri": normalized_uri,
                "md5": md5s[index] if index < len(md5s) else None,
                "bytes": bytes_values[index] if index < len(bytes_values) else None,
            })
        return fastqs

    def _parse_sra_library(self, experiment: ET.Element) -> dict:
        library = experiment.find(".//LIBRARY_DESCRIPTOR") if experiment is not None else None
        layout_node = library.find("./LIBRARY_LAYOUT") if library is not None else None
        layout = None
        if layout_node is not None and list(layout_node):
            layout = self._strip_ns(tag=list(layout_node)[0].tag).upper()

        return {
            "library_layout": self._clean_sdrf_text(layout),
            "library_selection": self._find_text(node=library, path="./LIBRARY_SELECTION"),
            "library_source": self._find_text(node=library, path="./LIBRARY_SOURCE"),
            "library_strategy": self._find_text(node=library, path="./LIBRARY_STRATEGY"),
        }

    def _parse_sra_sample_ids(self, sample: ET.Element) -> dict:
        ids = {
            "sample": self._element_accession(node=sample),
            "biosample": None,
            "geo_sample": None,
        }
        if sample is None:
            return ids

        for external_id in sample.findall(".//EXTERNAL_ID"):
            namespace = (external_id.get("namespace") or "").lower()
            if namespace == "biosample":
                ids["biosample"] = self._clean_sdrf_text(external_id.text)
            elif namespace == "geo":
                ids["geo_sample"] = self._clean_sdrf_text(external_id.text)

        return ids

    def _parse_sra_instrument_model(self, experiment: ET.Element):
        if experiment is None:
            return None
        for instrument_model in experiment.findall(".//PLATFORM/*/INSTRUMENT_MODEL"):
            if instrument_model.text:
                return self._clean_sdrf_text(instrument_model.text)
        return None

    def _parse_sra_fastqs(self, run: ET.Element) -> list:
        fastqs = []
        for sra_file in run.findall(".//SRAFile"):
            if (sra_file.get("semantic_name") or "").lower() != "fastq":
                continue
            filename = self._clean_sdrf_text(sra_file.get("filename"))
            uri = sra_file.get("url")
            alternative = sra_file.find("./Alternatives")
            if not uri and alternative is not None:
                uri = alternative.get("url")
            fastqs.append({
                "filename": filename,
                "uri": self._clean_sdrf_text(uri),
                "md5": self._clean_sdrf_text(sra_file.get("md5")),
            })
        return fastqs

    def _split_ena_file_field(self, value) -> list:
        value = self._clean_sdrf_text(value)
        if not value:
            return []
        return [
            part
            for part in (self._clean_sdrf_text(part) for part in value.split(";"))
            if part
        ]

    def _normalize_ena_ftp_uri(self, uri):
        uri = self._clean_sdrf_text(uri)
        if not uri:
            return None
        if uri.startswith(("ftp://", "http://", "https://")):
            return uri
        return f"ftp://{uri}"

    def _filename_from_uri(self, uri):
        uri = self._clean_sdrf_text(uri)
        if not uri:
            return None
        parsed = urlparse(uri)
        path = parsed.path or uri
        return path.rstrip("/").rsplit("/", 1)[-1] or None

    def _element_accession(self, node: ET.Element):
        if node is None:
            return None
        return self._clean_sdrf_text(node.get("accession") or self._find_text(node=node, path=".//PRIMARY_ID"))

    def _find_text(self, node: ET.Element, path: str):
        if node is None:
            return None
        return self._clean_sdrf_text(node.findtext(path))

    def _strip_ns(self, tag: str) -> str:
        return tag.split("}", 1)[-1] if "}" in tag else tag

    def _clean_sdrf_text(self, value):
        if value is None:
            return None
        return " ".join(str(value).replace("\t", " ").replace("\n", " ").split())
