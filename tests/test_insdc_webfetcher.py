# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
import os
import sys
import unittest
import xml.etree.ElementTree as ET
from unittest.mock import Mock

import requests


ROOT = os.path.dirname(os.path.dirname(__file__))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from meta_standards_converter.ae_handlers.ae_sdrf_handlers import SDRFConstructor  # noqa: E402
from meta_standards_converter.insdc_handlers.insdc_webfetcher import INSDCWebfetcher  # noqa: E402


class TestINSDCWebfetcher(unittest.TestCase):
    def test_extract_sra_finds_sra_style_accessions(self):
        self.assertEqual(
            ["SRX1", "ERR22", "drr333"],
            INSDCWebfetcher()._extract_sra("SRX1 https://example/ERR22 and drr333"),
        )

    def test_ncbi_nrx_uses_ncbi_requester(self):
        response = Mock(content=b"<EXPERIMENT_PACKAGE_SET />")
        response.raise_for_status = Mock()
        requester = Mock()
        requester.get.return_value = response

        root = INSDCWebfetcher(ncbi_requester=requester)._ncbi_nrx("SRX1")

        requester.get.assert_called_once_with(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id=SRX1&retmode=xml"
        )
        response.raise_for_status.assert_called_once()
        self.assertEqual("EXPERIMENT_PACKAGE_SET", root.tag)

    def test_fetch_ena_file_report_uses_ena_requester(self):
        response = Mock()
        response.raise_for_status = Mock()
        response.json.return_value = [{"run_accession": "SRR1"}]
        requester = Mock()
        requester.get.return_value = response

        rows = INSDCWebfetcher(ena_requester=requester).fetch_ena_file_report("SRX1")

        requester.get.assert_called_once_with(
            "https://www.ebi.ac.uk/ena/portal/api/filereport",
            params={
                "accession": "SRX1",
                "result": "read_run",
                "fields": "run_accession,fastq_ftp,fastq_md5,fastq_bytes",
                "format": "json",
            },
        )
        response.raise_for_status.assert_called_once()
        self.assertEqual([{"run_accession": "SRR1"}], rows)

    def test_fetch_sra_runs_parses_experiment_package_xml(self):
        xml = """<?xml version="1.0"?>
<EXPERIMENT_PACKAGE_SET>
  <EXPERIMENT_PACKAGE>
    <STUDY accession="ERP137216" />
    <EXPERIMENT accession="SRX1">
      <PLATFORM><ILLUMINA><INSTRUMENT_MODEL>NovaSeq 6000</INSTRUMENT_MODEL></ILLUMINA></PLATFORM>
      <DESIGN>
        <LIBRARY_DESCRIPTOR>
          <LIBRARY_LAYOUT><PAIRED /></LIBRARY_LAYOUT>
          <LIBRARY_SELECTION>cDNA</LIBRARY_SELECTION>
          <LIBRARY_SOURCE>TRANSCRIPTOMIC</LIBRARY_SOURCE>
          <LIBRARY_STRATEGY>RNA-Seq</LIBRARY_STRATEGY>
        </LIBRARY_DESCRIPTOR>
      </DESIGN>
    </EXPERIMENT>
    <SAMPLE accession="SRS1">
      <IDENTIFIERS>
        <EXTERNAL_ID namespace="BioSample">SAMN1</EXTERNAL_ID>
        <EXTERNAL_ID namespace="GEO">GSM1</EXTERNAL_ID>
      </IDENTIFIERS>
    </SAMPLE>
    <RUN_SET>
      <RUN accession="SRR1" alias="lane 1">
        <Statistics><Read average="151" /></Statistics>
        <SRAFiles>
          <SRAFile semantic_name="fastq" filename="r1.fastq.gz" md5="md5-r1" url="ftp://example/r1.fastq.gz" />
          <SRAFile semantic_name="fastq" filename="r2.fastq.gz" md5="md5-r2">
            <Alternatives url="ftp://example/r2.fastq.gz" />
          </SRAFile>
        </SRAFiles>
      </RUN>
    </RUN_SET>
  </EXPERIMENT_PACKAGE>
</EXPERIMENT_PACKAGE_SET>
"""
        fetcher = INSDCWebfetcher()
        fetcher._ncbi_nrx = Mock(return_value=ET.fromstring(xml))
        fetcher.fetch_ena_file_report = Mock(return_value=[
            {
                "run_accession": "SRR1",
                "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/SRR1/001/SRR1/ena-r1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/SRR1/001/SRR1/ena-r2.fastq.gz",
                "fastq_md5": "ena-md5-r1;ena-md5-r2",
                "fastq_bytes": "123;456",
            }
        ])

        runs = fetcher.fetch_sra_runs("SRX1")

        fetcher._ncbi_nrx.assert_called_once_with(nrx="SRX1")
        fetcher.fetch_ena_file_report.assert_called_once_with(accession="SRX1")
        self.assertEqual(1, len(runs))
        self.assertEqual(
            {
                "sample": "SRS1",
                "biosample": "SAMN1",
                "geo_sample": "GSM1",
                "library_layout": "PAIRED",
                "library_selection": "cDNA",
                "library_source": "TRANSCRIPTOMIC",
                "library_strategy": "RNA-Seq",
                "study": "ERP137216",
                "experiment": "SRX1",
                "run": "SRR1",
                "scan_name": "lane 1",
                "instrument_model": "NovaSeq 6000",
                "fastq_files": [
                    {
                        "filename": "ena-r1.fastq.gz",
                        "uri": "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR1/001/SRR1/ena-r1.fastq.gz",
                        "md5": "ena-md5-r1",
                        "bytes": "123",
                    },
                    {
                        "filename": "ena-r2.fastq.gz",
                        "uri": "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR1/001/SRR1/ena-r2.fastq.gz",
                        "md5": "ena-md5-r2",
                        "bytes": "456",
                    },
                ],
                "submitted_file_name": "ena-r1.fastq.gz",
                "md5": "ena-md5-r1",
                "read_lengths": ["151"],
            },
            runs[0],
        )

    def test_parse_ena_fastq_report_handles_empty_and_scheme_less_values(self):
        fetcher = INSDCWebfetcher()

        self.assertEqual([], fetcher._parse_ena_fastq_report({"fastq_ftp": ""}))
        self.assertEqual(
            [
                {
                    "filename": "SRR1_1.fastq.gz",
                    "uri": "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR1/001/SRR1/SRR1_1.fastq.gz",
                    "md5": "md5-r1",
                    "bytes": "123",
                },
                {
                    "filename": "SRR1_2.fastq.gz",
                    "uri": "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR1/001/SRR1/SRR1_2.fastq.gz",
                    "md5": "md5-r2",
                    "bytes": None,
                },
            ],
            fetcher._parse_ena_fastq_report({
                "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/SRR1/001/SRR1/SRR1_1.fastq.gz;ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR1/001/SRR1/SRR1_2.fastq.gz",
                "fastq_md5": "md5-r1;md5-r2",
                "fastq_bytes": "123",
            }),
        )

    def test_fetch_sra_runs_falls_back_to_ncbi_xml_urls_when_ena_has_no_fastqs(self):
        xml = """<?xml version="1.0"?>
<EXPERIMENT_PACKAGE_SET>
  <EXPERIMENT_PACKAGE>
    <EXPERIMENT accession="SRX1" />
    <SAMPLE accession="SRS1" />
    <RUN_SET>
      <RUN accession="SRR1784623">
        <SRAFiles>
          <SRAFile semantic_name="fastq" filename="SRR1784623_1.fastq.gz" md5="md5-r1" url="ftp://ncbi.example/SRR1784623_1.fastq.gz" />
          <SRAFile semantic_name="fastq" filename="SRR1784623_2.fastq.gz" md5="md5-r2">
            <Alternatives url="ftp://ncbi.example/SRR1784623_2.fastq.gz" />
          </SRAFile>
        </SRAFiles>
      </RUN>
    </RUN_SET>
  </EXPERIMENT_PACKAGE>
</EXPERIMENT_PACKAGE_SET>
"""
        fetcher = INSDCWebfetcher()
        fetcher._ncbi_nrx = Mock(return_value=ET.fromstring(xml))
        fetcher.fetch_ena_file_report = Mock(return_value=[])

        runs = fetcher.fetch_sra_runs("SRX1")

        self.assertEqual("ftp://ncbi.example/SRR1784623_1.fastq.gz", runs[0]["fastq_files"][0]["uri"])
        self.assertEqual("ftp://ncbi.example/SRR1784623_2.fastq.gz", runs[0]["fastq_files"][1]["uri"])
        self.assertIsNone(runs[0]["study"])

    def test_fetch_sra_runs_falls_back_to_ncbi_xml_urls_when_ena_request_fails(self):
        xml = """<?xml version="1.0"?>
<EXPERIMENT_PACKAGE_SET>
  <EXPERIMENT_PACKAGE>
    <EXPERIMENT accession="SRX1" />
    <SAMPLE accession="SRS1" />
    <RUN_SET>
      <RUN accession="SRR1">
        <SRAFiles>
          <SRAFile semantic_name="fastq" filename="r1.fastq.gz" md5="md5-r1" url="ftp://ncbi.example/r1.fastq.gz" />
        </SRAFiles>
      </RUN>
    </RUN_SET>
  </EXPERIMENT_PACKAGE>
</EXPERIMENT_PACKAGE_SET>
"""
        fetcher = INSDCWebfetcher()
        fetcher._ncbi_nrx = Mock(return_value=ET.fromstring(xml))
        fetcher.fetch_ena_file_report = Mock(side_effect=requests.RequestException("ena unavailable"))

        runs = fetcher.fetch_sra_runs("SRX1")

        self.assertEqual("ftp://ncbi.example/r1.fastq.gz", runs[0]["fastq_files"][0]["uri"])

    def test_sdrf_constructor_delegates_sra_lookup_to_fetcher(self):
        fetcher = Mock()
        fetcher.fetch_sra_runs.return_value = [{"run": "SRR1"}]

        result = SDRFConstructor(insdc_fetcher=fetcher)._lookup_sra("SRX1")

        fetcher.fetch_sra_runs.assert_called_once_with(accession="SRX1")
        self.assertEqual([{"run": "SRR1"}], result)

    def test_sdrf_constructor_returns_empty_runs_on_fetch_error(self):
        fetcher = Mock()
        fetcher.fetch_sra_runs.side_effect = requests.RequestException("network unavailable")

        self.assertEqual([], SDRFConstructor(insdc_fetcher=fetcher)._lookup_sra("SRX1"))


if __name__ == "__main__":
    unittest.main()
