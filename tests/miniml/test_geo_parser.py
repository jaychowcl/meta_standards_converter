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
from typing import get_type_hints
import unittest
from unittest.mock import patch


ROOT = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from meta_standards_converter.miniml.geo_parser import (  # noqa: E402
    GEOParser,
)
from meta_standards_converter.sources.geo import GEOSource, RelatedSeriesParseResult
from meta_standards_converter.miniml import MINiMLPackage  # noqa: E402


def miniml_body(body: str) -> str:
    return f"""<?xml version="1.0" encoding="UTF-8"?>
<MINiML xmlns="http://www.ncbi.nlm.nih.gov/geo/info/MINiML" version="0.5.4">
{body}
</MINiML>
"""


class TestGEOParser(unittest.TestCase):
    def test_xml_type_hints_resolve_at_runtime(self):
        hints = get_type_hints(GEOParser._top_level_nodes)

        self.assertEqual("Element", hints["root"].__name__)

    def test_parser_emits_canonical_versioned_package(self):
        package = GEOParser().parse(
            miniml_body(
                '<Sample iid="GSM1" /><Series iid="GSE1"><Sample-Ref ref="GSM1" /></Series>'
            )
        )[0]

        self.assertIsInstance(package, MINiMLPackage)
        self.assertEqual("3.0", package.miniml_schema_version)
        self.assertEqual("GEO MINiML", package.source.format)
        self.assertEqual("0.5.4", package.source.version)
        self.assertEqual(package, MINiMLPackage.from_mapping(package.to_mapping()))

    def test_parser_preserves_ambiguous_library_kit_protocol_as_authored(self):
        protocol = "Chromium Single Cell 3' Library &amp; Gel Bead Kit v2 or v3."
        package = GEOParser().parse(
            miniml_body(
                '<Sample iid="GSM1"><Channel-Count>1</Channel-Count><Channel>'
                f"<Extract-Protocol>{protocol}</Extract-Protocol>"
                '</Channel></Sample><Series iid="GSE1"><Sample-Ref ref="GSM1" />'
                "</Series>"
            )
        )[0]

        self.assertEqual(
            "Chromium Single Cell 3' Library & Gel Bead Kit v2 or v3.",
            package["sample"][0]["channel"][0]["extract_protocol"],
        )

    def test_single_series_resolves_relevant_records(self):
        xml = miniml_body(
            """
  <Organization iid="org1"><Name>Org One</Name><Address><City>London</City><Postal-Code>NW1</Postal-Code><Country>UK</Country></Address></Organization>
  <Contributor iid="contrib1"><Person><First>Ada</First><Last>Lovelace</Last></Person><Organization-Ref ref="org1" /></Contributor>
  <Contributor iid="contrib2"><Person><First>Grace</First><Last>Hopper</Last></Person></Contributor>
  <Database iid="GEO"><Name>GEO</Name><Public-ID>GEO</Public-ID><Organization-Ref ref="org1" /></Database>
  <Platform iid="GPL1"><Accession database="GEO">GPL1</Accession><Contact-Ref ref="contrib2" /></Platform>
  <Sample iid="GSM1"><Accession database="GEO">GSM1</Accession><Platform-Ref ref="GPL1" /><Contact-Ref ref="contrib1" /></Sample>
  <Sample iid="GSM2"><Accession database="GEO">GSM2</Accession></Sample>
  <Series iid="GSE1"><Title>Study</Title><Accession database="GEO">GSE1</Accession><Summary>Summary</Summary><Sample-Ref ref="GSM1" /></Series>
"""
        )

        parsed = GEOParser().parse(xml)

        self.assertEqual(1, len(parsed))
        package = parsed[0]
        self.assertEqual("0.5.4", package["source"]["version"])
        self.assertEqual("GSE1", package["series"]["iid"])
        self.assertEqual(["GSM1"], [sample["iid"] for sample in package["sample"]])
        self.assertEqual(["GPL1"], [platform["iid"] for platform in package["platform"]])
        self.assertEqual(
            ["contrib1", "contrib2"],
            sorted(contributor["iid"] for contributor in package["contributor"]),
        )
        self.assertEqual(["GEO"], [database["iid"] for database in package["database"]])
        self.assertEqual(["org1"], [organization["iid"] for organization in package["organization"]])

    def test_multiple_series_return_multiple_packages(self):
        xml = miniml_body(
            """
  <Sample iid="GSM1" />
  <Sample iid="GSM2" />
  <Series iid="GSE1"><Title>One</Title><Summary>One summary</Summary><Sample-Ref ref="GSM1" /></Series>
  <Series iid="GSE2"><Title>Two</Title><Summary>Two summary</Summary><Sample-Ref ref="GSM2" /></Series>
"""
        )

        parsed = GEOParser().parse(xml)

        self.assertEqual(["GSE1", "GSE2"], [package["series"]["iid"] for package in parsed])
        self.assertEqual(["GSM1"], [sample["iid"] for sample in parsed[0]["sample"]])
        self.assertEqual(["GSM2"], [sample["iid"] for sample in parsed[1]["sample"]])

    def test_repeated_and_singleton_fields_follow_xsd_cardinality(self):
        xml = miniml_body(
            """
  <Contributor iid="contrib1"><Person><First>Ada</First><Last>Lovelace</Last></Person><Address><City>London</City><Postal-Code>NW1</Postal-Code><Country>UK</Country></Address></Contributor>
  <Sample iid="GSM1"><Type>SRA</Type><Platform-Ref ref="GPL1" /><Extra>one</Extra><Extra>two</Extra></Sample>
  <Series iid="GSE1">
    <Title>Study</Title>
    <Accession database="GEO">GSE1</Accession>
    <Summary>Summary</Summary>
    <Type>Expression profiling by high throughput sequencing</Type>
    <Contact-Ref ref="contrib1" />
    <Sample-Ref ref="GSM1" />
    <Relation type="SuperSeries of" target="GSE2" />
  </Series>
"""
        )

        package = GEOParser().parse(xml)[0]

        self.assertIsInstance(package["series"]["accession"], list)
        self.assertEqual("GSE1", package["series"]["accession"][0]["value"])
        self.assertIsInstance(package["series"]["sample_ref"], list)
        self.assertIsInstance(package["series"]["relation"], list)
        self.assertEqual(
            [{"value": "Expression profiling by high throughput sequencing"}],
            package["series"]["type"],
        )
        self.assertEqual("SRA", package["sample"][0]["type"])
        self.assertEqual({"ref": "GPL1"}, package["sample"][0]["platform_ref"])
        self.assertEqual(["one", "two"], package["sample"][0]["extra"])
        self.assertNotIn("extras", package["sample"][0])
        self.assertIsInstance(package["contributor"][0]["address"], dict)
        self.assertEqual("London", package["contributor"][0]["address"]["city"])

    def test_attributes_and_empty_cleanup(self):
        xml = miniml_body(
            """
  <Sample iid="GSM1"><Title></Title><Channel position="1"><Source>CSF</Source><Characteristics tag="tissue">CSF</Characteristics></Channel></Sample>
  <Series iid="GSE1"><Title>Study</Title><Summary>Summary</Summary><Sample-Ref ref="GSM1" /></Series>
"""
        )

        package = GEOParser().parse(xml, remove_empty=True)[0]

        self.assertNotIn("title", package["sample"][0])
        channel = package["sample"][0]["channel"][0]
        self.assertNotIn("position", channel)
        self.assertEqual({"name": "tissue", "value": "CSF"}, channel["characteristics"][0])

    def test_real_fixture_parses_one_series_package(self):
        fixture = os.path.join(ROOT, "tests", "GSE328265_family.xml")
        if not os.path.exists(fixture):
            self.skipTest("tests/GSE328265_family.xml is not available")
        with open(fixture, encoding="utf-8") as handle:
            xml = handle.read()

        package = GEOParser().parse(xml)[0]

        self.assertEqual("GSE328265", package["series"]["accession"][0]["value"])
        self.assertEqual(31, len(package["sample"]))
        self.assertEqual(["GPL30173"], [platform["iid"] for platform in package["platform"]])
        self.assertEqual(["GEO"], [database["iid"] for database in package["database"]])

    def test_parse_logs_structural_counts_without_xml_payload(self):
        xml = miniml_body(
            '<Platform iid="GPL1"><Title>secret-title</Title></Platform>'
            '<Sample iid="GSM1"><Platform-Ref ref="GPL1" /></Sample>'
            '<Series iid="GSE1"><Sample-Ref ref="GSM1" /></Series>'
        )

        with self.assertLogs(
            "meta_standards_converter.miniml.geo_parser", level="INFO"
        ) as logs:
            GEOParser().parse(xml)

        output = "\n".join(logs.output)
        self.assertIn("MINiML parse stats packages=1 series=1 samples=1 platforms=1", output)
        self.assertNotIn("secret-title", output)

    @patch("meta_standards_converter.sources.geo.GEOWebFetcher")
    def test_related_series_are_fetched_recursively_and_deduplicated(self, fetcher_mock):
        root_xml = miniml_body(
            """
  <Series iid="GSE1"><Title>Root</Title><Accession database="GEO">GSE1</Accession><Summary>Root summary</Summary><Relation type="SuperSeries of" target="GSE2" /></Series>
"""
        )
        related_xml = miniml_body(
            """
  <Series iid="GSE2"><Title>Related</Title><Accession database="GEO">GSE2</Accession><Summary>Related summary</Summary><Relation type="SubSeries of" target="GSE1" /></Series>
"""
        )
        fetcher_mock.return_value.fetch_gse_miniml.return_value = related_xml

        parsed = GEOSource().parse(root_xml, related_series=True)

        self.assertEqual(["GSE1", "GSE2"], [package["series"]["iid"] for package in parsed])
        fetcher_mock.return_value.fetch_gse_miniml.assert_called_once_with(gse="GSE2")

    def test_permissive_related_series_returns_typed_partial_summary(self):
        root_xml = miniml_body(
            """
  <Series iid="GSE1"><Accession database="GEO">GSE1</Accession><Relation type="SuperSeries of" target="GSE2" /><Relation type="SuperSeries of" target="GSE3" /></Series>
"""
        )
        related_xml = miniml_body(
            """
  <Series iid="GSE3"><Accession database="GEO">GSE3</Accession></Series>
"""
        )

        class PartialFetcher:
            def fetch_gse_miniml(self, gse):
                if gse == "GSE2":
                    raise RuntimeError("private-provider-detail")
                return related_xml

        with self.assertLogs(
            "meta_standards_converter.sources.geo", level="WARNING"
        ) as logs:
            result = GEOSource(fetcher=PartialFetcher()).parse_related_series(
                root_xml,
                strict=False,
            )

        self.assertIsInstance(result, RelatedSeriesParseResult)
        self.assertEqual(["GSE3"], [item["series"]["iid"] for item in result])
        self.assertEqual("degraded", result.status.execution.value)
        self.assertEqual("partial", result.status.completeness.value)
        self.assertEqual("review_required", result.status.publication.value)
        self.assertEqual(("GSE2", "GSE3"), result.attempted_accessions)
        self.assertEqual(("GSE2",), result.failed_accessions)
        self.assertEqual("RuntimeError", result.status.errors[0].error_type)
        self.assertEqual("GSE2", result.status.errors[0].item_id)
        self.assertNotIn("private-provider-detail", "\n".join(logs.output))


if __name__ == "__main__":
    unittest.main()


def test_reference_fields_follow_document_order_regardless_of_key_set():
    parser = GEOParser()
    # Both orders must work in the same interpreter, independently of its hash seed.
    for first, second in [('contact_ref', 'contributor_ref'), ('contributor_ref', 'contact_ref')]:
        element = {first: [{'ref': 'first'}, {'ref': 'repeat'}],
                   second: [{'ref': 'second'}, {'ref': 'repeat'}],
                   'channel': [{'contact_ref': {'ref': 'nested'}}]}
        assert parser._reference_values(element, {'contact_ref', 'contributor_ref'}) == [
            'first', 'repeat', 'second', 'repeat', 'nested']
        contributors = parser._resolve_contributors(element,
            [{'contact_ref': {'ref': 'sample'}}], [{'contributor_ref': {'ref': 'platform'}}],
            {iid: {'iid': iid} for iid in ['first', 'second', 'repeat', 'nested', 'sample', 'platform']})
        assert [c['iid'] for c in contributors] == ['first', 'repeat', 'second', 'nested', 'sample', 'platform']


def test_gse60450_contributors_follow_real_xml_reference_order():
    from pathlib import Path
    from xml.etree import ElementTree as ET
    raw = (Path(__file__).parents[1] / 'fixtures/studies/GSE60450/inputs/geo.xml').read_text()
    root = ET.fromstring(raw)
    ns = {'g': 'http://www.ncbi.nlm.nih.gov/geo/info/MINiML'}
    series = root.find('g:Series', ns)
    refs = [node.attrib['ref'] for node in series if node.tag.rsplit('}', 1)[-1] in {'Contributor-Ref', 'Contact-Ref'}]
    package = GEOParser().parse(raw)[0]
    assert [c['iid'] for c in package['contributor']][:len(dict.fromkeys(refs))] == list(dict.fromkeys(refs))
    assert len(package['sample']) == 12
