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
from unittest.mock import patch


ROOT = os.path.dirname(os.path.dirname(__file__))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from meta_standards_converter.geo_handlers.geo_parser import GEOParser  # noqa: E402


def miniml_body(body: str) -> str:
    return f"""<?xml version="1.0" encoding="UTF-8"?>
<MINiML xmlns="http://www.ncbi.nlm.nih.gov/geo/info/MINiML" version="0.5.4">
{body}
</MINiML>
"""


class TestGEOParser(unittest.TestCase):
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
        self.assertEqual("0.5.4", package["version"])
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
        self.assertEqual(["Expression profiling by high throughput sequencing"], package["series"]["type"])
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
        self.assertEqual("1", channel["position"])
        self.assertEqual({"tag": "tissue", "value": "CSF"}, channel["characteristics"][0])

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

    @patch("meta_standards_converter.geo_handlers.geo_parser.GEOWebFetcher")
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

        parsed = GEOParser().parse(root_xml, related_series=True)

        self.assertEqual(["GSE1", "GSE2"], [package["series"]["iid"] for package in parsed])
        fetcher_mock.return_value.fetch_gse_miniml.assert_called_once_with(gse="GSE2")


if __name__ == "__main__":
    unittest.main()
