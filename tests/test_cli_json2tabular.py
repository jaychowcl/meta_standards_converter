# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from unittest.mock import patch

from meta_standards_converter.cli import json2csv, json2tsv


def test_json2tsv_cli_converts_inputs_in_order():
    with patch.object(json2tsv, "JSON2TSVConverter") as factory:
        factory.return_value.convert_source.return_value.partial = False

        status = json2tsv.main(["one.json", "two.json", "--out", "tables"])

    assert status == 0
    assert [call.args[0].name for call in factory.return_value.convert_source.call_args_list] == [
        "one.json",
        "two.json",
    ]


def test_json2csv_cli_returns_one_after_partial_result():
    with patch.object(json2csv, "JSON2CSVConverter") as factory:
        factory.return_value.convert_source.return_value.partial = True

        status = json2csv.main(["one.json"])

    assert status == 1
