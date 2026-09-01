# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from __future__ import annotations

from unittest.mock import Mock, patch

from meta_standards_converter.cli import ena2json as ena_cli
from meta_standards_converter.cli import sra2json as sra_cli


def test_cli_exposes_no_public_record_cap() -> None:
    for module in (sra_cli, ena_cli):
        parser = module._parser()
        destinations = {action.dest for action in parser._actions}
        assert "max_records" not in destinations
        assert {"enrich_geo", "enrich_ae", "out", "resource_profile"} <= destinations


def test_cli_continues_after_one_input_fails() -> None:
    converter = Mock()
    converter.convert.side_effect = [ValueError("bad"), []]
    with patch.object(sra_cli, "sra2json", return_value=converter):
        status = sra_cli.main(["SRP000001", "SRP000002", "--out", "."])

    assert status == 1
    assert [call.args[0] for call in converter.convert.call_args_list] == [
        "SRP000001",
        "SRP000002",
    ]
