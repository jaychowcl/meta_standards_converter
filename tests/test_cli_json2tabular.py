# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
import json
from unittest.mock import Mock, patch

from meta_standards_converter.cli import json2tsv


def test_json2tsv_cli_converts_inputs_in_order():
    summary = {"operation": "manifest", "status": "complete", "datasets": []}
    result = Mock(partial=False)
    result.to_dict.return_value = summary
    with patch.object(json2tsv, "JSONDataOutputOrchestrator") as factory:
        factory.return_value.export_manifest.return_value = result

        with patch("builtins.print") as emit:
            status = json2tsv.main(
                ["one.json", "two.json", "--outdir", "tables", "--format", "csv"]
            )

    assert status == 0
    assert [call.args[0] for call in factory.return_value.export_manifest.call_args_list] == [
        "one.json",
        "two.json",
    ]
    assert all(
        call.kwargs["output_format"] == "csv"
        for call in factory.return_value.export_manifest.call_args_list
    )
    assert json.loads(emit.call_args.args[0])["status"] == "complete"


def test_json2tsv_failure_is_sanitized_and_reported(capsys):
    canary = "private-manifest-detail"
    with patch.object(json2tsv, "JSONDataOutputOrchestrator") as factory:
        factory.return_value.export_manifest.side_effect = RuntimeError(
            f"manifest failed: {canary}"
        )
        status = json2tsv.main(["private/input.json"])

    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert status == 1
    assert canary not in captured.err
    assert canary not in captured.out
    assert payload["datasets"][0]["source"] == "input.json"
    assert payload["datasets"][0]["error"]["error_type"] == "RuntimeError"
