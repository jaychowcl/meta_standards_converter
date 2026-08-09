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

from meta_standards_converter.cli import json2obs


def test_json2obs_requires_outdir():
    with __import__("pytest").raises(SystemExit) as raised:
        json2obs.main(["input.json"])
    assert raised.value.code == 2


def test_json2obs_forwards_asset_and_component_options(capsys):
    result = Mock(partial=False)
    result.to_dict.return_value = {
        "operation": "anndata_metadata",
        "status": "complete",
    }
    with patch.object(json2obs, "JSONDataOutputOrchestrator") as factory:
        factory.return_value.export_anndata_metadata.return_value = result
        status = json2obs.main(
            [
                "input.json",
                "--outdir",
                "metadata",
                "--asset",
                "GSM1=source.h5ad",
                "--include-var",
                "--include-uns",
                "--allow-invalid",
                "--overwrite",
                "--resume",
                "--processed-checkpoint-dir",
                "checkpoints",
            ]
        )

    assert status == 0
    factory.return_value.export_anndata_metadata.assert_called_once_with(
        "input.json",
        outdir="metadata",
        include_var=True,
        include_uns=True,
        asset_manifest=None,
        asset_specs=["GSM1=source.h5ad"],
        force_reprocess=False,
        pipeline="auto",
        genome=None,
        fasta=None,
        gtf=None,
        gff=None,
        accept_inferred_reference=False,
        profile="docker",
        revision=None,
        params_file=None,
        nextflow_config=None,
        work_dir=None,
        resume=True,
        processed_checkpoint_dir="checkpoints",
        overwrite=True,
        allow_invalid=True,
        matrix_orientation="auto",
    )
    assert json.loads(capsys.readouterr().out)["status"] == "complete"


def test_json2obs_failure_is_sanitized_and_reported(capsys):
    canary = "private-observation-detail"
    with patch.object(json2obs, "JSONDataOutputOrchestrator") as factory:
        factory.return_value.export_anndata_metadata.side_effect = RuntimeError(
            f"export failed: {canary}"
        )
        status = json2obs.main(["private/input.json", "--outdir", "metadata"])

    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert status == 1
    assert canary not in captured.err
    assert canary not in captured.out
    assert payload["datasets"][0]["source"] == "input.json"
    assert payload["datasets"][0]["error"]["error_type"] == "RuntimeError"
