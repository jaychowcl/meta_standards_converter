# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""The bare MINiML payload and its digest-bound scientific disposition travel together."""
import json
from pathlib import Path
from copy import deepcopy
import pytest
from meta_standards_converter.converters import Converter, PreparationPolicy
from meta_standards_converter.sources.json import JSONPackageSource

OPTS = dict(preparation=PreparationPolicy(enrichment="off", expand_studies=False))


def atlas():
    path = Path(__file__).parents[1] / "fixtures/contracts/atlas-document-v1.json"
    data = json.loads(path.read_text())
    data["schema_version"] = "1.2"
    clean, reviewed = deepcopy(data["datasets"][0]), deepcopy(data["datasets"][0])
    for i, dataset in enumerate([clean, reviewed], 1):
        sid = f"GSE{i}"
        dataset.update(dataset_id=sid, source_ordinal=i, publication_ids=[], diagnostics=[], review=None)
        dataset["metadata"] = {"miniml_schema_version": "3.1", "source": {"format": "geo"},
            "series": {"iid": sid}, "sample": [{"iid": f"GSM{i}"}]}
        dataset["harmonization"].update(targets=[], degraded_stages=[], status={
            "contract_version": "2.0", "execution": "succeeded" if i == 1 else "degraded",
            "completeness": "complete" if i == 1 else "partial", "evidence_confidence": "high",
            "validation": "valid", "publication": "publishable" if i == 1 else "review_required",
            "terminal_reason": "harmonization_result", "errors": []})
    data.update(datasets=[clean, reviewed])
    data["summary"].update(dataset_count=2, completed_dataset_count=2, failed_dataset_count=0)
    return data


def test_aggregate_preserves_per_dataset_status_and_raw_list(tmp_path):
    path = tmp_path / "aggregate.json"
    item = Converter().convert(atlas(), out_type="json", outfile=path, output_options={"aggregate": True},
                               **OPTS).items[0]
    assert isinstance(json.loads(path.read_text()), list)
    assert item.publication == "review_required"
    assert path.with_name(path.name + ".msc-receipt").exists()
    item = Converter().convert(path, out_type="tsv", output_options={"aggregate": True}, **OPTS).items[0]
    table, = item.payload.values()
    assert [(r["msc.sample.accession"], r["msc.source.publication"]) for r in table.rows] == [
        ("GSM1", "publishable"), ("GSM2", "review_required")]
    assert item.publication == "review_required"
    assert item.completeness == "partial"


def test_manual_save_and_receipt_mismatch(tmp_path):
    projected = Converter().convert(atlas(), out_type="json", **OPTS).items[0]
    dest = tmp_path / "transfer.json"
    projected.save_json(dest)
    with pytest.raises(FileExistsError):
        projected.save_json(dest)
    contents = json.loads(dest.read_text())
    assert isinstance(contents, list)
    dest.write_text(json.dumps(contents[::-1]))
    with pytest.raises(ValueError, match="receipt"):
        JSONPackageSource().load(dest)
    dest.with_name(dest.name + ".msc-receipt").unlink()
    reopened = Converter().convert(dest, out_type="tsv", **OPTS).items[0]
    for table in reopened.payload.values():
        assert all(row["msc.source.publication"] == "unknown" for row in table.rows)


@pytest.mark.fake_process
def test_overwrite_interrupted_receipt_cannot_clear_review(tmp_path, monkeypatch):
    from meta_standards_converter.sources import receipts
    from meta_standards_converter.converters.json2tsv import JSON2DelimitedConverter
    path = tmp_path / "transfer.json"
    old = Converter().convert(atlas(), out_type="json", **OPTS).items[0]
    old.save_json(path)
    replace = receipts.os.replace

    def interrupted(src, dest):
        if Path(dest) == receipts.receipt_path(path):
            raise OSError("interrupted after payload publication")
        return replace(src, dest)

    with monkeypatch.context() as patch:
        patch.setattr(receipts.os, "replace", interrupted)
        with pytest.raises(OSError, match="interrupted"):
            old.save_json(path, overwrite=True)
    assert not receipts.receipt_path(path).exists()
    assert JSON2DelimitedConverter().convert_source(path, tmp_path / "unknown.tsv").dataset_outcomes[0][1]["publication"] == "unknown"
    old.save_json(path, overwrite=True)
    import os, subprocess, sys
    code = """import json, sys, socket
def offline(*args, **kwargs): raise AssertionError("Receipt readback must remain offline")
socket.socket.connect = offline
from meta_standards_converter.converters.json2tsv import JSON2DelimitedConverter
r = JSON2DelimitedConverter().convert_source(sys.argv[1], sys.argv[2])
print(json.dumps(r.to_dict()))
"""
    record = subprocess.run([sys.executable, "-c", code, str(path), str(tmp_path / "fresh.tsv")],
                            capture_output=True, text=True, check=True, timeout=45,
                            env={**os.environ, "PATH": str(tmp_path)})
    restored = json.loads(record.stdout.splitlines()[-1])
    assert restored["dataset_outcomes"]["GSE2"]["publication"] == "review_required"
    assert restored["status"] == "partial"


def test_invalid_material_collision_survives_a_saved_export(tmp_path):
    from meta_standards_converter.metadata.workflows import prepare_workflows
    from tests.magetab.test_metadata_integrity import geo
    data = prepare_workflows(geo("GSE13398")).to_mapping()
    paths = data["series"]["assay_paths"]
    paths[1]["steps"][0]["name"] = paths[0]["steps"][0]["name"]
    item = Converter().convert(data, out_type="json", **OPTS).items[0]
    assert item.validation == "invalid"
    assert item.publication == "blocked"
    dest = tmp_path / "invalid.json"
    item.save_json(dest)
    restored = Converter().convert(dest, out_type="tsv", **OPTS).items[0]
    assert restored.validation == "invalid"
    assert restored.publication == "blocked"


def test_resaving_a_bare_file_cannot_establish_its_prior_publication(tmp_path):
    original = tmp_path / "bare.json"
    original.write_text(json.dumps([atlas()["datasets"][0]["metadata"]]))
    output = tmp_path / "new.json"
    conversion = Converter().convert(original, out_type="json", **OPTS).items[0]
    conversion.save_json(output)
    item = Converter().convert(output, out_type="tsv", **OPTS).items[0]
    assert next(iter(item.payload.values())).rows[0]["msc.source.publication"] == "unknown"
