# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
import hashlib
from pathlib import Path
from meta_standards_converter.converters import JSON2H5ADConverter
from meta_standards_converter.expression.readers import ProcessedAssetReader
from tests.support.contracts import PBMC, prepare


def test_version_change_recomputes_and_preserves_old_checkpoint_bytes(workspace, monkeypatch):
    prepare(PBMC, workspace)
    calls = []
    class Reader(ProcessedAssetReader):
        def read(self, asset, *, orientation, localize):
            calls.append(asset.scope_id)
            return super().read(asset, orientation=orientation, localize=localize)
    options = dict(asset_specs=["PBMC3K_SAMPLE=counts.tsv"], matrix_orientation="genes-by-observations", resume=True, processed_checkpoint_dir="checkpoints")
    monkeypatch.setattr(JSON2H5ADConverter, "_package_version", lambda self: "7.0.0")
    JSON2H5ADConverter(reader=Reader()).convert("miniml.json", out="old", **options)
    old = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in Path("checkpoints").rglob("*") if p.is_file()}
    assert len(old) == 2
    monkeypatch.setattr(JSON2H5ADConverter, "_package_version", lambda self: "8.0.0")
    JSON2H5ADConverter(reader=Reader()).convert("miniml.json", out="new", **options)
    assert calls == ["PBMC3K_SAMPLE", "PBMC3K_SAMPLE"]
    assert all(hashlib.sha256(p.read_bytes()).hexdigest() == digest for p, digest in old.items())
    assert len(list(Path("checkpoints").rglob("*.h5ad"))) == 2
    JSON2H5ADConverter(reader=Reader()).convert("miniml.json", out="replay", **options)
    assert calls == ["PBMC3K_SAMPLE", "PBMC3K_SAMPLE"]
    assert (workspace / "replay/PBMC3K_SAMPLE.h5ad").is_file()


def test_profile_change_recomputes_even_when_sample_values_are_unchanged(workspace):
    prepare(PBMC, workspace)
    calls = []
    class Reader(ProcessedAssetReader):
        def read(self, asset, *, orientation, localize):
            calls.append(asset.scope_id)
            return super().read(asset, orientation=orientation, localize=localize)
    options = dict(asset_specs=["PBMC3K_SAMPLE=counts.tsv"], matrix_orientation="genes-by-observations", resume=True, processed_checkpoint_dir="checkpoints")
    profile = {"schema_version": "1.0", "replacements": {"disease": ["absent"]}}
    converter = JSON2H5ADConverter(reader=Reader())
    converter.convert("miniml.json", out="raw", **options)
    converter.convert("miniml.json", out="profile", replacement_profile=profile, **options)
    assert calls == ["PBMC3K_SAMPLE", "PBMC3K_SAMPLE"]
    converter.convert("miniml.json", out="same-profile", replacement_profile=dict(reversed(list(profile.items()))), **options)
    assert calls == ["PBMC3K_SAMPLE", "PBMC3K_SAMPLE"]
