# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from __future__ import annotations

import json


import hashlib


import os


import tempfile


from pathlib import Path

from typing import Any, Mapping, Sequence


from meta_standards_converter.expression.assets import Asset


from meta_standards_converter.expression.readers import scientific_modules

class ProcessedCheckpointStore:
    def __init__(self, package_version):
        self._package_version = package_version

    def key(
        self,
        root: Path | None,
        *,
        sample_id: str,
        source_json_sha256: str,
        sample: Mapping[str, Any],
        asset: Asset,
        orientation: str,
    ) -> tuple[Path, Path, str] | None:
        if root is None:
            return None
        payload = {
            "schema_version": "1",
            "converter_version": self._package_version(),
            "source_json_sha256": source_json_sha256,
            "sample_id": sample_id,
            "sample": sample,
            "asset": vars(asset),
            "orientation": orientation,
        }
        fingerprint = hashlib.sha256(
            json.dumps(payload, ensure_ascii=False, sort_keys=True).encode("utf-8")
        ).hexdigest()
        # Version-specific destinations preserve historical processed checkpoints.
        # The fingerprint payload and algorithm remain unchanged.
        version_key = hashlib.sha256(self._package_version().encode("utf-8")).hexdigest()[:20]
        root = root / f"msc-{version_key}"
        key = hashlib.sha256(sample_id.encode("utf-8")).hexdigest()[:20]
        return root / f"{key}.h5ad", root / f"{key}.json", fingerprint

    def _load_processed_checkpoint(self, checkpoint):
        metadata = self.metadata(checkpoint)
        if metadata is None:
            return None
        h5ad_path = checkpoint[0]
        try:
            anndata = scientific_modules()[0]
            return anndata.read_h5ad(h5ad_path), metadata
        except (OSError, TypeError, ValueError):
            return None

    def metadata(self, checkpoint):
        if checkpoint is None:
            return None
        h5ad_path, manifest_path, fingerprint = checkpoint
        if not h5ad_path.is_file() or not manifest_path.is_file():
            return None
        try:
            metadata = json.loads(manifest_path.read_text(encoding="utf-8"))
            if metadata.get("fingerprint") != fingerprint:
                return None
            return metadata
        except (OSError, TypeError, ValueError, json.JSONDecodeError):
            return None

    def observation_ids(self, path: Path) -> set[str] | None:
        """Read checkpoint observation metadata without materialising its matrix."""

        backed = None
        try:
            anndata = scientific_modules()[0]
            backed = anndata.read_h5ad(path, backed="r")
            return set(backed.obs_names.astype(str))
        except (OSError, TypeError, ValueError):
            return None
        finally:
            if backed is not None:
                backed.file.close()

    def write(
        self,
        checkpoint,
        adata,
        *,
        warnings: Sequence[str],
        errors: Sequence[str],
    ) -> None:
        if checkpoint is None:
            return
        h5ad_path, manifest_path, fingerprint = checkpoint
        h5ad_path.parent.mkdir(parents=True, exist_ok=True)
        self._write_h5ad(adata, h5ad_path, overwrite=True)
        with tempfile.NamedTemporaryFile(
            dir=manifest_path.parent,
            suffix=".json",
            mode="w",
            encoding="utf-8",
            delete=False,
        ) as handle:
            temporary = Path(handle.name)
            json.dump(
                {
                    "schema_version": "1",
                    "fingerprint": fingerprint,
                    "warnings": list(warnings),
                    "errors": list(errors),
                },
                handle,
                ensure_ascii=False,
                sort_keys=True,
            )
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        try:
            os.replace(temporary, manifest_path)
        finally:
            temporary.unlink(missing_ok=True)

    def _write_h5ad(self, adata, path: Path, overwrite: bool) -> None:
        if path.exists() and not overwrite:
            raise FileExistsError(f"Output already exists: {path}")
        path.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(dir=path.parent, suffix=".h5ad", delete=False) as handle:
            temporary = Path(handle.name)
        try:
            adata.write_h5ad(temporary, compression="gzip")
            os.replace(temporary, path)
            path.chmod(0o660)
        finally:
            if temporary.exists():
                temporary.unlink()
