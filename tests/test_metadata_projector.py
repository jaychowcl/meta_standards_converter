# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


ROOT = os.path.dirname(os.path.dirname(__file__))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from meta_standards_converter.metadata.projection.anndata import AnnDataMetadataProjection, AnnDataProjectionError
from meta_standards_converter.converters.json2h5ad import JSON2H5ADConverter


def test_asset_downloader_is_part_of_the_converter_public_api():
    from meta_standards_converter.retrieval import AssetDownloader

    assert AssetDownloader.__name__ == "AssetDownloader"


def _package(h5ad_path: str) -> dict:
    return {
        "miniml_schema_version": "2.0",
        "source": {"format": "test"},
        "series": {"accession": [{"value": "GSE1"}]},
        "sample": [
            {
                "iid": "Sample1",
                "accession": [{"value": "GSM1"}],
                "supplementary_data": [{"value": h5ad_path}],
            }
        ],
    }


class RecordingProjector:
    def __init__(self):
        self.sample_contexts = []
        self.combined_contexts = []

    def project_sample(self, *, adata, context):
        self.sample_contexts.append(context)
        return AnnDataMetadataProjection(
            obs={"private.sample": context.sample["iid"]},
            var={"private.genome_build": "GRCh38"},
            uns={"private": {"schema_version": "1.0"}},
            warnings=("private projection warning",),
        )

    def project_combined(self, *, adata, contexts):
        self.combined_contexts.append(contexts)
        return AnnDataMetadataProjection(
            uns={"private_combined": {"sample_count": len(contexts)}}
        )


class TestMetadataProjectorHook(unittest.TestCase):
    def setUp(self):
        import anndata
        import pandas
        from scipy import sparse

        self.anndata = anndata
        self.pandas = pandas
        self.sparse = sparse

    def _fixture(self, tmpdir: str, observations: int = 2) -> tuple[str, str]:
        h5ad_path = os.path.join(tmpdir, "source.h5ad")
        self.anndata.AnnData(
            X=self.sparse.csr_matrix([[1, 2]] * observations),
            obs=self.pandas.DataFrame(
                index=[f"cell-{index}" for index in range(observations)]
            ),
            var=self.pandas.DataFrame(index=["ENSG1", "ENSG2"]),
        ).write_h5ad(h5ad_path)
        json_path = os.path.join(tmpdir, "GSE1.json")
        Path(json_path).write_text(
            json.dumps([_package(h5ad_path)]),
            encoding="utf-8",
        )
        return h5ad_path, json_path

    def test_projector_adds_sample_metadata_without_combined_hook(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            _source, json_path = self._fixture(tmpdir)
            projector = RecordingProjector()

            result = JSON2H5ADConverter(
                metadata_projectors=[projector]
            ).convert(json_path=json_path, out=os.path.join(tmpdir, "out"))

            sample = self.anndata.read_h5ad(result.sample_h5ads["GSM1"])
            self.assertEqual(["Sample1", "Sample1"], sample.obs["private.sample"].tolist())
            self.assertEqual(
                ["GRCh38", "GRCh38"],
                sample.var["private.genome_build"].tolist(),
            )
            self.assertEqual("1.0", sample.uns["private"]["schema_version"])
            self.assertEqual(["private projection warning"], result.warnings)
            self.assertEqual("GSM1", projector.sample_contexts[0].sample_accession)
            self.assertEqual("GSE1", projector.sample_contexts[0].study_accession)
            self.assertEqual([], projector.combined_contexts)
            self.assertIsNone(result.combined_h5ad)

    def test_projector_cannot_overwrite_existing_metadata(self):
        class CollisionProjector:
            def project_sample(self, *, adata, context):
                return AnnDataMetadataProjection(
                    obs={"msc.sample.accession": "replacement"}
                )

        with tempfile.TemporaryDirectory() as tmpdir:
            _source, json_path = self._fixture(tmpdir)

            with self.assertRaisesRegex(
                ValueError, r"msc\.sample\.accession.*already exists"
            ):
                JSON2H5ADConverter(
                    metadata_projectors=[CollisionProjector()]
                ).convert(json_path=json_path, out=os.path.join(tmpdir, "out"))

    def test_projector_atomically_drops_and_renames_observation_metadata(self):
        class TransformingProjector:
            def project_sample(self, *, adata, context):
                return AnnDataMetadataProjection(
                    obs={"sample_name": "curated"},
                    obs_renames={"msc.sample.title": "author_sample_name"},
                    obs_drops=("msc.sample.accession",),
                )

        with tempfile.TemporaryDirectory() as tmpdir:
            _source, json_path = self._fixture(tmpdir)
            result = JSON2H5ADConverter(
                metadata_projectors=[TransformingProjector()]
            ).convert(json_path=json_path, out=os.path.join(tmpdir, "out"))

            sample = self.anndata.read_h5ad(result.sample_h5ads["GSM1"])
            self.assertNotIn("msc.sample.accession", sample.obs)
            self.assertNotIn("msc.sample.title", sample.obs)
            self.assertEqual(
                ["", ""], sample.obs["author_sample_name"].tolist()
            )
            self.assertEqual(["curated", "curated"], sample.obs["sample_name"].tolist())

    def test_projector_rejects_missing_rename_source(self):
        class MissingSourceProjector:
            def project_sample(self, *, adata, context):
                return AnnDataMetadataProjection(
                    obs_renames={"missing": "author_missing"}
                )

        with tempfile.TemporaryDirectory() as tmpdir:
            _source, json_path = self._fixture(tmpdir)
            with self.assertRaisesRegex(ValueError, "rename source 'missing'.*does not exist"):
                JSON2H5ADConverter(
                    metadata_projectors=[MissingSourceProjector()]
                ).convert(json_path=json_path, out=os.path.join(tmpdir, "out"))

    def test_projector_rejects_rename_target_collision(self):
        class CollisionProjector:
            def project_sample(self, *, adata, context):
                return AnnDataMetadataProjection(
                    obs_renames={"msc.sample.title": "msc.sample.accession"}
                )

        with tempfile.TemporaryDirectory() as tmpdir:
            _source, json_path = self._fixture(tmpdir)
            with self.assertRaisesRegex(
                ValueError, r"rename target 'msc\.sample\.accession'.*already exists"
            ):
                JSON2H5ADConverter(
                    metadata_projectors=[CollisionProjector()]
                ).convert(json_path=json_path, out=os.path.join(tmpdir, "out"))

    def test_projection_context_retains_declared_asset_path(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            source, json_path = self._fixture(tmpdir)
            projector = RecordingProjector()

            JSON2H5ADConverter(metadata_projectors=[projector]).convert(
                json_path=json_path,
                out=os.path.join(tmpdir, "out"),
            )

            self.assertEqual(source, projector.sample_contexts[0].asset.path)

    def test_projector_vector_lengths_must_match_axis(self):
        class InvalidVectorProjector:
            def project_sample(self, *, adata, context):
                return AnnDataMetadataProjection(obs={"private.vector": ["one"]})

        with tempfile.TemporaryDirectory() as tmpdir:
            _source, json_path = self._fixture(tmpdir, observations=2)

            with self.assertRaisesRegex(ValueError, "private.vector.*2"):
                JSON2H5ADConverter(
                    metadata_projectors=[InvalidVectorProjector()]
                ).convert(json_path=json_path, out=os.path.join(tmpdir, "out"))

    def test_projection_errors_fail_closed_without_final_artifacts(self):
        class InvalidProjector:
            def project_sample(self, *, adata, context):
                return AnnDataMetadataProjection(errors=("private validation failed",))

        with tempfile.TemporaryDirectory() as tmpdir:
            _source, json_path = self._fixture(tmpdir)
            out = Path(tmpdir) / "out"

            with self.assertRaisesRegex(
                AnnDataProjectionError, "private validation failed"
            ):
                JSON2H5ADConverter(
                    metadata_projectors=[InvalidProjector()]
                ).convert(json_path=json_path, out=str(out))

            self.assertEqual([], list(out.glob("*.h5ad")))
            self.assertEqual([], list(out.glob("*.json2h5ad.json")))

    def test_allow_invalid_writes_artifacts_and_manifest_diagnostics(self):
        class InvalidProjector:
            def project_sample(self, *, adata, context):
                return AnnDataMetadataProjection(errors=("private validation failed",))

        with tempfile.TemporaryDirectory() as tmpdir:
            _source, json_path = self._fixture(tmpdir)
            result = JSON2H5ADConverter(
                metadata_projectors=[InvalidProjector()]
            ).convert(
                json_path=json_path,
                out=os.path.join(tmpdir, "out"),
                allow_invalid=True,
            )

            self.assertEqual(["private validation failed"], result.errors)
            self.assertTrue(result.partial)
            self.assertTrue(Path(result.sample_h5ads["GSM1"]).is_file())
            self.assertIsNone(result.combined_h5ad)
            manifest = json.loads(Path(result.manifest_path).read_text(encoding="utf-8"))
            self.assertEqual(["private validation failed"], manifest["errors"])
            self.assertTrue(manifest["partial"])

    def test_structural_projection_errors_remain_unconditional_when_allowing_invalid(self):
        class WrongTypeProjector:
            def project_sample(self, *, adata, context):
                return {"errors": ["ignored"]}

        with tempfile.TemporaryDirectory() as tmpdir:
            _source, json_path = self._fixture(tmpdir)

            with self.assertRaisesRegex(TypeError, "AnnDataMetadataProjection"):
                JSON2H5ADConverter(
                    metadata_projectors=[WrongTypeProjector()]
                ).convert(
                    json_path=json_path,
                    out=os.path.join(tmpdir, "out"),
                    allow_invalid=True,
                )

    def test_combined_projection_hook_is_not_called_for_catalogue(self):
        class InvalidCombinedProjector:
            def project_sample(self, *, adata, context):
                return AnnDataMetadataProjection()

            def project_combined(self, *, adata, contexts):
                return AnnDataMetadataProjection(errors=("combined invalid",))

        with tempfile.TemporaryDirectory() as tmpdir:
            _source, json_path = self._fixture(tmpdir)
            out = Path(tmpdir) / "out"

            result = JSON2H5ADConverter(
                metadata_projectors=[InvalidCombinedProjector()]
            ).convert(json_path=json_path, out=str(out))

            self.assertTrue(Path(result.sample_h5ads["GSM1"]).is_file())
            self.assertIsNone(result.combined_h5ad)
            self.assertTrue(Path(result.manifest_path).is_file())

    def test_bundle_commit_failure_restores_previous_artifacts(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            _source, json_path = self._fixture(tmpdir)
            out = Path(tmpdir) / "out"
            converter = JSON2H5ADConverter()
            first = converter.convert(json_path=json_path, out=str(out))
            paths = [
                Path(first.sample_h5ads["GSM1"]),
                Path(first.manifest_path),
            ]
            originals = {path: path.read_bytes() for path in paths}
            real_replace = os.replace
            injected = False

            def fail_manifest_once(source, destination):
                nonlocal injected
                destination = Path(destination)
                if destination == Path(first.manifest_path) and not injected:
                    injected = True
                    raise OSError("injected bundle commit failure")
                return real_replace(source, destination)

            with patch(
                "meta_standards_converter.expression.catalogue.os.replace",
                side_effect=fail_manifest_once,
            ):
                with self.assertRaisesRegex(OSError, "injected bundle commit failure"):
                    converter.convert(
                        json_path=json_path,
                        out=str(out),
                        overwrite=True,
                    )

            self.assertTrue(injected)
            self.assertEqual(originals, {path: path.read_bytes() for path in paths})
            self.assertEqual([], list(out.glob(".json2h5ad-staging-*")))


if __name__ == "__main__":
    unittest.main()
