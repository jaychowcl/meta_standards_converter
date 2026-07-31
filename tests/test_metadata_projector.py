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


ROOT = os.path.dirname(os.path.dirname(__file__))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from meta_standards_converter.converters.json2h5ad import (  # noqa: E402
    AnnDataMetadataProjection,
    AnnDataProjectionError,
    JSON2H5ADConverter,
)


def _package(h5ad_path: str) -> dict:
    return {
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

    def test_projector_adds_sample_and_combined_metadata_and_warnings(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            _source, json_path = self._fixture(tmpdir)
            projector = RecordingProjector()

            result = JSON2H5ADConverter(
                metadata_projectors=[projector]
            ).convert(json_path=json_path, out=os.path.join(tmpdir, "out"))

            sample = self.anndata.read_h5ad(result.sample_h5ads["GSM1"])
            combined = self.anndata.read_h5ad(result.combined_h5ad)
            self.assertEqual(["Sample1", "Sample1"], sample.obs["private.sample"].tolist())
            self.assertEqual(
                ["GRCh38", "GRCh38"],
                sample.var["private.genome_build"].tolist(),
            )
            self.assertEqual("1.0", sample.uns["private"]["schema_version"])
            self.assertEqual(1, combined.uns["private_combined"]["sample_count"])
            self.assertEqual(["private projection warning"], result.warnings)
            self.assertEqual("GSM1", projector.sample_contexts[0].sample_accession)
            self.assertEqual("GSE1", projector.sample_contexts[0].study_accession)
            self.assertEqual(1, len(projector.combined_contexts[0]))

    def test_projector_cannot_overwrite_existing_metadata(self):
        class CollisionProjector:
            def project_sample(self, *, adata, context):
                return AnnDataMetadataProjection(obs={"msc_accession": "replacement"})

        with tempfile.TemporaryDirectory() as tmpdir:
            _source, json_path = self._fixture(tmpdir)

            with self.assertRaisesRegex(ValueError, "msc_accession.*already exists"):
                JSON2H5ADConverter(
                    metadata_projectors=[CollisionProjector()]
                ).convert(json_path=json_path, out=os.path.join(tmpdir, "out"))

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
            self.assertTrue(Path(result.combined_h5ad).is_file())
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


if __name__ == "__main__":
    unittest.main()
