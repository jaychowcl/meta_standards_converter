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
    Asset,
    AssetDownloader,
    AssetManifest,
    ConversionResult,
    RawProcessingResult,
    SourcePlanner,
    json2h5ad,
)


def package(*files, accession="GSM1"):
    supplementary_data = [{"value": path} for path in files]
    return {
        "series": {"accession": [{"value": "GSE1"}]},
        "sample": [
            {
                "iid": "Sample1",
                "accession": [{"value": accession}],
                "supplementary_data": supplementary_data,
                "sra_run": [
                    {
                        "run": "SRR1",
                        "fastq_files": [
                            {"uri": "https://example/R1.fastq.gz", "md5": "a"},
                            {"uri": "https://example/R2.fastq.gz", "md5": "b"},
                        ],
                    }
                ],
            }
        ],
    }


class TestSourcePlanner(unittest.TestCase):
    def test_discovers_and_prefers_h5ad_over_matrix_and_raw(self):
        planner = SourcePlanner()

        plan = planner.plan([package("counts.tsv.gz", "provided.h5ad")])

        self.assertEqual("h5ad", plan["GSM1"].kind)
        self.assertEqual("provided.h5ad", plan["GSM1"].path)
        self.assertEqual("json", plan["GSM1"].source)

    def test_explicit_asset_overrides_json_asset(self):
        planner = SourcePlanner()
        explicit = Asset(
            scope_id="GSM1",
            path="local/override.h5ad",
            kind="h5ad",
            source="cli",
        )

        plan = planner.plan([package("provided.h5ad")], explicit_assets=[explicit])

        self.assertEqual(explicit, plan["GSM1"])

    def test_force_reprocess_uses_raw_fastqs(self):
        planner = SourcePlanner()

        plan = planner.plan([package("provided.h5ad")], force_reprocess=True)

        self.assertEqual("raw", plan["GSM1"].kind)
        self.assertEqual(2, len(plan["GSM1"].members))

    def test_separate_raw_data_entries_are_grouped_for_one_sample(self):
        data = package()
        data["sample"][0]["sra_run"] = []
        data["sample"][0]["raw_data"] = [
            {"value": "GSM1_R1.fastq.gz"},
            {"value": "GSM1_R2.fastq.gz"},
        ]

        plan = SourcePlanner().plan([data], force_reprocess=True)

        self.assertEqual(2, len(plan["GSM1"].members))

    def test_separate_cli_fastqs_are_grouped_for_one_sample(self):
        manifest = AssetManifest()
        assets = [
            manifest.parse_spec("GSM1=GSM1_R1.fastq.gz"),
            manifest.parse_spec("GSM1=GSM1_R2.fastq.gz"),
        ]

        plan = SourcePlanner().plan([package()], explicit_assets=assets, force_reprocess=True)

        self.assertEqual(2, len(plan["GSM1"].members))

    def test_force_reprocess_requires_raw_fastqs(self):
        data = package("provided.h5ad")
        data["sample"][0]["sra_run"] = []

        with self.assertRaisesRegex(ValueError, "GSM1.*raw FASTQ"):
            SourcePlanner().plan([data], force_reprocess=True)

    def test_manifest_asset_has_highest_precedence(self):
        planner = SourcePlanner()
        cli_asset = Asset("GSM1", "cli.h5ad", "h5ad", source="cli")
        manifest_asset = Asset("GSM1", "manifest.h5ad", "h5ad", source="manifest")

        plan = planner.plan(
            [package("json.h5ad")],
            explicit_assets=[cli_asset, manifest_asset],
        )

        self.assertEqual("manifest.h5ad", plan["GSM1"].path)


class TestAssetInputs(unittest.TestCase):
    def test_manifest_loads_processed_and_grouped_raw_assets(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "assets.csv")
            with open(path, "w", encoding="utf-8") as handle:
                handle.write(
                    "scope_id,path,kind,role,orientation,read,lane\n"
                    "GSM1,/data/GSM1.h5ad,h5ad,primary,auto,,\n"
                    "GSM2,https://example/GSM2_R1.fastq.gz,raw,primary,auto,1,L001\n"
                    "GSM2,https://example/GSM2_R2.fastq.gz,raw,primary,auto,2,L001\n"
                )

            assets = AssetManifest().load(path)

            self.assertEqual(2, len(assets))
            self.assertEqual("manifest", assets[0].source)
            self.assertEqual("h5ad", assets[0].kind)
            self.assertEqual(2, len(assets[1].members))
            self.assertEqual("2", assets[1].members[1]["read"])

    def test_cli_asset_spec_infers_kind(self):
        asset = AssetManifest().parse_spec("GSM1=/data/GSM1.h5ad")

        self.assertEqual(Asset("GSM1", "/data/GSM1.h5ad", "h5ad", source="cli"), asset)

    def test_streaming_downloader_caches_and_verifies_md5(self):
        class Response:
            def raise_for_status(self):
                return None

            def iter_content(self, chunk_size):
                return iter([b"abc", b"123"])

        class Session:
            def __init__(self):
                self.calls = []

            def get(self, url, **kwargs):
                self.calls.append((url, kwargs))
                return Response()

        with tempfile.TemporaryDirectory() as tmpdir:
            session = Session()
            downloader = AssetDownloader(cache_dir=tmpdir, session=session)

            first = downloader.localize(
                "https://example/data.h5ad",
                md5="e99a18c428cb38d5f260853678922e03",
            )
            second = downloader.localize(
                "https://example/data.h5ad",
                md5="e99a18c428cb38d5f260853678922e03",
            )

            self.assertEqual(first, second)
            self.assertEqual(b"abc123", Path(first).read_bytes())
            self.assertEqual(1, len(session.calls))
            self.assertTrue(session.calls[0][1]["stream"])


class TestConversionContract(unittest.TestCase):
    def test_result_exposes_primary_combined_output(self):
        result = ConversionResult(
            study_accession="GSE1",
            combined_h5ad="out/GSE1.h5ad",
            sample_h5ads={"GSM1": "out/GSM1.h5ad"},
        )

        self.assertEqual("out/GSE1.h5ad", result.primary_h5ad)
        self.assertFalse(result.partial)

    def test_result_is_partial_when_combination_fails(self):
        result = ConversionResult(
            study_accession="GSE1",
            sample_h5ads={"GSM1": "out/GSM1.h5ad"},
            failures=["incompatible genome builds"],
        )

        self.assertEqual("out/GSM1.h5ad", result.primary_h5ad)
        self.assertTrue(result.partial)

    def test_missing_json_path_raises_file_not_found(self):
        with self.assertRaises(FileNotFoundError):
            json2h5ad().convert(json_path="/missing/GSE1.json", out=".")

    def test_empty_package_list_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            json_path = os.path.join(tmpdir, "GSE1.json")
            with open(json_path, "w", encoding="utf-8") as handle:
                json.dump([], handle)

            with self.assertRaisesRegex(ValueError, "non-empty list"):
                json2h5ad().convert(json_path=json_path, out=tmpdir)


class TestProcessedAssetConversion(unittest.TestCase):
    def setUp(self):
        import anndata
        import numpy
        import pandas
        from scipy import sparse

        self.anndata = anndata
        self.numpy = numpy
        self.pandas = pandas
        self.sparse = sparse

    def _write_json(self, tmpdir, data):
        path = os.path.join(tmpdir, "GSE1.json")
        with open(path, "w", encoding="utf-8") as handle:
            json.dump([data], handle)
        return path

    def test_normalizes_supplied_h5ad_without_mutating_source(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            source_path = os.path.join(tmpdir, "source.h5ad")
            source = self.anndata.AnnData(
                X=self.sparse.csr_matrix([[1, 2], [3, 4]]),
                obs=self.pandas.DataFrame(index=["cell1", "cell2"]),
                var=self.pandas.DataFrame(
                    {"gene_ids": ["ENSG1", "ENSG2"]},
                    index=["Gene1", "Gene2"],
                ),
            )
            source.write_h5ad(source_path)
            data = package(source_path)
            data["sample"][0]["title"] = "Control sample"
            data["sample"][0]["channel"] = [
                {
                    "source": "blood",
                    "organism": [{"taxid": "9606", "value": "Homo sapiens"}],
                    "characteristics": [{"tag": "disease", "value": "healthy"}],
                }
            ]
            json_path = self._write_json(tmpdir, data)
            out = os.path.join(tmpdir, "out")

            result = json2h5ad().convert(json_path=json_path, out=out)

            normalized = self.anndata.read_h5ad(result.sample_h5ads["GSM1"])
            original = self.anndata.read_h5ad(source_path)
            self.assertEqual(["cell1-GSM1", "cell2-GSM1"], list(normalized.obs_names))
            self.assertEqual(["GSM1", "GSM1"], list(normalized.obs["geo_accession"]))
            self.assertEqual(["Control sample"] * 2, list(normalized.obs["geo_title"]))
            self.assertEqual(["Homo sapiens"] * 2, list(normalized.obs["geo_organism"]))
            self.assertEqual(["healthy"] * 2, list(normalized.obs["geo_disease"]))
            self.assertEqual(["cell1", "cell2"], list(original.obs_names))
            self.assertEqual("h5ad", normalized.uns["meta_standards_converter"]["source_tier"])

    def test_reads_delimited_gene_by_observation_matrix_as_sparse_counts(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            matrix_path = os.path.join(tmpdir, "counts.tsv")
            with open(matrix_path, "w", encoding="utf-8") as handle:
                handle.write("gene\tcell1\tcell2\nGene1\t1\t0\nGene2\t2\t3\n")
            json_path = self._write_json(tmpdir, package(matrix_path))

            result = json2h5ad().convert(
                json_path=json_path,
                out=os.path.join(tmpdir, "out"),
                matrix_orientation="genes-by-observations",
            )

            adata = self.anndata.read_h5ad(result.sample_h5ads["GSM1"])
            self.assertTrue(self.sparse.issparse(adata.X))
            self.assertEqual((2, 2), adata.shape)
            self.assertEqual(["Gene1", "Gene2"], list(adata.var_names))
            self.assertEqual([[1, 2], [0, 3]], adata.X.toarray().tolist())

    def test_combines_compatible_samples_with_outer_sparse_join(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = []
            for sample_id, genes in (("GSM1", ["ENSG1", "ENSG2"]), ("GSM2", ["ENSG2", "ENSG3"])):
                path = os.path.join(tmpdir, f"{sample_id}.h5ad")
                adata = self.anndata.AnnData(
                    X=self.sparse.csr_matrix([[1, 2]]),
                    obs=self.pandas.DataFrame(index=["cell"]),
                    var=self.pandas.DataFrame({"gene_ids": genes}, index=genes),
                )
                adata.write_h5ad(path)
                paths.append(path)
            first = package(paths[0], accession="GSM1")
            second = package(paths[1], accession="GSM2")["sample"][0]
            first["sample"].append(second)
            json_path = self._write_json(tmpdir, first)

            result = json2h5ad().convert(json_path=json_path, out=os.path.join(tmpdir, "out"))

            combined = self.anndata.read_h5ad(result.combined_h5ad)
            self.assertEqual((2, 3), combined.shape)
            self.assertEqual({"GSM1", "GSM2"}, set(combined.obs["geo_accession"]))
            self.assertTrue(self.sparse.issparse(combined.X))

    def test_incompatible_organisms_keep_samples_and_mark_partial(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            samples = []
            for sample_id, organism in (("GSM1", "Homo sapiens"), ("GSM2", "Mus musculus")):
                path = os.path.join(tmpdir, f"{sample_id}.h5ad")
                self.anndata.AnnData(
                    X=self.sparse.csr_matrix([[1]]),
                    obs=self.pandas.DataFrame(index=["cell"]),
                    var=self.pandas.DataFrame({"gene_ids": ["GENE1"]}, index=["GENE1"]),
                ).write_h5ad(path)
                sample = package(path, accession=sample_id)["sample"][0]
                sample["channel"] = [{"organism": [{"value": organism}]}]
                samples.append(sample)
            data = package()
            data["sample"] = samples
            json_path = self._write_json(tmpdir, data)

            result = json2h5ad().convert(json_path=json_path, out=os.path.join(tmpdir, "out"))

            self.assertIsNone(result.combined_h5ad)
            self.assertEqual({"GSM1", "GSM2"}, set(result.sample_h5ads))
            self.assertTrue(result.partial)
            self.assertIn("organisms", result.failures[0])

    def test_incompatible_declared_references_keep_samples_and_mark_partial(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            samples = []
            for sample_id, genome in (("GSM1", "GRCh37"), ("GSM2", "GRCh38")):
                path = os.path.join(tmpdir, f"{sample_id}.h5ad")
                adata = self.anndata.AnnData(
                    X=self.sparse.csr_matrix([[1]]),
                    obs=self.pandas.DataFrame(index=["cell"]),
                    var=self.pandas.DataFrame({"gene_ids": ["ENSG1"]}, index=["ENSG1"]),
                )
                adata.uns["genome"] = genome
                adata.write_h5ad(path)
                sample = package(path, accession=sample_id)["sample"][0]
                sample["channel"] = [{"organism": [{"value": "Homo sapiens"}]}]
                samples.append(sample)
            data = package()
            data["sample"] = samples
            json_path = self._write_json(tmpdir, data)

            result = json2h5ad().convert(json_path=json_path, out=os.path.join(tmpdir, "out"))

            self.assertIsNone(result.combined_h5ad)
            self.assertEqual({"GSM1", "GSM2"}, set(result.sample_h5ads))
            self.assertIn("reference builds", result.failures[0])

    def test_incompatible_bulk_and_single_cell_modalities_are_not_combined(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            samples = []
            for sample_id, library_source in (
                ("GSM1", "single cell transcriptomic"),
                ("GSM2", "TRANSCRIPTOMIC"),
            ):
                path = os.path.join(tmpdir, f"{sample_id}.h5ad")
                self.anndata.AnnData(
                    X=self.sparse.csr_matrix([[1]]),
                    obs=self.pandas.DataFrame(index=["cell"]),
                    var=self.pandas.DataFrame({"gene_ids": ["ENSG1"]}, index=["ENSG1"]),
                ).write_h5ad(path)
                sample = package(path, accession=sample_id)["sample"][0]
                sample["library_source"] = library_source
                sample["channel"] = [{"organism": [{"value": "Homo sapiens"}]}]
                samples.append(sample)
            data = package()
            data["sample"] = samples
            json_path = self._write_json(tmpdir, data)

            result = json2h5ad().convert(json_path=json_path, out=os.path.join(tmpdir, "out"))

            self.assertIsNone(result.combined_h5ad)
            self.assertIn("modalities", result.failures[0])

    def test_rnaseq_counts_select_sample_column_and_add_tpm_layer(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            counts = os.path.join(tmpdir, "salmon.merged.gene_counts.tsv")
            tpm = os.path.join(tmpdir, "salmon.merged.gene_tpm.tsv")
            with open(counts, "w", encoding="utf-8") as handle:
                handle.write("gene\tGSM1\tGSM2\nENSG1\t1\t9\nENSG2\t2\t8\n")
            with open(tpm, "w", encoding="utf-8") as handle:
                handle.write("gene\tGSM1\tGSM2\nENSG1\t3\t7\nENSG2\t4\t6\n")
            json_path = self._write_json(tmpdir, package())
            asset = Asset(
                "GSM1",
                counts,
                "matrix",
                role="rnaseq_counts",
                source="nfcore",
                features_path=tpm,
                orientation="genes-by-observations",
            )

            result = json2h5ad().convert(
                json_path=json_path,
                out=os.path.join(tmpdir, "out"),
                explicit_assets=[asset],
            )

            adata = self.anndata.read_h5ad(result.sample_h5ads["GSM1"])
            self.assertEqual((1, 2), adata.shape)
            self.assertEqual([[1, 2]], adata.X.toarray().tolist())
            self.assertEqual([[3, 4]], adata.layers["tpm"].toarray().tolist())

    def test_raw_assets_are_replaced_by_nfcore_outputs(self):
        class FakeRunner:
            def __init__(self, output):
                self.output = output
                self.calls = []

            def process(self, assets, **kwargs):
                self.calls.append((assets, kwargs))
                return RawProcessingResult(
                    assets={"GSM1": Asset("GSM1", self.output, "h5ad", source="nfcore")},
                    retained_h5ads=[self.output],
                )

        with tempfile.TemporaryDirectory() as tmpdir:
            pipeline_h5ad = os.path.join(tmpdir, "pipeline.h5ad")
            self.anndata.AnnData(
                X=self.sparse.csr_matrix([[1]]),
                obs=self.pandas.DataFrame(index=["cell"]),
                var=self.pandas.DataFrame(index=["ENSG1"]),
            ).write_h5ad(pipeline_h5ad)
            data = package()
            json_path = self._write_json(tmpdir, data)
            runner = FakeRunner(pipeline_h5ad)

            result = json2h5ad(pipeline_runner=runner).convert(
                json_path=json_path,
                out=os.path.join(tmpdir, "out"),
                force_reprocess=True,
                pipeline="scrnaseq",
                genome="GRCh38",
            )

            self.assertEqual([pipeline_h5ad], result.retained_h5ads)
            self.assertEqual("raw", runner.calls[0][0]["GSM1"].kind)
            self.assertEqual("scrnaseq", runner.calls[0][1]["pipeline"])
            self.assertEqual("GRCh38", runner.calls[0][1]["genome"])

    def test_study_h5ad_is_split_by_sample_accession(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            source_path = os.path.join(tmpdir, "study.h5ad")
            self.anndata.AnnData(
                X=self.sparse.csr_matrix([[1], [2]]),
                obs=self.pandas.DataFrame(
                    {"sample_id": ["GSM1", "GSM2"]},
                    index=["cell1", "cell2"],
                ),
                var=self.pandas.DataFrame(index=["ENSG1"]),
            ).write_h5ad(source_path)
            data = package(accession="GSM1")
            data["sample"].append(package(accession="GSM2")["sample"][0])
            json_path = self._write_json(tmpdir, data)

            result = json2h5ad().convert(
                json_path=json_path,
                out=os.path.join(tmpdir, "out"),
                explicit_assets=[Asset("GSE1", source_path, "h5ad", source="manifest")],
            )

            first = self.anndata.read_h5ad(result.sample_h5ads["GSM1"])
            second = self.anndata.read_h5ad(result.sample_h5ads["GSM2"])
            self.assertEqual((1, 1), first.shape)
            self.assertEqual([[1]], first.X.toarray().tolist())
            self.assertEqual([[2]], second.X.toarray().tolist())


if __name__ == "__main__":
    unittest.main()
