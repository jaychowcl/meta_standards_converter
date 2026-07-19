# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
import json
import gzip
import hashlib
import os
import shutil
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
    JSON2H5ADConverter,
    PipelineRun,
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
    def test_portable_locations_preserve_remote_urls(self):
        converter = JSON2H5ADConverter()

        value, scope = converter._portable_location(
            "https://example.org/data/source.h5ad", Path("/tmp/output")
        )

        self.assertEqual("https://example.org/data/source.h5ad", value)
        self.assertEqual("remote", scope)

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
                obs=self.pandas.DataFrame(
                    {"pipeline_qc_pass": [True, False]},
                    index=["cell1", "cell2"],
                ),
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
            self.assertEqual(["GSM1", "GSM1"], list(normalized.obs["msc_accession"]))
            self.assertEqual(["Control sample"] * 2, list(normalized.obs["msc_title"]))
            self.assertEqual(["Homo sapiens"] * 2, list(normalized.obs["msc_organism"]))
            self.assertEqual(["healthy"] * 2, list(normalized.obs["msc_disease"]))
            sectioned_columns = {
                "msc_accession": "msc.sample.accession",
                "msc_series_accession": "msc.series.accession",
                "msc_title": "msc.sample.title",
                "msc_description": "msc.sample.description",
                "msc_organism": "msc.sample.channel.organism.value",
                "msc_organism_taxid": "msc.sample.channel.organism.taxid",
                "msc_organism_part": "msc.sample.channel.organism_part",
                "msc_developmental_stage": "msc.sample.channel.developmental_stage",
                "msc_disease": "msc.sample.channel.disease",
                "msc_genotype": "msc.sample.channel.genotype",
                "msc_source_name": "msc.sample.channel.source",
                "msc_biomaterial_provider": "msc.sample.channel.biomaterial_provider",
                "msc_material_type": "msc.sample.channel.material_type",
                "msc_molecule": "msc.sample.channel.molecule",
                "msc_platform_accession": "msc.platform.accession",
                "msc_sra_accession": "msc.archive.sra_accession",
                "msc_ena_accession": "msc.archive.ena_accession",
                "msc_biosample_accession": "msc.archive.biosample_accession",
                "msc_sra_run_accessions": "msc.archive.sra_run_accessions",
                "msc_library_strategy": "msc.library.strategy",
                "msc_library_source": "msc.library.source",
                "msc_library_selection": "msc.library.selection",
                "msc_library_layout": "msc.library.layout",
                "msc_instrument_model": "msc.instrument.model",
                "msc_protocol_types": "msc.protocol.types",
                "msc_protocol_term_source_refs": "msc.protocol.term_source_refs",
                "msc_protocol_term_accession_numbers": "msc.protocol.term_accession_numbers",
                "msc_metadata_source": "msc.database.identifier",
                "msc_metadata_source_name": "msc.database.name",
                "msc_metadata_source_uri": "msc.database.uri",
                "msc_source_tier": "msc.asset.tier",
                "msc_source_uri": "msc.asset.uri",
                "msc_source_uri_scope": "msc.asset.uri_scope",
                "msc_modality": "msc.expression.modality",
            }
            for legacy, sectioned in sectioned_columns.items():
                with self.subTest(legacy=legacy, sectioned=sectioned):
                    self.assertIn(sectioned, normalized.obs)
                    self.assertEqual(
                        normalized.obs[legacy].tolist(),
                        normalized.obs[sectioned].tolist(),
                    )
            self.assertEqual([True, False], normalized.obs["pipeline_qc_pass"].tolist())
            self.assertEqual(
                ["healthy"] * 2,
                normalized.obs["msc.characteristics.disease"].tolist(),
            )
            self.assertFalse(any(column.startswith("geo_") for column in normalized.obs))
            self.assertEqual(["cell1", "cell2"], list(original.obs_names))
            self.assertEqual("h5ad", normalized.uns["meta_standards_converter"]["source_tier"])
            self.assertEqual("artifact_parent", normalized.uns["meta_standards_converter"]["path_base"])
            self.assertEqual("../source.h5ad", normalized.uns["meta_standards_converter"]["source_uri"])
            self.assertEqual(["../source.h5ad"] * 2, list(normalized.obs["msc_source_uri"]))

            manifest = json.loads(Path(result.manifest_path).read_text())
            self.assertEqual("artifact_parent", manifest["path_base"])
            self.assertEqual("../GSE1.json", manifest["source_json"])
            self.assertEqual("GSE1.h5ad", manifest["combined_h5ad"])
            self.assertEqual({"GSM1": "GSM1.h5ad"}, manifest["sample_h5ads"])
            self.assertEqual("../source.h5ad", manifest["assets"]["GSM1"]["path"])
            self.assertEqual("external", manifest["assets"]["GSM1"]["path_scope"])

    def test_does_not_duplicate_sample_accession_observation_name(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            source_path = os.path.join(tmpdir, "source.h5ad")
            self.anndata.AnnData(
                X=self.sparse.csr_matrix([[1]]),
                obs=self.pandas.DataFrame(index=["GSM1"]),
                var=self.pandas.DataFrame(index=["ENSG1"]),
            ).write_h5ad(source_path)
            json_path = self._write_json(tmpdir, package(source_path))

            result = json2h5ad().convert(json_path=json_path, out=os.path.join(tmpdir, "out"))

            sample = self.anndata.read_h5ad(result.sample_h5ads["GSM1"])
            combined = self.anndata.read_h5ad(result.combined_h5ad)
            self.assertEqual(["GSM1"], list(sample.obs_names))
            self.assertEqual(["GSM1"], list(combined.obs_names))

    def test_derives_organism_from_harmonized_and_raw_channels(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            source_path = os.path.join(tmpdir, "source.h5ad")
            self.anndata.AnnData(
                X=self.sparse.csr_matrix([[1]]),
                obs=self.pandas.DataFrame(index=["cell"]),
                var=self.pandas.DataFrame(index=["ENSG1"]),
            ).write_h5ad(source_path)
            data = package(source_path)
            data["sample"][0]["channel"] = [
                {
                    "organism": [{"taxid": "9606", "value": "human"}],
                    "hz_organism": [
                        {
                            "value": "Homo sapiens",
                            "id": "NCBITaxon_9606",
                            "onto": "ncbitaxon",
                        }
                    ],
                },
                {
                    "organism": [{"taxid": "10090", "value": "Mus musculus"}],
                    "hz_organism": [],
                },
                {"organism": [{"taxid": "10090", "value": "mus musculus"}]},
            ]
            json_path = self._write_json(tmpdir, data)

            result = json2h5ad().convert(json_path=json_path, out=os.path.join(tmpdir, "out"))

            converted = self.anndata.read_h5ad(result.sample_h5ads["GSM1"])
            self.assertEqual(["Homo sapiens; Mus musculus"], converted.obs["msc_organism"].unique().tolist())
            self.assertEqual(
                converted.obs["msc_organism"].tolist(),
                converted.obs["msc.sample.channel.organism.value"].tolist(),
            )
            self.assertEqual(
                ["9606; 10090"],
                converted.obs["msc_organism_taxid"].unique().tolist(),
            )
            fields = converted.uns["msc_miniml"]["fields"]
            self.assertIn("channel[0].hz_organism[0].value", set(fields["path"]))

    def test_derives_organism_from_scalar_harmonization_and_preserves_missing_as_empty(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            for case, channel, expected in (
                (
                    "scalar",
                    {
                        "organism": [{"value": "human"}],
                        "hz_organism": "Homo sapiens",
                        "hz_organism_id": "NCBITaxon_9606",
                        "hz_organism_onto": "ncbitaxon",
                    },
                    "Homo sapiens",
                ),
                ("missing", {}, ""),
            ):
                with self.subTest(case=case):
                    source_path = os.path.join(tmpdir, f"{case}.h5ad")
                    self.anndata.AnnData(
                        X=self.sparse.csr_matrix([[1]]),
                        obs=self.pandas.DataFrame(index=["cell"]),
                        var=self.pandas.DataFrame(index=["ENSG1"]),
                    ).write_h5ad(source_path)
                    data = package(source_path)
                    data["sample"][0]["channel"] = [channel]
                    json_path = self._write_json(tmpdir, data)

                    result = json2h5ad().convert(
                        json_path=json_path,
                        out=os.path.join(tmpdir, f"out-{case}"),
                    )

                    converted = self.anndata.read_h5ad(result.sample_h5ads["GSM1"])
                    self.assertEqual([expected], converted.obs["msc_organism"].unique().tolist())
                    self.assertEqual(
                        [expected],
                        converted.obs["msc.sample.channel.organism.value"].unique().tolist(),
                    )

    def test_reads_gzip_compressed_h5ad(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            source_path = os.path.join(tmpdir, "source.h5ad")
            compressed_path = f"{source_path}.gz"
            self.anndata.AnnData(
                X=self.sparse.csr_matrix([[1, 2]]),
                obs=self.pandas.DataFrame(index=["cell1"]),
                var=self.pandas.DataFrame(index=["Gene1", "Gene2"]),
            ).write_h5ad(source_path)
            with open(source_path, "rb") as source, gzip.open(compressed_path, "wb") as target:
                shutil.copyfileobj(source, target)
            os.unlink(source_path)
            json_path = self._write_json(tmpdir, package(compressed_path))

            result = json2h5ad().convert(
                json_path=json_path,
                out=os.path.join(tmpdir, "out"),
            )

            converted = self.anndata.read_h5ad(result.sample_h5ads["GSM1"])
            self.assertEqual((1, 2), converted.shape)
            self.assertEqual([[1, 2]], converted.X.toarray().tolist())
            self.assertEqual(0o660, Path(result.sample_h5ads["GSM1"]).stat().st_mode & 0o777)

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
            self.assertEqual({"GSM1", "GSM2"}, set(combined.obs["msc_accession"]))
            self.assertIn("msc_batch", combined.obs)
            self.assertEqual(
                combined.obs["msc_batch"].tolist(),
                combined.obs["msc.combination.batch"].tolist(),
            )
            self.assertTrue(self.sparse.issparse(combined.X))

    def test_enriches_msc_metadata_and_flattens_relevant_miniml(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = {}
            for sample_id in ("GSM1", "GSM2"):
                path = os.path.join(tmpdir, f"{sample_id}.h5ad")
                self.anndata.AnnData(
                    X=self.sparse.csr_matrix([[1]]),
                    obs=self.pandas.DataFrame(index=["cell"]),
                    var=self.pandas.DataFrame(index=["ENSG1"]),
                ).write_h5ad(path)
                paths[sample_id] = path
            data = {
                "schema_location": "MINiML.xsd",
                "version": "1.0",
                "database": [
                    {
                        "iid": "GEO",
                        "public_id": "GEO",
                        "name": "Gene Expression Omnibus (GEO)",
                        "web_link": "https://www.ncbi.nlm.nih.gov/geo",
                    }
                ],
                "contributor": [
                    {"iid": "C1", "person": {"first": "Alice", "last": "Example"}},
                    {"iid": "C2", "person": {"first": "Unrelated"}},
                ],
                "platform": [
                    {
                        "iid": "P1",
                        "accession": [{"database": "GEO", "value": "GPL1"}],
                        "contributor_ref": [{"ref": "C1"}],
                    },
                    {
                        "iid": "P2",
                        "accession": [{"database": "GEO", "value": "GPL2"}],
                        "contributor_ref": [{"ref": "C2"}],
                    },
                ],
                "series": {
                    "accession": [{"database": "GEO", "value": "GSE1"}],
                    "title": "Study title",
                    "summary": "GEO experiment summary",
                    "contributor_ref": [{"ref": "C1"}],
                    "sample_ref": [{"ref": "S1"}, {"ref": "S2"}],
                    "pubmed_publication": [
                        {
                            "pubmed_id": "123",
                            "doi": "10.1/example",
                            "title": "Citation title",
                            "author_list": "A Example, B Example",
                            "status": "published",
                            "abstract": "must not be embedded",
                            "full_text": "must not be embedded",
                            "article_body": {"section": "must not be embedded"},
                        }
                    ],
                },
                "sample": [
                    {
                        "iid": "S1",
                        "accession": [{"database": "GEO", "value": "GSM1"}],
                        "title": "Sample one",
                        "description": "Sample description",
                        "supplementary_data": [{"value": paths["GSM1"]}],
                        "platform_ref": {"ref": "P1"},
                        "contact_ref": [{"ref": "C1"}],
                        "library_strategy": "RNA-Seq",
                        "library_source": ["transcriptomic", "TRANSCRIPTOMIC"],
                        "library_selection": "cDNA",
                        "instrument_model": {"predefined": "Illumina Test"},
                        "ena_accession": "SRS1",
                        "sra_accession": "SRX1",
                        "sra_run": [
                            {
                                "run": "SRR1",
                                "biosample": "SAMN1",
                                "library_layout": "PAIRED",
                                "fastq_files": [{"uri": "https://example/R1.fastq.gz"}],
                            },
                            {
                                "run": "SRR2",
                                "biosample": "SAMN1",
                                "library_layout": "PAIRED",
                            },
                        ],
                        "channel": [
                            {
                                "source": "blood",
                                "molecule": "total RNA",
                                "biomaterial_provider": ["Example Biobank"],
                                "organism": [{"taxid": "9606", "value": "Homo sapiens"}],
                                "characteristics": [
                                    {"tag": "cell type", "value": "Treg"},
                                    {"tag": "developmental stage", "value": "adult"},
                                    {"tag": "treatment", "value": "CPI-703"},
                                    {"tag": "hz_cell_type", "value": "regulatory T cell"},
                                    {"tag": "hz_cell_type_id", "value": "CL:0000815"},
                                    {"tag": "hz_cell_type_onto", "value": "cl"},
                                ],
                                "treatment_protocol": "Long treatment protocol",
                            },
                            {
                                "source": "blood",
                                "characteristics": [
                                    {"tag": "cell-type", "value": "Activated Treg"},
                                    {"tag": "treatment", "value": "CPI-703"},
                                ],
                            },
                        ],
                    },
                    {
                        "iid": "S2",
                        "accession": [{"database": "GEO", "value": "GSM2"}],
                        "title": "Sample two",
                        "supplementary_data": [{"value": paths["GSM2"]}],
                        "platform_ref": {"ref": "P2"},
                        "contact_ref": [{"ref": "C2"}],
                        "channel": [{"characteristics": [{"tag": "dose", "value": "5 uM"}]}],
                    },
                ],
            }
            json_path = self._write_json(tmpdir, data)

            result = json2h5ad().convert(json_path=json_path, out=os.path.join(tmpdir, "out"))

            first = self.anndata.read_h5ad(result.sample_h5ads["GSM1"])
            second = self.anndata.read_h5ad(result.sample_h5ads["GSM2"])
            combined = self.anndata.read_h5ad(result.combined_h5ad)
            self.assertEqual(["GEO"], first.obs["msc_metadata_source"].unique().tolist())
            self.assertEqual(
                ["Gene Expression Omnibus (GEO)"],
                first.obs["msc_metadata_source_name"].unique().tolist(),
            )
            self.assertEqual(["GPL1"], first.obs["msc_platform_accession"].unique().tolist())
            self.assertEqual(["GPL2"], second.obs["msc_platform_accession"].unique().tolist())
            self.assertEqual(["transcriptomic"], first.obs["msc_library_source"].unique().tolist())
            self.assertEqual(["SRR1; SRR2"], first.obs["msc_sra_run_accessions"].unique().tolist())
            self.assertEqual(["SAMN1"], first.obs["msc_biosample_accession"].unique().tolist())
            self.assertEqual(["Illumina Test"], first.obs["msc_instrument_model"].unique().tolist())
            self.assertEqual(["adult"], first.obs["msc_developmental_stage"].unique().tolist())
            self.assertEqual(["Example Biobank"], first.obs["msc_biomaterial_provider"].unique().tolist())
            self.assertEqual(["RNA"], first.obs["msc_material_type"].unique().tolist())
            self.assertEqual(
                ["sample treatment protocol"],
                first.obs["msc_protocol_types"].unique().tolist(),
            )
            self.assertEqual(["EFO"], first.obs["msc_protocol_term_source_refs"].unique().tolist())
            self.assertEqual(
                ["EFO_0003809"],
                first.obs["msc_protocol_term_accession_numbers"].unique().tolist(),
            )
            self.assertEqual(
                ["Treg; Activated Treg"],
                first.obs["msc_characteristic_cell_type"].unique().tolist(),
            )
            self.assertEqual(["CPI-703"], first.obs["msc_characteristic_treatment"].unique().tolist())
            self.assertEqual(
                ["regulatory T cell"],
                first.obs["msc.characteristics.hz_cell_type"].unique().tolist(),
            )
            self.assertEqual(
                ["CL:0000815"],
                first.obs["msc.characteristics.hz_cell_type_id"].unique().tolist(),
            )
            self.assertEqual(
                ["cl"],
                first.obs["msc.characteristics.hz_cell_type_onto"].unique().tolist(),
            )
            self.assertEqual([""], first.obs["msc_characteristic_dose"].unique().tolist())
            self.assertEqual([""], second.obs["msc_characteristic_cell_type"].unique().tolist())
            self.assertEqual([""], second.obs["msc.characteristics.hz_cell_type"].unique().tolist())
            self.assertFalse(any(column.startswith("geo_") for column in combined.obs))
            self.assertIn("msc_batch", combined.obs)

            miniml = first.uns["msc_miniml"]
            self.assertEqual("1.0", miniml["schema_version"])
            self.assertEqual("citation_metadata_only", miniml["publication_policy"])
            self.assertEqual(hashlib.sha256(Path(json_path).read_bytes()).hexdigest(), miniml["source_sha256"])
            fields = miniml["fields"]
            sample_entities = set(fields.loc[fields["entity_type"] == "sample", "entity_id"])
            contributor_entities = set(fields.loc[fields["entity_type"] == "contributor", "entity_id"])
            platform_entities = set(fields.loc[fields["entity_type"] == "platform", "entity_id"])
            self.assertEqual({"GSM1"}, sample_entities)
            self.assertEqual({"C1"}, contributor_entities)
            self.assertEqual({"GPL1"}, platform_entities)
            paths_in_sample = set(fields["path"])
            self.assertIn("pubmed_publication[0].title", paths_in_sample)
            self.assertIn("channel[0].treatment_protocol", paths_in_sample)
            self.assertNotIn("pubmed_publication[0].abstract", paths_in_sample)
            self.assertNotIn("pubmed_publication[0].full_text", paths_in_sample)
            self.assertFalse(any("article_body" in path for path in paths_in_sample))

            combined_fields = combined.uns["msc_miniml"]["fields"]
            self.assertEqual(
                {"GSM1", "GSM2"},
                set(combined_fields.loc[combined_fields["entity_type"] == "sample", "entity_id"]),
            )
            self.assertEqual(
                {"C1", "C2"},
                set(combined_fields.loc[combined_fields["entity_type"] == "contributor", "entity_id"]),
            )
            self.assertEqual(
                {"GPL1", "GPL2"},
                set(combined_fields.loc[combined_fields["entity_type"] == "platform", "entity_id"]),
            )

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

    def test_matching_harmonized_organisms_allow_combination(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            samples = []
            for sample_id, raw_organism in (("GSM1", "human"), ("GSM2", "Homo sapiens")):
                path = os.path.join(tmpdir, f"{sample_id}.h5ad")
                self.anndata.AnnData(
                    X=self.sparse.csr_matrix([[1]]),
                    obs=self.pandas.DataFrame(index=["cell"]),
                    var=self.pandas.DataFrame({"gene_ids": ["ENSG1"]}, index=["ENSG1"]),
                ).write_h5ad(path)
                sample = package(path, accession=sample_id)["sample"][0]
                sample["channel"] = [
                    {
                        "organism": [{"value": raw_organism}],
                        "hz_organism": [
                            {
                                "value": "Homo sapiens",
                                "id": "NCBITaxon_9606",
                                "onto": "ncbitaxon",
                            }
                        ],
                    }
                ]
                samples.append(sample)
            data = package()
            data["sample"] = samples
            json_path = self._write_json(tmpdir, data)

            result = json2h5ad().convert(json_path=json_path, out=os.path.join(tmpdir, "out"))

            self.assertIsNotNone(result.combined_h5ad)
            self.assertEqual([], result.failures)
            combined = self.anndata.read_h5ad(result.combined_h5ad)
            self.assertEqual(["Homo sapiens"], combined.obs["msc_organism"].unique().tolist())

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
                handle.write(
                    "gene_id\tgene_name\tGSM1\tGSM2\n"
                    "ENSG1\tGENE1\t1\t9\n"
                    "ENSG2\tGENE2\t2\t8\n"
                )
            with open(tpm, "w", encoding="utf-8") as handle:
                handle.write(
                    "gene_id\tgene_name\tGSM1\tGSM2\n"
                    "ENSG1\tGENE1\t3\t7\n"
                    "ENSG2\tGENE2\t4\t6\n"
                )
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
            self.assertEqual(["GENE1", "GENE2"], adata.var["gene_name"].tolist())

    def test_annotation_provenance_is_written_to_h5ad_and_manifest(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            source = os.path.join(tmpdir, "source.h5ad")
            annotation = os.path.join(tmpdir, "genes.gtf")
            with open(annotation, "w", encoding="utf-8") as handle:
                handle.write("chr1\ttest\texon\t1\t2\t.\t+\t.\tgene_id \"g1\";\n")
            digest = hashlib.sha256(Path(annotation).read_bytes()).hexdigest()
            self.anndata.AnnData(
                X=self.sparse.csr_matrix([[1]]),
                obs=self.pandas.DataFrame(index=["cell"]),
                var=self.pandas.DataFrame(index=["ENSG1"]),
            ).write_h5ad(source)
            json_path = self._write_json(tmpdir, package())
            asset = Asset(
                "GSM1", source, "h5ad", source="nfcore", reference="GRCh38",
                annotation_source=annotation, annotation_format="gtf",
                annotation_sha256=digest, effective_annotation=annotation,
            )

            result = json2h5ad().convert(
                json_path=json_path, out=os.path.join(tmpdir, "out"), explicit_assets=[asset]
            )

            adata = self.anndata.read_h5ad(result.sample_h5ads["GSM1"])
            provenance = adata.uns["meta_standards_converter"]
            self.assertEqual("../genes.gtf", provenance["annotation_source"])
            self.assertEqual("external", provenance["annotation_source_scope"])
            self.assertEqual("gtf", provenance["annotation_format"])
            self.assertEqual(digest, provenance["annotation_sha256"])
            manifest = json.loads(Path(result.manifest_path).read_text())
            manifest_asset = manifest["assets"]["GSM1"]
            self.assertEqual("../genes.gtf", manifest_asset["annotation_source"])
            self.assertEqual("external", manifest_asset["annotation_source_scope"])
            self.assertEqual(digest, manifest_asset["annotation_sha256"])

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
                    runs=[
                        PipelineRun(
                            pipeline="scrnaseq",
                            revision="4.2.0",
                            command=[
                                "nextflow",
                                "run",
                                "nf-core/scrnaseq",
                                "-work-dir",
                                os.path.dirname(self.output),
                                "https://example.org/reference.fa",
                            ],
                            work_dir=os.path.dirname(self.output),
                            out_dir=os.path.dirname(self.output),
                            warnings=["Unrecognized config option 'example'"],
                        )
                    ],
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
            self.assertEqual(["scrnaseq: Unrecognized config option 'example'"], result.warnings)
            manifest = json.loads(Path(result.manifest_path).read_text())
            self.assertEqual(result.warnings, manifest["warnings"])
            self.assertEqual(
                ["Unrecognized config option 'example'"],
                manifest["pipeline_runs"][0]["warnings"],
            )
            recorded_command = manifest["pipeline_runs"][0]["command"]
            self.assertFalse(os.path.isabs(recorded_command[4]))
            self.assertEqual("https://example.org/reference.fa", recorded_command[5])

    def test_study_h5ad_is_split_by_canonical_and_legacy_sample_accession(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            for accession_column in ("msc.sample.accession", "msc_accession", "geo_accession"):
                with self.subTest(accession_column=accession_column):
                    source_path = os.path.join(tmpdir, f"study-{accession_column}.h5ad")
                    self.anndata.AnnData(
                        X=self.sparse.csr_matrix([[1], [2]]),
                        obs=self.pandas.DataFrame(
                            {accession_column: ["GSM1", "GSM2"]},
                            index=["cell1", "cell2"],
                        ),
                        var=self.pandas.DataFrame(index=["ENSG1"]),
                    ).write_h5ad(source_path)
                    data = package(accession="GSM1")
                    data["sample"].append(package(accession="GSM2")["sample"][0])
                    json_path = self._write_json(tmpdir, data)

                    result = json2h5ad().convert(
                        json_path=json_path,
                        out=os.path.join(tmpdir, f"out-{accession_column}"),
                        explicit_assets=[Asset("GSE1", source_path, "h5ad", source="manifest")],
                    )

                    first = self.anndata.read_h5ad(result.sample_h5ads["GSM1"])
                    second = self.anndata.read_h5ad(result.sample_h5ads["GSM2"])
                    self.assertEqual((1, 1), first.shape)
                    self.assertEqual([[1]], first.X.toarray().tolist())
                    self.assertEqual([[2]], second.X.toarray().tolist())


if __name__ == "__main__":
    unittest.main()
