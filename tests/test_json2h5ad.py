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

import pytest


ROOT = os.path.dirname(os.path.dirname(__file__))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from meta_standards_converter.converters.json2h5ad import (  # noqa: E402
    Asset,
    AssetDownloader,
    AssetManifest,
    ConversionResult,
    DatasetBundleRecoveryError,
    JSON2H5ADConverter,
    PipelineRun,
    RawProcessingResult,
    SourcePlanner,
    json2h5ad,
)
from meta_standards_converter.converters.json_source import (
    DatasetPackageGroup,
    SourceLoadResult,
)
from meta_standards_converter.miniml import MINiMLCodec
from meta_standards_converter.retrieval import RetrievalPolicy


def package(*files, accession="GSM1"):
    supplementary_data = [{"value": path} for path in files]
    return {
        "miniml_schema_version": "2.0",
        "source": {"format": "test"},
        "series": {"accession": [{"value": "GSE1"}]},
        "sample": [
            {
                "iid": accession,
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


def test_converter_accepts_injected_dataset_combination_policy() -> None:
    class CombinationPolicy:
        @staticmethod
        def combine(adatas, **options):
            assert adatas == {"GSM1": "sample"}
            assert options["allow_unverified"] is True
            return "combined"

    converter = JSON2H5ADConverter(combination_policy=CombinationPolicy())

    assert converter._combine(
        {"GSM1": "sample"},
        allow_unverified=True,
    ) == "combined"


def atlas_v1_dataset(dataset_id, metadata):
    return {
        "dataset_id": dataset_id,
        "source_repository": "geo",
        "source_ordinal": 0,
        "status": "harmonized",
        "metadata": metadata,
        "publication_ids": [],
        "review": None,
        "harmonization": None,
        "diagnostics": [],
    }


def atlas_v1_payload(datasets):
    return {
        "schema_version": "1.0",
        "atlas": {"atlas_id": "atlas-test", "title": "Test", "theme": "test"},
        "run": {
            "run_id": "run-test",
            "created_at": "2026-07-31T00:00:00Z",
            "config": {"queries": [], "metadata_repositories": [], "options": {}},
            "status": "complete",
        },
        "datasets": datasets,
        "publications": [],
        "summary": {
            "dataset_count": len(datasets),
            "publication_count": 0,
            "completed_dataset_count": len(datasets),
            "failed_dataset_count": 0,
        },
    }


class TestSourcePlanner(unittest.TestCase):
    def test_indexes_assets_by_scope_without_changing_candidate_order(self):
        assets = [
            Asset("GSM2", "second.h5ad", "h5ad", source="json"),
            Asset("GSE1", "study.tsv", "matrix", source="json"),
            Asset("GSM2", "preferred.h5ad", "h5ad", source="manifest"),
        ]

        indexed = SourcePlanner()._index_assets_by_scope(assets)

        self.assertEqual([(0, assets[0]), (2, assets[2])], indexed["GSM2"])
        self.assertEqual([(1, assets[1])], indexed["GSE1"])

    def test_scope_index_preserves_global_tie_precedence(self):
        study = Asset("GSE1", "study.h5ad", "h5ad", source="json")
        sample = Asset("GSM1", "sample.h5ad", "h5ad", source="json")

        class FixedPlanner(SourcePlanner):
            def discover(self, packages):
                return [study, sample]

            def samples(self, packages):
                return ["GSM1"]

            def _study_by_sample(self, packages):
                return {"GSM1": "GSE1"}

        planned = FixedPlanner().plan([{}])

        self.assertEqual("study.h5ad", planned["GSM1"].path)
        self.assertEqual("GSE1", planned["GSM1"].study_scope)

    def test_groups_10x_matrix_barcode_and_gene_companions(self):
        data = package(
            "GSM1_brain.barcodes.tsv.gz",
            "GSM1_brain.genes.tsv.gz",
            "GSM1_brain.matrix.mtx.gz",
        )
        data["sample"][0]["sra_run"] = []

        plan = SourcePlanner().plan([data])

        self.assertEqual("10x_mtx", plan["GSM1"].role)
        self.assertEqual("GSM1_brain.matrix.mtx.gz", plan["GSM1"].path)
        self.assertEqual(
            "GSM1_brain.barcodes.tsv.gz", plan["GSM1"].barcodes_path
        )
        self.assertEqual("GSM1_brain.genes.tsv.gz", plan["GSM1"].features_path)

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


def test_convert_source_runs_each_atlas_study_independently(tmp_path):
    source = tmp_path / "atlas.json"
    source.write_text(
        json.dumps(
            atlas_v1_payload(
                [
                    atlas_v1_dataset(
                        "GSE1", package("one.h5ad", accession="GSM1")
                    ),
                    atlas_v1_dataset(
                        "GSE2",
                        {
                            **package("two.h5ad", accession="GSM2"),
                            "series": {"accession": [{"value": "GSE2"}]},
                        },
                    ),
                ]
            )
        ),
        encoding="utf-8",
    )
    converter = JSON2H5ADConverter()
    calls = []

    def convert_packages(packages, *, source_json, out, **options):
        study = packages[0]["series"]["accession"][0]["value"]
        calls.append((study, Path(out).name, Path(source_json).name))
        return ConversionResult(
            study_accession=study,
            sample_h5ads={f"GSM-{study}": f"{out}/sample.h5ad"},
        )

    converter._convert_packages = convert_packages

    result = converter.convert_source(str(source), out=str(tmp_path / "out"))

    assert ["GSE1", "GSE2"] == list(result.conversions)
    assert (
        [
            ("GSE1", "GSE1", "atlas.json"),
            ("GSE2", "GSE2", "atlas.json"),
        ]
        == calls
    )
    assert not result.partial


def test_convert_source_rejects_unsafe_dataset_id_without_escaping_output(tmp_path):
    class UnsafeSource:
        def load(self, _path):
            return SourceLoadResult(
                groups=(DatasetPackageGroup(
                    "../escape",
                    (MINiMLCodec().decode(package("one.h5ad")).package,),
                ),)
            )

    source = tmp_path / "input.json"
    source.write_text("{}", encoding="utf-8")
    converter = JSON2H5ADConverter(package_source=UnsafeSource())

    with pytest.raises(
        ValueError, match=r"^Unsafe dataset_id path component: '\.\./escape'$"
    ):
        converter.convert_source(str(source), out=str(tmp_path / "out"))

    assert not (tmp_path / "escape").exists()


def test_convert_source_persists_only_safe_group_failure_summaries(tmp_path):
    typed = MINiMLCodec().decode(package(), strict=True).package

    class FixedSource:
        def load(self, _path):
            return SourceLoadResult(
                groups=(
                    DatasetPackageGroup("GSE1", (typed,)),
                    DatasetPackageGroup("GSE2", (typed,)),
                )
            )

    class FailingConverter(JSON2H5ADConverter):
        def _convert_packages(self, *args, **kwargs):
            raise RuntimeError(
                "https://user:super-secret@example.org/data?token=super-secret"
            )

    source = tmp_path / "source.json"
    source.write_text("{}", encoding="utf-8")

    result = FailingConverter(package_source=FixedSource()).convert_source(
        str(source),
        out=str(tmp_path / "out"),
    )

    assert len(result.failures) == 2
    assert all("RuntimeError" in failure for failure in result.failures)
    assert all("correlation_id=" in failure for failure in result.failures)
    assert all("super-secret" not in failure for failure in result.failures)
    assert all("example.org" not in failure for failure in result.failures)

def test_convert_rejects_unsafe_single_dataset_id_before_conversion(tmp_path):
    class UnsafeSource:
        def load(self, _path):
            return SourceLoadResult(
                groups=(DatasetPackageGroup(
                    "../escape",
                    (MINiMLCodec().decode(package("one.h5ad")).package,),
                ),)
            )

    source = tmp_path / "input.json"
    source.write_text("{}", encoding="utf-8")
    converter = JSON2H5ADConverter(package_source=UnsafeSource())

    with pytest.raises(
        ValueError, match=r"^Unsafe dataset_id path component: '\.\./escape'$"
    ):
        converter.convert(str(source), out=str(tmp_path / "out"))

    assert not (tmp_path / "escape").exists()


@pytest.mark.parametrize("value", ["../escape", "..\\escape", ".", "..", "", "/tmp/x"])
def test_path_components_reject_unsafe_sample_and_study_ids(value):
    with pytest.raises(ValueError, match="Unsafe sample_id path component"):
        JSON2H5ADConverter._validate_path_component(value, "sample_id")


def test_dataset_bundle_restoration_failure_preserves_recovery_paths(
    tmp_path, monkeypatch
):
    output = tmp_path / "out"
    output.mkdir()
    old = output / "sample.h5ad"
    old.write_text("old", encoding="utf-8")
    staging = output / "staging"
    staging.mkdir()
    new = staging / "sample.h5ad"
    new.write_text("new", encoding="utf-8")
    real_replace = os.replace

    def fail_publish_and_restore(source, destination):
        source = Path(source)
        destination = Path(destination)
        if source == new:
            raise OSError("publish failed")
        if ".json2h5ad-recovery-" in source.parent.name and destination == old:
            raise OSError("restore failed")
        return real_replace(source, destination)

    monkeypatch.setattr(
        "meta_standards_converter.converters.json2h5ad.os.replace",
        fail_publish_and_restore,
    )

    with pytest.raises(DatasetBundleRecoveryError) as raised:
        JSON2H5ADConverter()._commit_dataset_bundle(
            [(new, old)], overwrite=True, staging=staging
        )

    assert "publish failed" in str(raised.value.publication_error)
    assert any("restore failed" in str(error) for error in raised.value.recovery_errors)
    assert raised.value.recovery_paths
    assert all(path.exists() for path in raised.value.recovery_paths)


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
            status_code = 200
            headers = {"Content-Length": "6"}

            def raise_for_status(self):
                return None

            def iter_content(self, chunk_size):
                return iter([b"abc", b"123"])

            def close(self):
                return None

        class Session:
            def __init__(self):
                self.calls = []

            def get(self, url, **kwargs):
                self.calls.append((url, kwargs))
                return Response()

        with tempfile.TemporaryDirectory() as tmpdir:
            session = Session()
            policy = RetrievalPolicy(
                allowed_hosts=frozenset({"example"}),
                resolver=lambda host, port, **kwargs: [
                    (2, 1, 6, "", ("93.184.216.34", port))
                ],
                disk_preflight=lambda *args, **kwargs: None,
            )
            downloader = AssetDownloader(
                cache_dir=tmpdir,
                session=session,
                policy=policy,
            )

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

    def test_result_exposes_primary_sample_catalogue_output(self):
        result = ConversionResult(
            study_accession="GSE1",
            combined_h5ad="out/GSE1.h5ad",
            sample_h5ads={"GSM1": "out/GSM1.h5ad"},
        )

        self.assertEqual("out/GSM1.h5ad", result.primary_h5ad)
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

            with self.assertRaisesRegex(
                ValueError, "JSON source contains no convertible package groups"
            ):
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

    def test_processed_conversion_resumes_from_atomic_sample_checkpoint(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            source_paths = []
            for sample_id in ("GSM1", "GSM2"):
                source_path = os.path.join(tmpdir, f"{sample_id}.h5ad")
                self.anndata.AnnData(
                    X=self.sparse.csr_matrix([[1]]),
                    obs=self.pandas.DataFrame(index=[f"cell-{sample_id}"]),
                    var=self.pandas.DataFrame(index=["ENSG1"]),
                ).write_h5ad(source_path)
                source_paths.append(source_path)
            data = package(source_paths[0], accession="GSM1")
            data["sample"].append(
                package(source_paths[1], accession="GSM2")["sample"][0]
            )
            json_path = self._write_json(tmpdir, data)
            checkpoint_dir = os.path.join(tmpdir, "checkpoints")

            class InterruptingConverter(JSON2H5ADConverter):
                calls = []

                def _read_processed_asset(self, asset, orientation="auto"):
                    self.calls.append(asset.scope_id)
                    if asset.scope_id == "GSM2":
                        raise RuntimeError("interrupted")
                    return super()._read_processed_asset(asset, orientation)

            with self.assertRaisesRegex(RuntimeError, "interrupted"):
                InterruptingConverter().convert(
                    json_path,
                    out=os.path.join(tmpdir, "first-out"),
                    resume=True,
                    processed_checkpoint_dir=checkpoint_dir,
                )
            self.assertEqual(1, len(list(Path(checkpoint_dir).rglob("*.h5ad"))))

            class CountingConverter(JSON2H5ADConverter):
                calls = []

                def _read_processed_asset(self, asset, orientation="auto"):
                    self.calls.append(asset.scope_id)
                    return super()._read_processed_asset(asset, orientation)

            result = CountingConverter().convert(
                json_path,
                out=os.path.join(tmpdir, "second-out"),
                resume=True,
                processed_checkpoint_dir=checkpoint_dir,
            )

            self.assertEqual(["GSM2"], CountingConverter.calls)
            self.assertEqual({"GSM1", "GSM2"}, set(result.sample_h5ads))
            self.assertIsNone(result.combined_h5ad)
            self.assertTrue(
                all(Path(path).is_file() for path in result.sample_h5ads.values())
            )

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
                    "characteristics": [{"name": "disease", "value": "healthy"}],
                }
            ]
            json_path = self._write_json(tmpdir, data)
            out = os.path.join(tmpdir, "out")

            result = json2h5ad().convert(json_path=json_path, out=out)

            normalized = self.anndata.read_h5ad(result.sample_h5ads["GSM1"])
            original = self.anndata.read_h5ad(source_path)
            self.assertEqual(["cell1-GSM1", "cell2-GSM1"], list(normalized.obs_names))
            self.assertEqual(["cell1", "cell2"], list(normalized.obs["msc.observation.original_id"]))
            self.assertEqual(["GSM1", "GSM1"], list(normalized.obs["msc.sample.accession"]))
            self.assertEqual(["Control sample"] * 2, list(normalized.obs["msc.sample.title"]))
            self.assertEqual(
                ["Homo sapiens"] * 2,
                list(normalized.obs["msc.sample.channel.organism.value"]),
            )
            self.assertEqual(["healthy"] * 2, list(normalized.obs["msc.sample.channel.disease"]))
            self.assertFalse(
                any(column.startswith("msc_") for column in normalized.obs),
                normalized.obs.columns,
            )
            self.assertEqual([True, False], normalized.obs["pipeline_qc_pass"].tolist())
            self.assertEqual(
                ["healthy"] * 2,
                normalized.obs["msc.characteristics.disease"].tolist(),
            )
            self.assertFalse(any(column.startswith("geo_") for column in normalized.obs))
            self.assertEqual(["cell1", "cell2"], list(original.obs_names))
            self.assertEqual("h5ad", normalized.uns["meta_standards_converter"]["source_tier"])
            self.assertEqual(
                "1.0",
                normalized.uns["meta_standards_converter"]["metadata_schema_version"],
            )
            self.assertEqual("artifact_parent", normalized.uns["meta_standards_converter"]["path_base"])
            self.assertEqual("../source.h5ad", normalized.uns["meta_standards_converter"]["source_uri"])
            self.assertEqual(["../source.h5ad"] * 2, list(normalized.obs["msc.asset.uri"]))

            sample_values = normalized.uns["msc_metadata"]["sample_values"]
            self.assertEqual("1.0", normalized.uns["msc_metadata"]["schema_version"])
            self.assertEqual(
                ["sample_accession", "field", "ordinal", "value", "value_type"],
                list(sample_values.columns),
            )
            self.assertEqual(
                ["healthy"],
                sample_values.loc[
                    sample_values["field"] == "msc.sample.channel.disease", "value"
                ].tolist(),
            )

            manifest = json.loads(Path(result.manifest_path).read_text())
            self.assertEqual("artifact_parent", manifest["path_base"])
            self.assertEqual("../GSE1.json", manifest["source_json"])
            self.assertIsNone(manifest["combined_h5ad"])
            self.assertEqual("per_sample_h5ad_catalogue", manifest["artifact_kind"])
            self.assertEqual("none", manifest["expression_integration"])
            self.assertEqual("1.0", manifest["h5ad_metadata_schema_version"])
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
            self.assertEqual(["GSM1"], list(sample.obs_names))
            self.assertIsNone(result.combined_h5ad)

    def test_preserves_delimiter_qualified_ids_and_qualifies_only_unqualified_ids(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            source_path = os.path.join(tmpdir, "source.h5ad")
            self.anndata.AnnData(
                X=self.sparse.csr_matrix([[1], [2], [3], [4]]),
                obs=self.pandas.DataFrame(
                    index=["GSM1-cell-1", "cell-2_GSM1", "notGSM10-cell", "barcode"]
                ),
                var=self.pandas.DataFrame(index=["ENSG1"]),
            ).write_h5ad(source_path)
            json_path = self._write_json(tmpdir, package(source_path))

            result = json2h5ad().convert(
                json_path=json_path, out=os.path.join(tmpdir, "out")
            )

            sample = self.anndata.read_h5ad(result.sample_h5ads["GSM1"])
            expected = [
                "GSM1-cell-1",
                "cell-2_GSM1",
                "notGSM10-cell-GSM1",
                "barcode-GSM1",
            ]
            self.assertEqual(expected, list(sample.obs_names))
            self.assertEqual(
                ["GSM1-cell-1", "cell-2_GSM1", "notGSM10-cell", "barcode"],
                sample.obs["msc.observation.original_id"].tolist(),
            )

    def test_resolves_duplicate_observation_ids_before_catalogue_writes(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = []
            for sample_id in ("GSM1", "GSM2"):
                source_path = os.path.join(tmpdir, f"{sample_id}.h5ad")
                adata = self.anndata.AnnData(
                    X=self.sparse.csr_matrix([[1], [2]]),
                    obs=self.pandas.DataFrame(index=["barcode", "barcode"]),
                    var=self.pandas.DataFrame(index=["ENSG1"]),
                )
                adata.write_h5ad(source_path)
                paths.append(source_path)
            data = package(paths[0], accession="GSM1")
            data["sample"].append(package(paths[1], accession="GSM2")["sample"][0])
            json_path = self._write_json(tmpdir, data)

            result = json2h5ad().convert(
                json_path=json_path, out=os.path.join(tmpdir, "out")
            )

            expected = {
                "GSM1": ["barcode-GSM1-1", "barcode-GSM1-2"],
                "GSM2": ["barcode-GSM2-1", "barcode-GSM2-2"],
            }
            for sample_id, names in expected.items():
                sample = self.anndata.read_h5ad(result.sample_h5ads[sample_id])
                self.assertEqual(names, list(sample.obs_names))
                self.assertEqual(["barcode", "barcode"], sample.obs["msc.observation.original_id"].tolist())
            self.assertIsNone(result.combined_h5ad)

    def test_preserves_legacy_underscore_columns_as_opaque_source_metadata(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            source_path = os.path.join(tmpdir, "source.h5ad")
            self.anndata.AnnData(
                X=self.sparse.csr_matrix([[1]]),
                obs=self.pandas.DataFrame(
                    {"msc_accession": ["opaque-source-value"]}, index=["cell"]
                ),
                var=self.pandas.DataFrame(index=["ENSG1"]),
            ).write_h5ad(source_path)
            json_path = self._write_json(tmpdir, package(source_path))

            result = json2h5ad().convert(
                json_path=json_path, out=os.path.join(tmpdir, "out")
            )

            sample = self.anndata.read_h5ad(result.sample_h5ads["GSM1"])
            self.assertEqual(["opaque-source-value"], sample.obs["msc_accession"].tolist())
            self.assertEqual(["GSM1"], sample.obs["msc.sample.accession"].tolist())

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
                    "organism": [{
                        "taxid": "9606", "value": "human",
                        "annotations": [{
                            "field": "organism", "value": "Homo sapiens",
                            "term_source_ref": "ncbitaxon",
                            "term_accession_number": "NCBITaxon_9606",
                        }],
                    }],
                },
                {
                    "organism": [{"taxid": "10090", "value": "Mus musculus"}],
                },
                {"organism": [{"taxid": "10090", "value": "mus musculus"}]},
            ]
            json_path = self._write_json(tmpdir, data)

            result = json2h5ad().convert(json_path=json_path, out=os.path.join(tmpdir, "out"))

            converted = self.anndata.read_h5ad(result.sample_h5ads["GSM1"])
            self.assertEqual(
                ["Homo sapiens; Mus musculus"],
                converted.obs["msc.sample.channel.organism.value"].unique().tolist(),
            )
            self.assertEqual(
                converted.obs["msc.sample.channel.organism.value"].tolist(),
                converted.obs["msc.sample.channel.organism.value"].tolist(),
            )
            self.assertEqual(
                ["9606; 10090"],
                converted.obs["msc.sample.channel.organism.taxid"].unique().tolist(),
            )
            fields = converted.uns["msc_miniml"]["fields"]
            self.assertIn("channel[0].organism[0].annotations[0].value", set(fields["path"]))

    def test_derives_organism_from_scalar_harmonization_and_preserves_missing_as_empty(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            for case, channel, expected in (
                (
                    "scalar",
                    {
                        "organism": [{
                            "value": "human",
                            "annotations": [{
                                "field": "organism", "value": "Homo sapiens",
                                "term_source_ref": "ncbitaxon",
                                "term_accession_number": "NCBITaxon_9606",
                            }],
                        }],
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
                    self.assertEqual(
                        [expected],
                        converted.obs["msc.sample.channel.organism.value"].unique().tolist(),
                    )
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

    def test_publishes_compatible_samples_as_catalogue_without_expression_join(self):
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

            self.assertIsNone(result.combined_h5ad)
            self.assertEqual({"GSM1", "GSM2"}, set(result.sample_h5ads))
            manifest = json.loads(Path(result.manifest_path).read_text())
            self.assertEqual("per_sample_h5ad_catalogue", manifest["artifact_kind"])
            self.assertEqual("none", manifest["expression_integration"])
            self.assertEqual(2, manifest["sample_count"])
            self.assertIsNone(manifest["combined_h5ad"])
            first_sample = self.anndata.read_h5ad(result.sample_h5ads["GSM1"])
            second_sample = self.anndata.read_h5ad(result.sample_h5ads["GSM2"])
            self.assertEqual((1, 2), first_sample.shape)
            self.assertEqual((1, 2), second_sample.shape)

    def test_all_unknown_samples_are_not_falsely_compatibility_verified(self):
        adatas = {
            sample_id: self.anndata.AnnData(
                X=self.sparse.csr_matrix([[1]]),
                obs=self.pandas.DataFrame(index=[f"{sample_id}-cell"]),
                var=self.pandas.DataFrame(index=["feature:1"]),
            )
            for sample_id in ("GSM1", "GSM2")
        }

        missing = JSON2H5ADConverter()._missing_combination_evidence(adatas)

        self.assertEqual(
            {
                "organism": ["GSM1", "GSM2"],
                "reference": ["GSM1", "GSM2"],
                "modality": ["GSM1", "GSM2"],
                "feature_namespace": ["GSM1", "GSM2"],
            },
            missing,
        )

    def test_feature_namespace_distinguishes_entrez_ids_from_gene_symbols(self):
        converter = JSON2H5ADConverter()
        entrez = self.anndata.AnnData(
            X=self.sparse.csr_matrix([[1, 2, 3]]),
            obs=self.pandas.DataFrame(index=["cell"]),
            var=self.pandas.DataFrame(index=["7157", "1956", "7422"]),
        )
        symbols = self.anndata.AnnData(
            X=self.sparse.csr_matrix([[1, 2, 3]]),
            obs=self.pandas.DataFrame(index=["cell"]),
            var=self.pandas.DataFrame(index=["TP53", "EGFR", "VEGFA"]),
        )

        self.assertEqual("entrez", converter._feature_namespace(entrez))
        self.assertEqual("symbol", converter._feature_namespace(symbols))

    def test_memory_preflight_skips_then_force_resume_bypasses_only_fixed_profile(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            source_path = os.path.join(tmpdir, "source.h5ad")
            self.anndata.AnnData(
                X=self.sparse.csr_matrix([[1]]),
                obs=self.pandas.DataFrame(index=["cell"]),
                var=self.pandas.DataFrame(index=["ENSG1"]),
            ).write_h5ad(source_path)
            json_path = self._write_json(tmpdir, package(source_path))
            out = os.path.join(tmpdir, "out")
            converter = JSON2H5ADConverter(
                resource_overrides={
                    "max_in_memory_matrix_bytes": 100,
                    "available_memory_fraction": 0.70,
                    "force_memory_fraction": 0.90,
                },
                available_memory=lambda: 1_000,
                memory_estimator=lambda _path, _asset: 800,
            )

            skipped = converter.convert(json_path=json_path, out=out)

            self.assertEqual({}, skipped.sample_h5ads)
            self.assertTrue(skipped.partial)
            self.assertEqual("skipped", skipped.memory_report[0]["decision"])
            self.assertEqual(100, skipped.memory_report[0]["limit_bytes"])
            first_manifest = json.loads(Path(skipped.manifest_path).read_text())
            self.assertEqual(skipped.memory_report, first_manifest["memory_report"])

            resumed = converter.convert(
                json_path=json_path,
                out=out,
                resume=True,
                force_memory=True,
                overwrite=True,
            )

            self.assertEqual({"GSM1"}, set(resumed.sample_h5ads))
            self.assertEqual("admitted", resumed.memory_report[0]["decision"])
            self.assertEqual(900, resumed.memory_report[0]["limit_bytes"])
            self.assertTrue(Path(resumed.sample_h5ads["GSM1"]).is_file())
            self.assertTrue(any((Path(out) / ".processed").rglob("*.h5ad")))

    def test_force_memory_never_exceeds_ninety_percent_available(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            source_path = os.path.join(tmpdir, "source.h5ad")
            self.anndata.AnnData(
                X=self.sparse.csr_matrix([[1]]),
                obs=self.pandas.DataFrame(index=["cell"]),
                var=self.pandas.DataFrame(index=["ENSG1"]),
            ).write_h5ad(source_path)
            json_path = self._write_json(tmpdir, package(source_path))
            converter = JSON2H5ADConverter(
                resource_overrides={"max_in_memory_matrix_bytes": 100},
                available_memory=lambda: 1_000,
                memory_estimator=lambda _path, _asset: 901,
            )

            result = converter.convert(
                json_path=json_path,
                out=os.path.join(tmpdir, "out"),
                resume=True,
                force_memory=True,
            )

            self.assertEqual({}, result.sample_h5ads)
            self.assertEqual("skipped", result.memory_report[0]["decision"])
            self.assertEqual(900, result.memory_report[0]["limit_bytes"])

    def test_force_memory_requires_resume(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            json_path = self._write_json(tmpdir, package("source.h5ad"))

            with self.assertRaisesRegex(ValueError, "force_memory requires resume"):
                JSON2H5ADConverter().convert(
                    json_path=json_path,
                    out=os.path.join(tmpdir, "out"),
                    force_memory=True,
                )

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
                "miniml_schema_version": "2.0",
                "source": {
                    "format": "GEO MINiML",
                    "version": "1.0",
                    "schema_location": "MINiML.xsd",
                },
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
                                    {
                                        "name": "cell type", "value": "Treg; memory",
                                        "annotations": [{
                                            "field": "cell_type", "value": "regulatory T cell",
                                            "term_source_ref": "cl",
                                            "term_accession_number": "CL:0000815",
                                        }],
                                    },
                                    {"name": "developmental stage", "value": "adult"},
                                    {"name": "treatment", "value": "CPI-703"},
                                ],
                                "treatment_protocol": "Long treatment protocol",
                            },
                            {
                                "source": "blood",
                                "characteristics": [
                                    {"name": "cell-type", "value": "Activated Treg"},
                                    {"name": "treatment", "value": "CPI-703"},
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
                        "library_strategy": "RNA-Seq",
                        "channel": [{
                            "organism": [{"taxid": "9606", "value": "Homo sapiens"}],
                            "characteristics": [{"name": "dose", "value": "5 uM"}],
                        }],
                    },
                ],
            }
            json_path = self._write_json(tmpdir, data)

            result = json2h5ad().convert(json_path=json_path, out=os.path.join(tmpdir, "out"))

            first = self.anndata.read_h5ad(result.sample_h5ads["GSM1"])
            second = self.anndata.read_h5ad(result.sample_h5ads["GSM2"])
            self.assertEqual(["GEO"], first.obs["msc.database.identifier"].unique().tolist())
            self.assertEqual(
                ["Gene Expression Omnibus (GEO)"],
                first.obs["msc.database.name"].unique().tolist(),
            )
            self.assertEqual(["GPL1"], first.obs["msc.platform.accession"].unique().tolist())
            self.assertEqual(["GPL2"], second.obs["msc.platform.accession"].unique().tolist())
            self.assertEqual(["transcriptomic"], first.obs["msc.library.source"].unique().tolist())
            self.assertEqual(["SRR1; SRR2"], first.obs["msc.archive.sra_run_accessions"].unique().tolist())
            self.assertEqual(["SAMN1"], first.obs["msc.archive.biosample_accession"].unique().tolist())
            self.assertEqual(["Illumina Test"], first.obs["msc.instrument.model"].unique().tolist())
            self.assertEqual(["adult"], first.obs["msc.sample.channel.developmental_stage"].unique().tolist())
            self.assertEqual(["Example Biobank"], first.obs["msc.sample.channel.biomaterial_provider"].unique().tolist())
            self.assertEqual(["RNA"], first.obs["msc.sample.channel.material_type"].unique().tolist())
            self.assertEqual(
                ["sample treatment protocol"],
                first.obs["msc.protocol.types"].unique().tolist(),
            )
            self.assertEqual(["EFO"], first.obs["msc.protocol.term_source_refs"].unique().tolist())
            self.assertEqual(
                ["EFO_0003809"],
                first.obs["msc.protocol.term_accession_numbers"].unique().tolist(),
            )
            self.assertEqual(
                ["Treg; memory; Activated Treg"],
                first.obs["msc.characteristics.cell_type"].unique().tolist(),
            )
            self.assertEqual(["CPI-703"], first.obs["msc.characteristics.treatment"].unique().tolist())
            self.assertEqual(
                ["regulatory T cell"],
                first.obs["msc.characteristics.harmonized_cell_type"].unique().tolist(),
            )
            self.assertEqual(
                ["CL:0000815"],
                first.obs["msc.characteristics.harmonized_cell_type_id"].unique().tolist(),
            )
            self.assertEqual(
                ["cl"],
                first.obs["msc.characteristics.harmonized_cell_type_onto"].unique().tolist(),
            )
            self.assertEqual([""], first.obs["msc.characteristics.dose"].unique().tolist())
            self.assertEqual([""], second.obs["msc.characteristics.cell_type"].unique().tolist())
            self.assertEqual([""], second.obs["msc.characteristics.harmonized_cell_type"].unique().tolist())
            self.assertIsNone(result.combined_h5ad)

            values = first.uns["msc_metadata"]["sample_values"]
            cell_types = values.loc[
                values["field"] == "msc.characteristics.cell_type"
            ]
            self.assertEqual([0, 1], cell_types["ordinal"].tolist())
            self.assertEqual(
                ["Treg; memory", "Activated Treg"], cell_types["value"].tolist()
            )
            self.assertEqual(["string", "string"], cell_types["value_type"].tolist())
            self.assertNotIn(
                "msc.sample.channel.disease", set(values["field"])
            )

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

            second_fields = second.uns["msc_miniml"]["fields"]
            self.assertEqual(
                {"GSM2"},
                set(second_fields.loc[second_fields["entity_type"] == "sample", "entity_id"]),
            )
            self.assertEqual(
                {"C1", "C2"},
                set(second_fields.loc[second_fields["entity_type"] == "contributor", "entity_id"]),
            )
            self.assertEqual(
                {"GPL2"},
                set(second_fields.loc[second_fields["entity_type"] == "platform", "entity_id"]),
            )

    def test_catalogue_keeps_samples_with_different_organisms_without_integration_claim(self):
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
            self.assertFalse(result.partial)
            self.assertEqual([], result.failures)

    def test_matching_harmonized_organisms_remain_separate_catalogue_members(self):
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
                        "organism": [{
                            "value": raw_organism,
                            "annotations": [{
                                "field": "organism", "value": "Homo sapiens",
                                "term_source_ref": "ncbitaxon",
                                "term_accession_number": "NCBITaxon_9606",
                            }],
                        }],
                    }
                ]
                samples.append(sample)
            data = package()
            data["sample"] = samples
            json_path = self._write_json(tmpdir, data)

            result = json2h5ad().convert(json_path=json_path, out=os.path.join(tmpdir, "out"))

            self.assertIsNone(result.combined_h5ad)
            self.assertEqual([], result.failures)
            for path in result.sample_h5ads.values():
                sample = self.anndata.read_h5ad(path)
                self.assertEqual(
                    ["Homo sapiens"],
                    sample.obs["msc.sample.channel.organism.value"].unique().tolist(),
                )

    def test_catalogue_keeps_samples_with_different_declared_references(self):
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
            self.assertEqual([], result.failures)

    def test_catalogue_never_requires_combination_acknowledgement(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            samples = []
            for sample_id, genome in (("GSM1", "GRCh38"), ("GSM2", None)):
                path = os.path.join(tmpdir, f"{sample_id}.h5ad")
                adata = self.anndata.AnnData(
                    X=self.sparse.csr_matrix([[1]]),
                    obs=self.pandas.DataFrame(index=["cell"]),
                    var=self.pandas.DataFrame(
                        {"gene_ids": ["ENSG1"]}, index=["ENSG1"]
                    ),
                )
                if genome is not None:
                    adata.uns["genome"] = genome
                adata.write_h5ad(path)
                sample = package(path, accession=sample_id)["sample"][0]
                sample["library_source"] = "single cell transcriptomic"
                sample["channel"] = [
                    {"organism": [{"value": "Homo sapiens"}]}
                ]
                samples.append(sample)
            data = package()
            data["sample"] = samples
            json_path = self._write_json(tmpdir, data)

            strict = json2h5ad().convert(
                json_path=json_path,
                out=os.path.join(tmpdir, "strict"),
            )
            acknowledged = json2h5ad().convert(
                json_path=json_path,
                out=os.path.join(tmpdir, "acknowledged"),
                allow_unverified_combination=True,
            )

            self.assertIsNone(strict.combined_h5ad)
            self.assertEqual([], strict.failures)
            self.assertIsNone(acknowledged.combined_h5ad)
            self.assertFalse(acknowledged.partial)
            self.assertIn("ignored", acknowledged.warnings[0])

    def test_catalogue_keeps_bulk_and_single_cell_modalities_separate(self):
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
            self.assertEqual([], result.failures)

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

    def test_study_h5ad_is_split_by_canonical_and_generic_sample_accession(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            for accession_column in ("msc.sample.accession", "geo_accession"):
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
