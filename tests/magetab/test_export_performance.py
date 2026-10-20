"""Export telemetry and work reuse must preserve the actual MAGE-TAB contract."""
import logging

import pytest

from meta_standards_converter.conversion_logging import operation
from meta_standards_converter.converters.json2ae import JSON2AEConverter
from meta_standards_converter.magetab.constructor import AEConstructor
from meta_standards_converter.miniml import MINiMLCodec
from meta_standards_converter.sources.json import DatasetPackageGroup, SourceLoadResult
from tests.test_archive_export_cleanup import legacy
from tests.test_archive_publications import NoRuns, PubMed
from tests.test_conversion_logging import events
from meta_standards_converter.miniml import MINiMLPackage
from copy import deepcopy
from dataclasses import replace


def test_export_stages_report_work_cpu_and_process_peak_without_changing_output(caplog, tmp_path):
    package = MINiMLCodec().decode(legacy()).package
    source = SourceLoadResult((DatasetPackageGroup("study", (package,)),))
    converter = JSON2AEConverter(
        ae_constructor=AEConstructor(pubmed_client=PubMed(), insdc_client=NoRuns())
    )
    expected = converter.convert_loaded(source, enrich=False)
    with caplog.at_level(logging.DEBUG):
        result = converter.convert_loaded(source, enrich=False, out=str(tmp_path))
    assert result == expected
    phases = [e for e in events(caplog)
              if e["event_type"] == "stage_finished" and
              str(e["stage"]).startswith("magetab_")]
    names = {e["stage"] for e in phases}
    assert {
        "magetab_package_resolution", "magetab_model_reconstruction",
        "magetab_normalization", "magetab_protocol_matching",
        "magetab_file_projection", "magetab_path_construction",
        "magetab_column_planning", "magetab_row_rendering",
        "magetab_semantic_overlay", "magetab_validation", "magetab_writing",
    } <= names
    for event in phases:
        metrics = event["attributes"]
        assert metrics["wall_seconds"] >= metrics["exclusive_instrumented_wall_seconds"] >= 0
        assert metrics["process_cpu_seconds"] >= metrics["exclusive_instrumented_cpu_seconds"] >= 0
        assert metrics["peak_rss_scope"] == "process_lifetime_high_water"
        assert metrics["peak_rss_bytes"] is None or metrics["peak_rss_bytes"] > 0
    assert len(list(tmp_path.glob("*.idf.txt"))) == 1
    assert len(list(tmp_path.glob("*.sdrf.txt"))) == 1


def test_export_phase_excludes_nested_instrumented_work_and_restores_failure_context(caplog, monkeypatch):
    from meta_standards_converter.magetab import _performance

    clock = [0.0]
    cpu = [0.0]
    monkeypatch.setattr(_performance.time, "perf_counter", lambda: clock[0])
    monkeypatch.setattr(_performance.time, "process_time", lambda: cpu[0])
    with caplog.at_level(logging.DEBUG), operation("test", heartbeat=False):
        with _performance.phase("outer"):
            clock[0], cpu[0] = 2, 1
            with pytest.raises(ValueError), _performance.phase("inner"):
                clock[0], cpu[0] = 5, 3
                raise ValueError("failed projection")
            clock[0], cpu[0] = 9, 5
        with _performance.phase("next"):
            clock[0], cpu[0] = 10, 6
    completed = {e["stage"]: e for e in events(caplog) if e["event_type"] == "stage_finished"}
    assert completed["inner"]["status"] == "failed"
    assert completed["outer"]["attributes"]["wall_seconds"] == 9
    assert completed["outer"]["attributes"]["exclusive_instrumented_wall_seconds"] == 6
    assert completed["outer"]["attributes"]["exclusive_instrumented_cpu_seconds"] == 3
    assert completed["next"]["attributes"]["exclusive_instrumented_wall_seconds"] == 1


def count_package_construction(monkeypatch):
    calls = []
    original = MINiMLPackage.from_mapping.__func__

    def construct(cls, value):
        calls.append(None)
        return original(cls, value)

    monkeypatch.setattr(MINiMLPackage, "from_mapping", classmethod(construct))
    return calls


def test_constructor_validates_once_and_owns_one_detached_export_view(monkeypatch):
    package = MINiMLCodec().decode(legacy()).package
    before = package.to_mapping()
    constructor = AEConstructor(pubmed_client=PubMed(), insdc_client=NoRuns())
    expected = constructor.miniml2magetab(package)
    calls = count_package_construction(monkeypatch)
    assert constructor.miniml2magetab(package) == expected
    assert package.to_mapping() == before
    assert len(calls) == 1


def test_default_resolution_validates_once_without_disabled_profile_roundtrip(monkeypatch):
    package = MINiMLCodec().decode(legacy()).package
    source = SourceLoadResult((DatasetPackageGroup("study", (package,)),))
    converter = JSON2AEConverter()
    calls = count_package_construction(monkeypatch)
    actual = converter._packages_from_loaded(source)
    assert actual == [package]
    assert actual[0] is not package  # The external typed input is revalidated.
    assert len(calls) == 1


def test_study_identity_lookup_does_not_materialize_the_whole_package(monkeypatch):
    package = MINiMLCodec().decode(legacy()).package
    expected = JSON2AEConverter()._usable_study_accession(package)

    def forbidden(*args):
        raise AssertionError("Identity lookup traversed all samples and source evidence")

    monkeypatch.setattr(MINiMLPackage, "to_mapping", forbidden)
    assert JSON2AEConverter()._usable_study_accession(package) == expected


def test_changed_external_typed_input_is_not_trusted_by_export():
    package = MINiMLCodec().decode(legacy()).package
    malformed = replace(package, samples=package.samples + package.samples)
    with pytest.raises(ValueError, match="duplicate sample iid"):
        AEConstructor(pubmed_client=PubMed(), insdc_client=NoRuns()).miniml2magetab(malformed)
    with pytest.raises(ValueError, match="duplicate sample iid"):
        JSON2AEConverter()._packages_from_loaded(
            SourceLoadResult((DatasetPackageGroup("bad", (malformed,)),)))


def test_custom_group_resolution_and_injected_constructor_keep_their_contract():
    package = MINiMLCodec().decode(legacy()).package
    supplied = []

    class CustomGroup(DatasetPackageGroup):
        def resolved(self, *, replacement_profile=None):
            supplied.append(replacement_profile)
            return DatasetPackageGroup("changed", (replace(package, series=replace(
                package.series, title="injected resolution")),))

    class Constructor:
        def miniml2magetab(self, data):
            assert isinstance(data, MINiMLPackage)
            return data.series.title

    result = JSON2AEConverter(ae_constructor=Constructor()).convert_loaded(
        SourceLoadResult((CustomGroup("original", (package,)),)), enrich=False)
    assert result == ["injected resolution"]
    assert supplied == [None]


def test_private_view_is_discarded_after_failure_and_next_export_is_clean(monkeypatch):
    package = MINiMLCodec().decode(legacy()).package
    before = deepcopy(package.to_mapping())
    constructor = AEConstructor(pubmed_client=PubMed(), insdc_client=NoRuns())
    expected = constructor.miniml2magetab(package)
    from meta_standards_converter.magetab import native_files
    original = native_files.project_native_files

    def fail(data):
        data["sample"].clear()
        raise RuntimeError("projection failure")

    with monkeypatch.context() as scoped:
        scoped.setattr(native_files, "project_native_files", fail)
        with pytest.raises(RuntimeError, match="projection failure"):
            constructor.miniml2magetab(package)
    assert package.to_mapping() == before
    assert native_files.project_native_files is original
    assert constructor.miniml2magetab(package) == expected


@pytest.mark.parametrize("forced", [None, "bulk_sequencing", "generic"])
def test_canonical_table_renders_once_and_retains_preparation_and_audit(monkeypatch, forced):
    from meta_standards_converter.magetab.sdrf.constructor import SDRFConstructor
    from meta_standards_converter.magetab.sdrf.renderer import SDRFRenderer

    class InjectedConstructor(SDRFConstructor):
        builds = 0

        def build(self, handler):
            self.builds += 1
            return super().build(handler)

    package = MINiMLCodec().decode(legacy()).package
    injected = InjectedConstructor()
    reference = AEConstructor(sdrf_constructor=injected, pubmed_client=PubMed(), insdc_client=NoRuns())
    expected = reference.miniml2magetab(package, platform_handler=forced)
    assert injected.builds == 1
    calls = []
    original = SDRFRenderer.render_paths

    def render(self, columns, paths):
        calls.append((len(columns), len(paths)))
        return original(self, columns, paths)

    monkeypatch.setattr(SDRFRenderer, "render_paths", render)
    candidate = AEConstructor(pubmed_client=PubMed(), insdc_client=NoRuns())
    assert candidate.miniml2magetab(package, platform_handler=forced) == expected
    assert candidate.sdrf_constructor.last_sdrf_audit == injected.last_sdrf_audit
    assert len(calls) == 1


def test_inputs_without_canonical_paths_keep_fallback_rendering(monkeypatch):
    from meta_standards_converter.magetab.sdrf.renderer import SDRFRenderer
    package = MINiMLCodec().decode({
        "miniml_schema_version": "3.0", "source": {"format": "test"},
        "series": {"iid": "GSE1", "title": "Sparse", "sample_ref": [{"ref": "S1"}]},
        "sample": [{"iid": "S1", "title": "Sample", "channel": [{"source": "source"}]}],
    }).package
    calls = []
    original = SDRFRenderer.render_paths

    def render(self, columns, paths):
        calls.append(len(paths))
        return original(self, columns, paths)

    monkeypatch.setattr(SDRFRenderer, "render_paths", render)
    constructor = AEConstructor(pubmed_client=PubMed(), insdc_client=NoRuns())
    rows = constructor.miniml2magetab(package)
    table = next(r[1] for r in rows if r[0] == "SDRF File")
    assert table[0][0] == "Source Name" and len(table) == 2
    assert len(calls) == 1
