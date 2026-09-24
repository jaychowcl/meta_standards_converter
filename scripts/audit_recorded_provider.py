# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================

"""Audit explicit recorded cohorts. No live fallback or golden regeneration.

Run with the existing environment: python scripts/audit_recorded_provider.py
tests/fixtures/recorded_provider/success.manifest.json. External corpora require
--evidence-root; relative file IDs always resolve against that one root.
"""
import argparse
import csv
from contextlib import contextmanager, ExitStack
from dataclasses import asdict, is_dataclass
from datetime import date
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile
import time

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path: sys.path.insert(0, str(REPO))
from tests.support.recorded_replay import HarnessError, ReplayTransport, ReplayBudgetExceeded, assert_occurrences
from tests.support.recorded_replay import occurrence_counter


def primitive(value):
    if hasattr(value, "to_mapping"): return primitive(value.to_mapping())
    if is_dataclass(value): return primitive(asdict(value))
    if isinstance(value, dict): return {key: primitive(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)): return [primitive(item) for item in value]
    return value


def pointer(value, path):
    for key in path.lstrip("/").split("/"):
        key = key.replace("~1", "/").replace("~0", "~")
        value = value[int(key)] if isinstance(value, list) else value[key]
    return value


def first_difference(actual, expected, path=""):
    """Bounded mismatch location; never dumps multi-megabyte golden payloads."""
    if type(actual) is not type(expected): return path or "/", actual, expected
    if isinstance(actual, dict):
        if actual.keys() != expected.keys(): return path + "/<keys>", sorted(actual), sorted(expected)
        for key in actual:
            found = first_difference(actual[key], expected[key], path + "/" + str(key).replace("~", "~0").replace("/", "~1"))
            if found: return found
    elif isinstance(actual, list):
        if len(actual) != len(expected): return path + "/<length>", len(actual), len(expected)
        for index, (left, right) in enumerate(zip(actual, expected)):
            found = first_difference(left, right, path + "/" + str(index))
            if found: return found
    elif actual != expected: return path or "/", actual, expected
    return None


def compare_golden(payload, expected):
    kind = expected["kind"]
    if kind == "pointers":
        for path, value in expected["pointers"].items():
            assert pointer(payload, path) == value, f"scientific value mismatch at {path}"
    elif kind == "sdrf-columns":
        rows = []
        for tables in payload:
            table = next(row[1] for row in tables if row[0] == "SDRF File")
            # Repeated columns must be addressed explicitly by ordinal.
            indices = []
            for column in expected["columns"]:
                name, ordinal = (column, 1) if isinstance(column, str) else column
                matches = [i for i, header in enumerate(table[0]) if header == name]
                indices.append(matches[ordinal - 1])
            rows.extend([[row[i] for i in indices] for row in table[1:]])
        assert_occurrences(expected["rows"], rows)
    elif kind == "tsv-columns":
        rows = [{column: row[column] for column in expected["columns"]}
                for table in payload.values() for row in table["rows"]]
        assert_occurrences(expected["rows"], rows)
    elif kind == "exact":
        mismatch = first_difference(payload, expected["value"])
        if mismatch:
            path, actual, wanted = mismatch
            raise AssertionError(f"exact scientific payload mismatch at {path}: actual={str(actual)[:240]!r}; expected={str(wanted)[:240]!r}")
    else:
        raise HarnessError(f"unknown golden comparator: {kind}")


def normalize_runtime_date(payload, expected, declaration):
    """Only declared generated IDF dates; every scientific cell is untouched."""
    if declaration is None: return payload, expected
    if set(declaration) != {"recorded_submission_date"}: raise HarnessError("unsupported golden normalization")
    def normalized(value, required):
        result = []
        for table in value:
            found = 0
            rows = []
            for row in table:
                if row and row[0] == "Comment[ArrayExpressSubmissionDate]":
                    assert row[1:] == [required], "unexpected generated submission date"
                    rows.append([row[0], "<AUDIT_RUN_DATE>"])
                    found += 1
                else: rows.append(row)
            assert found == 1, "missing or repeated generated submission date in table"
            result.append(rows)
        return result
    return normalized(payload, date.today().isoformat()), normalized(expected, declaration["recorded_submission_date"])


def build_services(replay, declarations):
    """Explicit parser/merge adapters; these never claim SRA HTTP discovery."""
    services = {}
    for name, spec in declarations.items():
        if name == "sra2json" and spec["kind"] == "frozen-source":
            if replay.manifest["evidence_kind"] != "parser-merge-replay":
                raise HarnessError("FrozenSource must be labelled parser-merge-replay")
            from meta_standards_converter.converters.sra2json import SRA2JSONConverter
            from meta_standards_converter.sources.archive_support import Resolution, StudySeed, StudyRecords
            from meta_standards_converter.xml_safety import parse_xml
            def snapshot(file_id, records=False):
                try:
                    value = json.loads(replay.read(file_id))
                    if value["schema_version"] != "1.0": raise ValueError("unknown snapshot version")
                    if not records:
                        return Resolution([StudySeed(**seed) for seed in value["studies"]], value["issues"])
                    return StudyRecords(seed=StudySeed(**value["seed"]),
                        xml=[parse_xml(xml, max_bytes=64 * 1024 * 1024) for xml in value["xml"]],
                        indexed=value["indexed"], linked=value["linked"], issues=value["issues"])
                except Exception as exc:
                    replay.fail(f"invalid manifested snapshot {file_id}: {exc}")
            class FrozenSource:
                def __init__(self, declaration): self.declaration = declaration
                def resolve(self, accession):
                    if accession not in self.declaration["requested_accessions"]:
                        replay.fail(f"unmanifested frozen resolution accession: {accession}")
                    return snapshot(self.declaration["resolution"])
                def fetch(self, seed):
                    if seed.study not in self.declaration["records"]: replay.fail(f"unmanifested frozen study: {seed.study}")
                    record = snapshot(self.declaration["records"][seed.study], records=True)
                    resolution = snapshot(self.declaration["resolution"])
                    if record.seed != seed or seed not in resolution.studies:
                        replay.fail("frozen record seed identity differs from requested/resolved seed")
                    return record
            services[name] = SRA2JSONConverter(source=FrozenSource(spec))
        elif name == "ae2json" and spec["kind"] == "local-ae":
            from meta_standards_converter.converters.ae2json import AE2JSONConverter
            from meta_standards_converter.sources.magetab import MAGETabInput, TextResource
            class LocalAE:
                def __init__(self, declaration): self.declaration = declaration
                def convert(self, accession, **kwargs):
                    try: entry = self.declaration["files"][accession]
                    except KeyError: replay.fail(f"unexpected local AE accession: {accession}")
                    idf_id = entry["idf"]
                    verify_local_companions(replay, idf_id, entry["sdrfs"])
                    def resource(file_id):
                        path = replay.root / replay.files[file_id]["path"]
                        return TextResource(path.name, replay.read(file_id).decode("utf-8"), str(path))
                    idf = resource(idf_id)
                    resolved = MAGETabInput(idf, tuple(resource(file_id) for file_id in entry["sdrfs"]), idf.origin, "path")
                    class Fetcher:
                        def resolve(self, source, **kwargs): return resolved
                    return AE2JSONConverter(fetcher=Fetcher()).convert(idf.origin)
            services[name] = LocalAE(spec)
        elif name == "enricher" and spec["kind"] == "no-publication":
            class NoPublication:
                publication_issues = []
                def enrich_selected(self, package, **kwargs): return package
            services[name] = NoPublication()
        else:
            raise HarnessError(f"unsupported explicit service: {name}")
    return services


def verify_local_companions(replay, idf_id, expected_ids=None):
    """Every implicit SDRF must be explicitly manifested and rehashed."""
    path = replay.root / replay.files[idf_id]["path"]
    references = [value for row in csv.reader(replay.read(idf_id).decode("utf-8").splitlines(), delimiter="\t")
                  if row and row[0] == "SDRF File" for value in row[1:] if value]
    known = {(replay.root / value["path"]).resolve(): key for key, value in replay.files.items()}
    actual = []
    for reference in references:
        candidate = (path.parent / reference).resolve()
        if candidate not in known: replay.fail(f"unmanifested SDRF companion: {reference}")
        file_id = known[candidate]
        replay.read(file_id)
        actual.append(file_id)
    if not actual: replay.fail("local IDF has no manifested SDRF companion")
    if expected_ids is not None and actual != expected_ids: replay.fail("explicit SDRF companion list differs from IDF")
    return actual


def validate_route_intent(scenario):
    intent = scenario["route_intent"]
    services = scenario["services"]
    preset = scenario["options"].get("enrichment", "standard")
    if intent == "linked-only":
        if services.get("enricher", {}).get("kind") != "no-publication" or services.get("ae2json", {}).get("kind") != "local-ae":
            raise HarnessError("linked-only requires explicit no-publication and local-ae services")
        options = scenario["options"].get("input_options", {})
        if options.get("include_peer") is not False or options.get("enrich_from_geo_ae") is not True:
            raise HarnessError("linked-only must explicitly disable peer and enable linked metadata")
    elif intent == "production-standard":
        if preset != "standard" or any(name in services for name in ("enricher", "ae2json", "linked_enricher")):
            raise HarnessError("production-standard requires production enrichment services and standard preset")
    elif intent in {"deposited-off", "parser-native-off", "loaded-export-off"}:
        if preset != "off" or any(name in services for name in ("enricher", "ae2json", "linked_enricher")):
            raise HarnessError("off route requires explicit off preset and no enrichment replacements")
    else:
        raise HarnessError("unknown or missing explicit route intent")
    formats = scenario["formats"]
    if not isinstance(formats, list) or not formats or formats[0] != "json" or len(set(formats)) != len(formats):
        raise HarnessError("formats must be a nonempty unique sequence beginning with json")
    if set(formats) != set(scenario["expected"]): raise HarnessError("format expectation coverage differs from declared calls")
    for options in (scenario["options"].get("input_options", {}), scenario["options"].get("options", {})):
        if "sdrf_sources" in options:
            raise HarnessError("SDRF overrides must be declared by a local-ae IDF/SDRF manifest, not option paths")


def validate_payload_count(payload, fmt, wanted):
    if type(wanted) is not int or wanted < 0: raise HarnessError("invalid required payload count")
    if payload is None: actual = 0
    elif fmt in {"json", "magetab"} and isinstance(payload, (list, tuple)): actual = len(payload)
    elif fmt in {"tsv", "csv"} and isinstance(payload, dict): actual = len(payload)
    else: raise HarnessError(f"unexpected {fmt} payload shape")
    if actual != wanted: raise HarnessError(f"{fmt} artifact/package count {actual} != {wanted}")


def preparation_snapshot(package):
    """Audit-only semantic fingerprint plus interpretable scientific counts."""
    value = primitive(package)
    samples = value.get("sample", [])
    paths = value.get("series", {}).get("assay_paths", [])
    runs = [run for sample in samples for run in sample.get("sra_run", [])]
    conflicts = [record for record in value.get("extensions", {}).get("insdc", {}).get("records", [])
                 if record.get("kind") == "source_conflict"]
    file_identities, run_identities = [], []
    for sample in samples:
        for run in sample.get("sra_run", []):
            scope = {"sample": sample.get("iid"), "run": run.get("run"), "library": run.get("experiment")}
            run_identities.append(scope | {key: run[key] for key in ("library_layout", "library_source", "library_strategy", "library_selection") if key in run})
            for file in run.get("fastq_files", []):
                file_identities.append(scope | {key: file[key] for key in ("uri", "filename", "format", "role", "md5", "bytes", "submitted_file_name") if key in file})
    def counted(records):
        return [{"identity": json.loads(key), "occurrences": count} for key, count in sorted(occurrence_counter(records).items())]
    return {"sha256": hashlib.sha256(json.dumps(value, sort_keys=True, ensure_ascii=False).encode()).hexdigest(),
            "sample_ids": [sample.get("iid") for sample in samples],
            "run_occurrences": len(runs), "assay_path_occurrences": len(paths),
            "native_file_occurrences": len(file_identities),
            "retained_archive_file_occurrences": sum(len(run.get("files", [])) for run in runs),
            "run_identities": counted(run_identities), "native_file_identities": counted(file_identities),
            "protocol_count": len(value.get("series", {}).get("protocols", [])),
            "conflicts": conflicts}


@contextmanager
def observe_preparation(replay, observations):
    """Read-only wrappers around real preparation services; restored on exit."""
    from meta_standards_converter.converters.unified.preparation import MetadataPreparation
    from functools import wraps
    from unittest.mock import patch
    with ExitStack() as stack:
        for name in ("_package", "_peer", "_linked", "_standard"):
            original = getattr(MetadataPreparation, name)
            def instrument(method, label):
                @wraps(method)
                def wrapped(self, package, *args, **kwargs):
                    observation = {"route": replay.route, "call": replay.call, "operation": label,
                                   "http_before": len(replay.events), "before": preparation_snapshot(package)}
                    observations.append(observation)
                    start = time.monotonic()
                    try:
                        result = method(self, package, *args, **kwargs)
                        observation["after"] = preparation_snapshot(result[0])
                        observation["status"] = "returned"
                        return result
                    except BaseException as exc:
                        observation.update(status="raised", error_type=type(exc).__name__)
                        raise
                    finally:
                        observation.update(http_after=len(replay.events), seconds=time.monotonic() - start)
                return wrapped
            stack.enter_context(patch.object(MetadataPreparation, name, instrument(original, name.lstrip("_"))))
        yield


def audit(manifest, *, evidence_root=None, artifacts_dir=None):
    report = {"classification": "harness-error", "scenarios": [], "errors": []}
    replay = None
    try:
        if sys.flags.optimize:
            raise HarnessError("optimized Python disables acceptance assertions; run without -O")
        declaration = json.loads(Path(manifest).read_text())
        if declaration.get("external_evidence") and evidence_root is None:
            raise HarnessError("unmet prerequisite: explicit --evidence-root is required")
        if declaration.get("unmet_prerequisites"):
            raise HarnessError(f"unmet prerequisites: {declaration['unmet_prerequisites']}")
        replay = ReplayTransport(manifest, evidence_root=evidence_root)
        doc = replay.manifest
        report.update(cohort=doc["cohort"], evidence_kind=doc["evidence_kind"],
                      source_revision=doc["source_revision"], captured_at=doc["captured_at"])
        scenarios = doc.get("scenarios", [])
        if not scenarios: raise HarnessError("no declared executable scenarios")
        if artifacts_dir is not None: Path(artifacts_dir).mkdir(parents=True, exist_ok=True)
        with replay.offline():
            # Bootstrap the DNS patch before importing MSC. Pacing is disabled
            # only for recorded replies, never for any live capture operation.
            from meta_standards_converter.converters import Converter
            from meta_standards_converter.miniml import MINiMLCodec
            from meta_standards_converter.helpers.request_helper import HostRequestGate
            from unittest.mock import patch
            observations = []
            report["preparation_observations"] = observations
            with tempfile.TemporaryDirectory(prefix="msc-replay-gates-") as gate_dir, \
                    patch.dict(os.environ, {"SCIENTIFIC_PROVIDER_GATE_DIR": gate_dir}), \
                    patch.object(HostRequestGate, "_default", None), \
                    patch.object(HostRequestGate, "wait", return_value=0.0), observe_preparation(replay, observations):
                for scenario in scenarios:
                    replay.route = scenario["id"]
                    validate_route_intent(scenario)
                    entry = {"id": scenario["id"], "route_intent": scenario["route_intent"], "classification": scenario["classification"], "calls": []}
                    report["scenarios"].append(entry)
                    if scenario["classification"] not in {"complete-captured-coverage", "expected-failure"}:
                        raise HarnessError("scenario is not declared executable")
                    services = build_services(replay, scenario["services"])
                    spec = scenario["source"]
                    if spec["kind"] == "file":
                        replay.read(spec["file"])
                        source = replay.root / replay.files[spec["file"]]["path"]
                        if source.name.endswith(".idf.txt"): verify_local_companions(replay, spec["file"])
                    elif spec["kind"] == "json": source = json.loads(replay.read(spec["file"]))
                    elif spec["kind"] == "accession": source = spec["value"]
                    else: raise HarnessError("unsupported source kind")
                    golden = json.loads(replay.read(scenario["golden"])) if "golden" in scenario else None
                    packages = None
                    for call_index, fmt in enumerate(scenario["formats"]):
                        replay.call, replay.stage = call_index + 1, fmt
                        start = time.monotonic()
                        input_evidence = primitive(packages) if packages is not None else spec
                        call = {"format": fmt, "http_before": len(replay.events), "preparation_before": input_evidence}
                        call["status"] = "running"
                        entry["calls"].append(call)
                        # Evidence uses bounded hashes/counts for preparation
                        # input rather than duplicating full scientific output.
                        before = json.dumps(call.pop("preparation_before"), sort_keys=True).encode()
                        call["input_sha256"] = hashlib.sha256(before).hexdigest()
                        wanted = scenario["expected"][fmt]
                        try:
                            result = Converter(services=services if packages is None else {}).convert(source if packages is None else packages, out_type=fmt,
                                **(scenario["options"] if packages is None else {"enrichment": "off", "options": {"expand_studies": False}}))
                        except ReplayBudgetExceeded:
                            replay.assert_no_errors()
                            if scenario["classification"] != "expected-failure" or wanted.get("exception") != "budget-exhausted":
                                raise HarnessError("undeclared budget exhaustion")
                            call.update(exception="budget-exhausted", seconds=time.monotonic() - start,
                                        http_after=len(replay.events), status="expected-failure")
                            continue
                        call.update(result=result.to_dict(), preparation=[primitive(item.preparation) for item in result.items])
                        replay.assert_no_errors()  # Mandatory even when MSC caught the sentinel.
                        expected_count = scenario["expected_item_count"]
                        if expected_count != 1:
                            raise HarnessError("this runner requires one explicitly declared logical outcome per scenario")
                        if len(result.items) != expected_count:
                            raise HarnessError(f"facade outcome count {len(result.items)} != {expected_count}")
                        for item in result.items:
                            for axis in ("execution", "validation", "completeness"):
                                assert getattr(item, axis) == wanted[axis], f"{fmt}/{axis}: {getattr(item, axis)} != {wanted[axis]}"
                            causes = [{"code": d.code, "stage": d.stage, "severity": d.severity, "message": d.message}
                                      for d in item.diagnostics if d.severity in {"warning", "error"}]
                            assert_occurrences(wanted["causes"], causes)
                            if wanted["execution"] == "succeeded":
                                assert item.payload is not None, f"{fmt}: required in-memory artifact missing"
                        payload = primitive(result.items[0].payload)
                        if artifacts_dir is not None:
                            label = hashlib.sha256(scenario["id"].encode()).hexdigest()[:16] + "." + fmt + ".json"
                            artifact = Path(artifacts_dir) / label
                            with artifact.open("x") as handle: json.dump(payload, handle, ensure_ascii=False)
                            call["audit_artifact"] = str(artifact.resolve())
                        validate_payload_count(payload, fmt, wanted["payload_count"])
                        expected_payload = golden[fmt] if golden is not None else {"kind": "exact", "value": json.loads(replay.read(scenario["goldens"][fmt]))}
                        normalization = scenario.get("runtime_date") if fmt == "magetab" else None
                        if normalization is not None:
                            if expected_payload["kind"] != "exact": raise HarnessError("date normalization requires exact MAGE-TAB golden")
                            compared, expected_payload["value"] = normalize_runtime_date(payload, expected_payload["value"], normalization)
                        else: compared = payload
                        compare_golden(compared, expected_payload)
                        if fmt == "json" and result.items[0].payload is not None:
                            packages = result.items[0].payload
                            for package in packages: MINiMLCodec().decode(package, strict=True)
                        call.update(seconds=time.monotonic() - start, http_after=len(replay.events),
                                    preparation=[primitive(item.preparation) for item in result.items],
                                    result=result.to_dict(), status="passed", output_sha256=hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest())
                    replay.assert_no_errors()
                replay.assert_complete()
        report["classification"] = "expected-failure" if any(s["classification"] == "expected-failure" for s in scenarios) else "complete-captured-coverage"
    except Exception as exc:
        report["errors"].append(f"{type(exc).__name__}: {exc}")
        if report["scenarios"] and report["scenarios"][-1]["calls"]:
            failed = report["scenarios"][-1]["calls"][-1]
            if failed["status"] == "running":
                failed.update(status="harness-error", seconds=time.monotonic() - start,
                              http_after=len(replay.events), error=report["errors"][-1])
    if replay is not None:
        report["transport"] = {"events": replay.events, "errors": replay.errors,
                               "unexpected_misses": len(replay.errors)}
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--evidence-root", type=Path)
    parser.add_argument("--report", type=Path)
    parser.add_argument("--artifacts-dir", type=Path, help="Optional fresh audit payload files; existing files are never overwritten")
    args = parser.parse_args(argv)
    report = audit(args.manifest, evidence_root=args.evidence_root, artifacts_dir=args.artifacts_dir)
    rendered = json.dumps(report, indent=2)
    if args.report: args.report.write_text(rendered + "\n")
    print(rendered)
    return 2 if report["classification"] == "harness-error" else 0


if __name__ == "__main__": raise SystemExit(main())
