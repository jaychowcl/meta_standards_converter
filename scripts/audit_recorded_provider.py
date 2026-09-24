#!/usr/bin/env python3
"""Audit explicit recorded cohorts. No live fallback or golden regeneration.

Run with the existing environment: python scripts/audit_recorded_provider.py
tests/fixtures/recorded_provider/success.manifest.json. External corpora require
--evidence-root; relative file IDs always resolve against that one root.
"""
import argparse
import csv
from dataclasses import asdict, is_dataclass
import json
from pathlib import Path
import sys
import time

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path: sys.path.insert(0, str(REPO))
from tests.support.recorded_replay import HarnessError, ReplayTransport, ReplayBudgetExceeded, assert_occurrences


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
        assert payload == expected["value"], "exact scientific payload mismatch"
    else:
        raise HarnessError(f"unknown golden comparator: {kind}")


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
                    return snapshot(self.declaration["resolution"])
                def fetch(self, seed):
                    if seed.study not in self.declaration["records"]: replay.fail(f"unmanifested frozen study: {seed.study}")
                    return snapshot(self.declaration["records"][seed.study], records=True)
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


def audit(manifest, *, evidence_root=None):
    report = {"classification": "harness-error", "scenarios": [], "errors": []}
    replay = None
    try:
        declaration = json.loads(Path(manifest).read_text())
        if declaration.get("external_evidence") and evidence_root is None:
            raise HarnessError("unmet prerequisite: explicit --evidence-root is required")
        replay = ReplayTransport(manifest, evidence_root=evidence_root)
        doc = replay.manifest
        report.update(cohort=doc["cohort"], evidence_kind=doc["evidence_kind"],
                      source_revision=doc["source_revision"], captured_at=doc["captured_at"])
        scenarios = doc.get("scenarios", [])
        if not scenarios: raise HarnessError("no declared executable scenarios")
        with replay.offline():
            # Bootstrap the DNS patch before importing MSC. Pacing is disabled
            # only for recorded replies, never for any live capture operation.
            from meta_standards_converter.converters import Converter
            from meta_standards_converter.miniml import MINiMLCodec
            from meta_standards_converter.helpers.request_helper import HostRequestGate
            from unittest.mock import patch
            with patch.object(HostRequestGate, "wait", return_value=0.0):
                for scenario in scenarios:
                    replay.route = scenario["id"]
                    entry = {"id": scenario["id"], "classification": scenario["classification"], "calls": []}
                    report["scenarios"].append(entry)
                    if scenario["classification"] not in {"complete-captured-coverage", "expected-failure"}:
                        raise HarnessError("scenario is not declared executable")
                    services = build_services(replay, scenario["services"])
                    if scenario["options"].get("enrichment") == "standard" and any(name in services for name in ("ae2json", "enricher")):
                        raise HarnessError("standard replay requires production enrichment services")
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
                        # Evidence uses bounded hashes/counts for preparation
                        # input rather than duplicating full scientific output.
                        before = json.dumps(call.pop("preparation_before"), sort_keys=True).encode()
                        import hashlib
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
                                        http_after=len(replay.events))
                            entry["calls"].append(call)
                            continue
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
                        expected_payload = golden[fmt] if golden is not None else {"kind": "exact", "value": json.loads(replay.read(scenario["goldens"][fmt]))}
                        compare_golden(payload, expected_payload)
                        if fmt == "json" and result.items[0].payload is not None:
                            packages = result.items[0].payload
                            for package in packages: MINiMLCodec().decode(package, strict=True)
                        call.update(seconds=time.monotonic() - start, http_after=len(replay.events),
                                    preparation=[primitive(item.preparation) for item in result.items],
                                    result=result.to_dict(), output_sha256=hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest())
                        entry["calls"].append(call)
                    replay.assert_no_errors()
                replay.assert_complete()
        report["classification"] = "expected-failure" if any(s["classification"] == "expected-failure" for s in scenarios) else "complete-captured-coverage"
    except Exception as exc:
        report["errors"].append(f"{type(exc).__name__}: {exc}")
    if replay is not None:
        report["transport"] = {"events": replay.events, "errors": replay.errors,
                               "unexpected_misses": len(replay.errors)}
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--evidence-root", type=Path)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args(argv)
    report = audit(args.manifest, evidence_root=args.evidence_root)
    rendered = json.dumps(report, indent=2)
    if args.report: args.report.write_text(rendered + "\n")
    print(rendered)
    return 2 if report["classification"] == "harness-error" else 0


if __name__ == "__main__": raise SystemExit(main())
