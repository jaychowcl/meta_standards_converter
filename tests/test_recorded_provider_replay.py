"""Fail-closed replay contracts, independent of converter exception handling."""
import copy
import hashlib
import json
import socket
import subprocess
import sys
from pathlib import Path

import pytest
import requests

from tests.support.recorded_replay import (
    HarnessError, ReplayTransport, assert_occurrences, request_identity,
)


def corpus(tmp_path, *, status=200, repeat=1):
    body = b'{"scientific": "literal value"}'
    (tmp_path / "response.body").write_bytes(body)
    request = requests.Request("GET", "https://www.ebi.ac.uk/ena/test?x=1").prepare()
    manifest = {
        "schema_version": "1.0", "cohort": "unit-explicit", "source_revision": None,
        "captured_at": None, "evidence_kind": "recorded-http",
        "files": {"response": {"path": "response.body", "sha256": hashlib.sha256(body).hexdigest()}},
        "requests": [{"identity": request_identity(request), "attempts": [{
            "response": "response", "status": status, "headers": {"Content-Type": "application/json"},
            "encoding": "utf-8", "repeat": repeat,
        }]}],
    }
    path = tmp_path / "manifest.json"
    path.write_text(json.dumps(manifest))
    return path, request, manifest


def test_selected_response_repeat_status_headers_and_restoration(tmp_path):
    path, request, _ = corpus(tmp_path, status=403, repeat=2)
    original = requests.Session.send
    replay = ReplayTransport(path)
    with replay.offline():
        for _ in range(2):
            response = requests.Session().send(request)
            assert response.status_code == 403
            assert response.encoding == "utf-8"
            assert response.json() == {"scientific": "literal value"}
        replay.assert_complete()
    assert requests.Session.send is original
    assert [event["status"] for event in replay.events] == [403, 403]


@pytest.mark.parametrize("fault", ["absent", "corrupt", "traversal", "duplicate", "unknown-version", "bad-repeat"])
def test_manifest_faults_are_harness_errors(tmp_path, fault):
    path, _, manifest = corpus(tmp_path)
    if fault == "absent": (tmp_path / "response.body").unlink()
    if fault == "corrupt": (tmp_path / "response.body").write_bytes(b"different")
    if fault == "traversal": manifest["files"]["response"]["path"] = "../response.body"
    if fault == "duplicate": manifest["requests"] *= 2
    if fault == "unknown-version": manifest["schema_version"] = "2.0"
    if fault == "bad-repeat": manifest["requests"][0]["attempts"][0]["repeat"] = 0
    path.write_text(json.dumps(manifest))
    with pytest.raises(HarnessError): ReplayTransport(path)


def test_swallowed_miss_and_excess_repeat_still_fail_postcall(tmp_path):
    path, request, _ = corpus(tmp_path)
    replay = ReplayTransport(path)
    with replay.offline():
        requests.Session().send(request)
        for wrong in [request, requests.Request("GET", "https://www.ebi.ac.uk/missing").prepare()]:
            try: requests.Session().send(wrong)
            except Exception: pass  # Production code may retain partial scientific output.
        with pytest.raises(HarnessError, match="unexpected"):
            replay.assert_complete()
    assert len(replay.errors) == 2


def test_unconsumed_expected_attempt_is_error(tmp_path):
    path, _, _ = corpus(tmp_path)
    with pytest.raises(HarnessError, match="unconsumed"):
        ReplayTransport(path).assert_complete()


@pytest.mark.parametrize("fault", ["timeout", "connection-error", "budget-exhausted"])
def test_explicit_transport_fault_followed_by_selected_response(tmp_path, fault):
    from tests.support.recorded_replay import TRANSPORT_EXCEPTIONS
    path, request, manifest = corpus(tmp_path)
    manifest["requests"][0]["attempts"].insert(0, {"exception": fault, "repeat": 1})
    path.write_text(json.dumps(manifest))
    replay = ReplayTransport(path)
    with replay.offline():
        with pytest.raises(TRANSPORT_EXCEPTIONS[fault]): requests.Session().send(request)
        assert requests.Session().send(request).status_code == 200
        replay.assert_complete()
    assert not replay.errors
    assert replay.events[0]["exception"] == fault


def test_response_is_rehashed_at_use(tmp_path):
    path, request, _ = corpus(tmp_path)
    replay = ReplayTransport(path)
    (tmp_path / "response.body").write_bytes(b"changed after preflight")
    with replay.offline(), pytest.raises(HarnessError):
        requests.Session().send(request)
    with pytest.raises(HarnessError): replay.assert_complete()


def test_http_adapter_seam_records_redirect_attempts_separately(tmp_path):
    path, request, manifest = corpus(tmp_path, status=302)
    manifest["transport"] = "http-adapter"
    manifest["requests"][0]["attempts"][0]["headers"]["Location"] = "https://www.ebi.ac.uk/second"
    second = requests.Request("GET", "https://www.ebi.ac.uk/second").prepare()
    attempt = copy.deepcopy(manifest["requests"][0]["attempts"][0])
    attempt.update(status=200, headers={"Content-Type": "application/json"})
    manifest["requests"].append({"identity": request_identity(second), "attempts": [attempt]})
    path.write_text(json.dumps(manifest))
    replay = ReplayTransport(path)
    with replay.offline():
        response = requests.Session().send(request)
        assert response.status_code == 200
        assert [r.status_code for r in response.history] == [302]
        replay.assert_complete()
    assert [event["status"] for event in replay.events] == [302, 200]


def test_offline_blocks_dns_network_and_subprocess_without_disabling_url_policy(tmp_path):
    path, _, _ = corpus(tmp_path)
    replay = ReplayTransport(path)
    original_dns = socket.getaddrinfo
    with replay.offline():
        from meta_standards_converter.retrieval import RetrievalPolicy, RetrievalSecurityError
        RetrievalPolicy().validate_url("https://www.ebi.ac.uk/ena/test")
        for url in ["https://127.0.0.1/private", "https://unlisted.invalid/", "https://user:pass@www.ebi.ac.uk/x"]:
            with pytest.raises(RetrievalSecurityError): RetrievalPolicy().validate_url(url)
        for operation in [lambda: socket.create_connection(("www.ebi.ac.uk", 443)), lambda: subprocess.run(["true"])]:
            with pytest.raises(HarnessError): operation()
    assert socket.getaddrinfo is original_dns


def test_request_identity_preserves_query_order_duplicates_and_hashes_body():
    a = requests.Request("POST", "https://www.ebi.ac.uk/test?a=1&a=2&api_key=secret&email=user", data=b"x=0").prepare()
    b = requests.Request("POST", "https://www.ebi.ac.uk/test?a=2&a=1&api_key=other", data=b"x=0").prepare()
    identity = request_identity(a)
    assert "secret" not in json.dumps(identity)
    assert identity["body_sha256"] == hashlib.sha256(b"x=0").hexdigest()
    assert identity != request_identity(b)


@pytest.mark.parametrize("name,args", [("gethostbyname", ("www.ebi.ac.uk",)),
                                     ("gethostbyname_ex", ("www.ebi.ac.uk",)),
                                     ("gethostbyaddr", ("93.184.216.34",)),
                                     ("getnameinfo", (("93.184.216.34", 443), 0))])
def test_alternative_dns_escapes_are_ledgered(tmp_path, name, args):
    path, _, _ = corpus(tmp_path)
    replay = ReplayTransport(path)
    with replay.offline():
        with pytest.raises(HarnessError): getattr(socket, name)(*args)
    with pytest.raises(HarnessError): replay.assert_no_errors()


def test_swallowed_unmanifested_dns_fails_postcall(tmp_path):
    path, _, _ = corpus(tmp_path)
    replay = ReplayTransport(path)
    with replay.offline():
        from meta_standards_converter.retrieval import RetrievalPolicy, RetrievalSecurityError
        try: RetrievalPolicy().validate_url("https://ftp.ebi.ac.uk/unmanifested")
        except RetrievalSecurityError: pass
    with pytest.raises(HarnessError, match="DNS"): replay.assert_no_errors()


@pytest.mark.fake_process
def test_first_import_dns_default_restored_in_fresh_python(tmp_path):
    path, _, _ = corpus(tmp_path)
    code = '''
import socket, sys
from tests.support.recorded_replay import ReplayTransport
assert 'meta_standards_converter.retrieval' not in sys.modules
original = socket.getaddrinfo
with ReplayTransport(sys.argv[1]).offline():
    from meta_standards_converter.retrieval import RetrievalPolicy
    assert RetrievalPolicy().resolver is socket.getaddrinfo
    RetrievalPolicy().validate_url('https://www.ebi.ac.uk/ena/test')
assert socket.getaddrinfo is original
assert RetrievalPolicy().resolver is original
assert RetrievalPolicy.__dataclass_fields__['resolver'].default is original
'''
    result = subprocess.run([sys.executable, "-c", code, str(path)], capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize("mutation", ["duplicate", "whitespace", "zero", "owner", "md5", "conflict"])
def test_occurrence_comparator_detects_scientific_loss(mutation):
    expected = [{"sample": "s1", "library": "l1", "run": "r1", "path": "p1", "step": 2,
                 "column": 3, "owner": "protocol", "value": " literal ", "bytes": 0,
                 "uri": "https://example/a", "md5": "abc", "conflict": ["a", "b"]}] * 2
    actual = copy.deepcopy(expected)
    if mutation == "duplicate": actual.pop()
    else:
        key, value = {"whitespace": ("value", "literal"), "zero": ("bytes", ""),
                      "owner": ("owner", "factor"), "md5": ("md5", "def"),
                      "conflict": ("conflict", ["a"])}[mutation]
        actual[0][key] = value
    with pytest.raises(AssertionError): assert_occurrences(expected, actual)
    assert_occurrences(expected, copy.deepcopy(expected))


def test_reduced_scientific_success_all_formats(tmp_path):
    from scripts.audit_recorded_provider import audit
    path = Path(__file__).parent / "fixtures/recorded_provider/success.manifest.json"
    report = audit(path)
    assert report["classification"] == "complete-captured-coverage", report
    assert [call["format"] for call in report["scenarios"][0]["calls"]] == ["json", "magetab", "tsv"]
    assert all(call["preparation"] for call in report["scenarios"][0]["calls"])


def test_reduced_conflict_is_partial_valid_in_all_formats_and_reimports():
    from scripts.audit_recorded_provider import audit
    from meta_standards_converter.converters import Converter
    from tests.metadata.test_portable_conflicts import import_tables, records
    root = Path(__file__).parent / "fixtures/recorded_provider"
    report = audit(root / "conflict.manifest.json")
    assert report["classification"] == "complete-captured-coverage", report
    assert all(call["result"]["items"][0]["completeness"] == "partial" for call in report["scenarios"][0]["calls"])
    data = json.loads((root / "E-MTAB-16253.reduced.json").read_text())
    expected = records(data)
    for _ in range(2):
        result = Converter().convert(data, out_type="magetab", enrichment="off", options={"expand_studies": False})
        assert result.items[0].execution == "succeeded"
        assert result.items[0].validation == "valid"
        assert result.items[0].completeness == "partial"
        data, _ = import_tables(result.items[0].payload[0])
        assert_occurrences(expected, records(data))


def test_full_row_evidence_manifest_hashes_and_exact_source_counts():
    root = Path(__file__).parent / "fixtures/recorded_provider"
    replay = ReplayTransport(root / "row-evidence.manifest.json")
    for study, count in {"E-MTAB-16253": 16, "E-MTAB-14566": 1558, "E-MTAB-14416": 104,
                         "E-GEOD-52564": 51, "E-GEOD-60450": 24}.items():
        ledger = json.loads(replay.read(study + ".rows"))
        assert len(ledger["rows"]) == ledger["source_rows"] == count


def test_external_corpus_requires_explicit_evidence_root(tmp_path):
    from scripts.audit_recorded_provider import audit
    path, _, manifest = corpus(tmp_path)
    manifest["external_evidence"] = True
    path.write_text(json.dumps(manifest))
    report = audit(path)
    assert report["classification"] == "harness-error"
    assert "explicit --evidence-root" in report["errors"][0]


def test_explicit_local_ae_companions_are_verified_and_no_implicit_reads(tmp_path):
    from scripts.audit_recorded_provider import build_services
    root = Path(__file__).parent / "fixtures/recorded_provider"
    replay = ReplayTransport(root / "success.manifest.json")
    services = build_services(replay, {"ae2json": {"kind": "local-ae", "files": {
        "E-MTAB-16847": {"idf": "idf", "sdrfs": ["sdrf"]}}}})
    with replay.offline():
        assert services["ae2json"].convert("E-MTAB-16847")[0]["series"]["iid"] == "E-MTAB-16847"
    replay.files.pop("sdrf")
    with pytest.raises(HarnessError, match="unmanifested SDRF"):
        services["ae2json"].convert("E-MTAB-16847")
    with pytest.raises(HarnessError): replay.assert_no_errors()


def test_frozen_source_reconstructs_only_json_and_safe_xml(tmp_path):
    from scripts.audit_recorded_provider import build_services
    from meta_standards_converter.sources.archive_support import StudySeed
    path, _, manifest = corpus(tmp_path)
    manifest["evidence_kind"] = "parser-merge-replay"
    snapshots = {"resolution": {"schema_version": "1.0", "studies": [{"study": "SRP1", "primary": "SRP1"}], "issues": []},
                 "records": {"schema_version": "1.0", "seed": {"study": "SRP1", "primary": "SRP1"},
                             "xml": ['<STUDY accession="SRP1"/>'], "indexed": {}, "linked": [], "issues": []}}
    for name, value in snapshots.items():
        file = tmp_path / (name + ".json")
        file.write_text(json.dumps(value))
        manifest["files"][name] = {"path": file.name, "sha256": hashlib.sha256(file.read_bytes()).hexdigest()}
    path.write_text(json.dumps(manifest))
    replay = ReplayTransport(path)
    service = build_services(replay, {"sra2json": {"kind": "frozen-source", "resolution": "resolution", "records": {"SRP1": "records"}}})["sra2json"]
    source = service.source
    seed = source.resolve("SRP1").studies[0]
    assert seed == StudySeed("SRP1", "SRP1")
    assert source.fetch(seed).xml[0].get("accession") == "SRP1"
    # A valid-hash file is still parsed as data; it cannot become executable pickle.
    file = tmp_path / "records.json"
    file.write_bytes(b"cos\nsystem\n(S'touch should-not-exist'\ntR.")
    replay.files["records"]["sha256"] = hashlib.sha256(file.read_bytes()).hexdigest()
    with pytest.raises((ValueError, HarnessError)): source.fetch(seed)
    assert not (tmp_path / "should-not-exist").exists()


def test_runner_reports_harness_error_and_nonzero_for_missing_evidence(tmp_path, capsys):
    from scripts.audit_recorded_provider import main
    path, _, _ = corpus(tmp_path)
    (tmp_path / "response.body").unlink()
    assert main([str(path)]) == 2
    assert json.loads(capsys.readouterr().out)["classification"] == "harness-error"
