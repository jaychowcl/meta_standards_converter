"""Explicit, hash-verified provider replay; never searches development folders.

Importing this module imports no MSC code. Enter ``offline`` before importing
MSC in standalone runners so dataclass-bound DNS defaults are also offline.
"""
from collections import Counter
from contextlib import ExitStack, contextmanager
import hashlib
import ipaddress
import json
import os
from pathlib import Path
import socket
import subprocess
import sys
from urllib.parse import parse_qsl, urlencode, urlsplit, urlunsplit
from unittest.mock import patch

import requests


class HarnessError(RuntimeError):
    """Missing/corrupt evidence or an undeclared side effect, never science."""


class ReplayBudgetExceeded(BaseException):
    """An explicitly captured study budget exhaustion, not a replay miss."""


TRANSPORT_EXCEPTIONS = {"timeout": requests.Timeout, "connection-error": requests.ConnectionError,
                        "budget-exhausted": ReplayBudgetExceeded}


def _canonical(value):
    return json.dumps(value, sort_keys=True, ensure_ascii=False, separators=(",", ":"), allow_nan=False)


def request_identity(request):
    """Identity of a prepared request, redacting only declared credentials.

    The prepared query's order and repeated keys are significant. Body contents
    are never logged; their hash distinguishes POST requests without disclosure.
    """
    parts = urlsplit(request.url)
    if parts.username is not None or parts.password is not None:
        raise HarnessError("prepared request contains URL credentials")
    query = [(key, value) for key, value in parse_qsl(parts.query, keep_blank_values=True)
             if key.casefold() not in {"api_key", "email", "tool", "token", "access_token"}]
    body = request.body or b""
    if isinstance(body, str): body = body.encode("utf-8")
    if not isinstance(body, bytes): raise HarnessError("streaming request body cannot be replayed")
    return {"method": request.method.upper(), "url": urlunsplit(parts._replace(query=urlencode(query), fragment="")),
            "body_sha256": hashlib.sha256(body).hexdigest()}


def occurrence_counter(records, identity_fields=None):
    """Count exact typed occurrences; no trimming, truthiness or set collapse.

    Omit identity_fields to compare every field, including owner/path/step,
    file URI/name/role/md5/bytes, factor cells and conflict identities.
    """
    return Counter(_canonical(record if identity_fields is None else {key: record[key] for key in identity_fields})
                   for record in records)


def assert_occurrences(expected, actual, identity_fields=None):
    left, right = occurrence_counter(expected, identity_fields), occurrence_counter(actual, identity_fields)
    assert left == right, f"lost occurrences: {dict(left - right)}; extra occurrences: {dict(right - left)}"


def table_occurrences(table):
    """Lossless row projection with repeated-column ordinal identities."""
    if not table: raise HarnessError("empty table")
    headers = table[0]
    counts = Counter()
    columns = []
    for header in headers:
        counts[header] += 1
        columns.append((header, counts[header]))
    result = []
    for row in table[1:]:
        if len(row) != len(headers): raise HarnessError("nonrectangular table")
        result.append([[name, ordinal, value] for (name, ordinal), value in zip(columns, row)])
    return result


class ReplayTransport:
    """Strict version-1 cohort with explicit ordered attempts and repetitions."""

    def __init__(self, manifest, *, evidence_root=None):
        self.manifest_path = Path(manifest).resolve()
        self.root = Path(evidence_root).resolve() if evidence_root is not None else self.manifest_path.parent
        self.errors, self.events = [], []
        self.stage = "setup"
        self.route = None
        self.call = None
        try:
            self.manifest = json.loads(self.manifest_path.read_text())
            doc = self.manifest
            if doc["schema_version"] != "1.0": raise HarnessError("unsupported manifest version")
            for field in ("cohort", "source_revision", "captured_at", "evidence_kind", "files", "requests"):
                if field not in doc: raise HarnessError(f"missing manifest field: {field}")
            if doc["evidence_kind"] not in {"recorded-http", "parser-merge-replay", "reduced-source"}:
                raise HarnessError("unknown evidence kind")
            if doc.get("transport", "requests-session") not in {"requests-session", "http-adapter"}:
                raise HarnessError("unsupported transport interception seam")
            self.files = doc["files"]
            for file_id in self.files: self.read(file_id)
            self.attempts = {}
            for item in doc["requests"]:
                identity = item["identity"]
                if set(identity) != {"method", "url", "body_sha256"}: raise HarnessError("invalid request identity")
                key = _canonical(identity)
                if key in self.attempts: raise HarnessError("duplicate request identity")
                if request_identity(requests.Request(identity["method"], identity["url"]).prepare())["url"] != identity["url"]:
                    raise HarnessError("request URL is not prepared and redacted")
                if len(identity["body_sha256"]) != 64: raise HarnessError("invalid request body hash")
                selected = []
                for attempt in item["attempts"]:
                    repeat = attempt["repeat"]
                    if type(repeat) is not int or repeat < 1: raise HarnessError("invalid expected repetition")
                    if "exception" in attempt:
                        if attempt["exception"] not in TRANSPORT_EXCEPTIONS or "response" in attempt:
                            raise HarnessError("invalid selected transport exception")
                    else:
                        if attempt["response"] not in self.files: raise HarnessError("unmanifested response")
                        if type(attempt["status"]) is not int or not 100 <= attempt["status"] <= 599:
                            raise HarnessError("invalid HTTP status")
                        if not isinstance(attempt["headers"], dict): raise HarnessError("invalid response headers")
                        if attempt["encoding"] is not None and not isinstance(attempt["encoding"], str):
                            raise HarnessError("invalid encoding")
                    selected.extend([attempt] * repeat)
                if not selected: raise HarnessError("missing selected response attempts")
                self.attempts[key] = selected
            self.consumed = Counter()
        except HarnessError:
            raise
        except (OSError, ValueError, TypeError, KeyError) as exc:
            raise HarnessError(f"invalid manifest: {exc}") from exc

    def read(self, file_id):
        try:
            item = self.files[file_id]
            relative = Path(item["path"])
            if relative.is_absolute() or ".." in relative.parts: raise HarnessError("unsafe evidence path")
            path = self.root / relative
            if not path.resolve().is_relative_to(self.root): raise HarnessError("evidence symlink leaves root")
            content = path.read_bytes()
            if hashlib.sha256(content).hexdigest() != item["sha256"]: raise HarnessError(f"corrupt evidence: {file_id}")
            return content
        except HarnessError as exc:
            self.fail(str(exc))
        except (OSError, KeyError, TypeError) as exc:
            self.fail(f"missing evidence: {file_id}")

    def fail(self, reason):
        self.errors.append({"route": self.route, "call": self.call, "stage": self.stage, "reason": reason})
        raise HarnessError(reason)

    def send(self, session, request, **kwargs):
        try:
            identity = request_identity(request)
            key = _canonical(identity)
            attempt_index = self.consumed[key]
            if key not in self.attempts: self.fail(f"unexpected request: {identity}")
            if attempt_index >= len(self.attempts[key]): self.fail(f"unexpected repetition: {identity}")
            attempt = self.attempts[key][attempt_index]
            if "exception" in attempt:
                self.consumed[key] += 1
                self.events.append({"route": self.route, "call": self.call, "stage": self.stage,
                                    "identity": identity, "attempt": attempt_index + 1,
                                    "exception": attempt["exception"]})
                raise TRANSPORT_EXCEPTIONS[attempt["exception"]]("Recorded " + attempt["exception"])
            content = self.read(attempt["response"])
            response = requests.Response()
            response.status_code = attempt["status"]
            response.headers.update(attempt["headers"])
            response.encoding = attempt["encoding"]
            response._content = content
            response._content_consumed = True
            response.url, response.request = request.url, request
            self.consumed[key] += 1
            self.events.append({"route": self.route, "call": self.call, "stage": self.stage,
                                "identity": identity, "attempt": attempt_index + 1,
                                "status": response.status_code, "response": attempt["response"],
                                "sha256": hashlib.sha256(content).hexdigest()})
            return response
        except HarnessError as exc:
            if not self.errors or self.errors[-1]["reason"] != str(exc):
                self.errors.append({"route": self.route, "call": self.call, "stage": self.stage, "reason": str(exc)})
            raise

    def assert_no_errors(self):
        if self.errors: raise HarnessError(f"unexpected replay errors: {self.errors}")

    def assert_complete(self):
        self.assert_no_errors()
        remaining = {key: len(attempts) - self.consumed[key] for key, attempts in self.attempts.items()
                     if self.consumed[key] != len(attempts)}
        if remaining: raise HarnessError(f"unconsumed expected attempts: {remaining}")

    @contextmanager
    def offline(self):
        """Isolate DNS/network/processes and restore every patch on exit.

        Only manifest hosts receive a deterministic public DNS answer. URL
        validation still rejects credentials, private IPs and disallowed hosts.
        Already-imported RetrievalPolicy bound defaults are patched reversibly
        for pytest; standalone callers enter before the first MSC import.
        """
        hosts = {urlsplit(json.loads(key)["url"]).hostname for key in self.attempts}
        original_dns = socket.getaddrinfo
        imported_before = "meta_standards_converter.retrieval" in sys.modules
        def dns(host, port, *args, **kwargs):
            try: address = str(ipaddress.ip_address(host))
            except ValueError:
                if host not in hosts:
                    self.errors.append({"route": self.route, "call": self.call, "stage": self.stage,
                                        "reason": f"unexpected offline DNS host: {host}"})
                    raise OSError("unmanifested offline DNS host")
                address = "93.184.216.34"
            return [(socket.AF_INET, socket.SOCK_STREAM, 6, "", (address, port))]
        def blocked(*args, **kwargs): self.fail("unexpected network or subprocess escape")
        with ExitStack() as stack:
            stack.enter_context(patch.object(socket, "getaddrinfo", dns))
            for name in ("gethostbyname", "gethostbyname_ex", "gethostbyaddr", "getnameinfo"):
                stack.enter_context(patch.object(socket, name, blocked))
            for target, names in [(socket.socket, ("connect", "connect_ex", "sendto")),
                                  (socket, ("create_connection",)),
                                  (subprocess, ("Popen", "run", "call", "check_call", "check_output")),
                                  (os, ("system", "popen", "posix_spawn", "posix_spawnp", "spawnl", "spawnlp", "spawnv", "spawnvp"))]:
                for name in names:
                    if hasattr(target, name): stack.enter_context(patch.object(target, name, blocked))
            target = requests.adapters.HTTPAdapter if self.manifest.get("transport") == "http-adapter" else requests.Session
            stack.enter_context(patch.object(target, "send", lambda session, request, **kw: self.send(session, request, **kw)))
            # Also handle pytest modules which imported MSC during collection.
            module = sys.modules.get("meta_standards_converter.retrieval")
            if module is not None:
                policy = module.RetrievalPolicy
                original = policy.__dataclass_fields__["resolver"].default
                defaults = policy.__init__.__defaults__
                stack.enter_context(patch.object(policy.__init__, "__defaults__", tuple(dns if value is original else value for value in defaults)))
                stack.enter_context(patch.object(policy, "resolver", dns))
            try:
                yield self
            finally:
                # A class imported inside the context captured dns in generated
                # defaults. ExitStack cannot restore an object absent at entry.
                if not imported_before and "meta_standards_converter.retrieval" in sys.modules:
                    policy = sys.modules["meta_standards_converter.retrieval"].RetrievalPolicy
                    policy.__init__.__defaults__ = tuple(original_dns if value is dns else value
                                                        for value in policy.__init__.__defaults__)
                    policy.__dataclass_fields__["resolver"].default = original_dns
                    policy.resolver = original_dns
